import numpy as np
import pandas as pd
import pysam
import re


def em_k_cluster_methylation(
    df,
    alpha_init,    # array-like of shape (K,) for initial mixing proportions
    p_init,        # array-like of shape (K,) or (K, num_sites) for init methylation probabilities
    max_iter=100,
    tol=1e-6,
    random_state=42,
    fix_alpha=False,
    fix_p=False,
    drop_na=2,
):
    """
    EM for a K-component mixture of 'Bernoulli-like' distributions on [0,1] methylation calls.
    
    This version:
      - Drops rows that have > drop_na NaNs (if drop_na>0) for the EM calculations,
        but keeps the original DataFrame intact, only adding a 'gamma' column at the end.
      - Stores a single posterior probability (gamma) for each read in df['gamma'].
        By default, we store gamma for cluster 0 (gamma_cluster=0).
        If there's more than 2 clusters, you can pick which cluster's gamma to store.
      - Rows that were dropped (or never used) will have gamma=NaN.

    Parameters
    ----------
    df : pd.DataFrame
        Rows = reads, columns = CpG sites, values in [0,1] or NaN.
    alpha_init : array-like, shape (K,)
        Initial mixing proportions for K cell types/clusters.
        If fix_alpha=False, these can be updated by EM.
    p_init : array-like, shape (K,) or (K, num_sites)
        Initial methylation probabilities per cluster. If shape = (K,),
        we replicate across all CpG sites. If shape = (K, num_sites),
        each cluster k has site-specific probabilities.
    max_iter : int
        Maximum EM iterations.
    tol : float
        Convergence threshold for the change in log-likelihood.
    random_state : int
        For reproducibility.
    fix_alpha : bool
        If True, do NOT update alpha in the M-step.
    fix_p : bool
        If True, do NOT update the per-site probabilities in the M-step.
    drop_na : int or None
        If an integer d>0, drop any read that has more than d NaNs for EM.
        If 0 or None, do not drop any rows.
    gamma_cluster : int
        Which cluster's gamma to store in the final DataFrame column 'gamma'.
        E.g., 0 if you want the posterior probability for cluster #0.

    Returns
    -------
    alpha : np.ndarray, shape (K,)
        The final (or fixed) mixing proportions.
    p : np.ndarray, shape (K, num_sites)
        The final (or fixed) site-specific probabilities for each cluster.
    df_out : pd.DataFrame
        A **copy** of the original df with an added 'gamma' column containing
        the posterior for the chosen cluster. Rows dropped from EM (if any)
        get gamma=NaN.
    """
    np.random.seed(random_state)

    # Make a copy so we don't modify the original in-place
    df_out = df.copy()

    # ---------------------------
    # (A) CREATE EM WORKING DF
    # ---------------------------
    # If drop_na>0, we drop rows that exceed that many NaNs for the EM step only
    if drop_na and drop_na > 0:
        num_sites_all = df_out.shape[1]
        # threshold: require at least num_sites_all - drop_na valid (non-NaN) entries
        thresh = num_sites_all - drop_na
        df_em = df_out.dropna(thresh=thresh, axis=0)
    else:
        # if drop_na=0 or None, don't drop any rows
        df_em = df_out

    # Convert df_em -> numeric arrays
    data = df_em.values  # shape = (num_reads_em, num_sites)
    mask = ~np.isnan(data)
    num_reads_em, num_sites = data.shape

    # Number of clusters
    alpha = np.array(alpha_init, dtype=float)
    K = alpha.shape[0]

    # Initialize p array
    p_init = np.array(p_init, dtype=float)
    if p_init.ndim == 1:
        # shape=(K,) => replicate across all sites
        if p_init.shape[0] != K:
            raise ValueError("p_init must have length K or shape (K, num_sites).")
        p = np.tile(p_init.reshape(K, 1), (1, num_sites))
    else:
        # shape=(K, num_sites)
        if p_init.shape != (K, num_sites):
            raise ValueError(f"p_init must be shape (K,) or (K, {num_sites}).")
        p = p_init.copy()

    # ---------------------------
    # (B) LOG-LIKELIHOOD FUNCTION
    # ---------------------------
    def log_likelihood(i, k_):
        valid_cols = mask[i, :]
        x = data[i, valid_cols]
        pk = p[k_, valid_cols]
        eps = 1e-12
        pk = np.clip(pk, eps, 1 - eps)
        # Bernoulli-like formula for fractional x
        return np.sum(x * np.log(pk) + (1 - x) * np.log(1 - pk))

    # ---------------------------
    # (C) EM LOOP
    # ---------------------------
    prev_log_lik = None
    for iteration in range(max_iter):
        # ===== E STEP =====
        log_resp = np.zeros((num_reads_em, K))
        for i in range(num_reads_em):
            for k_ in range(K):
                log_resp[i, k_] = np.log(alpha[k_]) + log_likelihood(i, k_)

        # Convert to gamma
        max_log = np.max(log_resp, axis=1, keepdims=True)
        resp = np.exp(log_resp - max_log)
        row_sums = np.sum(resp, axis=1, keepdims=True)
        gamma = resp / row_sums  # shape=(num_reads_em, K)

        # ===== M STEP =====
        # Update alpha
        if not fix_alpha:
            N_k = np.sum(gamma, axis=0)
            alpha = N_k / num_reads_em

        # Update p
        if not fix_p:
            for k_ in range(K):
                for s in range(num_sites):
                    valid_rows = mask[:, s]
                    if valid_rows.any():
                        numerator = np.sum(gamma[valid_rows, k_] * data[valid_rows, s])
                        denominator = np.sum(gamma[valid_rows, k_])
                        if denominator > 0:
                            p[k_, s] = numerator / denominator

        # Check log-likelihood for convergence
        total_ll = 0.0
        for i in range(num_reads_em):
            for k_ in range(K):
                total_ll += gamma[i, k_] * (np.log(alpha[k_]) + log_likelihood(i, k_))

        if (prev_log_lik is not None) and (abs(total_ll - prev_log_lik) < tol):
            break
        prev_log_lik = total_ll

    # -----------------------------
    # 5) Store the entire gamma vector in one column
    # -----------------------------
    # Build a Series whose values are lists of length K, one for each row in df_em
    gamma_lists = [ gamma[i, :].tolist() for i in range(num_reads_em) ]
    gamma_series_em = pd.Series(gamma_lists, index=df_em.index, name="gamma")

    # Reindex to df_out, so dropped rows become NaN
    gamma_series_out = gamma_series_em.reindex(df_out.index)

    # Attach as a new column in df_out
    df_out["gamma"] = gamma_series_out

    return alpha, p, df_out



def em_k_cluster_methylation_avg(
    df,
    alpha_init,  # array-like of shape (K,) for initial mixing proportions
    p_init,      # array-like of shape (K,) for initial methylation probabilities
    max_iter=100,
    tol=1e-6,
    random_state=42,
    fix_alpha=False,
    fix_p=False,
    drop_na=2,
):
    """
    EM for a K-component mixture of beta-like distributions on average methylation values.
    
    Instead of treating each CpG site separately, this version **averages methylation per read**.
    This reduces dimensionality and assumes a single probability per read rather than per CpG site.
    
    Parameters
    ----------
    df : pd.DataFrame
        Rows = reads, columns = CpG sites, values in [0,1] or NaN.
    alpha_init : array-like, shape (K,)
        Initial mixing proportions for K cell types/clusters.
    p_init : array-like, shape (K,)
        Initial methylation probabilities per cluster.
    max_iter : int
        Maximum EM iterations.
    tol : float
        Convergence threshold for log-likelihood.
    random_state : int
        For reproducibility.
    fix_alpha : bool
        If True, do NOT update alpha in the M-step.
    fix_p : bool
        If True, do NOT update the cluster probabilities in the M-step.
    drop_na : int or None
        If an integer d>0, drop any read that has more than d NaNs for EM.
        If 0 or None, do not drop any rows.
    
    Returns
    -------
    alpha : np.ndarray, shape (K,)
        The final (or fixed) mixing proportions.
    p : np.ndarray, shape (K,)
        The final (or fixed) probabilities for each cluster.
    df_out : pd.DataFrame
        A **copy** of the original df with an added 'gamma' column containing
        the posterior probability of the most probable cluster for each read.
    """
    np.random.seed(random_state)
    df_out = df.copy()

    # Drop rows with excessive NaNs if specified
    if drop_na and drop_na > 0:
        df_em = df_out.dropna(thresh=df.shape[1] - drop_na, axis=0)
    else:
        df_em = df_out
    
    # Compute mean methylation per read
    avg_meth = df_em.mean(axis=1, skipna=True).values  # shape=(num_reads,)
    num_reads = avg_meth.shape[0]
    
    # Convert inputs to numpy arrays
    alpha = np.array(alpha_init, dtype=float)
    p = np.array(p_init, dtype=float)
    K = alpha.shape[0]
    
    if p.shape[0] != K:
        raise ValueError("p_init must have length K.")
    
    # Log-likelihood function for average methylation
    def log_likelihood(i, k_):
        eps = 1e-12  # Avoid log(0)
        pk = np.clip(p[k_], eps, 1 - eps)
        x = avg_meth[i]
        return x * np.log(pk) + (1 - x) * np.log(1 - pk)
    
    prev_log_lik = None
    for iteration in range(max_iter):
        # ===== E STEP =====
        log_resp = np.zeros((num_reads, K))
        for i in range(num_reads):
            for k_ in range(K):
                log_resp[i, k_] = np.log(alpha[k_]) + log_likelihood(i, k_)

        max_log = np.max(log_resp, axis=1, keepdims=True)
        resp = np.exp(log_resp - max_log)
        row_sums = np.sum(resp, axis=1, keepdims=True)
        gamma = resp / row_sums

        # ===== M STEP =====
        if not fix_alpha:
            N_k = np.sum(gamma, axis=0)
            alpha = N_k / num_reads
        
        if not fix_p:
            for k_ in range(K):
                numerator = np.sum(gamma[:, k_] * avg_meth)
                denominator = np.sum(gamma[:, k_])
                if denominator > 0:
                    p[k_] = numerator / denominator
        
        # Check for convergence
        total_ll = np.sum(gamma * (np.log(alpha) + log_resp))
        if prev_log_lik is not None and abs(total_ll - prev_log_lik) < tol:
            break
        prev_log_lik = total_ll
    
    # Store gamma in df_out
    df_out['gamma'] = pd.Series(np.max(gamma, axis=1), index=df_em.index)
    
    return alpha, p, df_out
