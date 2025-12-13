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
    fix_alpha=True,
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


def em_on_haplotypes_and_harmonize(df, alpha_init, p_init, **em_kwargs):
    """
    Run EM on each haplotype (HP) cluster in a MultiIndex DataFrame, then harmonize gamma values across haplotypes.
    If only haplotype 0 is present, run standard EM on all data.
    If haplotypes 1 and 2 are present, run EM only on 1 and 2, harmonize, and assign to all (including 0 if present).
    Returns alpha, p, and a DataFrame with a single 'gamma' column (harmonized, as a list per row),
    indexed by (read_name, haplotype), matching the output format of em_k_cluster_methylation.
    """
    # Ensure MultiIndex with 'read_name' and 'haplotype'
    if not isinstance(df.index, pd.MultiIndex) or df.index.names != ["read_name", "haplotype"]:
        raise ValueError("Input DataFrame must have MultiIndex with ['read_name', 'haplotype'].")

    hps = set(df.index.get_level_values('haplotype').unique())
    hps_int = set([int(h) if str(h).isdigit() else h for h in hps])

    # If only haplotype 0, run standard EM
    if hps_int == {0}:
        alpha, p, df_em = em_k_cluster_methylation(df, alpha_init, p_init, **em_kwargs)
        return alpha, p, df_em[["gamma"]]

    # Otherwise, run EM on haplotypes 1 and 2 only
    em_results = {}
    alpha_last = None
    p_last = None
    for hp in [1, 2]:
        if hp in hps_int:
            df_hp = df.xs(hp, level='haplotype', drop_level=False)
            alpha, p, df_em = em_k_cluster_methylation(df_hp.drop(columns=['gamma'], errors='ignore'), alpha_init, p_init, **em_kwargs)
            em_results[hp] = df_em[['gamma']]
            alpha_last = alpha
            p_last = p

    # Concatenate results for 1 and 2
    gamma_concat = pd.concat(em_results.values(), keys=em_results.keys(), names=['haplotype'])
    gamma_concat = gamma_concat.reset_index(level=0).sort_index()

    # Harmonize: average gamma for each read_name across haplotypes 1 and 2
    gamma_lists = gamma_concat.groupby('read_name')['gamma'].apply(list)
    # Average each cluster's gamma across haplotypes
    def avg_gamma_lists(gamma_list_of_lists):
        # gamma_list_of_lists: list of lists (one per haplotype)
        arr = np.array(gamma_list_of_lists)
        return arr.mean(axis=0).tolist()
    gamma_harmonized = gamma_lists.apply(avg_gamma_lists)

    # Assign harmonized gamma back to all (read_name, haplotype) pairs in a new DataFrame
    gamma_df = pd.DataFrame(index=df.index)
    gamma_df['gamma'] = df.index.get_level_values('read_name').map(gamma_harmonized)

    return alpha_last, p_last, gamma_df


def em_haplotype_and_combined(df, alpha_init, p_init, **em_kwargs):
    """
    Run the Expectation-Maximization (EM) algorithm on haplotype-specific and combined data.

    This function performs the EM algorithm three times:
    1. On data corresponding to haplotype 1.
    2. On data corresponding to haplotype 2.
    3. On all data combined, ignoring haplotype information.

    For each read, it stores a list of gamma vectors from each run (or a single value if only one run applies).

    Args:
        df (pd.DataFrame): Input DataFrame with a MultiIndex of ['read_name', 'haplotype'].
                           The DataFrame should contain methylation data for clustering.
        alpha_init (np.ndarray): Initial alpha values for the EM algorithm.
        p_init (np.ndarray): Initial p values for the EM algorithm.
        **em_kwargs: Additional keyword arguments to pass to the `em_k_cluster_methylation` function.

    Returns:
        tuple:
            - alphas (list): A list of up to 3 alpha arrays (one for each EM run: haplotype 1, haplotype 2, and combined).
                             Each alpha array is of shape (K,), where K is the number of clusters.
            - ps (list): A list of up to 3 p arrays (one for each EM run: haplotype 1, haplotype 2, and combined).
                         Each p array is of shape (K, num_sites) for site-specific probabilities.
            - gamma_df (pd.DataFrame): A DataFrame with the same index as the input `df` and a column 'gamma_list',
                                       which contains a list of gamma vectors for each read (or a single value if only one run applies).
                                       The structure of the 'gamma_list' column is as follows:
                                       - If only one EM run applies to a read, the column contains a single gamma vector (list of floats).
                                       - If multiple EM runs apply, the column contains a list of gamma vectors (one per run).

    Raises:
        ValueError: If the input DataFrame does not have a MultiIndex with ['read_name', 'haplotype'].

    Notes:
        - The function checks for the presence of haplotype 1 and haplotype 2 in the input data.
          If a haplotype is not present, the corresponding alpha and p values will be `None`.
        - The `gamma_list` column in the output `gamma_df` contains gamma values from each EM run.
          If only one non-None gamma value exists for a read, it is stored directly; otherwise, a list of gamma values is stored.
    """
    if not isinstance(df.index, pd.MultiIndex) or df.index.names != ["read_name", "haplotype"]:
        raise ValueError("Input DataFrame must have MultiIndex with ['read_name', 'haplotype'].")

    # Prepare outputs
    alphas = []
    ps = []
    gamma_dicts = [None, None, None]  # To store gamma dictionaries for each haplotype and combined
    # 1. EM on haplotype 1
    phased_flag = False
    if 1 in set([int(h) if str(h).isdigit() else h for h in df.index.get_level_values('haplotype').unique()]):
        phased_flag = True
        df_h1 = df.xs(1, level='haplotype', drop_level=False)
        alpha1, p1, df_em1 = em_k_cluster_methylation(df_h1.drop(columns=['gamma_list', 'gamma'], errors='ignore'), alpha_init, p_init, **em_kwargs)
        alphas.append(alpha1)
        ps.append(p1)
        gamma_dicts[0] = df_em1['gamma'].to_dict()
    else:
        alphas.append(None)
        ps.append(None)
        gamma_dicts[0] = {}

    # 2. EM on haplotype 2
    if 2 in set([int(h) if str(h).isdigit() else h for h in df.index.get_level_values('haplotype').unique()]):
        phased_flag = True
        df_h2 = df.xs(2, level='haplotype', drop_level=False)
        alpha2, p2, df_em2 = em_k_cluster_methylation(df_h2.drop(columns=['gamma_list', 'gamma'], errors='ignore'), alpha_init, p_init, **em_kwargs)
        alphas.append(alpha2)
        ps.append(p2)
        gamma_dicts[1] = df_em2['gamma'].to_dict()
    else:
        alphas.append(None)
        ps.append(None)
        gamma_dicts[1] = {}

    # 3. EM on all data (ignore haplotype)
    df_all = df.copy()
    alpha_all, p_all, df_em_all = em_k_cluster_methylation(df_all.drop(columns=['gamma_list', 'gamma'], errors='ignore'), alpha_init, p_init, **em_kwargs)
    alphas.append(alpha_all)
    ps.append(p_all)
    gamma_dicts[2] = df_em_all['gamma'].to_dict()



    # For each row in the original df, collect gamma from each run (or None if not present)
    gamma_list_col = []
    for idx in df.index:
        gammas = []
        for d in gamma_dicts:
            gammas.append(d.get(idx, None))
        gamma_list_col.append(gammas)

    gamma_df = pd.DataFrame(index=df.index)
    gamma_df['gamma_list'] = gamma_list_col
    return alphas, ps, gamma_df,  # return sizes of each haplotype and combined data
