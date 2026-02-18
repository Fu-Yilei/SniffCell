from scipy.stats import ranksums
import numpy as np
import scipy
import warnings

def calculate_ranksum(data_dict):
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = values[1]
        # Perform ranksum test but remove the warning 
        if len(list1) >= 2 and len(list2) >= 2:
            _, p_value = ranksums(list1, list2)
        else:
            p_value = np.nan 
        # Store the p-value in the result dictionary
        result[i] = p_value
    return result

def calculate_ranksum_basic(data_dict, data_dict2):
    """
        for single list dict. 
    """
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = data_dict2[i][0]
        # Perform ranksum test but remove the warning 
        if len(list1) >= 2 and len(list2) >= 2:
            _, p_value = ranksums(list1, list2)
        else:
            p_value = np.nan 
        # Store the p-value in the result dictionary
        result[i] = p_value
    return result


def calculate_ttest(data_dict):
    result = {}

    for i, values in data_dict.items():
        list1 = np.array(values[0])
        list2 = np.array(values[1])

        # Ensure there are at least 2 values in each list
        if len(list1) >= 2 and len(list2) >= 2:
            # Compute variance to check for near-identical data
            var1, var2 = np.var(list1), np.var(list2)

            if var1 < 1e-10 and var2 < 1e-10:
                result[i] = 1  # No meaningful difference, assign p=1
                continue

            try:
                # Suppress warnings inside the try block
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    
                    # Run Welch’s t-test (handles unequal variance)
                    _, p_value = scipy.stats.ttest_ind(list1, list2, equal_var=False, nan_policy='omit')

                # If precision loss is suspected (very low p-value or NaN), switch to Mann-Whitney U test
                if np.isnan(p_value) or p_value < 1e-300:
                    _, p_value = scipy.stats.mannwhitneyu(list1, list2, alternative='two-sided')

            except RuntimeWarning:
                # If catastrophic cancellation occurs, add small noise (jitter) and retry t-test
                list1 += np.random.normal(0, 1e-10, size=list1.shape)
                list2 += np.random.normal(0, 1e-10, size=list2.shape)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    _, p_value = scipy.stats.ttest_ind(list1, list2, equal_var=False, nan_policy='omit')

        else:
            p_value = np.nan  # Not enough data points
        
        result[i] = p_value

    return result

def calculate_ttest_basic(data_dict, data_dict2):
    """
        for single list dict. 
    """
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = data_dict2[i][0]
        if len(list1) >= 2 and len(list2) >= 2:
            _, p_value = scipy.stats.ttest_ind(list1, list2, equal_var=False)
        else:
            p_value = np.nan 
        result[i] = p_value
    return result

def calculate_fisher(data_dict):
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = values[1]
        if len(list1) >= 2 and len(list2) >= 2:
            methylated_list1 = sum(1 for x in list1 if x > 0.66 * 255)
            unmethylated_list1 = sum(1 for x in list1 if x < 0.33 * 255)
            methylated_list2 = sum(1 for x in list2 if x > 0.66 * 255)
            unmethylated_list2 = sum(1 for x in list2 if x < 0.33 * 255)
            contingency_table = [[methylated_list1, unmethylated_list1], [methylated_list2, unmethylated_list2]]
            _, p_value = scipy.stats.fisher_exact(contingency_table)
        else:
            p_value = np.nan 
        result[i] = p_value
    return result

def calculate_fisher_basic(data_dict, data_dict2):
    """
        for single list dict. 
    """
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = data_dict2[i][0]
        if len(list1) >= 2 and len(list2) >= 2:
            methylated_list1 = sum(1 for x in list1 if x > 0.66 * 255)
            unmethylated_list1 = sum(1 for x in list1 if x < 0.33 * 255)
            methylated_list2 = sum(1 for x in list2 if x > 0.66 * 255)
            unmethylated_list2 = sum(1 for x in list2 if x < 0.33 * 255)
            contingency_table = [[methylated_list1, unmethylated_list1], [methylated_list2, unmethylated_list2]]
            _, p_value = scipy.stats.fisher_exact(contingency_table)
        else:
            p_value = np.nan 
        result[i] = p_value
    return result


def calculate_chi2(data_dict):
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = values[1]
        if len(list1) >= 2 and len(list2) >= 2:
            methylated_list1 = sum(1 for x in list1 if x > 0.66 * 255)
            unmethylated_list1 = sum(1 for x in list1 if x < 0.33 * 255)
            methylated_list2 = sum(1 for x in list2 if x > 0.66 * 255)
            unmethylated_list2 = sum(1 for x in list2 if x < 0.33 * 255)
            contingency_table = [[methylated_list1, unmethylated_list1], [methylated_list2, unmethylated_list2]]
            if all(all(cell > 0 for cell in row) for row in contingency_table):
                _, p_value, _, _ = scipy.stats.chi2_contingency(contingency_table)
            else:
                p_value = np.nan
        else:
            p_value = np.nan 
        result[i] = p_value
    return result

def calculate_chi2_basic(data_dict, data_dict2):
    """
        for single list dict. 
    """
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = data_dict2[i][0]
        if len(list1) >= 2 and len(list2) >= 2:
            methylated_list1 = sum(1 for x in list1 if x > 0.66 * 255)
            unmethylated_list1 = sum(1 for x in list1 if x < 0.33 * 255)
            methylated_list2 = sum(1 for x in list2 if x > 0.66 * 255)
            unmethylated_list2 = sum(1 for x in list2 if x < 0.33 * 255)
            contingency_table = [[methylated_list1, unmethylated_list1], [methylated_list2, unmethylated_list2]]
            if all(all(cell > 0 for cell in row) for row in contingency_table):
                _, p_value, _, _ = scipy.stats.chi2_contingency(contingency_table)
            else:
                p_value = np.nan        
        else:
            p_value = np.nan 
        result[i] = p_value
    return result

def calculate_mannwhitneyu(data_dict):
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = values[1]
        if len(list1) >= 2 and len(list2) >= 2:
            _, p_value = scipy.stats.mannwhitneyu(list1, list2)
        else:
            p_value = np.nan 
        result[i] = p_value
    return result

def calculate_mannwhitneyu_basic(data_dict, data_dict2):
    """
        for single list dict. 
    """
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = data_dict2[i][0]
        if len(list1) >= 2 and len(list2) >= 2:
            _, p_value = scipy.stats.mannwhitneyu(list1, list2)
        else:
            p_value = np.nan 
        result[i] = p_value
    return result