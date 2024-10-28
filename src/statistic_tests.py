from scipy.stats import ranksums
import numpy as np
import scipy

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

def calculate_ranksum_basic(data_dict):
    """
        for single list dict. 
    """
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = values[2]
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
        list1 = values[0]
        list2 = values[1]
        if len(list1) >= 2 and len(list2) >= 2:
            _, p_value = scipy.stats.ttest_ind(list1, list2, equal_var=False)
        else:
            p_value = np.nan 
        result[i] = p_value
    return result

def calculate_ttest_basic(data_dict):
    """
        for single list dict. 
    """
    result = {}
    for i, values in data_dict.items():
        list1 = values[0]
        list2 = values[2]
        if len(list1) >= 2 and len(list2) >= 2:
            _, p_value = scipy.stats.ttest_ind(list1, list2, equal_var=False)
        else:
            p_value = np.nan 
        result[i] = p_value
    return result

