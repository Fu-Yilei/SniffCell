def smoothing_basic(data_dict, key_diff_threshold=10):
    merged_dict = {}
    sorted_keys = sorted(data_dict.keys())  # Sort keys to compare neighbors
    if not sorted_keys:
        return {}
    merged_key = sorted_keys[0]  # Start with the first key
    merged_value = data_dict[merged_key]
    for key in sorted_keys[1:]:
        if abs(key - merged_key) < key_diff_threshold:
            # Merge lists and add numbers
            merged_value = [merged_value[0] + data_dict[key][0], merged_value[1] + data_dict[key][1]]
        else:
            # Store the merged result and reset for the new key
            merged_dict[merged_key] = merged_value
            merged_key = key
            merged_value = data_dict[merged_key]
    merged_dict[merged_key] = merged_value
    return merged_dict
    
def smoothing_comp(data_dict, key_diff_threshold=10):
    if len(data_dict) < key_diff_threshold:
        return data_dict
    elif len(data_dict) == 0:
        return {}
    merged_dict = {}
    sorted_keys = sorted(data_dict.keys())  # Sort keys to compare neighbors
    merged_key = sorted_keys[0]  # Start with the first key
    merged_value = data_dict[merged_key]
    for key in sorted_keys[1:]:
        if abs(key - merged_key) < key_diff_threshold:
            # Merge lists and sum numbers
            merged_value = [
                merged_value[0] + data_dict[key][0],  # Concatenate first list
                merged_value[1] + data_dict[key][1],  # Sum num1
            ]
        else:
            merged_dict[merged_key] = merged_value
            merged_key = key
            merged_value = data_dict[merged_key]
    merged_dict[merged_key] = merged_value
    return merged_dict
