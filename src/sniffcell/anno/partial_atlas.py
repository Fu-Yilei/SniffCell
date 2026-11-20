"""Consensus for a shared schema with explicitly unprofiled assay-specific classes."""
import pandas as pd


def consensus_with_partial_profiles(assignment, threshold, complete_consensus):
    if 'unprofiled_celltypes' not in assignment:
        result = complete_consensus(assignment, threshold)
        result['assay_consensus_rule'] = 'all_assays_fully_profiled'
        return result
    unknown = assignment.unprofiled_celltypes.fillna('').astype(str).str.strip()
    if not unknown.ne('').any():
        result = complete_consensus(assignment, threshold)
        result['assay_consensus_rule'] = 'all_assays_fully_profiled'
        return result
    assert set(unknown[unknown.ne('')]) == {'Granulocyte'}
    schemas = assignment.code_order.dropna().unique()
    assert len(schemas) == 1, 'Partial-profile consensus requires one shared class schema'
    schema = schemas[0]
    labels = schema.split('|')
    gi = labels.index('Granulocyte')
    # Deconv drops the schema prefix when every row shares the same order.
    # Restore it locally so partial-profile rules behave identically in deconv
    # and annotation, while leaving the caller's evidence table unchanged.
    assignment = assignment.copy(deep=False)
    assignment['code_token'] = schema + '::' + assignment.code.astype(str).str.zfill(len(labels))
    names = assignment.index.astype(str)
    complete = assignment.loc[unknown.eq('')]
    partial = assignment.loc[unknown.ne('')]
    assert set(partial.modification) <= {'5mC', '5hmC'}
    assert complete.empty or set(complete.modification) == {'modifiedC'}
    # Never interpret a missing granulocyte m/h mean as a negative profile.
    if complete.empty:
        base = pd.DataFrame(columns=['readname', 'read_consensus_code', 'read_is_mixed'])
    else:
        base = complete_consensus(complete, threshold)
    def excludes_granulocyte(token):
        if not isinstance(token, str) or '::' not in token:
            return False
        bits = token.split('::', 1)[1]
        return len(bits) == len(labels) and bits[gi] == '0' and '1' in bits
    eligible = base.loc[~base.read_is_mixed.astype(bool) & base.read_consensus_code.map(excludes_granulocyte), 'readname']
    results = []
    retained = base.loc[~base.readname.isin(eligible)].copy()
    if not retained.empty:
        retained['assay_consensus_rule'] = 'modifiedC_only_granulocyte_unresolved_or_primary_mixed'
        results.append(retained)
    if len(eligible):
        combined = complete_consensus(assignment.loc[names.isin(eligible)], threshold)
        combined['assay_consensus_rule'] = 'combined_assays_after_modifiedC_granulocyte_exclusion'
        results.append(combined)
    # Reads with m/h evidence alone retain Granulocyte as an unresolved option.
    # First resolve the measured PBMC classes; conflicting PBMC evidence cannot
    # turn into a false singleton Granulocyte call by intersection.
    partial_only = partial.loc[~partial.index.astype(str).isin(base.readname)]
    if not partial_only.empty:
        result = complete_consensus(partial_only, threshold)
        def add_unknown(token):
            if not isinstance(token, str) or '::' not in token:
                return token
            prefix, bits = token.split('::', 1)
            assert prefix == schema and bits[gi] == '0'
            if '1' not in bits:
                return None
            bits = list(bits); bits[gi] = '1'
            return prefix + '::' + ''.join(bits)
        result['read_consensus_code'] = result.read_consensus_code.map(add_unknown)
        result['assay_consensus_rule'] = 'mh_only_granulocyte_unresolved'
        results.append(result)
    result = pd.concat(results, ignore_index=True)
    assert result.readname.is_unique
    assert set(result.readname) == set(names)
    return result
