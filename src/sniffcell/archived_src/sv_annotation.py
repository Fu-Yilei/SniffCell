from collections import defaultdict


def check_clusters_by_haplotype(read_dict):
    haplotype_clusters = defaultdict(set)

    # collect clusters per haplotype
    for read_id, info in read_dict.items():
        hap = info['haplotype']
        cluster = info['cluster']
        haplotype_clusters[hap].add(cluster)

    # check if each haplotype is homogeneous in cluster assignment
    results = {}
    for hap, clusters in haplotype_clusters.items():
        results[hap] = len(clusters) == 1
    return results

def sv_annotation_from_intermediate(intermediate_df):
    intermediate_df['supporting_reads_in_same_cluster_status'] = intermediate_df['supporting_reads_in_same_cluster'].apply(check_clusters_by_haplotype)
    return intermediate_df