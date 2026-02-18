import re, pysam, os
os.environ["HTS_LOG_LEVEL"] = "error"
import pandas as pd
def get_base_modification_dictionary(
    bam_file, ref_seq, chromosome, phase_region, sv_supporting_reads
):
    """
        Return value: a dictionary that contains cpg location and its haplotype related base modification score
    """
    methylation_identifier_0 = ('C', 0, 'm') # We only care about 5mc!!!
    methylation_identifier_1 = ('C', 1, 'm')
    phase_region_start = phase_region[0]
    phase_region_end = phase_region[1]

    phased_block_ref = ref_seq.fetch(
        chromosome, phase_region_start, phase_region_end)
    cg_loc = [
        m.start(0) for m in re.finditer("CG", str(phased_block_ref))
    ]  # Use regular expression to find all CpG locations on the reference
    # record G locatioon of 'CG's
    cg_loc = [x + phase_region_start + 1 for x in cg_loc]
    hp_myth_dict = dict()
    """
        Data structure:
        {i:[[list1], [list2], num_1, num_2]} 
            i: CpG locations
            list1: ML score for SV-supporting reads
            list2: ML score for non SV-supporting reads
            num1: coverage for SV-supporting reads
            num2: coverage for non SV-supporting reads
        Use dictionary so that when querying CpG locations the time complexity is O(1)
    """
    for i in cg_loc:
        # build the dictionary
        hp_myth_dict.update({i: [[],[],0,0]})
    phased_block_alignment = bam_file.fetch(
        chromosome, phase_region[0], phase_region[1], multiple_iterators=True
    )

    for reads in phased_block_alignment:
        read_base_ref_loc = reads.get_reference_positions(full_length=True)  
        # use full_length=True or the positions won't match
        mm = (reads.modified_bases)
        # mm is a dictionary that contains {score type: [(location, score)]}. score is 255-based
        if (mm != -1) and (mm != {}):  # update base modification scores list
            if methylation_identifier_0 in list(mm.keys()):
                methylation_identifier = methylation_identifier_0
            elif methylation_identifier_1 in list(mm.keys()):
                methylation_identifier = methylation_identifier_1
            else:
                continue
            for i in mm[methylation_identifier]:  # Remora only output one type of score: c 1 m/c 0 m, but this part can be improved for other methlyation callers
                if read_base_ref_loc[i[0]]:  # i format: (loc, score)
                    if reads.is_forward:  # cg/gc on forward and reverse reads
                        mm_ref_loc = read_base_ref_loc[i[0]] + 1
                    else:
                        mm_ref_loc = read_base_ref_loc[i[0]]
                    if mm_ref_loc in hp_myth_dict.keys():                       
                        modification_chance = i[1]  # 0 - 255 base d
                        if reads.query_name in sv_supporting_reads:
                            hp_myth_dict[mm_ref_loc][0].append(modification_chance)
                            hp_myth_dict[mm_ref_loc][2] += 1 
                        else:
                            hp_myth_dict[mm_ref_loc][1].append(modification_chance)
                            hp_myth_dict[mm_ref_loc][3] += 1 

    return hp_myth_dict




def get_base_modification_dictionary_new_bam(
    bam_file, ref_seq, chromosome, phase_region, sv_supporting_reads, sv_id, output_bam_folder, output_bam
):
    """
        output BAM file with SV supporting reads as a new read group (Now optional)
         Return value: a dictionary that contains cpg location and its haplotype related base modification score
    """
    methylation_identifier_0 = ('C', 0, 'm') # We only care about 5mc!!!
    methylation_identifier_1 = ('C', 1, 'm')
    phase_region_start = phase_region[0]
    phase_region_end = phase_region[1]
    output_bam_file = os.path.join(output_bam_folder, f"{sv_id}.bam")
    phased_block_ref = ref_seq.fetch(
        chromosome, phase_region_start, phase_region_end)
    cg_loc = [
        m.start(0) for m in re.finditer("CG", str(phased_block_ref))
    ]  # Use regular expression to find all CpG locations on the reference
    # record G locatioon of 'CG's
    cg_loc = [x + phase_region_start + 1 for x in cg_loc]
    hp_myth_dict = dict()
    """
        Data structure:
        {i:[[list1], [list2], num_1, num_2]} 
            i: CpG locations
            list1: ML score for SV-supporting reads
            list2: ML score for non SV-supporting reads
            num1: coverage for SV-supporting reads
            num2: coverage for non SV-supporting reads
        Use dictionary so that when querying CpG locations the time complexity is O(1)
    """
    for i in cg_loc:
        # build the dictionary
        hp_myth_dict.update({i: [[],[],0,0]})
    phased_block_alignment = bam_file.fetch(
        chromosome, phase_region[0], phase_region[1], multiple_iterators=True
    )
    if output_bam:
        with pysam.AlignmentFile(output_bam_file, "wb", header=bam_file.header) as outfile:
            for reads in phased_block_alignment:
                if reads.query_name in sv_supporting_reads:
                    reads.set_tag("RG", "sv_reads")
                else:
                    reads.set_tag("RG", "non_sv_reads")
                outfile.write(reads)
        pysam.index(output_bam_file)
    # print(output_bam)
    phased_block_alignment = bam_file.fetch(
        chromosome, phase_region[0], phase_region[1], multiple_iterators=True
    )
    for reads in phased_block_alignment:
        if not reads.is_secondary and not reads.is_supplementary:
            read_base_ref_loc = reads.get_reference_positions(full_length=True)  
            # use full_length=True or the positions won't match
            mm = (reads.modified_bases)
            # mm is a dictionary that contains {score type: [(location, score)]}. score is 255-based
            if (mm != -1) and (mm != {}):  # update base modification scores list
                # print(reads.query_name)
                if methylation_identifier_0 in list(mm.keys()):
                    methylation_identifier = methylation_identifier_0
                elif methylation_identifier_1 in list(mm.keys()):
                    methylation_identifier = methylation_identifier_1
                else:
                    continue
                for i in mm[methylation_identifier]:  # Remora only output one type of score: c 1 m/c 0 m, but this part can be improved for other methlyation callers
                    if read_base_ref_loc[i[0]]:  # i format: (loc, score)
                        if reads.is_forward:  # cg/gc on forward and reverse reads
                            mm_ref_loc = read_base_ref_loc[i[0]] + 1
                        else:
                            mm_ref_loc = read_base_ref_loc[i[0]]
                        if mm_ref_loc in hp_myth_dict.keys():                       
                            modification_chance = i[1]  # 0 - 255 base d
                            if reads.query_name in sv_supporting_reads:
                                # print(reads.query_name)
                                hp_myth_dict[mm_ref_loc][0].append(modification_chance)
                                hp_myth_dict[mm_ref_loc][2] += 1 
                            else:
                                hp_myth_dict[mm_ref_loc][1].append(modification_chance)
                                hp_myth_dict[mm_ref_loc][3] += 1 
                        
    return hp_myth_dict



def get_base_modification_dictionary_basic_supporting_reads(
    bam_file, ref_seq, sv_supporting_reads, chromosome, phase_region, sv_id, output_bam_folder, output_bam
):
    """
        _basic means for single list data structure
        Return value: a dictionary that contains cpg location and its haplotype related base modification score
    """
    methylation_identifier_0 = ('C', 0, 'm') # We only care about 5mc!!!
    methylation_identifier_1 = ('C', 1, 'm')
    phase_region_start = phase_region[0]
    phase_region_end = phase_region[1]

    phased_block_ref = ref_seq.fetch(
        chromosome, phase_region_start, phase_region_end)
    cg_loc = [
        m.start(0) for m in re.finditer("CG", str(phased_block_ref))
    ]  # Use regular expression to find all CpG locations on the reference
    # record G locatioon of 'CG's
    cg_loc = [x + phase_region_start + 1 for x in cg_loc]
    hp_myth_dict = dict()
    """
        Data structure:
        {i:[[], num_1]} 
            i: CpG locations
            list1: ML score for SV-supporting reads
            num1: coverage for SV-supporting reads
        Use dictionary so that when querying CpG locations the time complexity is O(1)
    """
    for i in cg_loc:
        # build the dictionary
        hp_myth_dict.update({i: [[],0]})
    phased_block_alignment = bam_file.fetch(
        chromosome, phase_region[0], phase_region[1], multiple_iterators=True
    )
    output_bam_file = os.path.join(output_bam_folder, f"{sv_id}.bam")
    if output_bam:
        with pysam.AlignmentFile(output_bam_file, "wb", header=bam_file.header) as outfile:
            for reads in phased_block_alignment:
                if reads.query_name in sv_supporting_reads:
                    reads.set_tag("RG", "sv_reads")
                else:
                    reads.set_tag("RG", "non_sv_reads")
                outfile.write(reads)
        pysam.index(output_bam_file)
    phased_block_alignment = bam_file.fetch(
        chromosome, phase_region[0], phase_region[1], multiple_iterators=True
    )
    for reads in phased_block_alignment:
        if not reads.is_secondary and not reads.is_supplementary:
            read_base_ref_loc = reads.get_reference_positions(full_length=True)  
            # use full_length=True or the positions won't match
            mm = (reads.modified_bases)
            # mm is a dictionary that contains {score type: [(location, score)]}. score is 255-based
            if (mm != -1) and (mm != {}):  # update base modification scores list
                if methylation_identifier_0 in list(mm.keys()):
                    methylation_identifier = methylation_identifier_0
                elif methylation_identifier_1 in list(mm.keys()):
                    methylation_identifier = methylation_identifier_1
                else:
                    continue
                for i in mm[methylation_identifier]:  # Remora only output one type of score: c 1 m/c 0 m, but this part can be improved for other methlyation callers
                    if read_base_ref_loc[i[0]]:  # i format: (loc, score)
                        if reads.is_forward:  # cg/gc on forward and reverse reads
                            mm_ref_loc = read_base_ref_loc[i[0]] + 1
                        else:
                            mm_ref_loc = read_base_ref_loc[i[0]]
                        if mm_ref_loc in hp_myth_dict.keys():                       
                            modification_chance = i[1]  # 0 - 255 base d
                            if sv_supporting_reads is None or reads.query_name in sv_supporting_reads:
                                hp_myth_dict[mm_ref_loc][0].append(modification_chance)
                                hp_myth_dict[mm_ref_loc][1] += 1 
    return hp_myth_dict



def get_base_modification_dictionary_basic_supporting_reads_df(
    bam_file: pysam.AlignmentFile,
    ref_seq: pysam.FastaFile,
    chromosome: str,
    phase_region: tuple[int, int],
    min_read_length: int = 0  # filter out reads shorter than this
) -> pd.DataFrame:
    """
    Build a matrix (rows = reads, columns = CpG genomic positions) holding the
    per–read 5‑methyl‑cytosine probabilities within a phased region,
    excluding reads with supplementary alignments or shorter than min_read_length.

    Parameters
    ----------
    bam_file : pysam.AlignmentFile
        BAM/CRAM opened in read mode with `modified_bases` tags.
    ref_seq : pysam.FastaFile
        Reference FASTA opened with `fetch`.
    chromosome : str
        Chromosome name (must match both BAM and FASTA).
    phase_region : (int, int)
        (start, end) 0‑based half‑open coordinates of the phased block.
    min_read_length : int
        Minimum read length to include (default = 0 = no filtering).

    Returns
    -------
    pd.DataFrame
        Index: MultiIndex ``(read_name, haplotype)``  
        Columns: genomic positions of the **C** in every CpG in the region  
        Values: probability of 5‑mC (0–1 float); missing → NaN
    """
    import re
    import numpy as np
    import pandas as pd

    start, end = phase_region

    # --------------------------------------------------------------------- #
    # 1. Identify every CpG *C* position in the reference segment
    # --------------------------------------------------------------------- #
    block_ref = ref_seq.fetch(chromosome, start, end)
    cpg_c_sites = [
        m.start() + start
        for m in re.finditer(r"CG", block_ref)
    ]
    if not cpg_c_sites:
        return pd.DataFrame()

    # --------------------------------------------------------------------- #
    # 2a. Collect reads with supplementary alignments
    # --------------------------------------------------------------------- #
    supplementary_reads = set(
        read.query_name
        for read in bam_file.fetch(chromosome, start, end, multiple_iterators=True)
        if read.is_supplementary
    )

    # --------------------------------------------------------------------- #
    # 2b. Collect 5‑mC calls for reads not in supplementary set and long enough
    # --------------------------------------------------------------------- #
    rows, row_names, hps = [], [], []
    col_template = {site: np.nan for site in cpg_c_sites}
    wanted_bases = {('C', 0, 'm'), ('C', 1, 'm')}

    for read in bam_file.fetch(chromosome, start, end, multiple_iterators=True):
        if (
            read.query_name in supplementary_reads or
            read.is_secondary or
            read.is_supplementary or
            read.query_length < min_read_length
        ):
            continue

        mm_dict = read.modified_bases
        if not mm_dict:
            continue

        row = col_template.copy()
        ref_positions = read.get_reference_positions(full_length=True)

        for m_key, mods in mm_dict.items():
            if m_key not in wanted_bases:
                continue
            for read_pos, mod_score in mods:
                ref_pos = ref_positions[read_pos]
                if ref_pos is None:
                    continue
                if ref_pos in row:
                    row[ref_pos] = mod_score / 255.0

        if not np.all(pd.isna(list(row.values()))):
            try:
                hp_tag = read.get_tag("HP")
            except (KeyError, AttributeError):
                hp_tag = -1
            rows.append(row)
            row_names.append(read.query_name)
            hps.append(hp_tag)

    # --------------------------------------------------------------------- #
    # 3. Assemble DataFrame
    # --------------------------------------------------------------------- #
    if rows:
        df = pd.DataFrame(rows, columns=cpg_c_sites)
        df.index = pd.MultiIndex.from_arrays(
            [row_names, hps],
            names=["read_name", "haplotype"],
        )
        return df
    else:
        return pd.DataFrame(columns=cpg_c_sites)
