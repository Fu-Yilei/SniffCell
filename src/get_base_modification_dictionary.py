import re, pysam, os
os.environ["HTS_LOG_LEVEL"] = "error"

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
        # for i in read_base_ref_loc:
        #     if i in hp_myth_dict.keys():  # O(1) search
        #         hp_myth_dict[i][2] += 1 
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