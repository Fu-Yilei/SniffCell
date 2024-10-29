# SniffMeth
### Sniffing DNA methylation changes aound (Mosaic) structural variations. 


Usage:

```
usage: SniffMeth [-h] -b BAM -v VCF -r REFERENCE -o OUTPUT [-t THREADS] [-ob OUTPUT_BAM] [-s SMOOTHING] [-i INTERVAL] [-m MIN_SUPPORTING] [-th THRESHHOLD]
                 [-b2 SECOND_BAM] [-v2 SECOND_VCF]

Sniffing CpG methylaiton changed around a (Mosaic) SV

options:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads, default 1
  -ob OUTPUT_BAM, --output_bam OUTPUT_BAM
                        Output SV-related BAM file for each processed SVs. Default: True
  -s SMOOTHING, --smoothing SMOOTHING
                        Enable smoothing, which consider neighboring [input] bps' CpG as an unit. Default: 0 (no smoothing)
  -i INTERVAL, --interval INTERVAL
                        Inverval for checking methylation changes around a SV.
  -m MIN_SUPPORTING, --min_supporting MIN_SUPPORTING
                        Minimum supporting reads requirement for a SV.
  -th THRESHHOLD, --threshhold THRESHHOLD
                        Threshold for determining nearby region is differently methylated or not
  -tf TEST_FUNCTION, --test_function TEST_FUNCTION
                        Statistical test function for determining methylation difference. Default: ttest. Other options: ranksum, fisher
  -b2 SECOND_BAM, --second_bam SECOND_BAM
                        Another BAM file that used for benchmarking.
  -v2 SECOND_VCF, --second_vcf SECOND_VCF
                        Another SV VCF file that used for benchmarking. No SV-supporting reads needed.
  

Required arguments:
  -b BAM, --bam BAM     Input BAM file.
  -v VCF, --vcf VCF     Input SV VCF file (supporting reads needed).
  -r REFERENCE, --reference REFERENCE
                        Reference genome
  -o OUTPUT, --output OUTPUT
                        Output directory

Version 0.1
```

Now only works in single thread. 

Output:
- SV BAM file: Each BAM file contains a SV candidate with supporting reads and other reads in different read groups. This can help the user to visulize DNA methylation differences in IGV.
- CSV file: A CSV file that shows how different the DNA methylation in SV-supporting reads VS other reads. 