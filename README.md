# SniffCell

![SniffCell Workflow](./img/workflow.png)

### Identifying cell type specific SV from long-read bulk only data

SniffCell is a tool designed to analyze DNA methylation changes associated with structural variations (SVs), including mosaic SVs. It processes primary alignments from BAM files and provides detailed outputs for visualization and analysis.


---

## Usage

Run SniffCell using the following command:

```bash
sniffcell [-h] [-t THREADS] -b BAM -v VCF -r REFERENCE -o OUTPUT [-d] [-a ATLAS] [-c TISSUE] [-n REGION_NUMBER] [-me METHOD] [-vb] [-ot OUTLIER_THRESHOLD] [-conf CONFIDENCE] [-nc N_CLOSEST] [-wuoff] [-wgbs WGBS_PATH] [-uxm UXM_PATH] [-ua UXM_ATLAS]
      [-dis HP_DISTANCE] [-rf ASSIGNED_READ_FRACTION] [-ob] [-p] [-s SMOOTHING] [-i INTERVAL] [-m MIN_SUPPORTING] [-th THRESHOLD] [-tf TEST_FUNCTION] [-b2 SECOND_BAM]
```

### Required Arguments:
- `-b, --bam`: Input BAM file.
- `-v, --vcf`: Input SV VCF file (requires supporting reads).
- `-r, --reference`: Reference genome.
- `-o, --output`: Output directory (default: ).

### Optional Arguments:
- `-t, --threads`: Number of threads (default: 1).
- `-d, --deconv`: Enable deconvolution for cell type-specific methylation analysis.

#### Deconvolution-Specific Options:
- `-a, --atlas`: Path to the cell type atlas file (default: ).
- `-c, --tissue`: Tissue type (default: `brain_cereb`). Required if `-d` is specified.
- `-n, --region_number`: Number of regions to select (default: 300).
- `-me, --method`: Region selection method (`std` or `diff`, default: `diff`).
- `-ot, --outlier_threshold`: Threshold for filtering out unreliable regions (default: 0.8).
- `-conf, --confidence`: Minimum confidence threshold for the EM algorithm (default: 0.9).
- `-nc, --n_closest`: Number of closest methylation-informative regions for proportion estimation.

#### General Options:
- `-ob, --output_bam`: Output SV-related BAM file.
- `-p, --primary_only`: Retain only primary alignments in the BAM file.
- `-s, --smoothing`: Enable smoothing (default: 0, no smoothing).
- `-i, --interval`: Interval for checking methylation changes (default: ±1000).
- `-m, --min_supporting`: Minimum supporting reads for SV (default: 3).
- `-th, --threshold`: Threshold for differentially methylated regions (default: 0.2).
- `-tf, --test_function`: Statistical test function (default: `ttest`).
- `-b2, --second_bam`: Benchmark BAM file for comparison.

---

## Statistical Test Recommendations
- **t-test**: General-purpose test.
- **Mann-Whitney U**: Conservative approach.
- **Ranksum, chi2, Fisher's exact**: Use with caution when comparison lists differ in size.

---

## Outputs
1. **SV BAM File**: Contains SV candidates with supporting and other reads in separate read groups. Useful for visualizing methylation differences in IGV.
2. **CSV File**: Summarizes methylation differences between SV-supporting reads and other reads.

---

SniffCell provides a comprehensive framework for analyzing DNA methylation changes associated with structural variations, enabling researchers to gain insights into epigenetic modifications at the cellular level.
