# SniffCell - Identifying cell type specific SV from long-read bulk sequenced tissue only

SniffCell is a tool designed to analyze DNA methylation changes associated with structural variations (SVs), including mosaic SVs. It processes primary alignments from BAM files and provides detailed outputs for visualization and analysis.


    positional arguments:
    {find,deconv,anno,svanno,dmsv}
        find                Find cell type-specific DMRs.
        deconv              Deconvolve cell-type composition from methylation data.
        anno                Annotate variants with cell-type-specific methylation.
        svanno              Use pre-annotated reads csv to annotate variants' cell types
        dmsv                Find out which SV's supporting reads have differential methylation compared to non-supporting reads.

    options:
    -h, --help            show this help message and exit
    -v, --version         show program's version number and exit