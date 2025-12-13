# SniffCell - Identifying cell type specific SV from long-read bulk sequenced tissue only
[![PyPI version](https://img.shields.io/pypi/v/sniffcell.svg)](https://pypi.org/project/sniffcell/)
[![Install](https://img.shields.io/badge/Install-PyPI-3776AB?logo=pypi&logoColor=white)](https://pypi.org/project/sniffcell/)
[![Docs](https://img.shields.io/badge/Docs-GitHub-181717?logo=github)](https://github.com/Fu-Yilei/SniffCell/wiki)
[![Issues](https://img.shields.io/badge/Issues-GitHub-181717?logo=github)](https://github.com/Fu-Yilei/SniffCell/issues)

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
