# SH3RF3 dual SV/TR example

This public fixture is a 45-kb native-GRCh38 slice around the AATGG repeat at
`chr2:109199301-109199876`. It runs the complete SniffCell workflow from a
small atlas through deconvolution, SV/TR discovery, annotation, visualization,
and reporting.

Run it with an indexed GRCh38 no-alt reference:

```bash
./tests/data/sh3rf3_dual_sv_tr/run_example.sh /path/to/GRCh38_no_alt.fa
```

The fixture contains 976 alignment records from 914 reads, nine
atlas samples (six Neuron and three Oligodendrocyte), 164 atlas regions, and 24
nearby TR catalog entries. The BAM retains the normal GRCh38 sequence
dictionary and coordinates; it does not use a miniature reference or require
Sniffles `--all-contigs`.

Expected results are in `expected/`. The TR branch reports a Neuron-only
2,502-bp expansion with 5 versus 0 alternate-supporting reads. The SV branch
independently reports a Neuron-only 726-bp insertion supported by the 596-,
726-, and 942-bp insertion alignments.

See the [full wiki tutorial](../../../wiki/SH3RF3-Dual-SV-TR-Example.md) for
requirements, commands, output layout, and interpretation.
