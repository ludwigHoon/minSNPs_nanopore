## Krocus comparison
#
# MLST file required are downloaded with Krocus before the start of the following script with:
# krocus_database_downloader  --species "Staphylococcus aureus"
#

SLURM_TEMPLATE <- "slurm_bc_template.tmpl"

library(BiocParallel)
library(data.table)

mlst_database <- "~/mlst_files"
bp <- BatchtoolsParam(workers = 391, cluster="slurm", template=SLURM_TEMPLATE,
    resources=list(walltime=60*60*24*5, ncpus=1))

