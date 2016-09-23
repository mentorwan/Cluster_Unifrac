# Cluster_Unifrac

Run R package GUnifrac (https://cran.r-project.org/web/packages/GUniFrac/GUniFrac.pdf) in parallel

## Description

This script uses the GUniFrac package created by Dr. Jun Chen to separate the calculations of UniFrac distances into multiple jobs. These jobs can be run on a local cluster in parallel. Then the results can be merged back into the full distance matrix.

## Usage:

* Obtain the number of rows (n) from the OTU table using script:

Job 1: Rscript run.R 2
....
Job n-1: Rscript run.R n.

* Submit n-1 jobs to cluster.

* After you get n-1 Rdata files then run Load.R to get the full matrix.

 



