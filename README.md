Cluster_Unifrac

This script is addition to available GUnifrac package by Dr Chen Jun. It will separate the calculation of unifrac distances to multiple jobs. 
Those jobs can be run on local cluster in parallel.  Then the results can be merged back into the big distance matrix.

Usage:

1. Get number of row (n) from OTU table. Using script

Job 1: Rscript run.R 2
....
Job n-1: Rscript run.R n-1.

2. Submit n-1 jobs to cluster.

3. After you got n-1 Rdata files. Run Load.R to get a full matrix.

 



