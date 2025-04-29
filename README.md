# RepNano
This repository contains a Snakemake-based pipeline and stand-alone scripts for analyzing DNA replication using BrdU and EdU detection in nanopore sequencing reads. 

# load modules
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1

# check file path (dry run)
snakemake -npr

# submit snakejob
sbatch --wrap="snakemake -j 999 --use-envmodules --latency-wait 30 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task} --mem {cluster.mem} --output {cluster.output}'"
