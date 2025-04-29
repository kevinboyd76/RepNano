# RepNano
This repository contains a Snakemake-based pipeline and stand-alone scripts for analyzing DNA replication using BrdU and EdU detection in nanopore sequencing reads. 


# results location
/archive/sansamc/NanoporeSequencing/CS_20250428/20250428/20250428_1853_P2S-02862-A_PBA15434_3150224e/pod5


# load modules to run snakemake
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1


# run snakejob
sbatch --wrap="snakemake -j 999 \
  --use-envmodules \
  --latency-wait 60 \
  --cluster-config config/cluster_config.yml \
  --cluster 'sbatch -A {cluster.account} \
                   -p {cluster.partition} \
                   --cpus-per-task={cluster.cpus-per-task} \
                   --mem={cluster.mem} \
                   --output={cluster.output} \
                   --gres={cluster.gres} \
                   --constraint=\"{cluster.constraint}\"'"

