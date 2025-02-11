# Test brdu_edu_sliding_windows_from_bam.py 

If using a system running on slurm, use the following command to run a test. You should run this from the root directory of this repository. The test will output two bed files to the "results" directory called "BrdU.bedgraph" and "EdU.bedgraph". This runs the script on all 50 aligned reads in the file creating 32,715 lines in each bedgraph file.

```         
module load slurm/20.02 python/3.10

sbatch -c 2 --wrap="\
python scripts/brdu_edu_sliding_windows_from_bam.py \
-b tests/data/test_data.bam \
-p 2 \
-w 100 \
-s 10"
```

To limit the run to just the first two reads, use the "-n" option. This should output bedgraph files with only 7,754 lines.

```         
module load slurm/20.02 python/3.10

sbatch -c 2 --wrap="\
python scripts/brdu_edu_sliding_windows_from_bam.py \
-b tests/data/test_data.bam \
-p 2 \
-w 100 \
-s 10 \
-n 2"
```

To limit the run to a specific read, use the "-r" option. For the specified read, this should output bedgraph files with 12,035 lines.

```         
module load slurm/20.02 python/3.10

sbatch -c 2 --wrap="\
python scripts/brdu_edu_sliding_windows_from_bam.py \
-b tests/data/test_data.bam \
-p 2 \
-w 100 \
-s 10 \
-r 394b52e6-d5ce-436c-85a0-d301ab2a129e"
```
