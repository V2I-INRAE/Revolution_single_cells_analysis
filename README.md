# REVOLUTION

Analysis of the paired lungs lobes single cells experiment.

## Setting up

Our data and code is on the cluster @ `/work/project/revo-pig-sc`. A set of library have been installed. The script used for that is on `scripts` folder. When you `ssh` in the cluster:

```bash
srun -c 8 --mem=64G --time=08:00:00 --cpu=6 --pty bash -i
```

Move to `/work/project/revo-pig-sc` and run:

```bash
module load statistics/R/4.6.1
```

## Getting data locally

If you want to work locally you can download the data using the scirpt `scripts/download_reads.sh`

>Use the sequential downloads if you download sequentially so that rsync opens just one ssh connection. For parallel download it is better to setup a ssh key with the cluster, otehrwise it will ask for inputing password each time it opens a ssh connection
 
