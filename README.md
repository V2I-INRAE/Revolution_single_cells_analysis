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

