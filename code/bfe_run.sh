###========================================
#!/bin/bash
# Your job will use 1 node, 16 cores, and 62gb of memory total.
#PBS -q "el_high_pri" 
#PBS -l select=2:ncpus=16:mem=62gb:pcmem=4gb
### Specify a name for the job
#PBS -N bfc-c_test
### Specify the group name
#PBS -W group_list=kmkshare
### Used if job requires partial node only
#PBS -l place=free:shared
### CPUtime required in hhh:mm:ss.
### Leading 0's can be omitted e.g 48:0:0 sets 48 hours
#PBS -l cput=4000:00:00
### Walltime is how long your job will run
#PBS -l walltime=200:00:00
### Joins standard error and standard out
#PBS -j oe


in_path=/rsgrps/gbeslastudents/nicolas/MWLMC_sims/current_snaps_ascii/MW/
snap_name=MW_81_MWLMC3_100M_new_b0
out_path=/home/u9/jngaravitoc/codes/BFE-c/code/
out_name=test.txt
nmax=20
lmax=20
npart=81619663
n_samples=0
n_sampling=0
init_snap=90
final_snap=90
rs=40.85

OMP_NUM_THREADS=20 /home/u9/jngaravitoc/codes/BFE-c/code/bfe-c $nmax $lmax $npart $n_samples $in_path $snap_name $out_path $out_name $n_sampling $init_snap $final_snap $rs

###end of script
