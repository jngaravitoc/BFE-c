###========================================
#!/bin/bash
# Your job will use 1 node, 16 cores, and 62gb of memory total.
#PBS -q "standard" 
#PBS -l select=1:ncpus=2:mem=62gb:pcmem=4gb
### Specify a name for the job
#PBS -N sample_halo 
### Specify the group name
#PBS -W group_list=gbesla
### Used if job requires partial node only
#PBS -l place=free:shared
### CPUtime required in hhh:mm:ss.
### Leading 0's can be omitted e.g 48:0:0 sets 48 hours
#PBS -l cput=800:00:00
### Walltime is how long your job will run
#PBS -l walltime=40:00:00
### Joins standard error and standard out
#PBS -j oe




#module load openmpi
python3 "/home/u9/jngaravitoc/codes/SCF_tools/code/"gadget_to_scf_halo_com.py

###end of script
