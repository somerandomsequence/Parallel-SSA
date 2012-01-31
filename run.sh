#!/bin/bash

# Set a walltime for the job. The time format is HH:MM:SS. This program just 
# prints out information to STDOUT and therefore requires little walltime.

#PBS -l walltime=10:00:00

# Name the job
#PBS -N parallel_sa50

# Requesting the janus-debug queue (Limited to one hour of walltime)
#PBS -q janus-normal

# We are requesting one node per host and four hosts. Without the naccesspolicy
# set to singletask we could get four processors on a single host. Since we 
# want to use all the processors on a single host this isn't ideal.

#PBS -l nodes=31:ppn=12

# # Join the Output and Errors in one file. you can also set the path for this output.
# (see "man qsub" for details.)

#PBS -j oe

#
# Source the Dotkit init so you can pull in extra software via the "reuse" command:

. /curc/tools/utils/dkinit

#
# Your code goes here, typical MPI usage for example is:
#
# First we "reuse" to make sure nothing from the submitting environment confuses things.

reuse OpenMPI-1.4.3

#
# cd to the jobs working directory, which you can set above with a #PBS directive,
# (see "man qsub" for details)

cd $PBS_O_WORKDIR

# The OMP_NUM_THREADS environmental variable if set will limit the number of threads for
# each program. Since we want to use all available, we will leave it unset and the program
# will use all available cores. If you request fewer than 12 cores only use that number of threads.

#OMP_NUM_THREADS=4 


# If you get a segmentation fault when using OpenMP it might be due to an insufficient stack size.
# To increase the stack size uncomment the following line and adjust required memory.
#OMP_STACKSIZE=4G

#
# Submit the job.  Note: mpiexec, mpirun and orterun are synonymous. Note the "-pernode" option,
# this option will get mpi to only start one job per host allowing OpenMP to parallelize on that host.
# Since we only requested one core per host (ppn=1), this is superfluous but was left as an example.  

mpirun -np 31 -pernode /home/phillict/thud/parallel_sa4/master_slave
