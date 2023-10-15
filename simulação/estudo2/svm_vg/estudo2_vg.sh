#PBS -N Es2_vg
#PBS -l ncpus=1
#PBS -l walltime=96:00:00
#PBS -m abe
#PBS -M brunoeholtz@gmail.com
cd /home/beholtz/estudo2/svm_vg
module load gcc/4.9.2
module load R/4.2.2
R CMD BATCH estudo2_vg.R
