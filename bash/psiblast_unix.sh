


#  This bash script when runs executes PSI-BLAST using swissprot database* 
#  (swissprot database is in ../data file, but needs to be extracted first and then dbformat)
#  on every *.fasta file to get the corresponding .align and .pssm files
#  Locate to the fastafile folder before running PSI-BLAST locally (Default: '../data/fasta/')
#  Output of psiblast stores in pssm folder (Default: '../data/pssm/')and align folder (Default: '../data/align/')
#  In terms of align folder, it will not be pushed upto GitHub because of its large size
#  Also, when runs remembers to create a log file to record the standard output and error
#  The running command should be: bash psiblast_unix.sh >> psiblast.log 2>&1 &


#Locate to the folder contains fasta files for psiblast
cd ../data/fasta/


echo "Program is running"
echo "Starting time:$(date)"


for i in *.fasta
do
	echo "Start $i on PSI-BLAST at $(date +"%x %T")..."
	time psiblast -db /scratch/swissprot/uniprot.fasta -query $i -evalue 0.01 -num_iterations 3 -word_size 3 -num_threads 4 -out ../data/align/$i.align -out_ascii_pssm ../data/pssm/$i.pssm 
	echo "Finish $i at $(date +"%x %T")..."
done


echo "Ending time: $(date)"
echo "Mission Completed"
