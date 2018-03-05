

#  This bash script when runs executes PSI-BLAST using swissprot database*
#  on every *.fasta file to get the corresponding .align and .pssm files
#  Locate to the fastafile folder before running PSI-BLAST locally
#  fasta files are separate into two folders(fasta1,fasta2) for psiblast
#  Output of psiblast stores in pssm folders(pssm1,pssm2) and align folders(align1,align2)
#  This bash script takes fastafiles in one folder(fasta1) as an example
#  Also, when runs remembers to create a log file to record the standard output and error
#  The running command should be: bash psiblast_unix.sh >> psiblast.log 2>&1 &


#Locate to the folder contains fasta files for psiblast
cd ~/temp/fasta/fasta1


echo "Program is running"
echo "Starting time:$(date)"


for i in *.fasta
do
	echo "Start $i on PSI-BLAST at $(date +"%x %T")..."
	time psiblast -db /scratch/swissprot/uniprot.fasta -query $i -evalue 0.01 -num_iterations 3 -word_size 3 -num_threads 4 -out ~/temp/align/align1/$i.align -out_ascii_pssm ~/temp/pssm/pssm1/$i.pssm 
	echo "Finish $i at $(date +"%x %T")..."
done


echo "Ending time: $(date)"
echo "Mission Completed"
