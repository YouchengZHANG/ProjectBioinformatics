


#  Module: Predict_Preprocessor 
#  This module when runs takes each .fasta files in fasta folder (Default: dir_fasta) to PSIBLAST 
#  And outputs .pssm files to pssm folder (Default: dir_pssm), .align files to align folder (Default: dir_align)
#  dir_fasta = './' + str(unknown_fasta) + '_result/fasta/'
#  dir_pssm = './' + str(unknown_fasta) + '_result/pssm/'
#  dir_align = './' + str(unknown_fasta) + '_result/align/'
#  Here use Swiss-Prot database saved in '~/project/data/swissprot/uniprot_sprot.fasta'


def Predict_Preprocessor(dir_db,dir_fasta):


    #  STEP 1: Import the library and state the directory of database for PSIBLAST
    import os
    

    #  STEP 2: Run PSIBLAST and output .pssm files for each protein sequence
    #  dir_fasta = './' + str(unknown_fasta) + '_result/fasta/'
    #  dir_db = '~/project/data/swissprot/uniprot_sprot.fasta'
    Run_PSIBLAST = '''#!/bin/bash
    cd {dir_fasta}
    echo "PSI-BLAST is running"
    echo "Starting time:$(date)"
    for i in *.fasta
    do
        psiblast -db {dir_db} -query $i -evalue 0.01 -num_iterations 3 -word_size 3 -num_threads 4 -out ../align/$i.align -out_ascii_pssm ../pssm/$i.pssm 
    done
    echo "Ending time: $(date)"
    echo "Mission Completed"
    cd ../../
    '''.format(dir_db=dir_db, dir_fasta=dir_fasta)
    os.system(Run_PSIBLAST)


