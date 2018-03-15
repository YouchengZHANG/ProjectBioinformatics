


#  Module: Predict_Separater
#  This modules when runs takes each unknown sequence in the fasta file 
#  and outputs each protein and its sequence to each separated single .fasta files
#  This .fasta output is for further running PSI-BLAST to get homologs of each protein/sequence for prediction*
#  The output fasta files are saved in folder *_result/fasta
#  dir_fasta = './' + str(unknown_fasta) + '_result/fasta/'


def Predict_Separater(unknown_fasta, dir_fasta):


    print('Separating each protein and its sequence to individual fasta files...')
    f = open(unknown_fasta,'r+')
    fr = f.readlines()

    for i in range(0,len(fr),2):
        if '>' in fr[i]:
            fasta_name = fr[i].lstrip('>').rstrip('\n')
            path_fasta = dir_fasta + str(fasta_name) + '.fasta'
            o = open(path_fasta,'w')
            o.write(fr[i]+fr[i+1])
            o.close()
    f.close()
    print('Fasta files has been processed...')


