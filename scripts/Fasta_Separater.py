


#  This python script when runs takes the dataset and outputs every protein and its sequence to one individual file
#  The outputs are .fasta format
#  This script is for further running PSI-BLAST to get homologs of each protein/sequence
#  Input the dataset file below (Default: '../data/300train_ts_red_70.txt')
#  Change the path you want to save the output fasta files (Default: same path/folder of the input dataset file)


def Fasta_Separater(filename):
    f = open(filename,'r+')
    fr = f.readlines()
    for i in range(0,len(fr),3):
       fasta_name = fr[i].lstrip('>').rstrip('\n')
       #path = 'PATH' + str(fasta_name) + '.fasta'
       o = open(fasta_name + '.fasta','w')
       o.write(fr[i]+fr[i+1])
       o.close()
    f.close()


if __name__ == '__main__':
    from datetime import datetime
    start_time = datetime.now()
    print('Program is running...')
    
    #filename = input('Input filename: ')
    Fasta_Separater('../data/300train_ts_red_70.txt')

    end_time = datetime.now()
    print('Starting from',start_time,'to',end_time)
    print('Running Time: {}'.format(end_time - start_time))

