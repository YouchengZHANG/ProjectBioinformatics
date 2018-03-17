



#  PSSM_Based_Predictor
#  This is the PREDICTOR when runs takes all the sequences in a fasta file 
#  and creates an input to sklearn and predict features from each sequence
#  and finally writes all the predicted result on one output files
#  The predicted result will be saved in a new directory 
#  which contains one Predict_Result.txt and two subfolders (*_result/fasta/ and *_result/pssm/)


#  STEP 1: Import python library
import os
import sys
sys.path.append('../scripts/')
from datetime import datetime
import numpy as np
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
import Table_Creater
import Predict_Separater 
import Predict_Preprocessor
import Predict_PSSM_Processor


start_time = datetime.now()
print('Predicting Program is running...')


#  STEP 2: Check Swissprot database path (Default: '~/project/data/swissprot/uniprot_sprot.fasta')
dir_db = '~/project/data/swissprot/uniprot_sprot.fasta'


#  STEP 3: Load the trained model
#  The trained model applied Random Forest Classifier
print('Loading the model...')
model = joblib.load('pssm_based_model.pkl')
print(model)


#  STEP 4: Input the fastafile
fasta_file = sys.argv[1]
fr = open(fasta_file,'r+')
f = fr.read().splitlines()
test_pid_all, test_seq_all, test_stc_all = [], [], []
for i in range(0,len(f),2):
    if '>' in f[i]:
        test_pid_all.append(f[i].lstrip('>'))
        test_seq_all.append(f[i+1])
fr.close()


#  STEP 4: Build AA_table that convert each amino acid in each sequence into integer form
#  Create a string with twenty '0' and convert into list as [0,0,0,...,0,0,0]
#  Twentry common amino acid 'ARNDCQEGHILKMFPSTWYV'
print('Creating tables...')
AA_table = Table_Creater.Table_Creater()[0]


#  STEP 5: Build int_stc_table that convert each integer form to the structure
#  Example:  int_stc_table[1] = S; stc_table[0] = .
int_stc_table = Table_Creater.Table_Creater()[3]


#  STEP x: Create *_result folder to save predicted output (* refers to the name of unknown fasta file )
#  './*_result/fasta/' : stores individual .fasta files separated from original unknown fasta file
#  './*_result/pssm/' : saves the .pssm files searching from PSI-BLAST
print('Creating _result folder')
fasta_file = sys.argv[1]
unknown_fasta = os.path.splitext(fasta_file)[0]
dir_result = './' + str(unknown_fasta) + '_result/'
dir_fasta = './' + str(unknown_fasta) + '_result/fasta/'
os.makedirs(dir_fasta,exist_ok=True)
dir_pssm = './' + str(unknown_fasta) + '_result/pssm/'
os.makedirs(dir_pssm,exist_ok=True)
dir_align = './' + str(unknown_fasta) + '_result/align/'
os.makedirs(dir_align,exist_ok=True)


#  STEP x: Run Module: Predict_Separater
#  This modules when runs takes each unknown sequence in the fasta file 
#  and outputs each protein and its sequence to each separated single .fasta files
print('Separating input fasta file...')
print('Creating single sequence fasta file...')
Predict_Separater.Predict_Separater(fasta_file, dir_fasta)


#  STEP x: Run Module: Predict_Preprocessor 
#  This module when runs takes each .fasta files in fasta folder (Default: dir_fasta) to PSIBLAST 
#  And outputs .pssm files to pssm folder (Default: dir_pssm)
#  dir_fasta = './' + str(unknown_fasta) + '_result/fasta/'
#  dir_pssm = './' + str(unknown_fasta) + '_result/pssm/'
#  Here use Swiss-Prot database saved in '~/project/data/swissprot/uniprot_sprot.fasta'
Predict_Preprocessor.Predict_Preprocessor(dir_db,dir_fasta)


#  STEP 6: Run Module: Predict_PSSM_Processor 
#  Predict the sequence and convert integer output to character form using int_stc_table
print('Start processing the window...')
print('Window size: [j-18 ... j .. j+2]')
print('Predicting the structure...')

for names in test_pid_all:
    test_window_input = Predict_PSSM_Processor.Predict_PSSM_Processor(names,test_pid_all,test_seq_all,dir_pssm)
    test_int = model.predict(test_window_input)
    test_stc = ''
    for s in list(test_int):
        test_stc += str(int_stc_table[int(s)])
    if len(test_stc) != len(test_seq_all[test_pid_all.index(names)]):
        dif_len = int (len(test_seq_all[test_pid_all.index(names)]) - len(test_stc))
        test_stc += '.' * dif_len
    test_stc_all.append(test_stc)
#print(test_stc_all)


#  STEP 9: Save the predicted output to file
fo = open(dir_result + 'Predict_Result.txt','w')
for i in range(len(test_pid_all)):
    fo.write('>' + test_pid_all[i] + '\n')
    fo.write(test_seq_all[i] + '\n')
    fo.write(test_stc_all[i] + '\n')
fo.close()
print('Predicted ouput file has been created.')
print('Prediction program is completed.')
end_time = datetime.now()
print('Total Running Time: {}'.format(end_time - start_time))
#input("Press any keys to exit the program...")  '''


