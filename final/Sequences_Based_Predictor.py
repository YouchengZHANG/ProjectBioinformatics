


#  Single Sequence Predictor (Updated)
#  This is the PREDICTOR when runs takes all the sequences in a fasta file 
#  and creates an input to sklearn and predict features from each sequence
#  and finally writes all the predicted result on one output files
#  The predicted result will be saved in a new directory (/*_result/)
#  which contains one Predict_Result.txt (/*_result/Predict_Result.txt)


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


start_time = datetime.now()
print('Predicting Program is running...')


#  STEP 2: Load the trained model
#  The trained model applied Random Forest Classifier
print('Loading the model...')
model = joblib.load('sequences_based_model.pkl')
print(model)


#  STEP 3: Input the fastafile
fasta_file = sys.argv[1]
fr = open(fasta_file,'r+')
f = fr.read().splitlines()
test_pid_all, test_seq_all, test_stc_all = [], [], []
for i in range(0,len(f),2):
    if '>' in f[i]:
        test_pid_all.append(f[i].lstrip('>'))
        test_seq_all.append(f[i+1])
fr.close()
fr.close()


#  STEP 4: Build AA_table that convert each amino acid in each sequence into integer form
#  Twentry common amino acid 'ARNDCQEGHILKMFPSTWYV'
print('Creating tables...')
AA_table = Table_Creater.Table_Creater()[0]


#  STEP 5: Build int_stc_table that convert each integer form to the structure
#  Example:  int_stc_table[1] = S; stc_table[0] = .
int_stc_table = Table_Creater.Table_Creater()[3]


#  STEP 6: Create *_result folder to save predicted output (* refers to the name of unknown fasta file )
print('Creating _result folder')
fasta_file = sys.argv[1]
unknown_fasta = os.path.splitext(fasta_file)[0]
dir_result = './' + str(unknown_fasta) + '_result/'
os.makedirs(dir_result,exist_ok=True)


#  STEP 7: Process the window size
print('Start processing the window...')
print('Window size: [j-18 ... j .. j+2]')
print('Start predicting the structure...')
ws_left = 18
ws_right = 2
ws = ws_left + 1 + ws_right
#print(ws)
for test_seq in test_seq_all:
    test = []
    for j in range(0,len(test_seq)):
        if j < ws_left :
            test.append('X'*int(ws_left-j) + test_seq[0:int(j + ws_right + 1)])
        elif j >= len(test_seq) - ws_right:
            test.append(test_seq[int(j - ws_left):int(j + ws_right + 1)] + 'X'*int(ws_right + 1 - (len(test_seq) - j)))
        else:
            test.append(test_seq[int(j - ws_left):int(j + ws_right + 1)])
    #print(test)
    test_window_input = []
    for seq_each in test:
        test_input = []
        for AA in seq_each:
            test_input += AA_table[AA]
        test_window_input.append(test_input)
    #print(test_window_input)


#  STEP 8: Convert input into proper format for model training
    test_window_input = np.array(test_window_input)


#  STEP 9: Predict the sequence and convert integer output to character form using int_stc_table
    #print('Predicting the structure...')
    test_int = model.predict(test_window_input)
    test_stc = ''
    for s in list(test_int):
        test_stc += str(int_stc_table[int(s)])
    test_stc_all.append(test_stc)


#  STEP 10: Save the predicted output to file
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
#input("Press any keys to exit the program...") 


