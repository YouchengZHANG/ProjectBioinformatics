

#  Single Sequence Predictor
#  This is the python script when runs takes a sequence
#  and create an input to sklearn and predict your feature from a single sequence using one sequence


#  STEP 1: Import python library
from datetime import datetime
import numpy as np
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib


start_time = datetime.now()
print('Predicting Program is running...')


#  STEP 2: Load the trained model
#  The trained model applied Random Forest Classifier
print('Loading the model...')
clf5 = joblib.load('RFC_01.pkl')


#  STEP 3: Input the fastafile
#fastafile = input('Input fastafile:')
fr = open('unknown_seq.fasta','r+')
f = fr.read().splitlines()
for i in range(0,len(f),2):
    if '>' in f[i]:
        test_pid = f[i].lstrip('>')
        test_seq = f[i+1]
#print('The number of entry:',len(pid),len(seq),len(stc))
#print(pid,'\n',seq,'\n',stc,'\n')
fr.close()


#  STEP 4: Build AA_table that convert each amino acid in each sequence into integer form
#  Create a string with twenty '0' and convert into list as [0,0,0,...,0,0,0]
#  Twentry common amino acid 'ARNDCQEGHILKMFPSTWYV'
#print('Creating amino acid table...')
AA_num = 20
X = [ int(z) for z in '0'*AA_num ]
AA_table = {}
AA_table['X'] = X 
for m,n in enumerate('ARNDCQEGHILKMFPSTWYV'):
    X[m] = int(1)
    AA_table[n] = X.copy()
    X[m] = int(0)
#print(AA_table)
#test_seq_input = [AA_table[AA] for seq_each in seq for AA in seq_each]


#  STEP 5: Build int_stc_table that convert each predicted integer output to each structure
#  Example:  stc_table[1] = S; stc_table[2] = T
#print('Creating structure table...')
int_stc_table = dict(zip(list(range(4)),list('.STt')))
print(int_stc_table)


#  STEP 6: Choose the Window-size [j-ws_left,...,j,...j+ws_right], 
#  where ws_left/ws_right refers to the number of elements before/after central residue j 
#  Example: If ws_left,ws_right = 3,1 , then the window looks like [j-3,j-2,j-1,j,j+1] and the size of window is 5
#  For those position with no corresponding residues, 'X' is used
print('Start processing the window...')
ws_left = 1
ws_right = 1
ws = ws_left + 1 + ws_right
#print(ws)
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


#  STEP 7: Convert input into proper format for model training (without using OneHotEncoder() command)
test_window_input = np.array(test_window_input)
#print(test_window_input.shape)


#  STEP 8: Predict the sequence and convert integer output to character form using int_stc_table
print('Predicting the structure...')
test_int = clf5.predict(test_window_input)
#print(test_int)
test_stc = ''
for s in list(test_int):
    test_stc += str(int_stc_table[int(s)])
#print(test_stc)


#  STEP 9: Save the predicted output to file
fo = open('unknown_seq_output.txt','w')
fo.write('>' + test_pid + '\n')
fo.write(test_seq + '\n')
fo.write(test_stc + '\n')
fo.close()
print('Predicted ouput file has been created.')
print('Prediction program is completed.')
end_time = datetime.now()
print('Total Running Time: {}'.format(end_time - start_time))
#input("Press any keys to exit the program...")
