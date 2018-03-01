


#  Single Sequence Feature Extractor
#  Single Sequence Model Trainer
#  This is the python script when runs takes a dataset, creates sequence and feature inputs to sklearn 
#  and use classifier to train a model


#  STEP 1: Import python library
from datetime import datetime
import numpy as np
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib


start_time = datetime.now()
print('Model Training Program is running...')


#  STEP 2: Extract feature from input file with 3 lines Name,Sequence,Structure
#dataset = input('Input dataset:')
fr = open('train_set.txt','r+')
f = fr.read().splitlines()
pid,seq,stc = [],[],[]
for i in range(0,len(f),3):
    if '>' in f[i]:
        pid.append(f[i].lstrip('>'))
        seq.append(f[i+1])
        stc.append(f[i+2])
#print('The number of entry:',len(pid),len(seq),len(stc))
fr.close()


#  STEP 3: Build AA_table that convert each amino acid in each sequence into integer form
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
#seq_input = [AA_table[AA] for seq_each in seq for AA in seq_each]
#print(len(seq_input))


#  STEP 4: Build stc_table that convert each structure in each sequence into integer form
#  Example:  stc_table[S] = 1; stc_table[T] = 2
#print('Creating structure table...')
stc_table = dict(zip(list('.STt'),list(range(4))))
#print(stc_table)
stc_input = [stc_table[SS] for stc_each in stc for SS in stc_each]
#print(len(stc_input))


#  STEP 5: Choose the Window-size [j-ws_left,...,j,...j+ws_right], 
#  where ws_left/ws_right refers to the number of elements before/after central residue j 
#  Example: If ws_left,ws_right = 3,1 , then the window looks like [j-3,j-2,j-1,j,j+1] and the size of window is 5
#  For those position with no corresponding residues, 'X' is used
print('Start processing window-size...')
ws_left = 1
ws_right = 1
ws = ws_left + 1 + ws_right
#print(ws)
tmp = []
for sub_seq in seq:
    for j in range(0,len(sub_seq)):
        if j < ws_left :
            tmp.append('X'*int(ws_left-j) + sub_seq[0:int(j + ws_right + 1)])
        elif j >= len(sub_seq) - ws_right:
            tmp.append(sub_seq[int(j - ws_left):int(j + ws_right + 1)] + 'X'*int(ws_right + 1 - (len(sub_seq) - j)))
        else:
            tmp.append(sub_seq[int(j - ws_left):int(j + ws_right + 1)])           
#print(tmp)
seq_window_input = []
for seq_each in tmp:
    tmp_input = []
    for AA in seq_each:
        tmp_input += AA_table[AA]
    seq_window_input.append(tmp_input)
#print(seq_window_input)


#  STEP 6: Convert input into proper format for model training
#  Convert the seq_window_input to map_window_input (without using OneHotEncoder() command)
#  Convert the stc_input to stc_input
print('Start converting inputs to proper format for model training...')
map_window_input = np.array(seq_window_input)
np.set_printoptions(edgeitems =33)
#print(map_window_input.shape)
stc_input = np.array(stc_input)
#print(stc_input.shape)


#  STEP 7: Check the shape of label and feature input
print('Checking the shape of labels and features...')
if map_window_input.shape[0] != stc_input.shape[0]:
    print('Error: The shapes of label and feature are not corresponded.')


#  STEP 8: Train SVM model 
#  **Notes: predict input set shape should be same as training set, e.g. if test_input.shape[1] = 1 not 20*ws, it should be reshaped first
#  Use Random Forest Classifier first to roughly test because it is much more faster than SVM
print('The model is using SVM...')
clf1 = SVC(C=1,kernel='linear')
#clf5 = RandomForestClassifier(min_samples_leaf=20)
print('Start model training at '+ str(datetime.now()))
clf1.fit(map_window_input,stc_input)
#clf5.fit(map_window_input,stc_input)
print('Finish model training at '+ str(datetime.now()))
#roughly calucation: 100 samples linear kernel with ws = 1 takes 2 mins ; with ws = 3 takes 6 mins ; with ws = 7 takes 15 mins 
#roughly calucation: 9000 samples with ws = 14 takes about 2 days


#  STEP 9: Save the trained model for further reloading in the Single Sequence Predictor Program
print('Saving trained model as RFC_01.pkl...')
joblib.dump(clf5, 'RFC_01.pkl')
print('Model training program is completed.')
end_time = datetime.now()
print('Total Running Time: {}'.format(end_time - start_time))
input("Press any keys to exit the program...")

