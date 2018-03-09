


#  Single Sequence Feature Extractor (Updated)
#  This is the python script when runs takes a dataset, creates sequence and feature inputs to sklearn 
#  and use classifier to train a model
#  The updates are recorded under each step description
#  Model parameters are optimized and updated by mannually trying different values/options 
#  or by running the optimizer python script which uses the GridSearchCV command for automatically optimizing


#  STEP 1: Import python library
from datetime import datetime
import numpy as np
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib


start_time = datetime.now()
print('Model Training Program is running...')


#  STEP 2: Extract feature from input file with 3 lines Name,Sequence,Structure
#  Updates: original dataset has been manually randomly splitted into half train set and half test set
#  Updates: Before splitting, the dataset has been performed homology reduction
#  Updates: According to previous research, only the first 70 N-terminal residues will be used to ensure the model quality
#dataset = input('Input dataset:')
fr = open('300train_ts_red_70.txt','r+')
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
#  Notes: Use modules and import function to reduce code redundency 
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


#  STEP 4: Build stc_table that convert each structure in each sequence into integer form
#  Example:  stc_table[S] = 1; stc_table[T] = 0
#  Updates: Two different stc_table(below) are compared. With only 'S' assign the value 1 and other structure assign the value 0, the model 
#           has a higher accuracy as well as better MCC measurement
#  stc_table = dict(zip(list('.STt'),list(range(4))))
stc_table = dict(zip(list('.STt'),[0,1,0,0]))
stc_input = [stc_table[SS] for stc_each in stc for SS in stc_each]


#  STEP 5: Choose the Window-size [j-ws_left,...,j,...j+ws_right], 
#  Example: If ws_left,ws_right = 3,1 , then the window looks like [j-3,j-2,j-1,j,j+1] and the size of window is 5
#  For those position with no corresponding residues, 'X' is used
#  Updates: window-size has been optimized, however there's almost on difference between symmetric and asymmetric window size
#  Updates: It shows relatively best accuracy and MCC when the half window size equals to 9
#  Notes: Use modules and import function 
print('Start processing window-size...')
ws_left = 9
ws_right = 9
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
#  Updates: actually no need to convert the stc_input to array format
print('Start converting inputs to proper format for model training...')
map_window_input = np.array(seq_window_input)
#np.set_printoptions(edgeitems =33)
#print(map_window_input.shape)
stc_input = np.array(stc_input)
#print(stc_input.shape)


#  STEP 7: Check the shape of label and feature input
print('Checking the shape of labels and features...')
if map_window_input.shape[0] != stc_input.shape[0]:
    print('Error: The shapes of label and feature are not corresponded.')


#  STEP 8: Train SVM/RFC model 
#  **Notes: always remember to check the shape of input array!
#  Updates: uploaded three optimizers on GitHub that are used to GridSearch the parameters
#           For SVM model: linear works better enough so far while class_weight needs to be balanced for the imbalanced dataset
#           For RFC model: n_estimators larger works better, much faster than SVM, more suitable for imbalaced dataset
#  Updates: Currently, GridSearchCV is used to optimize some other parameters e.g. gamma,C, thus still waiting for better results
#  Updates: Models that apply PSSM are also being optimized, not sure whether the outcome would be largely different
#  Updates: The test set accuracy is about 0.88 - 0.90, but MCC is also calucalated with only 0.48 - 0.50, overfitting problem 
print('The model is using Random Forest Classifier...')
#clf1 = SVC(C=200,kernel='linear',class_weight='balanced',cache_size=8000,random_state=0,probability=True)
clf5 = RandomForestClassifier(n_estimators=500,min_samples_leaf=20,max_features='auto',criterion='gini')
print('Start model training at '+ str(datetime.now()))
clf5.fit(map_window_input,stc_input)
print('Finish model training at '+ str(datetime.now()))


#  STEP 9: Save the trained model for further reloading in the Single Sequence Predictor Program
print('Saving trained model as RFC_01.pkl...')
joblib.dump(clf5, 'RFC_01.pkl')
print('Model training program is completed.')
end_time = datetime.now()
print('Total Running Time: {}'.format(end_time - start_time))
#input("Press any keys to exit the program...")

