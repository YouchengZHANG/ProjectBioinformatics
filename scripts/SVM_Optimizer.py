


#  This script when runs use GridSearchCV() to optimize SVM parameters on single sequence model
#  and outputs .cv_results_ a dictionary which has saved evaluation outcome of every parameters combinations
#  and also save all the standard output to a log file
#  Input the training set used to optimize parameter first (Default: The training set used here is 30train_ts_red_SVM.txt)
#  Custom the window-size (Default: ws_left,ws_right = 9,9)


def SVM_Optimizer(filename):


    #  Open and process dataset file
    foo = open(filename,'r+')
    f = foo.read().splitlines()
    pid,seq,stc = [],[],[]
    for i in range(0,len(f),3):
        if '>' in f[i]:
            pid.append(f[i].lstrip('>'))
            seq.append(f[i+1])
            stc.append(f[i+2])
    foo.close()
    

    #  Build AA_table and stc_table
    AA_num = 20
    X = [ int(z) for z in '0'*AA_num ]
    AA_table = {}
    AA_table['X'] = X 
    for m,n in enumerate('ARNDCQEGHILKMFPSTWYV'):
        X[m] = int(1)
        AA_table[n] = X.copy()
        X[m] = int(0)
    
    stc_table = dict(zip(list('.STt'),[0,1,0,0]))
    stc_input = [stc_table[SS] for stc_each in stc for SS in stc_each]


    #  Set the SVM model and parameters 
    clf1 = SVC(class_weight='balanced',cache_size=7000,random_state=0,probability=True)
    clf1_para = {'kernel':['linear','rbf','sigmoid','poly'],'gamma':[0.01,0.001,0.0001],'C':[0.1,0.5,1,5,10,20]}
    #GS_SVM_2.dict: clf1_para = {'kernel':['linear','rbf','sigmoid','poly'],'gamma':[1,0.1,0.01,0.001,0.0001],'C':[0.1,0.5,1,5,10,20,50]}
    kf = KFold(n_splits=3,shuffle=True,random_state=0)
    GS_SVM = GridSearchCV(estimator=clf1, param_grid=clf1_para,cv=kf)
    classifier = [GS_SVM]


    #  Start process window-size
    print('Start processing window-size...',file=stot)
    for ws_left,ws_right in [[9,9]]:
        ws = ws_left + 1 + ws_right
        print('Window size: '+ str(ws) + ' -> ' + str(ws_left) + '+' + str(ws_right),file=stot)
        tmp = []
        for sub_seq in seq:
            for j in range(0,len(sub_seq)):
                if j < ws_left :
                    tmp.append('X'*int(ws_left-j) + sub_seq[0:int(j + ws_right + 1)])
                elif j >= len(sub_seq) - ws_right:
                    tmp.append(sub_seq[int(j - ws_left):int(j + ws_right + 1)] + 'X'*int(ws_right + 1 - (len(sub_seq) - j)))
                else:
                    tmp.append(sub_seq[int(j - ws_left):int(j + ws_right + 1)])           

        seq_window_input = []
        for seq_each in tmp:
            tmp_input = []
            for AA in seq_each:
                tmp_input += AA_table[AA]
            seq_window_input.append(tmp_input)
        

        #  Convert input to proper format for training
        print('Start converting seq_window_input to proper format for model training...',file=stot)
        map_window_input = np.array(seq_window_input)
        stc_input = np.array(stc_input)


        #  Train the GridSearchCV(SVM) model and save the output to .cv_results_ a dictionary file
        print('***Start training a new model***',file=stot)
        for modelnum,model in enumerate(classifier): 
            print(model,file=stot)
            time_s = datetime.now()
            print('Start model training at '+ str(time_s),file=stot)
            model.fit(map_window_input,stc_input)
            time_e = datetime.now()
            print('Finish model training at '+ str(time_e),file=stot)
            print('Model Training Time: {}'.format(time_e - time_s),file=stot)
            
            print('Saving .cv_results_ ')
            joblib.dump(model.cv_results_,'GS_SVM.dict')
            print('Saved .cv_results_ ',file=stot)



if __name__ == '__main__':
    import numpy as np
    from sklearn.svm import SVC
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.externals import joblib
    from sklearn.model_selection import KFold
    from sklearn import datasets
    from datetime import datetime
    from sklearn.model_selection import GridSearchCV

    stot = open('SVM_para.log','w')
    start_time = datetime.now()
    print('Program is running...',file=stot)

    #filename = input('Input the training set: ')
    SVM_Optimizer('../data/30train_ts_red_SVM.txt')

    end_time = datetime.now()
    print('Total Running Time: {}'.format(end_time - start_time),file=stot)
    stot.close()
