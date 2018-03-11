


#  This script when runs use GridSearchCV() to optimize SVM parameters on PSSM model
#  and outputs .cv_results_ a dictionary which has saved evaluation outcome of every parameters combinations
#  and also save all the standard output to a log file
#  Input the training set used to optimize parameter first (Default: The training set used here is '../data/30train_ts_red_SVM.txt')
#  Custom the window-size (Default: ws_left,ws_right = 9,9)
#  Locate to the folder that stores all the pssm files (Default: '../data/pssm/')
#  Notes: 6 Check -> name.dict, name.log, filename path, pssm path, window-size, model name


def SVM_PSSM_Optimizer(filename):
    
    
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
    #GS_PSSM_SVM_2.dict: clf1_para = {'kernel':['linear','rbf','sigmoid','poly'],'gamma':[1,0.1,0.01,0.001,0.0001],'C':[0.1,0.5,1,5,10,20,50]}
    kf = KFold(n_splits=3,shuffle=True,random_state=0)
    GS_PSSM_SVM = GridSearchCV(estimator=clf1, param_grid=clf1_para,cv=kf)
    classifier = [GS_PSSM_SVM]


    #  Start process window-size
    print('Start processing window-size...',file=stot)
    for ws_left,ws_right in [[9,9]]:

        pssm_input = []
        tmp = []
        ws = ws_left + 1 + ws_right
        print('Window size: '+ str(ws) + ' -> ' + str(ws_left) + '+' + str(ws_right),file=stot)

        for names in pid:
            path = '../data/pssm/' + str(names) + '.fasta.pssm'
            with open(path, 'r+') as fpssm:
                pssm = fpssm.read().splitlines()
                if len(stc[pid.index(names)]) >= 70:
                    pssm2 = pssm[3:73]
                else:
                    pssm2 = pssm[3:(3+len(stc[pid.index(names)]))]

                for s in pssm2:
                    tmp_str = s.split()[2:22]
                    tmp_int2 = list(map(float,tmp_str))
                    tmp_int = [round(1 / (1 + math.exp(-each_num)),2) for each_num in tmp_int2]
                    pssm_input.append(tmp_int)
                fpssm.close()
   
                for j in range(0,len(pssm_input)):
                    if j < ws_left :
                        tmp.append( X * int(ws_left-j) + sum([pos for pos in pssm_input[0:int(j + ws_right + 1)]], []) )
                    elif j >= len(pssm_input) - ws_right:
                        tmp.append( sum([pos for pos in pssm_input[int(j - ws_left):int(j + ws_right + 1)]], []) + X * int(ws_right + 1 - (len(pssm_input) - j)) )
                    else:
                        tmp.append( sum([pos for pos in pssm_input[int(j - ws_left):int(j + ws_right + 1)]], []) )     
                pssm_input = []
        
        pssm_window_input = tmp


        #  Convert input to proper format for training
        print('Start converting seq_window_input to proper format for model training...',file=stot)
        map_window_input = np.array(pssm_window_input)
        stc_input = np.array(stc_input)


        #  Train the GridSearchCV(SVM) model and save the output to .cv_results_ a dictionary file
        print('***Start training a new model***',file=stot)
        for modelnum,model in enumerate(classifier):
            print(model,file=stot)
            time_s = datetime.now()
            print('Start model training at '+ str(time_s),file=stot)
            model_save = model.fit(map_window_input,stc_input)
            time_e = datetime.now()
            print('Finish model training at '+ str(time_e),file=stot)
            print('Model Training Time: {}'.format(time_e - time_s),file=stot)
            
            print('Saving .cv_results_ ')
            joblib.dump(model.cv_results_,'GS_PSSM_SVM.dict')
            print('Saved .cv_results_ ',file=stot)


   
if __name__ == '__main__':
    import math
    import numpy as np
    from sklearn.svm import SVC
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.externals import joblib
    from sklearn.model_selection import KFold
    from sklearn import datasets
    from datetime import datetime
    from sklearn.model_selection import GridSearchCV

    stot = open('SVM_PSSM_para.log','w')
    start_time = datetime.now()
    print('Program is running...',file=stot)

    #filename = input('Input the training set: ')
    SVM_PSSM_Optimizer('../data/30train_ts_red_SVM.txt')

    end_time = datetime.now()
    print('Total Running Time: {}'.format(end_time - start_time),file=stot)
    stot.close()



