


#  This script when runs use GridSearchCV() to optimize RFC parameters on PSSM model
#  and outputs .cv_results_ a dictionary which has saved evaluation outcome of every parameters combinations
#  and also save all the standard output to a log file
#  Input the training set used to optimize parameter first (Default: The training set used here is '../data/300train_ts_red_70.txt')
#  Custom the window-size (Default: ws_left,ws_right = 9,9)
#  Locate to the folder that stores all the pssm files (Default: '../data/pssm/')
#  Notes: 6 Check -> name.dict, name.log, filename path, pssm path, window-size, model name


def RFC_PSSM_Optimizer(filename):


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


    #  Set the RFC model and parameters 
    clf5 = RandomForestClassifier(class_weight='balanced',min_samples_leaf=20)
    #clf5_para = {'n_estimators':[10,20,60,100,200,500,1000],'max_depth':[None,10,20,100,200,500],'class_weight':[None,'balanced'],'max_features':['auto','sqrt','log2'],'min_samples_leaf':[5,10,15,20,30]}
    clf5_para = {'n_estimators':[10,100,500,1000],'max_depth':[None,10,100,500],'max_features':['auto','sqrt','log2']}
    kf = KFold(n_splits=3,shuffle=True,random_state=0)
    GS_PSSM_RFC = GridSearchCV(estimator=clf5, param_grid=clf5_para,cv=kf)
    classifier = [GS_PSSM_RFC]


    #  Start process window-size
    print('Start processing window-size...',file=stot)
    for ws_left,ws_right in [[9,9]]:

        pssm_input = []
        tmp = []
        ws = ws_left + 1 + ws_right
        print('Window size: '+ str(ws) + ' -> ' + str(ws_left) + '+' + str(ws_right),file=stot)

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
            joblib.dump(model.cv_results_,'GS_RFC_PSSM.dict')
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

    stot = open('RFC_PSSM_para.log','w')
    start_time = datetime.now()
    print('Program is running...',file=stot)

    #filename = input('Input the training set: ')
    RFC_PSSM_Optimizer('../data/300train_ts_red_70.txt')

    end_time = datetime.now()
    print('Total Running Time: {}'.format(end_time - start_time),file=stot)
    stot.close()



