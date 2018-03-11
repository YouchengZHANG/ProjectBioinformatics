


#  This script when runs extract feature and structure from all PSSM files, and to convert results into proper inputs for model training
#  Input the training dateset first (Default: The training set used here is '300train_ts_red_70.txt') 
#  Custom the window-size ws_left,ws_right separated by ',' (Default: ws_left,ws_right = 9,9)
#  Locate the folder where stores all the .pssm files of each single sequence (Default: '../data/pssm/')
#  The function used here is sigmoid function
#  The output of this script returns the window-processed sequence input in array format 


def PSSM_Editor(filename):


    #  Open training dataset
    f = open(filename,'r+')
    fr = f.read().splitlines()
    pid,seq,stc = [],[],[]
    for i in range(0,len(fr),3):
        pid.append(fr[i].lstrip('>'))
        stc.append(fr[i+2])
    f.close()

  
    #  Create AA_table and stc_table
    AA_num = 20
    X = [ int(z) for z in '0'*AA_num ]
    stc_table = dict(zip(list('.STt'),[0,1,0,0]))
    stc_input = [stc_table[SS] for stc_each in stc for SS in stc_each]


    #  Start process window-size
    print('Start processing window-size...')    
    #ws_left,ws_right = input('Input ws_left,ws_right = (separate by ',')') 
    
    pssm_input = []
    tmp = []
    ws_left,ws_right = 9,9
    ws = ws_left + 1 + ws_right
    print('Window size: '+ str(ws) + ' -> ' + str(ws_left) + '+' + str(ws_right))

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
                tmp_int_raw = list(map(float,tmp_str))
                tmp_int = [round(1 / (1 + math.exp(-each_num)),2) for each_num in tmp_int_raw]
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


        #  Convert input to proper format for model training
        print('Start converting seq_window_input to proper format for model training...')
        print('Start converting stc_input to proper format for model training...')
        map_window_input = np.array(pssm_window_input)
        stc_input = np.array(stc_input)
        #print(map_window_input.shape)


        #  Return map_window_input
        return map_window_input



if __name__ == '__main__':
    import math
    import numpy as np
    from sklearn.externals import joblib
    from sklearn import datasets
    from datetime import datetime
    
    #stot = open('PSSM_Editor.log','w')
    start_time = datetime.now()
    print('Program is running...')
    
    #filename = input('Input the training dataset: ')
    PSSM_Editor('../data/300train_ts_red_70.txt')

    end_time = datetime.now()
    print('Total Running Time: {}'.format(end_time - start_time))
    #stot.close()

