


#  Module: Predict_PFM_Processor
#  This module when runs processes window-size for PFM model 
#  and convert pfm_window_input that contains window-size processed input (integer format) to map_window_input (array format)
#  then return variable map_window_input
#  When imported, the function needs to state four variables: 
#  names, pid(protein id list), seq(sequence list), dir_pssm(pssm folder)


def Predict_PFM_Processor(names,pid,seq,dir_pssm):


    #  STEP 1: Import the library and modules/functions
    import sys
    sys.path.append('./')
    import Table_Creater
    import numpy as np


    #  STEP 2: Import X 
    #  X residue uses [0,0,...,0,0] that contains 20 zeros to represent
    X = Table_Creater.Table_Creater()[2]


    #  STEP 3: Choose the Window-size [j-ws_left,...,j,...j+ws_right], 
    #  where ws_left/ws_right refers to the number of elements before/after central residue j 
    #  Example: If ws_left,ws_right = 18,2 , then the window looks like [j-18,...,j,j+2] and the size of window is 21
    #  For those position with no corresponding residues, 'X' is used
    #  Here uses frequency function to normalize each figure in the pfm
    
    ws_left, ws_right = 18, 2
    pfm_input = []
    tmp = []

    path = str(dir_pssm) + str(names) + '.fasta.pssm'
    with open(path, 'r+') as fpfm:
        pfm = fpfm.read().splitlines()
        if len(seq[pid.index(names)]) >= 70:
            pfm2 = pfm[3:73]
        else:
            pfm2 = pfm[3:(3+len(seq[pid.index(names)]))]

        for s in pfm2:
            tmp_str = s.split()[22:42]
            tmp_int_raw = list(map(float,tmp_str))
            tmp_int = [round( each_num / 100, 2) for each_num in tmp_int_raw]
            pfm_input.append(tmp_int)
        fpfm.close()

        for j in range(0,len(pfm_input)):
            if j < ws_left :
                tmp.append( X * int(ws_left-j) + sum([pos for pos in pfm_input[0:int(j + ws_right + 1)]], []) )
            elif j >= len(pfm_input) - ws_right:
                tmp.append( sum([pos for pos in pfm_input[int(j - ws_left):int(j + ws_right + 1)]], []) + X * int(ws_right + 1 - (len(pfm_input) - j)) )
            else:
                tmp.append( sum([pos for pos in pfm_input[int(j - ws_left):int(j + ws_right + 1)]], []) )     
        
        pfm_window_input = tmp

    
    #  STEP 4: Convert input into proper format for model training
    #  Convert the window-processed pssm_window_input to map_window_input array format
        #print('Start converting pssm_window_input to proper format for model training...')
        map_window_input = np.array(pfm_window_input)


    #  STEP 5: Return map_window_input
    return map_window_input


