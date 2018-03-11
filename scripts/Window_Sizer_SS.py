


#  Module: Window_Sizer_SS
#  This module when runs processes window-size for single sequence model 
#  and convert seq_window_input that contains window-size processed input (integer format) to map_window_input (array format)
#  then return variable map_window_input
#  When imported, the function needs to state three variables: ws_left(window size left), ws_right(window size right), seq(sequence list)


def Window_Sizer_SS(ws_left,ws_right,seq):


    #  STEP 1: Import the library and modules/functions
    import sys
    sys.path.append('./')
    import Table_Creater
    import numpy as np
    

    #  STEP 2: Import the AA_table
    AA_table = Table_Creater.Table_Creater()[0]


    #  STEP 3: Choose the Window-size [j-ws_left,...,j,...j+ws_right], 
    #  where ws_left/ws_right refers to the number of elements before/after central residue j 
    #  Example: If ws_left,ws_right = 3,1 , then the window looks like [j-3,j-2,j-1,j,j+1] and the size of window is 5
    #  For those position with no corresponding residues, 'X' is used
    print('Start processing window-size...')
    #ws_left = 9
    #ws_right = 9
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
    

    #  STEP 4: Convert input into proper format for model training
    #  Convert the window-processed seq_window_input to map_window_input array format
    print('Start converting seq_window_input to proper format for model training...')
    map_window_input = np.array(seq_window_input)


    #  STEP 5: Return map_window_input
    return map_window_input


