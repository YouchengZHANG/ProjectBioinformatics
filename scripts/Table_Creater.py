


#  Module: Table_Creater
#  This module when runs creates and returns amino acid table, structure table, and residue X
#  When imported, Table_Creater()[0] refers to AA_table, Table_Creater()[1] refers to stc_table, Table_Creater()[2] refers to X


def Table_Creater():


    # STEP 1: Build AA_table that convert each amino acid in each sequence into integer form
    # Create a string with twenty '0' and convert into list as [0,0,0,...,0,0,0]
    # Twentry common amino acid 'ARNDCQEGHILKMFPSTWYV'
    AA_num = 20
    X = [ int(z) for z in '0'*AA_num ]
    AA_table = {}
    AA_table['X'] = X 
    for m,n in enumerate('ARNDCQEGHILKMFPSTWYV'):
        X[m] = int(1)
        AA_table[n] = X.copy()
        X[m] = int(0)
    

    # STEP 2: Build stc_table that convert each structure in each sequence into integer form
    # Structure 'S' is assigned value 1, while other structures are assigned value 0
    stc_table = dict(zip(list('.STt'),[0,1,0,0]))


    # STEP 3: Return two tables and residue X 
    return AA_table,stc_table,X


