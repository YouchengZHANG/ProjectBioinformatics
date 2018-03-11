

#  This python script when runs detects the uncommon amino acid
#  This function returns a sorted list that contains all different amino acid
#  and print out the uncommon amino acid if there is any, otherwise print 'No uncommon amino acid.'
#  Input the filename below (Default: '../data/300train_ts_red_70.txt')


def AA_Detector(filename):
    f = open(filename,'r+')
    fr = f.read().splitlines()
    setAA = set()
    for i in range(0,len(fr),3):
        if '>' in fr[i]:
            setAA.update(fr[i+1])
    sl = sorted(list(setAA))
    f.close()
    return sl


if __name__ == '__main__':
    from datetime import datetime
    start_time = datetime.now()

    #filename = input('Input filename: ')
    print(AA_Detector('../data/300train_ts_red_70.txt'))
    for i in AA_Detector('300train_ts_red_70.txt'):
        if not i in 'ACDEFGHIKLMNPQRSTVWY':
            print(i)
    else:
        print('No uncommon amino acid.')
    
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

