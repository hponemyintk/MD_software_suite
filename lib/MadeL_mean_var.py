try:
    np.array()
except NameError:
    # print("importing numpy")
    import numpy as np
try:
    sys
except NameError:
    # print("importing sys")
    import sys
np.set_printoptions(threshold=sys.maxsize)
try:
    from lib.read_coord import get_Coul
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def doMadelungMeanVar(CoulDump,Tstep,itype,charge):
    ### get mean first ###
    mean = 0
    ct = 0
    coulfile0 = open(CoulDump,"r")
    cur_step = 0
    while cur_step < Tstep:
        cur_coul = get_Coul(coulfile0)
        cur_coul_i = cur_coul[ (cur_coul[:,0]==itype) ]
        # print(cur_coul_i[:,1])
        mean += np.sum(cur_coul_i[:,1]/(charge/2))
        ct += len(cur_coul_i)
        cur_step += 1
    coulfile0.close()
    mean /= ct


    var = 0
    ct = 0
    coulfile1 = open(CoulDump,"r")
    cur_step = 0
    while cur_step < Tstep:
        cur_coul = get_Coul(coulfile1)
        cur_coul_i = cur_coul[ (cur_coul[:,0]==itype) ]
        cur_coul_i = (cur_coul_i[:,1]/(charge/2)-mean)
        var += np.sum(cur_coul_i*cur_coul_i)
        ct += len(cur_coul_i)
        cur_step += 1
    var /= ct

    print("mean is :::",mean,"\n variance is :::",var)
