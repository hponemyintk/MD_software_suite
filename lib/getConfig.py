import numpy as np
import pandas as pd
from lib.math_lib import unwrapBond
from lib.read_coord import get_wrapped_coord, manual_parsing

# IDs = np.array([[0,290,6,12,459,3210,5307,7104,10069,10496,198,1556,1662,2375,3210,3782,4054,4445,4850,5013,5192,5307,5834,5907,6202,6350,6682,8479,8860,9490,9564,9693,9903,10069,10496,10696,459,813,1020,1517,1937,2024,2073,2129,2890,3025,3210,3664,4348,5187,5307,5483,6242,6760,7021,7104,8171,9106,9706,10051,10069,10133,10980],
# [1,45,0,2,503,3809,4515,4888,8778,10912,503,677,919,1928,2590,2873,3809,3918,4110,4539,4732,4739,4751,5765,6828,7283,7490,7609,7939,8116,8326,8732,8778,9783,23,225,361,503,528,1430,2452,2473,2768,2777,2989,3473,3809,3872,4043,4766,4774,4888,7957,8059,8116,8241,8248,8748,9473,9644,10209,10912,]])


def makeConfig(lammpstrj,RefFile,skipline):
    # IDs = pd.read_table(RefFile,header=None,sep='\s+').values
    IDs = manual_parsing(RefFile," ",skipline,dtype=np.float64)
    # get rid of distance/angles infos
    IDs = np.hstack((IDs[:,0:4],IDs[:,10:]))
    print(IDs,np.shape(IDs))

    tmpfname = open(lammpstrj,"r")
    cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(tmpfname)
    WaterCt = len(cur_pos[(cur_pos[:,0]==1)])
    ionCt = len(cur_pos[(cur_pos[:,0]==3)|(cur_pos[:,0]==4)])
    FerriCt = len(cur_pos[(cur_pos[:,0]==5)])*13
    print(WaterCt,ionCt,FerriCt)
    tmpfname.close()

    xyz_in = open(lammpstrj,"r")
    cur_t = -1
    for ct, config in enumerate(IDs,1):
        while config[0] != cur_t:
            # print("here***",cur_t,config[0])
            cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(xyz_in)
            cur_t += 1

        print("Tstep is :::",tcur_Tstep)
        print("now doing::",ct,cur_t)
        cat = np.array([cur_pos[cur_pos[:,0]==3][int(config[1])]])
        print(cat,np.shape(cat))
        Ferri = cur_pos[(cur_pos[:,0]==5)|(cur_pos[:,0]==6)|(cur_pos[:,0]==7)]
        Ferri = Ferri[(int(config[2])*13):(int(config[2])*13+13)]
        # Now make sure there is no bonds that are loop around the box
        Ferri[1:,1:4],Ferri[0,1:4] = unwrapBond(Ferri[1:,1:4],Ferri[0,1:4],boxL)
        # print(Ferri)

        Ferro = cur_pos[(cur_pos[:,0]==8)|(cur_pos[:,0]==9)|(cur_pos[:,0]==10)]
        Ferro = Ferro[(int(config[3])*13):(int(config[3])*13+13)]
        Ferro[1:,1:4],Ferro[0,1:4] = unwrapBond(Ferro[1:,1:4],Ferro[0,1:4],boxL)

        # Now unwrap the cation Ferri/Ferro distances
        Ferri[:,1:4],cat[0,1:4] = unwrapBond(Ferri[:,1:4],cat[0,1:4],boxL)
        Ferro[:,1:4],cat[0,1:4] = unwrapBond(Ferro[:,1:4],cat[0,1:4],boxL)

        Twater = cur_pos[(cur_pos[:,0]==1)|(cur_pos[:,0]==2)]
        for ctw, ii in enumerate(config[4:],0):
            if not np.isnan(ii):
                tmp = Twater[int(ii)*3:int(ii)*3+3]
                tmp[1:,1:4], tmp[0,1:4] = unwrapBond(tmp[1:,1:4],tmp[0,1:4],boxL)
                if ctw == 0:
                    Water = tmp
                else:
                    Water = np.vstack((Water,tmp))
        Water[:,1:4],cat[0,1:4] = unwrapBond(Water[:,1:4],cat[0,1:4],boxL)


        # move the cation to the origin
        Ferri[:,1:4] -= cat[0,1:4]
        Ferro[:,1:4] -= cat[0,1:4]
        Water[:,1:4] -= cat[0,1:4]
        cat[:,1:4] -= cat[0,1:4]

        outfile = open("MD2DFT"+str(ct)+".xyz","w")
        outfile.write(str(cat[0,0])+" "+str(cat[0,1])+" "+str(cat[0,2])+" "+str(cat[0,3])+"\n")
        for ii in range(13):
            outfile.write(str(Ferri[ii,0])+" "+str(Ferri[ii,1])+" "+str(Ferri[ii,2])+" "+str(Ferri[ii,3])+"\n")
        for ii in range(13):
            outfile.write(str(Ferro[ii,0])+" "+str(Ferro[ii,1])+" "+str(Ferro[ii,2])+" "+str(Ferro[ii,3])+"\n")
        for ii in range(len(Water)):
            outfile.write(str(Water[ii,0])+" "+str(Water[ii,1])+" "+str(Water[ii,2])+" "+str(Water[ii,3])+"\n")
        outfile.close()

