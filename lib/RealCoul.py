import numpy as np
import sys, scipy.special
np.set_printoptions(threshold=np.inf,suppress=True)				#force print entire array
np.seterr('raise')


try:
    from lib.read_coord import get_wrapped_coord, r_sqrd
    from lib.math_lib import ID2MolData, nearest_neigh_rsqrd
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def DoMeanVar(lammpstrj,K_lim):
    fname = open(lammpstrj,"r")
    cur_xyz, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(fname)

    # select redox
    i = 5
    Ferri = cur_xyz[( cur_xyz[:,0] ==i )|( cur_xyz[:,0] ==i+1 )|( cur_xyz[:,0] ==i+2 )]

    preambient = cur_xyz[( cur_xyz[:,0] !=i )&( cur_xyz[:,0] !=i+1 )&( cur_xyz[:,0] !=i+2 )]

    # replace water and redox charges
    Ferri[Ferri[:,0]==i,0] = -0.834
    Ferri[Ferri[:,0]==i+1,0] = 0.402
    Ferri[Ferri[:,0]==i+2,0] = -0.763

    j = 8
    preambient[preambient[:,0]==1,0] = -0.8476
    preambient[preambient[:,0]==2,0] = 0.4238
    preambient[preambient[:,0]==j,0] = -1.15
    preambient[preambient[:,0]==j+1,0] = 0.431
    preambient[preambient[:,0]==j+2,0] = -0.906

    sel_Ferri = Ferri[0:13,:]
    preambient = np.vstack((Ferri[13:,:],preambient))
    # print(sel_Ferri,preambient)

    ambient = np.copy(preambient)
    ### Now add periodic boxes ###
    for xx in [-1,0,1]:
        for yy in [-1,0,1]:
            for zz in [-1,0,1]:
                ambient1 = np.copy(preambient)
                ambient1[:,1:4] = preambient[:,1:4] + np.array([xx,yy,zz])*boxL
                ambient = np.vstack((ambient,ambient1))

                # print(preambient[0,1:4],np.array([xx,yy,zz])*boxL,ambient1[0,1:4],preambient[0,1:4] + np.array([xx,yy,zz])*boxL)

    # print(ambient)

    coul = np.zeros(13)
    lst = []
    for ii in range(len(ambient)):
        diff = sel_Ferri[:,1:4] - ambient[ii,1:4]
        dis = np.sqrt(nearest_neigh_rsqrd(diff,3*boxL))

    #     lst.append(dis.tolist())
    # # plot data for sanity check
    # print(np.sqrt(boxL[0]**2+2*boxL[0]**2)/2)
    # tohist = np.array(lst)
    # a,b = np.shape(tohist)
    # print(np.shape(tohist))
    # tohist = np.reshape(tohist,(a*b,1))
    # import matplotlib.pyplot as plt
    # plt.hist(tohist,bins=100)
    # plt.show()

        if dis[0]<K_lim:
            if ii == 5:
                pass
                # print(coul, sel_Ferri[:,0] *  ambient[ii,0]/dis)
            coul += sel_Ferri[:,0] *  ambient[ii,0]/dis
    print(coul)
