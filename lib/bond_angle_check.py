try:
    np.array()
except NameError:
    print("importing numpy")
    import numpy as np
try:
    sys
except NameError:
    print("importing sys")
    import sys
np.set_printoptions(threshold=sys.maxsize)
try:
    from lib.read_coord import get_wrapped_coord, truncate,r_sqrd
    from lib.math_lib import adjustL_conv2UV
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def doBondAnlgeHist(Tstep,lammpstrj,itype,jtype,scaling,R_lim):
    # get i,j count from the files for initializing arrays
    posfile0 = open(lammpstrj,"r")
    cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord(posfile0)
    posfile0.close()
    ict = len(cur_xyz[cur_xyz[:,0]==itype])
    jct = len(cur_xyz[cur_xyz[:,0]==jtype])
    jiRatio = int(jct/ict)
    # print(ict,jct,jiRatio)

    angle_hist = np.zeros((int(181*scaling),1))
    bond_hist = np.zeros((int(R_lim)*scaling,1))
    posfile = open(lammpstrj,"r")
    for cur_t in range(Tstep):
        cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord(posfile)
        cur_xyz_ij = cur_xyz[ (cur_xyz[:,0]==itype) |(cur_xyz[:,0]==jtype) ]
        # print(cur_xyz_ij)

        ji_uv = np.zeros((jiRatio,3))
        for ii in range(ict):
            ii_ind = ii*(jiRatio+1)
            i_pos = cur_xyz_ij[ii_ind][1:4]
            for jj in range(jiRatio):
                j_pos = cur_xyz_ij[ii_ind+jj+1][1:4]

                # do bond hist
                bondL = np.sqrt(r_sqrd(j_pos,i_pos, boxL))
                ind = int(truncate(bondL,np.log10(scaling))*scaling)
                bond_hist[ind] += 1

                # store unit vector between j,i
                ji_uv[jj] = adjustL_conv2UV(j_pos-i_pos,boxL)

            # do angle hist
            for ii in range(jiRatio):
                for jj in range(ii+1,jiRatio):
                    dot_pdt = np.sum(ji_uv[ii] * ji_uv[jj])

                    # to prevent round-off error in handling dot_pdt
                    error_threshold = 1e-16
                    if abs(dot_pdt)-1 > error_threshold:
                        print("current Tstep is:::",cur_t)
                        print("*** correcting for round-off errors, converting ",dot_pdt," to ",dot_pdt/abs(dot_pdt),"***")
                        dot_pdt = dot_pdt/abs(dot_pdt)

                    angle = np.arccos(dot_pdt)*180/np.pi
                    ind = int(truncate(angle,np.log10(scaling))*scaling)
                    angle_hist[ind] += 1

    tmpstr = "i_j"+str(itype)+"_"+str(jtype)+"_scaling"+str(scaling)+"_Rlim"+str(R_lim)
    np.savetxt(tmpstr+"BondHist.out",bond_hist)
    np.savetxt(tmpstr+"AngleHist.out",angle_hist)

    import matplotlib
    from matplotlib import pyplot as plt
    tmpstr = "i,j ="+str(itype)+","+str(jtype)
    plt.figure()
    x=np.array([ii for ii in range(int(R_lim)*scaling)])/scaling
    plt.plot(x,bond_hist,"r-o",ms=2,label=tmpstr)
    plt.xlabel(r'Bond Length [$\AA$]')
    plt.ylabel("Count")
    plt.xlim(1,R_lim)
    plt.legend(frameon=False)
    # plt.show()

    plt.figure()
    x=np.array([ii for ii in range(int(181*scaling))])/scaling
    plt.plot(x,angle_hist,"r-o",ms=2,label=tmpstr)
    plt.xlabel(r'Angle [Degree]')
    plt.ylabel("Count")
    plt.xlim(85,185)
    plt.legend(frameon=False)
    plt.show()
