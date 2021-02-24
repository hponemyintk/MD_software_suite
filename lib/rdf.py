try:
    np.array()
except NameError:
    print("importing numpy")
    import numpy as np
    import time
try:
    sys
except NameError:
    print("importing sys")
    import sys
np.set_printoptions(threshold=sys.maxsize)
try:
    from lib.read_coord import get_wrapped_coord, get_wrapped_xyz,r_sqrd
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def doRDF(lammpstrj,boxsize,Tstep,TskipB,itype,jtype,r_cut,scaling):
    r_cut = int(r_cut)
    # get i,j count from the files for initializing arrays
    posfile0 = open(lammpstrj,"r")
    cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord(posfile0) if isinstance(boxsize,bool) else get_wrapped_xyz(posfile0,boxsize)
    posfile0.close()
    ict = len(cur_xyz[cur_xyz[:,0]==itype])
    jct = len(cur_xyz[cur_xyz[:,0]==jtype])

    # keep track of time
    t0 = time.time()
    cur_prcnt = 0

    rdf_hist = np.zeros((r_cut*scaling+1))	# storing hist count
    posfile = open(lammpstrj,"r")
    ###################
    ### do rdf here ###
    ###################
    cur_step = 0
    while cur_step < Tstep:
        cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord(posfile) if isinstance(boxsize,bool) else get_wrapped_xyz(posfile,boxsize)
        cur_xyz_i = cur_xyz[ (cur_xyz[:,0]==itype) ]
        cur_xyz_j = cur_xyz[ (cur_xyz[:,0]==jtype) ]

        for i_rdf in range(0,ict):
            print(cur_step,i_rdf)
            for j_rdf in range(0,jct):
                dis = np.sqrt(r_sqrd(cur_xyz_i[i_rdf][1:4],cur_xyz_j[j_rdf][1:4], boxL))
                if dis < r_cut:
                    cur_ind = int(round(dis*scaling))
                    rdf_hist[cur_ind] += 1

        if (TskipB!=0) :
            for i in range(0,TskipB):
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord(posfile) if isinstance(boxsize,bool) else get_wrapped_xyz(posfile,boxsize)
            cur_step = cur_step + TskipB
        cur_step += 1

        # check whether to output time info
        TskipI = 0
        if cur_step*100/(Tstep-TskipI) >= cur_prcnt:
            if cur_prcnt == 0:
                print("*** Outputting time data. ***\nTstep,TskipI,TskipB,lammpstrj,itype,jtype,boxsize,scaling",
                      Tstep,TskipI,TskipB,lammpstrj,itype,jtype,boxsize,scaling,
                      "\n#######################################\nPercent[%] Time[s]")
            print("{0:3.2f}       {1:3.2f}".format(cur_step*100/(Tstep-TskipI), time.time()-t0))
            cur_prcnt += 10
            t0 = time.time()

    tmpstr = "i_"+str(itype)+"j_"+str(jtype)+"T_"+str(Tstep)+"r_cut_"+str(r_cut)+"scaling_"+str(scaling)+"TskipB_"+str(TskipB)+".rdf"
    outfile=open(tmpstr,"w")
    tmpstr = "# R rdf hist_ct\n"
    outfile.write(tmpstr)
    ### Histogramming done, let's do rdf
    for i in range(1,r_cut*scaling+1):
        tmpstr = str(i/float(scaling)) + " " + str( float(rdf_hist[i])/(ict*int(Tstep/(TskipB+1))*4*np.pi*(i/float(scaling))**2*(1/float(scaling))*(jct/(boxL[0]*boxL[1]*boxL[2]))) ) + " " + str(rdf_hist[i]) + "\n"
        outfile.write(tmpstr)
    outfile.close()

