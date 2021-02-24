# python dis_corr.py --lammpstrj equ_nve_woWater.lammpstrj --Tstep 100 --T0ct 3 --corr_t 10 --tau 70
try:
    np.array()
except NameError:
    import numpy as np
try:
    sys
except NameError:
    import sys
np.set_printoptions(threshold=sys.maxsize)
try:
    from read_coord import get_wrapped_coord
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")


def dis_corr(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau):
    # <S1*S2>-<S1><S2>
    corr = np.zeros(tau+1)                            # store sum of corr values for each t0 which will be averaged over by number of t0 in the end
    corr = np.zeros((tau+1,2))                            # store sum of corr values for each t0 which will be averaged over by number of t0 in the end
    for cur_t0 in range(T0ct):
        xyz_in = open(lammpstrj,"r")
        print("***** doing t0#",cur_t0,xyz_in.tell())
        for ii in range(cur_t0*corr_t):                # skip this many timestep before the new t0
            # print(">>>> Printing ii #",ii,xyz_in.tell())
            get_wrapped_coord(xyz_in)
        for cur_t in range(tau+1):
            # print("#### printing cur_t",cur_t,xyz_in.tell())
            cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(xyz_in)
            cur_pos=cur_pos[:,1:4]
            if cur_t==0:
                pos0 = np.copy(cur_pos)
            # do dot product here
            dpt_means = np.sum(np.sum(pos0/Ntotal,axis=0)*np.sum(cur_pos/Ntotal,axis=0),axis=0)     # axis goes Row, Column, Depth, and fill them in that order. 1D array only has rows
            m_S1S2 = np.sum(cur_pos*pos0/Ntotal)
            # corr[cur_t] += m_S1S2#-dpt_means
            # corr[cur_t] += m_S1S2-dpt_means

            #debug
            corr[cur_t][0] += m_S1S2
            corr[cur_t][1] += dpt_means
        xyz_in.close()
    corr = corr/float(T0ct)

    # outfile="dis_corr_T"+str(Tstep)+"TskipI"+str(TskipI)+"T0ct"+str(T0ct)+"corr_t"+str(corr_t)+"tau"+str(tau)+".corr"
    # out_fname=open(outfile,"w")
    # for i in range(tau+1):
    #     tmpstr=str(i)+" "+str(corr[i])+" "+str(corr[i]/corr[0])+"\n"
    #     out_fname.write(tmpstr)
    # out_fname.close()

    # debug
    outfile="dis_corr_T"+str(Tstep)+"TskipI"+str(TskipI)+"T0ct"+str(T0ct)+"corr_t"+str(corr_t)+"tau"+str(tau)+".corr"
    out_fname=open(outfile,"w")
    for i in range(tau+1):
        tmpstr=str(i)+" "+str((corr[i][0]-corr[i][1])/(corr[0][0]-corr[0][1]))+" "+str(corr[i][0])+" "+str(corr[i][1])+"\n"
        out_fname.write(tmpstr)
    out_fname.close()
