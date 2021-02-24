# formula for claulating total = (T0ct-1)*corr_t+501
# python3 main.py --routine 7 --lammpstrj equ_nve_woWater.lammpstrj --Tstep 10001 --stype 3 --itype 8 --R_lim 5.5 --threshold 20 --scaling 10
### doubled checked with excel and passed! Look at 20191001_MSD_FECN6.xlsx
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
    from lib.read_coord import get_wrapped_coord, truncate
    from lib.math_lib import r_iN, adjustL_conv2UV
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def get_topology(s_pos,full_redox,boxL,R_lim,config_hist,scaling):
    v_redox = np.zeros((6,3))
    v_Sredox = np.zeros((1,3))
    for cur_redox in range(int(len(full_redox)/13)):
        ### find the cation in the shell
        diff = r_iN(s_pos,full_redox[cur_redox*13],boxL)
        ind_list = np.where(diff<R_lim**2)

        #### create unit vector between Fe-C
        for CC in range(1,7):
            v_redox[CC-1] = adjustL_conv2UV(full_redox[cur_redox*13+CC]-full_redox[cur_redox*13],boxL)

        for ii in range(len(ind_list[-1])):                 # here ind_list[-1] has the same effect as ind_list[0] as we are just trying to remove extra [] around the data
            ### do spherical sum between Fe-C and Fe-stype unit vectors
            v_Sredox = adjustL_conv2UV(s_pos[ind_list[0][ii]]- (full_redox[cur_redox*13]),boxL)
            spherical_sum = np.sum(abs(np.sum(v_redox*v_Sredox,axis=1)))
            ind = int(truncate(spherical_sum, np.log10(scaling))*scaling)
            config_hist[ind] += 1
    return config_hist

# for t* use 2ps and tau = 100ps should be more than enough
def DoHist(lammpstrj,Tstep,stype,itype,R_lim,threshold,scaling):
    cur_step = 0
    xyz_in = open(lammpstrj,"r")
    config_hist = np.zeros((4*scaling))
    while cur_step < Tstep:
        cur_step += 1
        cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(xyz_in)
        s_pos=cur_pos[cur_pos[:,0]==stype][:,1:4]
        full_redox = cur_pos[(cur_pos[:,0]==itype) | (cur_pos[:,0]==itype+1) | (cur_pos[:,0]==itype+2) ][:,1:4]
        get_topology(s_pos,full_redox,boxL,R_lim,config_hist,scaling)
    outfile="topo_Hist_sitype"+str(stype)+"_"+str(itype)+"_"+"R_lim"+str(R_lim)+"Tthreshold"+str(threshold)+"T"+str(Tstep)+"scaling"+str(scaling)+".out"
    with open(outfile,"w") as fout:
        for ii in range(len(config_hist)):
            fout.write(str(ii/scaling)+" "+str(config_hist[ii])+"\n")
