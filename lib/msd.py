# formula for claulating total = (T0ct-1)*corr_t+501
# python dis_corr.py --lammpstrj equ_nve_woWater.lammpstrj --Tstep 100 --T0ct 3 --corr_t 10 --tau 70
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
    from lib.read_coord import get_unwrapped_coord
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")


def doMSD(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau,stype):
    if (T0ct-1)*corr_t+tau+1 > Tstep:
        print("The T0ct, corr_t and tau you have chose in larger than Ttep. Choose again!!!")
    else:
        corr = np.zeros(tau+1)                              # store sum of corr values for each t0 which will be averaged over by number of t0 in the end
        for cur_t0 in range(T0ct):
            xyz_in = open(lammpstrj,"r")
            print("***** doing t0#",cur_t0,xyz_in.tell())
            for ii in range(cur_t0*corr_t):                 # skip this many timestep before the new t0
                get_unwrapped_coord(xyz_in)
            for cur_t in range(tau+1):
                cur_pos, tcur_Tstep, Ntotal, boxL = get_unwrapped_coord(xyz_in)
                if not stype:                               # check whether we have selected to look at only the specific type of atoms 
                    cur_pos=cur_pos[:,1:4]
                else:
                    cur_pos=cur_pos[cur_pos[:,0]==stype][:,1:4]
                if cur_t==0:
                    pos0 = np.copy(cur_pos)
                    atom_ct = len(cur_pos)
                diff = cur_pos - pos0
                msd = diff*diff
                msd = np.sum(msd)/atom_ct

                corr[cur_t] += msd
            xyz_in.close()
        corr = corr/float(T0ct)

        # write to file
        outfile="msd_stype"+str(stype)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"T0ct"+str(T0ct)+"corr_t"+str(corr_t)+"tau"+str(tau)+".out"
        out_fname=open(outfile,"w")
        for i in range(tau+1):
            tmpstr=str(i)+" "+str(corr[i])+"\n"
            out_fname.write(tmpstr)
        out_fname.close()
