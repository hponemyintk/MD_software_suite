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
    from lib.read_coord import get_unwrapped_coord, get_wrapped_coord
    from lib.math_lib import r_iN
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def flag_solvation_ions(xyzsol_in,corr_t,tau,stype,itype,jtype,R_lim):
    shell_ions = 0          # store sum of shell_ions values for each t0 which will be normalized by the length of the time window (tau+1) in the end
    for cur_t in range(tau+1):
        cur_pos1, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(xyzsol_in)
        redox_pos=cur_pos1[(cur_pos1[:,0]==itype) | (cur_pos1[:,0]==jtype)][:,1:4]
        s_pos=cur_pos1[cur_pos1[:,0]==stype][:,1:4]
        if isinstance(shell_ions,int):
            shell_ions = np.zeros(len(s_pos))

        for cur_redox in range(len(redox_pos)):
            diff = r_iN(s_pos,redox_pos[cur_redox],boxL)
            shell_ions[diff<R_lim**2] +=1
            # print("cur_t:::\n",cur_t,"\ncur_redox",cur_redox,"\nredox_pos:::\n",redox_pos,"\ns_pos:::\n",s_pos,"\ndiff:::\n",diff,"\nshell_ions:::\n",shell_ions,"\n")
    shell_ions /= (tau+1.0)
    # print("\nshell_ions:::\n",shell_ions)
    return shell_ions


def doMSD_pop(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau,stype,itype,jtype,R_lim,threshold):
    if itype==False or jtype==False or stype==False:
        print("No s,i,j type define!!! Terminating analysis...")
    elif (T0ct-1)*corr_t+tau+1 > Tstep:
        print("The T0ct, corr_t and tau you have chose in larger than Ttep. Choose again!!!")
    else:
        ct_shell_skip = 0               # to accout for t0 that doesn't contribute to averaging
        ct_bulk_skip = 0
        shell_corr = np.zeros(tau+1)                              # store sum of corr values for each t0 which will be averaged over by number of t0 in the end
        bulk_corr = np.zeros(tau+1)                              # store sum of corr values for each t0 which will be averaged over by number of t0 in the end
        for cur_t0 in range(T0ct):
            xyz_in = open(lammpstrj,"r")
            print("***** doing t0#",cur_t0,xyz_in.tell())
            for ii in range(cur_t0*corr_t):                 # skip this many timestep before the new t0
                get_unwrapped_coord(xyz_in)

            # open a second file pointer for historgramming solvation shell ions
            try:
                xyzsol_in
            except NameError:
                xyzsol_in = open(lammpstrj,"r")
            xyzsol_in.seek(xyz_in.tell())
            shell_ions = flag_solvation_ions(xyzsol_in,corr_t,tau,stype,itype,jtype,R_lim)
            if len(shell_ions>threshold)==0:       # exclude the t0 that doesn't contains either any shell or bulk (stype) atom
                ct_shell_skip += 1
            if len(shell_ions<threshold)==0:
                ct_bulk_skip += 1

            for cur_t in range(tau+1):
                cur_pos, tcur_Tstep, Ntotal, boxL = get_unwrapped_coord(xyz_in)
                cur_pos = cur_pos[cur_pos[:,0]==stype][:,1:4]
                if cur_t==0:
                    pos0 = np.copy(cur_pos)
                    atom_ct = len(cur_pos)
                diff = cur_pos - pos0
                msd = diff*diff
                # print("\ndiff:::\n",diff,"\nmsd\n",msd,"\nin out num::\n",len(msd[shell_ions>threshold]),len(msd[shell_ions<threshold]))
                add_it1 = np.sum(msd[shell_ions>threshold])/float(len(msd[shell_ions>threshold])) if len(msd[shell_ions>threshold])>0 else 0
                shell_corr[cur_t] += add_it1
                add_it2 = np.sum(msd[shell_ions<threshold])/float(len(msd[shell_ions<threshold])) if len(msd[shell_ions<threshold])>0 else 0
                bulk_corr[cur_t] += add_it2
            xyz_in.close()
        xyzsol_in.close()
        shell_corr = shell_corr/float(T0ct-ct_shell_skip)
        bulk_corr = bulk_corr/float(T0ct-ct_bulk_skip)
        # print("shell MSD",shell_corr)

        # write to file
        outfile="shell_msd_sijtype"+str(stype)+"_"+str(itype)+"_"+str(jtype)+"R_lim"+str(R_lim)+"Tthreshold"+str(threshold)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"T0ct"+str(T0ct)+"corr_t"+str(corr_t)+"tau"+str(tau)+".out"
        outfile2="bulk_msd_sijtype"+str(stype)+"_"+str(itype)+"_"+str(jtype)+"R_lim"+str(R_lim)+"Tthreshold"+str(threshold)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"T0ct"+str(T0ct)+"corr_t"+str(corr_t)+"tau"+str(tau)+".out"
        out_fname=open(outfile,"w")
        out_fname2=open(outfile2,"w")
        for i in range(tau+1):
            tmpstr=str(i)+" "+str(shell_corr[i])+"\n"
            out_fname.write(tmpstr)
            tmpstr=str(i)+" "+str(bulk_corr[i])+"\n"
            out_fname2.write(tmpstr)
        out_fname.close()
        out_fname2.close()

