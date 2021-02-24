# formula for claulating total = (T0ct-1)*corr_t+501
# python3 main.py --routine 8 --lammpstrj equ_nve_wC.lammpstrj --RefFile Ref_O_dist --Tstep 10001 --stype 4 --itype 1 --scaling 10
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
    from lib.read_coord import get_wrapped_coord_Madel
    from lib.math_lib import r_iN, adjustL_conv2UV
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def DoLambda(lammpstrj,RefFile,Tstep,T0ct,corr_t,tau,stype,itype,scaling):
    if (T0ct-1)*corr_t+tau+1 > Tstep:
        print("The T0ct, corr_t and tau you have chosen is larger than Ttep. Choose again!!!")
    else:
        # read in boxsize from traj file
        xyz_0 = open(lammpstrj,"r")
        cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord_Madel(xyz_0)
        lambda_mean_var = np.zeros((int(round(boxL[0]/2*scaling))+1,T0ct))      # 0 for mean; 1 for variance
        xyz_0.close()

        Ref_Cdist = np.loadtxt(RefFile)
        # print(Ref_Cdist[530])               # will print probability of having Madelung potential value of -30
        Ref_Cdist = Ref_Cdist[:,1]          # madelung to coulomb conversion formula => ind = Coul*2/q+500

        for cur_t0 in range(T0ct):
            cur_lambda_hist = np.zeros((int(round(boxL[0]/2*scaling))+1,2))     # 0 for -log(P); 1 for ct
            xyz_in = open(lammpstrj,"r")
            print("***** doing t0#",cur_t0,xyz_in.tell())
            for ii in range(cur_t0*corr_t):                 # skip this many timestep before the new t0
                get_wrapped_coord_Madel(xyz_in)
            for cur_t in range(tau+1):
                cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord_Madel(xyz_in)
                if stype == -1:                 # -1 for pure water box
                    s_pos = np.array([boxL[0]/2,boxL[1]/2,boxL[2]/2,0])
                else:
                    s_pos = cur_pos[cur_pos[:,0]==stype][:,1:5]         # index 4 here is the Columb PE
                i_pos = cur_pos[cur_pos[:,0]==itype][:,1:5]
                # print(s_pos,i_pos[:5],i_pos[-5:])

                r_sqrd = r_iN(i_pos[:,0:3],s_pos[:,0:3],boxL)
                ind = np.around(np.sqrt(r_sqrd)*scaling)
                ind = np.reshape(ind, (len(ind),1))
                Madel = np.reshape(i_pos[:,3], (len(i_pos[:,3]),1))
                ind_mad = np.concatenate((ind,Madel), axis = 1)
                # print(ind_mad)
                # print(ind_mad[ind_mad[:,0]/scaling<5])
                # print(np.argwhere(ind_mad[:,0]/scaling<5))

                for ii in range(len(ind_mad)):
                    if int(ind_mad[ii,0]) > boxL[0]/2*scaling:
                        pass
                    else:
                        coul_ind = int(round(ind_mad[ii,1]))+500
                        # print(coul_ind,Ref_Cdist[coul_ind])
                        prob = Ref_Cdist[coul_ind] if Ref_Cdist[coul_ind]>0 else 1
                        # print(-np.log(prob),ind_mad[ii,0])
                        cur_lambda_hist[int(ind_mad[ii,0]),0] += -np.log(prob)
                        cur_lambda_hist[int(ind_mad[ii,0]),1] += 1
            xyz_in.close()
            cur_lambda_hist[cur_lambda_hist[:,1]==0,1] = 1
            # print(cur_lambda_hist)
            lambda_mean_var[:,cur_t0] = cur_lambda_hist[:,0]/cur_lambda_hist[:,1]


    outfile="lambda_sitype"+str(stype)+"_"+str(itype)+"T"+str(Tstep)+"T0ct"+str(T0ct)+"corr_t"+str(corr_t)+"tau"+str(tau)+"scaling"+str(scaling)+".out"
    with open(outfile,"w") as fout:
        for ii in range(len(lambda_mean_var)):
            fout.write(str(ii/scaling)+" "+str(np.mean(lambda_mean_var[ii,:]))+" "+str(np.std(lambda_mean_var[ii,:]))+"\n")
