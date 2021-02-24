import numpy as np
import sys, scipy.special
np.set_printoptions(threshold=np.inf,suppress=True)				#force print entire array
np.seterr('raise')


try:
    from lib.read_coord import get_Coul_coord, r_sqrd
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def DoMeanVar(lammpstrj,CoulDump,Tstep,TskipI,TskipB,x_1sideL,nEmoment):
    nEmoment = int(nEmoment)

    ### CONSTANTS ###
    E_conv_const = 332.06371			# copied from LAMMPS update.cpp [convert qq/r to kcal/mol]
    echarge = 1 # 1.60217662e-19                        # [Coulomb]
    kb = 0.001987204 									# [kcal/mol/K]
    temp = 300											# [K]

    ### INPUT PARAMS ###
    itype = 5
    jtype = 8                               ### this code assume j has higher index than i in the lammps output file
    icharge = -3 * echarge
    jcharge = -4 * echarge
    dq = icharge - jcharge                  ### *** the way this code is set up, dq = oxidized - reduced, always!!! So, make sure the j is more negatively charge than I when selecting i,j atom type!!!
    coulfile = open(CoulDump, "r")
    posfile = open(lammpstrj,"r")
    file_name_tail = str(Tstep)+"_TskipB" +str(TskipB)+"_upTonthMoment_MDv"

    # get i,j count from the files for initializing arrays
    coulfile0 = open(CoulDump, "r")
    posfile0 = open(lammpstrj,"r")
    cur_xyz, cur_Coul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile0,posfile0)
    coulfile0.close()
    posfile0.close()
    ict = len(cur_Coul[cur_Coul[:,0]==itype])
    jct = len(cur_Coul[cur_Coul[:,0]==jtype])
    mean_dE_A = 0
    mean_dE_B = 0
    dE_ct =0

    spherical_lim = boxL[0]/2	# 1e6 # boxL[0]/2

    ### For r distribution for each dE ###
    print(spherical_lim)
    rDist_per_dE_A = np.zeros((2*x_1sideL+1,int(round(spherical_lim))+1))           # r from 0 to spherical_lim (as dis always less than spehical_lim)
    rDist_per_dE_B = np.zeros((2*x_1sideL+1,int(round(spherical_lim))+1))
    rDist_mean_var_dE_A = np.zeros((int(round(spherical_lim))+1,nEmoment+2))           # index 0 is mean 1..n is 2..2n moment
    rDist_mean_var_dE_B = np.zeros((int(round(spherical_lim))+1,nEmoment+2))
    ### correlation Hist
    rDist_mean_var_oxd = np.zeros((int(round(spherical_lim))+1,2))           # index 0 is mean 1..n is 2..2n moment
    rDist_mean_var_red = np.zeros((int(round(spherical_lim))+1,2))


    ### Histogram to get dE distribution for state A & B ###
    A_hist = np.zeros((2*x_1sideL+1,1))
    B_hist = np.zeros((2*x_1sideL+1,1))

    cur_step=0
    while cur_step < Tstep:
        if(cur_step<TskipI):
            cur_xyz, cur_Coul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile,posfile)
        else:
            cur_xyz, cur_Coul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile,posfile)

            cur_xyz_ferri = cur_xyz[ (cur_xyz[:,0]==itype) ]
            cur_xyz_ferro = cur_xyz[ (cur_xyz[:,0]==jtype) ]
            cur_Coul_ferri = cur_Coul[ (cur_Coul[:,0]==itype) ]
            cur_Coul_ferro = cur_Coul[ (cur_Coul[:,0]==jtype) ]
            # print cur_xyz_ferri, cur_xyz_ferro, cur_Coul_ferri, cur_Coul_ferro, timestamp, Ntotal, dq
            # print cur_xyz, cur_Coul, timestamp, Ntotal

            for i_AB in range(0,ict):
                for j_AB in range(0,jct):
                    dis = np.sqrt(r_sqrd(cur_xyz_ferri[i_AB][1:4],cur_xyz_ferro[j_AB][1:4], boxL))
                    ferri_ML = cur_Coul_ferri[i_AB][1]/icharge*2.
                    ferro_ML = cur_Coul_ferro[j_AB][1]/jcharge*2.

                    # calculate and store dE_A and dE_B
                    cur_dE_A = dq *(E_conv_const * dq/dis + ferri_ML - ferro_ML)
                    cur_dE_B = dq *(-E_conv_const * dq/dis + ferro_ML - ferri_ML)
                    # print dis, E_conv_const * dq**2/dis, dq*ferri_ML-dq*ferro_ML, dq*ferri_ML, dq*ferro_ML

                    if dis < spherical_lim:          # only consider the dE with dis less than half a box length away from the central atom
                        mean_dE_A += cur_dE_A
                        mean_dE_B += cur_dE_B

                        # do for state B
                        dE_B_ind = int(round(cur_dE_B+x_1sideL))
                        B_hist[dE_B_ind] += 1
                        rDist_per_dE_B[dE_B_ind][int(round(dis))] += 1                  # for Histogramming rDist_per_dE_B
                        rDist_mean_var_dE_B[int(round(dis))][0] += cur_dE_B

                        # do for state A
                        dE_A_ind = int(round(cur_dE_A+x_1sideL))
                        A_hist[dE_A_ind] += 1
                        rDist_per_dE_A[dE_A_ind][int(round(dis))] += 1                  # for Histogramming rDist_per_dE_A
                        rDist_mean_var_dE_A[int(round(dis))][0] += cur_dE_A

                        rDist_mean_var_oxd[int(round(dis))][0] += ferri_ML
                        rDist_mean_var_red[int(round(dis))][0] += ferro_ML
                        dE_ct += 1

            if (TskipB!=0) :
                for i in range(0,TskipB):
                    cur_xyz, cur_Coul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile,posfile)
                cur_step = cur_step + TskipB
            cur_step += 1
            # print "cur_step is", cur_step
    coulfile.close()
    posfile.close()



    #################################
    ### second pass to get moment ###
    #################################
    var_dE_A=0
    var_dE_B=0

    ### Final Means ###
    mean_dE_A/=dE_ct
    mean_dE_B/=dE_ct
    for i_mean in range(int(round(spherical_lim))+1):                           # do avg for each 
        if sum(rDist_per_dE_A[:,i_mean]) != 0 and sum(rDist_per_dE_B[:,i_mean]) !=0:
            rDist_mean_var_dE_A[i_mean][0] /= sum(rDist_per_dE_A[:,i_mean])
            rDist_mean_var_dE_B[i_mean][0] /= sum(rDist_per_dE_B[:,i_mean])
            rDist_mean_var_oxd[i_mean][0] /= sum(rDist_per_dE_A[:,i_mean])
            rDist_mean_var_red[i_mean][0] /= sum(rDist_per_dE_A[:,i_mean])

    coulfile1 = open(CoulDump, "r")
    posfile1 = open(lammpstrj,"r")
    cur_step=0
    while cur_step < Tstep:
        if(cur_step<TskipI):
            cur_xyz, cur_Coul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile1,posfile1)
        else:
            cur_xyz, cur_Coul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile1,posfile1)

            cur_xyz_ferri = cur_xyz[ (cur_xyz[:,0]==itype) ]
            cur_xyz_ferro = cur_xyz[ (cur_xyz[:,0]==jtype) ]
            cur_Coul_ferri = cur_Coul[ (cur_Coul[:,0]==itype) ]
            cur_Coul_ferro = cur_Coul[ (cur_Coul[:,0]==jtype) ]

            for i_AB in range(0,ict):
                for j_AB in range(0,jct):
                    dis = np.sqrt(r_sqrd(cur_xyz_ferri[i_AB][1:4],cur_xyz_ferro[j_AB][1:4], boxL))
                    ferri_ML = cur_Coul_ferri[i_AB][1]/icharge*2.
                    ferro_ML = cur_Coul_ferro[j_AB][1]/jcharge*2.

                    # calculate and store dE_A and dE_B
                    cur_dE_A = dq *(E_conv_const * dq/dis + ferri_ML - ferro_ML)
                    cur_dE_B = dq *(-E_conv_const * dq/dis + ferro_ML - ferri_ML)

                    if dis < spherical_lim:          # only consider the dE half a box length away from the central atom
                        var_dE_A += (cur_dE_A - mean_dE_A)**2
                        var_dE_B += (cur_dE_B - mean_dE_B)**2

                        for i_moment in range(1,nEmoment+1):
                            # do for state A
                            rDist_mean_var_dE_A[int(round(dis))][i_moment] += (cur_dE_A - rDist_mean_var_dE_A[int(round(dis))][0])**(i_moment)
                            # do for state B
                            rDist_mean_var_dE_B[int(round(dis))][i_moment] += (cur_dE_B - rDist_mean_var_dE_B[int(round(dis))][0])**(i_moment)

                        # do correlation
                        rDist_mean_var_oxd[int(round(dis))][1] += (ferri_ML - rDist_mean_var_oxd[int(round(dis))][0])**2
                        rDist_mean_var_red[int(round(dis))][1] += (ferro_ML - rDist_mean_var_red[int(round(dis))][0])**2
                        rDist_mean_var_dE_A[int(round(dis))][nEmoment+1] += (ferri_ML - rDist_mean_var_oxd[int(round(dis))][0]) * (ferro_ML - rDist_mean_var_red[int(round(dis))][0])

            if (TskipB!=0) :
                for i in range(0,TskipB):
                    cur_xyz, cur_Coul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile1,posfile1)
                cur_step = cur_step + TskipB
            cur_step += 1
    coulfile1.close()
    posfile1.close()

    ### Final vairances ###
    var_dE_A/= dE_ct
    var_dE_B/= dE_ct
    for i_var in range(int(round(spherical_lim))+1):
        if sum(rDist_per_dE_A[:,i_var]) != 0 and sum(rDist_per_dE_B[:,i_var]) !=0:
            rDist_mean_var_dE_A[i_var][1:nEmoment+2]/=sum(rDist_per_dE_A[:,i_var])
            rDist_mean_var_dE_B[i_var][1:nEmoment+2]/=sum(rDist_per_dE_B[:,i_var])
            rDist_mean_var_oxd[i_var][1]/=sum(rDist_per_dE_A[:,i_var])
            rDist_mean_var_red[i_var][1]/=sum(rDist_per_dE_A[:,i_var])
        # do correlation
        tmp_SD_oxd = np.sqrt(rDist_mean_var_oxd[i_var][1])
        tmp_SD_red = np.sqrt(rDist_mean_var_red[i_var][1])
        if tmp_SD_oxd != 0 and tmp_SD_red != 0:
            rDist_mean_var_dE_A[i_var][nEmoment+1]/=(tmp_SD_oxd*tmp_SD_red)


    ### PRINT R_DIST PER dEs ###
    tmpstr = "rDists_per_dE"+ file_name_tail + ".Hist"
    R_out = open(tmpstr,"w")
    tmpstr = "# dEs"
    for ii_R in range(int(round(spherical_lim)+1)):
        tmpstr += " A_"+str(ii_R)
    for ii_R in range(int(round(spherical_lim)+1)):
        tmpstr += " B_"+str(ii_R)
    tmpstr += " \n"
    R_out.write(tmpstr)
    for EE in range(0,2*x_1sideL+1):
        tmpstr = str(EE-x_1sideL)
        for cur_r in range(int(round(spherical_lim)+1)):
            tmpstr += " " + str(rDist_per_dE_A[EE][cur_r])
        for cur_r in range(int(round(spherical_lim)+1)):
            tmpstr += " " + str(rDist_per_dE_B[EE][cur_r])
        tmpstr += " \n"
        R_out.write(tmpstr)
    R_out.close()

    ### WRITE TOTAL THE HIST TO FILES ###
    tmpstr = "dE_AB_T"+ file_name_tail +".Hist"
    AB_out = open(tmpstr,"w")
    AB_out.write("# A,B_means::: "+str(mean_dE_A)+" "+str(mean_dE_B)+"\n")
    AB_out.write("# A,B_variances::: "+str(var_dE_A)+" "+str(var_dE_B)+"\n")
    AB_out.write("# dE stateA_prob stateA_ct stateB_prob stateB_ct \n")
    for i in range(2*x_1sideL+1):
            tmpstr = str(i-x_1sideL)+" "+str(A_hist[i][0]/np.sum(A_hist[:]))+" "+str(A_hist[i][0])+" "+str(B_hist[i][0]/np.sum(B_hist[:]))+" "+str(B_hist[i][0])+"\n"
            AB_out.write(tmpstr)
    AB_out.close()

    ### WRITE R_DIST_MEAN_VAR TO FILE ###
    tmpstr = "rDists_mean_var"+ file_name_tail + ".Hist"
    R_mv_out = open(tmpstr,"w")
    tmpstr = "# moments"
    for ii_R in range(int(round(spherical_lim)+1)):
        tmpstr += " A_"+str(ii_R)
    for ii_R in range(int(round(spherical_lim)+1)):
        tmpstr += " B_"+str(ii_R)
    tmpstr += " \n"
    R_mv_out.write(tmpstr)
    for Rmoment in range(0,nEmoment+2):
        tmpstr = str(Rmoment)
        for cur_r in range(int(round(spherical_lim)+1)):
            tmpstr += " " + str(rDist_mean_var_dE_A[cur_r][Rmoment])
        for cur_r in range(int(round(spherical_lim)+1)):
            tmpstr += " " + str(rDist_mean_var_dE_B[cur_r][Rmoment])
        tmpstr += " \n"
        R_mv_out.write(tmpstr)
    R_mv_out.close()
