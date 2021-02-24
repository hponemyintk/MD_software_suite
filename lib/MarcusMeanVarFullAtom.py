import numpy as np
import sys, scipy.special
np.set_printoptions(threshold=np.inf,suppress=True)				#force print entire array
np.seterr('raise')


try:
    from lib.read_coord import get_Coul_coord, r_sqrd
    from lib.math_lib import ID2MolData, nearest_neigh_rsqrd
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def calcInterRedPE(xyz_i,xyz_j,dq_i,dq_j,boxL):
    sum = 0
    for ii in range(13):
        cur_xyz_i = xyz_i[ii]
        cur_dq_i = dq_i[ii]
        cur_dis = np.sqrt(nearest_neigh_rsqrd( (xyz_j-cur_xyz_i) ,boxL))
        sum += np.sum( (cur_dq_i * dq_j)/cur_dis )
        # if ii ==0:
            # print("ii,sum,(cur_dq_i * dq_j)/cur_dis,xyz_i,xyz_j,dq_i,dq_j,np.sum((cur_dq_i * dq_j)/cur_dis)\n",ii,sum,(cur_dq_i * dq_j)/cur_dis,xyz_i,xyz_j,dq_i,dq_j,np.sum((cur_dq_i * dq_j)/cur_dis),boxL)
        # print(np.sum( (cur_dq_i * dq_j)/cur_dis ))
    return sum


def DoMeanVar(lammpstrj,CoulDump,Tstep,TskipI,TskipB,x_1sideL,nEmoment):
    nEmoment = int(nEmoment)

    # to check dE dist and inerRedox term dist
    distCheck = np.zeros((2*x_1sideL+1,8))

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
    file_name_tail = str(Tstep)+"_TskipB" +str(TskipB)+"_upTonthMoment_MDvFullAtom"

    ### Define Ferri/Ferro Charges ###
    FerriC = np.array([-0.834,0.402,0.402,0.402,0.402,0.402,0.402,-0.763,-0.763,-0.763,-0.763,-0.763,-0.763])
    FerroC = np.array([-1.15,0.431,0.431,0.431,0.431,0.431,0.431,-0.906,-0.906,-0.906,-0.906,-0.906,-0.906])
    Ci2o = FerroC-FerriC
    Co2i = FerriC-FerroC
    print(FerriC,"\n",FerroC,"\n",Ci2o,"\n",Co2i)


    ### calculate changes in intra molecular energy ###
    FerriK = np.array([77.61701123,1286.327551])            # 0 is bond coeff for Fe-C, 1 is for C-N
    FerroK = np.array([41.31587596,1189.598196])            # 0 is bond coeff for Fe-C, 1 is for C-N
    FerriL = np.array([1.997,1.176])                        # 0 is equi len for Fe-C, 1 is for C-N
    FerroL = np.array([2.035,1.186])                        # 0 is equi len for Fe-C, 1 is for C-N
    mag_dE_intra = 6*( FerriK[0]*(FerroL[0]-FerriL[0])**2 + FerriK[1]*(FerroL[1]-FerriL[1])**2 + FerroK[0]*(FerriL[0]-FerroL[0])**2 + FerroK[1]*(FerriL[1]-FerroL[1])**2 )
    print("mag_dE_intra",mag_dE_intra)


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

    ### get 1-5 offset term ###
    tmp_xyz_ferri = cur_xyz[ (cur_xyz[:,0]==itype) | (cur_xyz[:,0]==itype+1) | (cur_xyz[:,0]==itype+2) ]
    tmp_xyz_ferro = cur_xyz[ (cur_xyz[:,0]==jtype) | (cur_xyz[:,0]==jtype+1) | (cur_xyz[:,0]==jtype+2) ]
    txyz_ferri = ID2MolData(0,tmp_xyz_ferri,13)[:,1:4]
    txyz_ferro = ID2MolData(0,tmp_xyz_ferro,13)[:,1:4]
    # get 1-5 for Ferri
    tInd = 7
    intra_15 = 0
    tC = FerriC
    dq = Ci2o
    for ii_i in range(6):
        for jj_i in range(6):
            if ii_i != jj_i:
                dis = np.sqrt(r_sqrd(txyz_ferri[ii_i+tInd],txyz_ferri[jj_i+tInd], boxL))
                intra_15 += ( -dq[ii_i+tInd]*dq[jj_i+tInd] + dq[ii_i+tInd]*tC[jj_i+tInd] - tC[ii_i+tInd]*dq[jj_i+tInd] )/dis
    # get 1-5 for ferro
    tC = FerroC
    dq = Co2i
    for ii_i in range(6):
        for jj_i in range(6):
            if ii_i != jj_i:
                dis = np.sqrt(r_sqrd(txyz_ferro[ii_i+tInd],txyz_ferro[jj_i+tInd], boxL))
                intra_15 += ( -dq[ii_i+tInd]*dq[jj_i+tInd] + dq[ii_i+tInd]*tC[jj_i+tInd] - tC[ii_i+tInd]*dq[jj_i+tInd] )/dis
    intra_15 = (intra_15 * E_conv_const)/2
    print("intra_15 is :::",intra_15)

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
            Txyz, TCoul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile,posfile)
        else:
            Txyz, TCoul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile,posfile)

            Txyz_ferri = Txyz[ (Txyz[:,0]==itype) | (Txyz[:,0]==itype+1) | (Txyz[:,0]==itype+2) ]
            Txyz_ferro = Txyz[ (Txyz[:,0]==jtype) | (Txyz[:,0]==jtype+1) | (Txyz[:,0]==jtype+2) ]
            TCoul_ferri = TCoul[ (TCoul[:,0]==itype) | (TCoul[:,0]==itype+1) | (TCoul[:,0]==itype+2) ]
            TCoul_ferro = TCoul[ (TCoul[:,0]==jtype) | (TCoul[:,0]==jtype+1) | (TCoul[:,0]==jtype+2) ]

            ### test ID2mol func ###
            # indTmp = 11
            # print("Ferri,Ferro,FerriC,FerroC\n",ID2MolData(indTmp,Txyz_ferri,13),"\n",ID2MolData(indTmp,Txyz_ferro,13),"\n",ID2MolData(indTmp,TCoul_ferri,13),"\n",ID2MolData(indTmp,TCoul_ferro,13),boxL)

            for i_AB in range(0,ict):
                for j_AB in range(0,jct):
                    ### recalc dE gaps again here ###
                    xyz_ferri = ID2MolData(i_AB,Txyz_ferri,13)[:,1:4]
                    xyz_ferro = ID2MolData(j_AB,Txyz_ferro,13)[:,1:4]
                    Madl_ferri = ID2MolData(i_AB,TCoul_ferri,13)[:,1]/FerriC*2
                    Madl_ferro = ID2MolData(j_AB,TCoul_ferro,13)[:,1]/FerroC*2
                    # print("Madl_ferri,Madl_ferro\n",Madl_ferri,Madl_ferro)

                    ### calculate each term in dE expression ###
                    FerriAmbient = np.sum(Madl_ferri * Ci2o)
                    FerroAmbient = np.sum(Madl_ferro * Co2i)
                    interRed_A = calcInterRedPE(xyz_ferro,xyz_ferri,Co2i,Ci2o,boxL) * E_conv_const
                    interRed_B = calcInterRedPE(xyz_ferri,xyz_ferro,Ci2o,Co2i,boxL) * E_conv_const
                    cur_dE_A = -FerriAmbient -FerroAmbient -interRed_A - mag_dE_intra + intra_15
                    cur_dE_B = +FerriAmbient +FerroAmbient +interRed_B + mag_dE_intra - intra_15
                    # if i_AB == indTmp and j_AB == indTmp:
                    #     print("interRed_A,interRed_B,cur_dE_A,cur_dE_B::: ",interRed_A,interRed_B,cur_dE_A,cur_dE_B)

                    dis = np.sqrt(r_sqrd(xyz_ferri[0],xyz_ferro[0], boxL))      # distance for binning per shell data

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

                        rDist_mean_var_oxd[int(round(dis))][0] += np.sum(Madl_ferri)
                        rDist_mean_var_red[int(round(dis))][0] += np.sum(Madl_ferro)
                        dE_ct += 1

                        distCheck[int(round(cur_dE_A+x_1sideL)),0] += 1
                        distCheck[int(round(interRed_A+x_1sideL)),1] += 1
                        distCheck[int(round(Madl_ferri[0]+x_1sideL)),2] += 1
                        for tmpii in range(6):
                            distCheck[int(round(Madl_ferri[tmpii+1]+x_1sideL)),3] += 1
                        for tmpii in range(6):
                            distCheck[int(round(Madl_ferri[tmpii+7]+x_1sideL)),4] += 1
                        distCheck[int(round(Madl_ferro[0]+x_1sideL)),5] += 1
                        for tmpii in range(6):
                            distCheck[int(round(Madl_ferro[tmpii+1]+x_1sideL)),6] += 1
                        for tmpii in range(6):
                            distCheck[int(round(Madl_ferro[tmpii+7]+x_1sideL)),7] += 1

            if (TskipB!=0) :
                for i in range(0,TskipB):
                    cur_xyz, cur_Coul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile,posfile)
                cur_step = cur_step + TskipB
            cur_step += 1
            # print "cur_step is", cur_step
    coulfile.close()
    posfile.close()
    np.savetxt("distCheck.out",distCheck)



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
            Txyz, TCoul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile1,posfile1)
        else:
            Txyz, TCoul, timestamp, Ntotal, boxL = get_Coul_coord(coulfile1,posfile1)

            Txyz_ferri = Txyz[ (Txyz[:,0]==itype) | (Txyz[:,0]==itype+1) | (Txyz[:,0]==itype+2) ]
            Txyz_ferro = Txyz[ (Txyz[:,0]==jtype) | (Txyz[:,0]==jtype+1) | (Txyz[:,0]==jtype+2) ]
            TCoul_ferri = TCoul[ (TCoul[:,0]==itype) | (TCoul[:,0]==itype+1) | (TCoul[:,0]==itype+2) ]
            TCoul_ferro = TCoul[ (TCoul[:,0]==jtype) | (TCoul[:,0]==jtype+1) | (TCoul[:,0]==jtype+2) ]

            for i_AB in range(0,ict):
                for j_AB in range(0,jct):
                    ### recalc dE gaps again here ###
                    xyz_ferri = ID2MolData(i_AB,Txyz_ferri,13)[:,1:4]
                    xyz_ferro = ID2MolData(j_AB,Txyz_ferro,13)[:,1:4]
                    Madl_ferri = ID2MolData(i_AB,TCoul_ferri,13)[:,1]/FerriC*2
                    Madl_ferro = ID2MolData(j_AB,TCoul_ferro,13)[:,1]/FerroC*2

                    ### calculate each term in dE expression ###
                    FerriAmbient = np.sum(Madl_ferri * Ci2o)
                    FerroAmbient = np.sum(Madl_ferro * Co2i)
                    interRed_A = calcInterRedPE(xyz_ferro,xyz_ferri,Co2i,Ci2o,boxL) * E_conv_const
                    interRed_B = calcInterRedPE(xyz_ferri,xyz_ferro,Ci2o,Co2i,boxL) * E_conv_const
                    cur_dE_A = -FerriAmbient -FerroAmbient -interRed_A - mag_dE_intra + intra_15
                    cur_dE_B = +FerriAmbient +FerroAmbient +interRed_B + mag_dE_intra - intra_15

                    dis = np.sqrt(r_sqrd(xyz_ferri[0],xyz_ferro[0], boxL))      # distance for binning per shell data

                    if dis < spherical_lim:          # only consider the dE half a box length away from the central atom
                        var_dE_A += (cur_dE_A - mean_dE_A)**2
                        var_dE_B += (cur_dE_B - mean_dE_B)**2

                        for i_moment in range(1,nEmoment+1):
                            # do for state A
                            rDist_mean_var_dE_A[int(round(dis))][i_moment] += (cur_dE_A - rDist_mean_var_dE_A[int(round(dis))][0])**(i_moment)
                            # do for state B
                            rDist_mean_var_dE_B[int(round(dis))][i_moment] += (cur_dE_B - rDist_mean_var_dE_B[int(round(dis))][0])**(i_moment)

                        # do correlation
                        rDist_mean_var_oxd[int(round(dis))][1] += (np.sum(Madl_ferri) - rDist_mean_var_oxd[int(round(dis))][0])**2
                        rDist_mean_var_red[int(round(dis))][1] += (np.sum(Madl_ferro) - rDist_mean_var_red[int(round(dis))][0])**2
                        rDist_mean_var_dE_A[int(round(dis))][nEmoment+1] += (np.sum(Madl_ferri) - rDist_mean_var_oxd[int(round(dis))][0]) * (np.sum(Madl_ferro) - rDist_mean_var_red[int(round(dis))][0])

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
