import numpy as np
import sys, scipy.special
np.set_printoptions(threshold=np.inf,suppress=True)             #force print entire array
np.seterr('raise')


try:
    from lib.read_coord import get_Coul_coord, r_sqrd
    from lib.math_lib import ID2MolData, nearest_neigh_rsqrd
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")


# use recurrance relation until we get big enough xx as argument in psi for approximation
# double checked to be working with analytical values
def psi(xx):
    if xx < 8.5:
        return psi(xx+1) - 1/(xx)
    else:
        #using natural log here (compared wiki and wolfram def of digamma function, and it seems it should be natural log.)
        return np.log(xx) - 1/(2*xx) - 1/(12*xx**2) + 1/(120*xx**4) - 1/(252*xx**6) + 1/(240*xx**8) - 5/(660*xx**10) + 691/(32760*xx**12) - 1/(12*xx**14)


# Need to put in adjusted Lz (the distance between two electrode for this)
def resum_f(dd1,dd2,gamma, Lz):
    fArray = np.zeros(13)
    for ii in range(13):
        d1 = dd1[ii]
        d2 = dd2[ii]
        fArray[ii] = 1/(2*Lz) * ( -1/d1 -1/d2 +2*gamma +psi(1+d1) +psi(1+d2) )
        # print("d1,d2,1/(2*Lz),-1/d1 -1/d2,+2*gamma,psi(1+d1),psi(1+d2):::",d1,d2,1/(2*Lz),-1/d1 -1/d2,+2*gamma,psi(1+d1),psi(1+d2))
    return fArray


### make a function that returns dE_A or dE_B along with adjusted z position, given madelung potential and coord of either oxidized or reduced, and the state (A/B) of the system
def calc_dE_x(state,redInd,redC,redTxyz,redTCoul,dCs,mag_dE_intra,intra_15,Lz,E_conv_const,gamma):        # state = A or B
    xyz = ID2MolData(redInd,redTxyz,13)[:,1:4]
    Madl = ID2MolData(redInd,redTCoul,13)[:,1]/redC*2
    redAmbient = np.sum(Madl * dCs)            # Ambient means not Image or Replica of alpha
    d1 = (xyz[:,2]+Lz/2)/Lz
    d2 = 1-d1
    # print("d1,d2,xyz[:,2],int(abs(round(xyz[2]/2)))",d1,d2,xyz[:,2])
    resum_term = np.sum( (dCs**2 * ( resum_f(d1,d2,gamma,Lz) ) )/2 ) * E_conv_const
    if state == "A":
        return int(abs(round(xyz[0,2]))), redAmbient + resum_term + mag_dE_intra - intra_15 * dCs[-1]        # red => oxd
    if state == "B":
        return int(abs(round(xyz[0,2]))), redAmbient - resum_term - mag_dE_intra - intra_15 * dCs[-1]       # odx => red


def doHeterodE(lammpstrj,CoulDump,Tstep,TskipI,TskipB,x_1sideL,nEmoment):
    nEmoment = int(nEmoment)

    ### CONSTANTS ###
    E_conv_const = 332.06371            # copied from LAMMPS update.cpp [convert qq/r to kcal/mol]
    gamma = 0.577215664901532860606512  # Euler–Mascheroni constant

    ### Changable Vars ###
    Lz = 80                                 # hardcoded to be 80 for Ferri/Ferro electrode system ***change this if z distance between 2 electrode changes!!!

    ### INPUT PARAMS ###
    itype = 5
    jtype = 8                               ### this code assume j has higher index than i in the lammps output file
    coulfile = open(CoulDump, "r")
    posfile = open(lammpstrj,"r")
    file_name_tail = str(Tstep)+"_TskipB" +str(TskipB)+"_upTonthMoment_ElectrodeReorg"

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

    ### these two might be flipped!!!!***
    mag_dE_intraB = 6*( FerroK[0]*(FerriL[0]-FerroL[0])**2 + FerroK[1]*(FerriL[1]-FerroL[1])**2 )           # [kcal/mol]
    mag_dE_intraA = 6*( FerriK[0]*(FerroL[0]-FerriL[0])**2 + FerriK[1]*(FerroL[1]-FerriL[1])**2 )           # [kcal/mol]
    print("mag_dE_intraI, mag_dE_intraO",mag_dE_intraA, mag_dE_intraB)


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
    # grab one ferri and one ferro for calculating 1-5 terms
    txyz_ferri = ID2MolData(0,tmp_xyz_ferri,13)[:,1:4]
    txyz_ferro = ID2MolData(0,tmp_xyz_ferro,13)[:,1:4]
    print("txyz_ferro:", txyz_ferro)
    # get 1-5 for Ferro
    tInd = 7
    tC = FerroC
    intra_15A = 0
    for ii_i in range(6):
        for jj_i in range(6):
            if ii_i != jj_i:
                dis = np.sqrt(r_sqrd(txyz_ferro[ii_i+tInd],txyz_ferro[jj_i+tInd], boxL))
                intra_15A += tC[ii_i+tInd]*tC[jj_i+tInd]/dis
    intra_15A = (intra_15A * E_conv_const)          #/2 # half-factor removed as following equation 20 in EnergyGapwFullAtoms.pdf
    # get 1-5 for ferri
    tC = FerriC
    intra_15B = 0
    for ii_i in range(6):
        for jj_i in range(6):
            if ii_i != jj_i:
                dis = np.sqrt(r_sqrd(txyz_ferri[ii_i+tInd],txyz_ferri[jj_i+tInd], boxL))
                intra_15B += tC[ii_i+tInd]*tC[jj_i+tInd]/dis
    intra_15B = (intra_15B * E_conv_const)          #/2 # half-factor removed as following equation 20 in EnergyGapwFullAtoms.pdf
    print("intra_15A, intra_15B :::",intra_15A, intra_15B)

    z_lim = int(round(boxL[2]/2))   # 1e6 # boxL[0]/2

    ### For r distribution for each dE ###
    print("z_lim is:::",z_lim)
    rDist_per_dE_A = np.zeros((2*x_1sideL+1,z_lim+1))           # r from 0 to z_lim (as dis always less than spehical_lim)
    rDist_per_dE_B = np.zeros((2*x_1sideL+1,z_lim+1))
    rDist_mean_var_dE_A = np.zeros((z_lim+1,nEmoment+2))           # index 0 is mean 1..n is 2..2n moment
    rDist_mean_var_dE_B = np.zeros((z_lim+1,nEmoment+2))


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

            ### for checking with analytical expression in 20210120_DigammaApproximation.xlsx ###
            # print("printing psi:::",[psi(1/ii) for ii in [1,2,3,4,6,8]])
            # print("printing psi:::",[psi(ii) for ii in [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]])
            # print("printing psi:::",[psi(ii) for ii in [2, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1]])
            ### for checking resummation function, f ###
            # print("printing resum_f:::",[ (np.sum( (Ci2o**2 * ( resum_f((ii*10)/Lz,(Lz-ii*10)/Lz,gamma,Lz) ) )/2 ) * E_conv_const) for ii in [1,2,3,4,6]])
            # tmpArray = np.ones((13))
            # tmpArray = tmpArray-.05
            # for ii in range(1,10,1):
            #     tmpArray[ii-1] = ii/10
            # print("tmpArray is:::",tmpArray)
            # print("printing resum_f:::", resum_f(tmpArray,1-tmpArray,gamma,80))

            for AA in range(0,jct):
                zInd_A, cur_dE_A = calc_dE_x("A",AA,FerroC,Txyz_ferro,TCoul_ferro,Co2i,mag_dE_intraA,intra_15A,Lz,E_conv_const,gamma)

                mean_dE_A += cur_dE_A
                # do for state A
                dE_A_ind = int(round(cur_dE_A+x_1sideL))
                A_hist[dE_A_ind] += 1
                rDist_per_dE_A[dE_A_ind][zInd_A] += 1                  # for Histogramming rDist_per_dE_A
                rDist_mean_var_dE_A[zInd_A][0] += cur_dE_A

                # ### check xyz and Coul in ###
                # if AA==0 and cur_step==8:
                #     print("AA Txyz_ferro,TCoul_ferro",Txyz_ferro,TCoul_ferro)

            for BB in range(0,ict):
                zInd_B, cur_dE_B = calc_dE_x("B",BB,FerriC,Txyz_ferri,TCoul_ferri,Ci2o,mag_dE_intraB,intra_15B,Lz,E_conv_const,gamma)

                mean_dE_B += cur_dE_B
                # do for state B
                dE_B_ind = int(round(cur_dE_B+x_1sideL))
                B_hist[dE_B_ind] += 1
                rDist_per_dE_B[dE_B_ind][zInd_B] += 1                  # for Histogramming rDist_per_dE_B
                rDist_mean_var_dE_B[zInd_B][0] += cur_dE_B

                dE_ct += 1

                # ### check xyz and Coul in ###
                # if BB==0 and cur_step==8:
                #     print("BB Txyz_ferri,TCoul_ferri",Txyz_ferri,TCoul_ferri)

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
    for i_mean in range(int(round(z_lim))+1):                           # do avg for each 
        if sum(rDist_per_dE_A[:,i_mean]) != 0 and sum(rDist_per_dE_B[:,i_mean]) !=0:
            rDist_mean_var_dE_A[i_mean][0] /= sum(rDist_per_dE_A[:,i_mean])
            rDist_mean_var_dE_B[i_mean][0] /= sum(rDist_per_dE_B[:,i_mean])

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

            for AA in range(0,jct):
                zInd_A, cur_dE_A = calc_dE_x("A",AA,FerroC,Txyz_ferro,TCoul_ferro,Co2i,mag_dE_intraA,intra_15A,Lz,E_conv_const,gamma)
                var_dE_A += (cur_dE_A - mean_dE_A)**2

                for i_moment in range(1,nEmoment+1):
                    # do for state A
                    rDist_mean_var_dE_A[zInd_A][i_moment] += (cur_dE_A - rDist_mean_var_dE_A[zInd_A][0])**(i_moment)

            for BB in range(0,ict):
                zInd_B, cur_dE_B = calc_dE_x("B",BB,FerriC,Txyz_ferri,TCoul_ferri,Ci2o,mag_dE_intraB,intra_15B,Lz,E_conv_const,gamma)
                var_dE_B += (cur_dE_B - mean_dE_B)**2

                for i_moment in range(1,nEmoment+1):
                    # do for state B
                    rDist_mean_var_dE_B[zInd_B][i_moment] += (cur_dE_B - rDist_mean_var_dE_B[zInd_B][0])**(i_moment)

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
    for i_var in range(int(round(z_lim))+1):
        if sum(rDist_per_dE_A[:,i_var]) != 0 and sum(rDist_per_dE_B[:,i_var]) !=0:
            rDist_mean_var_dE_A[i_var][1:nEmoment+2]/=sum(rDist_per_dE_A[:,i_var])
            rDist_mean_var_dE_B[i_var][1:nEmoment+2]/=sum(rDist_per_dE_B[:,i_var])


    ### PRINT R_DIST PER dEs ###
    tmpstr = "rDists_per_dE"+ file_name_tail + ".Hist"
    R_out = open(tmpstr,"w")
    tmpstr = "# dEs"
    for ii_R in range(int(round(z_lim)+1)):
        tmpstr += " A_"+str(ii_R)
    for ii_R in range(int(round(z_lim)+1)):
        tmpstr += " B_"+str(ii_R)
    tmpstr += " \n"
    R_out.write(tmpstr)
    for EE in range(0,2*x_1sideL+1):
        tmpstr = str(EE-x_1sideL)
        for cur_r in range(int(round(z_lim)+1)):
            tmpstr += " " + str(rDist_per_dE_A[EE][cur_r])
        for cur_r in range(int(round(z_lim)+1)):
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
    for ii_R in range(int(round(z_lim)+1)):
        tmpstr += " A_"+str(ii_R)
    for ii_R in range(int(round(z_lim)+1)):
        tmpstr += " B_"+str(ii_R)
    tmpstr += " \n"
    R_mv_out.write(tmpstr)
    for Rmoment in range(0,nEmoment+2):
        tmpstr = str(Rmoment)
        for cur_r in range(int(round(z_lim)+1)):
            tmpstr += " " + str(rDist_mean_var_dE_A[cur_r][Rmoment])
        for cur_r in range(int(round(z_lim)+1)):
            tmpstr += " " + str(rDist_mean_var_dE_B[cur_r][Rmoment])
        tmpstr += " \n"
        R_mv_out.write(tmpstr)
    R_mv_out.close()
