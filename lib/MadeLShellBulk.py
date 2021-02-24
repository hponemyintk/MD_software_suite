import numpy as np
import sys
import time

np.set_printoptions(threshold=sys.maxsize)
try:
    from lib.read_coord import get_TqxyzCoul, r_sqrd
    from lib.math_lib import nearest_neigh_rsqrd
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def find_saroundi(cur_ipos,s_xyz,boxL,Avg_saroundi,colind,Kcut):
    is_dist = s_xyz - cur_ipos
    is_dist = nearest_neigh_rsqrd(is_dist,boxL)
    Avg_saroundi[0,colind] += len(is_dist[is_dist<=Kcut**2])
    Avg_saroundi[1,colind] += 1
    Avg_saroundi[2,colind] += np.sum(is_dist[is_dist<=Kcut**2])
    return Avg_saroundi

def doMadeLHist(lammpstrj,Tstep,itype,jtype,stype,x_1sideL,Rcut,Kcut):
    MadeLHist = np.zeros((int(2*x_1sideL)+1,2))         # col 0 for shell, col 1 for bulk
    Avg_saroundi = np.zeros((3,3))                      # col 0 for shell, col 1 for bulk, row 0 for sum, row 1 for ct, row 2 for distance

    posfile = open(lammpstrj,"r")
    # do histogramming here
    curT = 0

    # keep track of time
    t0 = time.time()
    cur_prcnt = 0

    while curT < Tstep:
        TqxyzCoul, tcur_Tstep, Ntotal, boxL = get_TqxyzCoul(posfile)
        i_qxyzC = TqxyzCoul[TqxyzCoul[:,0]==itype]
        j_qxyzC = TqxyzCoul[TqxyzCoul[:,0]==jtype]
        s_qxyzC = TqxyzCoul[TqxyzCoul[:,0]==stype]

        # convert CoulPE to MadeLung Potential
        i_qxyzC[:,5] = i_qxyzC[:,5]/(i_qxyzC[:,1]/2)

        # keep only xyzCoul columns
        i_xyz = i_qxyzC[:,2:5]
        j_xyz = j_qxyzC[:,2:5]
        s_xyz = s_qxyzC[:,2:5]
        i_Coul = i_qxyzC[:,5]

        # print(i_qxyzC[-3:],"\n\n",j_qxyzC[-3:])

        ij_thred = 0
        # chek if s & j are of the same type
        if itype == jtype:
            ij_thred = 1

        # find distance and histgram the subpopulation
        for ii_i in range(len(i_qxyzC)):
            madel_ind = int(round(i_Coul[ii_i]))+x_1sideL
            ij_dist = j_xyz - i_xyz[ii_i]
            ij_dist = nearest_neigh_rsqrd(ij_dist,boxL)

            if len(ij_dist[ij_dist<=Rcut**2]) > ij_thred:
                MadeLHist[madel_ind,0] += 1
                Avg_saroundi = find_saroundi(i_xyz[ii_i],s_xyz,boxL,Avg_saroundi,0,Kcut)
            else:
                MadeLHist[madel_ind,1] += 1
                Avg_saroundi = find_saroundi(i_xyz[ii_i],s_xyz,boxL,Avg_saroundi,1,Kcut)
        curT += 1

        # print the progress of the code
        TskipI = 0
        if curT*100/(Tstep-TskipI) >= cur_prcnt:
            if cur_prcnt == 0:
                print("*** Outputting time data. ***\nTstep,TskipI,lammpstrj,itype,jtype,x_1sideL,Rcut",
                      Tstep,TskipI,lammpstrj,itype,jtype,x_1sideL,Rcut,
                      "\n#######################################\nPercent[%] Time[s]")
            print("{0:3.2f}       {1:3.2f}".format(curT*100/(Tstep-TskipI), time.time()-t0))
            cur_prcnt += 10
            t0 = time.time()

    tmpstr = "MadeLShellBulk_ijs"+str(itype)+str(jtype)+str(stype)+"T"+str(Tstep)+"x_1sideL"+str(x_1sideL)+"R_K_lim"+str(Rcut)+"_"+str(Kcut)+"sjthred"+str(ij_thred)+".out"
    with open(tmpstr,'w') as outfile:
        tmpstr1 = ("# MadelPotential[kcal/mol/e] P(shell) P(bulk). ShellCt "+str(np.sum(MadeLHist[:,0]))+" BulkCt "+str(np.sum(MadeLHist[:,1]))+"\n")
        outfile.write(tmpstr1)
        tmpstr1 = ("# avg s around i in i's Shell "+str(Avg_saroundi[0,0]/Avg_saroundi[1,0])+" Dis "+str(Avg_saroundi[2,0]/Avg_saroundi[1,0])+" Bulk "+str(Avg_saroundi[0,1]/Avg_saroundi[1,1])+" Dis "+str(Avg_saroundi[2,1]/Avg_saroundi[1,1])+"\n")
        outfile.write(tmpstr1)
        # do normalization
        MadeLHist[:,0] /= np.sum(MadeLHist[:,0])
        MadeLHist[:,1] /= np.sum(MadeLHist[:,1])

    # print the data
    with open(tmpstr,'a') as outfile:
        for ii in range(int(2*x_1sideL+1)):
            tmpstr = str(ii-x_1sideL)+" "+str(MadeLHist[ii,0])+" "+str(MadeLHist[ii,1])+"\n"
            outfile.write(tmpstr)


