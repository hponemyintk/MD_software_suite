import numpy as np
import pandas as pd
from lib.math_lib import unwrapBond, nearest_neigh_rsqrd
from lib.read_coord import get_wrapped_coord, get_Coul, r_sqrd, manual_parsing

def getIDs(lammpstrj,CoulDump,Tstep,stype,itype,jtype,uv,K_lim):
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
    file_name_tail = "T"+str(Tstep)+"_uv" +np.array2string(uv,separator='_')[1:-1].replace(" ","")+"K_lim"+str(K_lim)+".out"
    # print(uv,file_name_tail)

    # check if a snapshot has cations in cut-off range from Ferri/Ferro
    # if so, check dE in given range or not
    # separate configs w or w/o briding cation and output them in 2 diff files
    xyzFile = open(lammpstrj,"r")
    CoulFile = open(CoulDump,"r")

    # out files
    wcat_out = open("wcat_"+file_name_tail,"w")
    wcat_out.write( "# timestamp catID FerriID FerroID dis cur_dE_A cur_dE_B\n")
    wOcat_out = open("wocat_"+file_name_tail,"w")
    wOcat_out.write( "# timestamp catID FerriID FerroID dis cur_dE_A cur_dE_B\n")

    for cur_t in range(Tstep):
        cur_Coul = get_Coul(CoulFile)
        cur_xyz, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(xyzFile)
        # print(curCoul,cur_xyz,tcur_Tstep,Ntotal,boxL)

        ### get xyz and coul data for cat, ferri and ferro
        s_xyz = cur_xyz[ (cur_xyz[:,0]==stype) ]
        i_xyz = cur_xyz[ (cur_xyz[:,0]==itype) ]
        j_xyz = cur_xyz[ (cur_xyz[:,0]==jtype) ]
        i_coul = cur_Coul[ (cur_Coul[:,0]==itype) ]
        j_coul = cur_Coul[ (cur_Coul[:,0]==jtype) ]


        # # find the closest cation around Ferri and Ferro
        # is_dis = np.sqrt(nearest_neigh_rsqrd(s_xyz[:,1:4]-i_xyz[0,1:4],boxL))
        # js_dis = np.sqrt(nearest_neigh_rsqrd(s_xyz[:,1:4]-j_xyz[0,1:4],boxL))
        # i_ind, j_ind = np.argmin(is_dis),np.argmin(js_dis)
        # if is_dis[i_ind]< K_lim and js_dis[j_ind]<K_lim:
        #     haveBridgedCation = True

        ### this choose the cation bridged state with smallest Ferri-cat, Ferro-cat distance
        haveBridgedCation = False                   # check if a given snapshot contains bridging cations
        s_ind, i_ind, j_ind, shortestTdis = np.nan, np.nan, np.nan, np.inf
        for cur_s in range(len(s_xyz)):
            # print(s_xyz,i_xyz,j_xyz)
            is_dis = np.sqrt(nearest_neigh_rsqrd(i_xyz[:,1:4]-s_xyz[cur_s][1:4],boxL))
            js_dis = np.sqrt(nearest_neigh_rsqrd(j_xyz[:,1:4]-s_xyz[cur_s][1:4],boxL))
            i_indTMP, j_indTMP = np.argmin(is_dis),np.argmin(js_dis)
            if is_dis[i_indTMP]< K_lim and js_dis[j_indTMP]<K_lim and (is_dis[i_indTMP] + js_dis[j_indTMP])<shortestTdis:
                s_ind, i_ind, j_ind, shortestTdis = cur_s, i_indTMP, j_indTMP, (is_dis[i_indTMP] + js_dis[j_indTMP])
                haveBridgedCation = True
                # print(s_ind, i_ind, j_ind, shortestTdis)
                # break

        i_ind, j_ind = 0, 0                    # this is hard-coded for 1Ferri-1Ferro system
        dis = np.sqrt(r_sqrd(i_xyz[i_ind][1:4],j_xyz[j_ind][1:4], boxL))
        ferri_ML = i_coul[i_ind][1]/icharge*2.
        ferro_ML = j_coul[j_ind][1]/jcharge*2.

        # calculate and store dE_A and dE_B
        cur_dE_A = dq *(E_conv_const * dq/dis + ferri_ML - ferro_ML)
        cur_dE_B = dq *(-E_conv_const * dq/dis + ferro_ML - ferri_ML)
        for ii in range(int(len(uv)/2)):                 # using side A (state A) dE to do constraint here
            # print(np.abs(uv[ii*2])< np.abs(cur_dE_A) < np.abs(uv[ii*2+1]),(np.abs(uv[ii*2]) > np.abs(cur_dE_A) > np.abs(uv[ii*2+1])),haveBridgedCation)
            if ( (np.abs(uv[ii*2])< np.abs(cur_dE_A) < np.abs(uv[ii*2+1])) or (np.abs(uv[ii*2]) > np.abs(cur_dE_A) > np.abs(uv[ii*2+1])) ) and haveBridgedCation:
                # print("found briding cations:::",cur_t,s_ind,i_ind,j_ind,dis,cur_dE_A,cur_dE_B)
                wcat_out.write(str(cur_t)+" "+str(s_ind)+" "+str(i_ind)+" "+str(j_ind)+" "+str(dis)+" "+str(cur_dE_A)+" "+str(cur_dE_B)+"\n")
            if ( (np.abs(uv[ii*2])< np.abs(cur_dE_A) < np.abs(uv[ii*2+1])) or (np.abs(uv[ii*2]) > np.abs(cur_dE_A) > np.abs(uv[ii*2+1])) ) and not haveBridgedCation:
                # print("NoBridge:",cur_t,s_ind,i_ind,j_ind,dis,cur_dE_A,cur_dE_B)
                wOcat_out.write(str(cur_t)+" "+str(s_ind)+" "+str(i_ind)+" "+str(j_ind)+" "+str(dis)+" "+str(cur_dE_A)+" "+str(cur_dE_B)+"\n")
    wcat_out.close()
    wOcat_out.close()

