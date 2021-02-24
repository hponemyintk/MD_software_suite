# awk 'BEGIN{min=100}{if($1==1 && $5<min){min=$5}}END{print min}' tests/equ_nve.lammpstrj
# awk 'BEGIN{max=0}{if($1==1 && $5>max){max=$5}}END{print max}' tests/equ_nve.lammpstrj

import numpy as np
import time
import sys
np.set_printoptions(threshold=sys.maxsize)

try:
    from lib.read_coord import get_wrapped_coord_Tqxyz
    from lib.math_lib import adjustL_conv2UV, nearest_neigh_rsqrd
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def histZ(Z_array,Hist_array,col_ind,boxL,scaling):
    for ii_ion in range(len(Z_array)):
        z_ind = int(round((Z_array[ii_ion]+boxL[2]/2)*scaling))
        Hist_array[z_ind,col_ind] += 1
    return Hist_array

### return a boolean list for shell atoms
def flag_shell_water(Txyz_pos,itype,jtype,z_lim,boxL):
    cat_xyz = Txyz_pos[Txyz_pos[:,0]==itype+2][:,1:4]                # *** make sure you are gonna take the right type ---confirmed 3 is the pos ions
    ani_xyz = Txyz_pos[Txyz_pos[:,0]==jtype+2][:,1:4]                # *** make sure you are gonna take the right type ---confirmed 4 is the neg ions
    Oxy_xyz = Txyz_pos[Txyz_pos[:,0]==itype][:,1:4]
    all_atom_flagged = np.zeros(len(Txyz_pos),dtype=bool)
    # print(cat_xyz[0:3],"\n",Oxy_xyz[0:3])
    # print(cat_xyz[-3:],"\n",Oxy_xyz[-3:])

    flagged_oxy = np.zeros(len(Oxy_xyz),dtype=bool)
    flagged_ani = np.zeros(len(ani_xyz),dtype=bool)
    for ii_cat in range(len(cat_xyz)):
        cur_cat = cat_xyz[ii_cat]
        cur_oxy = nearest_neigh_rsqrd(Oxy_xyz-cur_cat,boxL)
        flagged_oxy[cur_oxy<z_lim**2] = True
        cur_ani = nearest_neigh_rsqrd(ani_xyz-cur_cat,boxL)
        flagged_ani[cur_ani<z_lim**2] = True

    ## now flag all the atoms
    # flag the waters
    tmp_lst = np.where(flagged_oxy)[0]
    for ii in range(3):         # add flags for O,H,H
        all_atom_flagged[tmp_lst*3+ii] = True

    # flag anions & cations
    water_Act = len(Oxy_xyz)*3
    all_atom_flagged[water_Act:water_Act+len(cat_xyz)] = True
    tmp_lst = np.where(flagged_ani)[0]
    all_atom_flagged[water_Act+len(cat_xyz):water_Act+len(cat_xyz)+len(ani_xyz)] = flagged_ani

    return all_atom_flagged

def doAngZhist(Tstep,TskipI,TskipB,lammpstrj,itype,jtype,uv,scaling,Qscaling,subroutine,z_lim):
    if (Tstep-TskipI)%(TskipB+1) != 0:
        print("Fatal Error: The time duration chosen (Tstep-TskipI) is not divisible by TskipB+1. Chose again to make sure the code will not run into error in the end.")
    else:

        ###########################################################################################
        ### will calc & sort by population (shell vs bulk) of water,cation,anion densities, and output theta distribution at each z-pos ###
        ###########################################################################################
        if subroutine == 3:
            # get boxL in z direction from the files for initializing arrays
            with open(lammpstrj,"r") as posfile0:
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile0)
            shell_cosQ_Hist = np.zeros((int(round(boxL[2]*scaling+1)),int(round(2*Qscaling+1))))
            shell_ions_Hist = np.zeros((int(round(boxL[2]*scaling+1)),4))         # ind::: 0 O, 1 H, 2 cat, 3 anion
            bulk_cosQ_Hist = np.zeros((int(round(boxL[2]*scaling+1)),int(round(2*Qscaling+1))))
            bulk_ions_Hist = np.zeros((int(round(boxL[2]*scaling+1)),4))         # ind::: 0 O, 1 H, 2 cat, 3 anion

            posfile = open(lammpstrj,"r")
            # skip this many time step before analyzing the data
            for ii in range(TskipI):
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile)

            # keep track of time
            t0 = time.time()
            cur_prcnt = 0

            cur_t = 0
            while cur_t < Tstep-TskipI:
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile)
                og_xyz = np.copy(cur_xyz)
                shell_flag = flag_shell_water(cur_xyz,itype,jtype,z_lim,boxL)

                ##### ***** do it for shell atoms ***** #####
                cur_xyz = og_xyz[shell_flag]                        # selecting only the shell atoms

                # histogram shell atom densities
                for i_hist in range(4):
                    Ion_Z = cur_xyz[cur_xyz[:,0]==(i_hist+1)][:,3]
                    shell_ions_Hist = histZ(Ion_Z,shell_ions_Hist,i_hist,boxL,scaling)

                # extract water pos
                cur_xyz = cur_xyz[(cur_xyz[:,0]==itype) | (cur_xyz[:,0]==jtype)]
                ict = len(cur_xyz[cur_xyz[:,0]==itype])
                cur_xyz = cur_xyz[:,1:4]

                for ii_O in range(ict):
                    cur_O = cur_xyz[ii_O*3]
                    for jj_H in range(1,3):     # find unit vectors for each Hydrogen (2 for each water)
                        cur_H = cur_xyz[ii_O*3+jj_H]
                        O2H_uv = adjustL_conv2UV(cur_H-cur_O,boxL)
                        cosQ = np.sum(uv*O2H_uv)
                        z_ind = int(round((cur_O[2]+boxL[2]/2)*scaling))
                        Q_ind = int(round((cosQ+1)*Qscaling))
                        shell_cosQ_Hist[z_ind,Q_ind] += 1

                ##### ***** do it for shell atoms ***** #####
                bulk_flag = np.logical_not(shell_flag)              # flip all True to False and vice versa
                cur_xyz = og_xyz[bulk_flag]                        # selecting only the shell atoms

                # histogram shell atom densities
                for i_hist in range(4):
                    Ion_Z = cur_xyz[cur_xyz[:,0]==(i_hist+1)][:,3]
                    bulk_ions_Hist = histZ(Ion_Z,bulk_ions_Hist,i_hist,boxL,scaling)

                # extract water pos
                cur_xyz = cur_xyz[(cur_xyz[:,0]==itype) | (cur_xyz[:,0]==jtype)]
                ict = len(cur_xyz[cur_xyz[:,0]==itype])
                cur_xyz = cur_xyz[:,1:4]

                for ii_O in range(ict):
                    cur_O = cur_xyz[ii_O*3]
                    for jj_H in range(1,3):     # find unit vectors for each Hydrogen (2 for each water)
                        cur_H = cur_xyz[ii_O*3+jj_H]
                        O2H_uv = adjustL_conv2UV(cur_H-cur_O,boxL)
                        cosQ = np.sum(uv*O2H_uv)
                        z_ind = int(round((cur_O[2]+boxL[2]/2)*scaling))
                        Q_ind = int(round((cosQ+1)*Qscaling))
                        bulk_cosQ_Hist[z_ind,Q_ind] += 1

                # skip this many timestep in between
                for ii in range(TskipB):
                    cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile)
                cur_t += TskipB + 1

                # check whether to output time info
                if cur_t*100/(Tstep-TskipI) >= cur_prcnt:
                    if cur_prcnt == 0:
                        print("*** Outputting time data. ***\nTstep,TskipI,TskipB,lammpstrj,itype,jtype,uv,scaling,Qscaling",
                              Tstep,TskipI,TskipB,lammpstrj,itype,jtype,uv,scaling,Qscaling,
                              "\n#######################################\nPercent[%] Time[s]")
                    print("{0:3.2f}       {1:3.2f}".format(cur_t*100/(Tstep-TskipI), time.time()-t0))
                    cur_prcnt += 10
                    t0 = time.time()

            ### print cosQ matrix to file
            outname = "shell_Qmat_ij_"+str(itype)+"_"+str(jtype)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"TskipB"+str(TskipB)+"uv"+np.array2string(uv,separator='')[1:-1]+"scaling"+str(scaling)+"Qscaling"+str(Qscaling)+"Rcut"+str(z_lim)+".out"
            tmpstr = "# boxL = "+str(boxL)+"\n"
            with open(outname, "w") as f:
                f.write(tmpstr)
            with open(outname, "ab") as f:
                np.savetxt(f,shell_cosQ_Hist)

            outname2 = "shell_zdist_ij_"+str(itype)+"_"+str(jtype)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"TskipB"+str(TskipB)+"uv"+np.array2string(uv,separator='')[1:-1]+"scaling"+str(scaling)+"Qscaling"+str(Qscaling)+"Rcut"+str(z_lim)+".out"
            with open(outname2,"w") as outfile2:
                bulkOCt = shell_ions_Hist[int(round(boxL[2]/2*scaling)),0]
                bulkHCt = shell_ions_Hist[int(round(boxL[2]/2*scaling)),1]
                bulkCatCt = shell_ions_Hist[int(round(boxL[2]/2*scaling)),2]
                bulkAnCt = shell_ions_Hist[int(round(boxL[2]/2*scaling)),3]
                outfile2.write("# z_pos O H Cat An")  # i is O from H2O, i+2 is cation, j+2 is anion
                for ii in range(int(round(boxL[2]*scaling+1))):
                    tmpstr = str(ii/scaling-boxL[2]/2)+" "+str(shell_ions_Hist[ii,0]/bulkOCt)+" "+str(shell_ions_Hist[ii,1]/bulkHCt)+" "+str(shell_ions_Hist[ii,2]/bulkCatCt)+" "+str(shell_ions_Hist[ii,3]/bulkAnCt)+"\n"
                    outfile2.write(tmpstr)

            outname = "bulk_Qmat_ij_"+str(itype)+"_"+str(jtype)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"TskipB"+str(TskipB)+"uv"+np.array2string(uv,separator='')[1:-1]+"scaling"+str(scaling)+"Qscaling"+str(Qscaling)+"Rcut"+str(z_lim)+".out"
            tmpstr = "# boxL = "+str(boxL)+"\n"
            with open(outname, "w") as f:
                f.write(tmpstr)
            with open(outname, "ab") as f:
                np.savetxt(f,bulk_cosQ_Hist)

            outname2 = "bulk_zdist_ij_"+str(itype)+"_"+str(jtype)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"TskipB"+str(TskipB)+"uv"+np.array2string(uv,separator='')[1:-1]+"scaling"+str(scaling)+"Qscaling"+str(Qscaling)+"Rcut"+str(z_lim)+".out"
            with open(outname2,"w") as outfile2:
                bulkOCt = bulk_ions_Hist[int(round(boxL[2]/2*scaling)),0]
                bulkHCt = bulk_ions_Hist[int(round(boxL[2]/2*scaling)),1]
                bulkCatCt = bulk_ions_Hist[int(round(boxL[2]/2*scaling)),2]
                bulkAnCt = bulk_ions_Hist[int(round(boxL[2]/2*scaling)),3]
                outfile2.write("# z_pos O H Cat An")  # i is O from H2O, i+2 is cation, j+2 is anion
                for ii in range(int(round(boxL[2]*scaling+1))):
                    tmpstr = str(ii/scaling-boxL[2]/2)+" "+str(bulk_ions_Hist[ii,0]/bulkOCt)+" "+str(bulk_ions_Hist[ii,1]/bulkHCt)+" "+str(bulk_ions_Hist[ii,2]/bulkCatCt)+" "+str(bulk_ions_Hist[ii,3]/bulkAnCt)+"\n"
                    outfile2.write(tmpstr)


        ##########################################################
        ### will output theta & phi distribution at each z-pos ###
        ##########################################################
        if subroutine == 2:
            # get boxL in z direction from the files for initializing arrays
            with open(lammpstrj,"r") as posfile0:
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile0)
            cosQP_Hist = np.zeros((int(round((z_lim+boxL[2]/2)*scaling+1)),int(round(2*Qscaling+1)),int(round(2*Qscaling+1))))

            posfile = open(lammpstrj,"r")
            # skip this many time step before analyzing the data
            for ii in range(TskipI):
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile)

            # keep track of time
            t0 = time.time()
            cur_prcnt = 0

            cur_t = 0
            while cur_t < Tstep-TskipI:
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile)
                # extract water pos
                cur_xyz = cur_xyz[(cur_xyz[:,0]==itype) | (cur_xyz[:,0]==jtype)]

                # flag water oxy below the z_lim
                oxy_xyz = cur_xyz[(cur_xyz[:,0]==itype)]
                tmp_lst = oxy_xyz[:,3]<z_lim
                tmp_lst = np.where(tmp_lst)[0]        # [0] as np.where() usually return index of each dimension in an array so, R3 will be [a,b,c] and R1 will be [a] where a,b,c are numpy arrays

                WaterInd2Take = np.zeros(len(cur_xyz),dtype=bool)
                WaterInd2Take[tmp_lst*3] = True           # indices of the Oxygen
                WaterInd2Take[tmp_lst*3+1] = True         # indices of 1st Hydrogen
                WaterInd2Take[tmp_lst*3+2] = True         # indices of 2nd Hydrogen

                # select only water under z_lim
                cur_xyz = cur_xyz[WaterInd2Take]
                ict = len(cur_xyz[cur_xyz[:,0]==itype])
                cur_xyz = cur_xyz[:,1:4]

                for ii_O in range(ict):
                    cur_O = cur_xyz[ii_O*3]
                    z_ind = int(round((cur_O[2]+boxL[2]/2)*scaling))

                    # get Q index
                    jj_H = 1
                    cur_H = cur_xyz[ii_O*3+jj_H]
                    O2H_uv = adjustL_conv2UV(cur_H-cur_O,boxL)
                    cosQ = np.sum(uv*O2H_uv)
                    Q_ind = int(round((cosQ+1)*Qscaling))

                    # get Phi index
                    jj_H = 2
                    cur_H = cur_xyz[ii_O*3+jj_H]
                    O2H_uv = adjustL_conv2UV(cur_H-cur_O,boxL)
                    cosP = np.sum(uv*O2H_uv)
                    P_ind = int(round((cosP+1)*Qscaling))

                    # increment the histogram count
                    cosQP_Hist[z_ind,Q_ind,P_ind] += 1

                # skip this many timestep in between
                for ii in range(TskipB):
                    cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile)
                cur_t += TskipB + 1

                # check whether to output time info
                if cur_t*100/(Tstep-TskipI) >= cur_prcnt:
                    if cur_prcnt == 0:
                        print("*** Outputting time data. ***\nTstep,TskipI,TskipB,lammpstrj,itype,jtype,uv,scaling,Qscaling",
                              Tstep,TskipI,TskipB,lammpstrj,itype,jtype,uv,scaling,Qscaling,
                              "\n#######################################\nPercent[%] Time[s]")
                    print("{0:3.2f}       {1:3.2f}".format(cur_t*100/(Tstep-TskipI), time.time()-t0))
                    cur_prcnt += 10
                    t0 = time.time()

            ### print cosQ matrix to file
            outname = "QPhiMat_ij_"+str(itype)+"_"+str(jtype)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"TskipB"+str(TskipB)+"uv"+np.array2string(uv,separator='')[1:-1]+"scaling"+str(scaling)+"Qscaling"+str(Qscaling)
            np.save(outname,cosQP_Hist)



        ###########################################################################################
        ### will calc water,cation,anion densities, and output theta distribution at each z-pos ###
        ###########################################################################################
        if subroutine == 1:
            # get boxL in z direction from the files for initializing arrays
            with open(lammpstrj,"r") as posfile0:
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile0)
            cosQ_Hist = np.zeros((int(round(boxL[2]*scaling+1)),int(round(2*Qscaling+1))))
            ions_Hist = np.zeros((int(round(boxL[2]*scaling+1)),4))         # 0 for i+2 ion (cation), 1 for j+2 ion, 2 for i, 3 for j

            posfile = open(lammpstrj,"r")
            # skip this many time step before analyzing the data
            for ii in range(TskipI):
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile)

            # keep track of time
            t0 = time.time()
            cur_prcnt = 0

            cur_t = 0
            while cur_t < Tstep-TskipI:
                cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile)
                # extract ions' z-coord
                ion_xyz = np.copy(cur_xyz)
                iIon_Z = ion_xyz[ion_xyz[:,0]==itype+2][:,3]
                jIon_Z = ion_xyz[ion_xyz[:,0]==jtype+2][:,3]
                ions_Hist = histZ(iIon_Z,ions_Hist,0,boxL,scaling)
                ions_Hist = histZ(jIon_Z,ions_Hist,1,boxL,scaling)
                # now do water oxygen and hydrogen
                iIon_Z = ion_xyz[ion_xyz[:,0]==itype][:,3]
                jIon_Z = ion_xyz[ion_xyz[:,0]==jtype][:,3]
                ions_Hist = histZ(iIon_Z,ions_Hist,2,boxL,scaling)
                ions_Hist = histZ(jIon_Z,ions_Hist,3,boxL,scaling)

                # extract water pos
                cur_xyz = cur_xyz[(cur_xyz[:,0]==itype) | (cur_xyz[:,0]==jtype)]
                ict = len(cur_xyz[cur_xyz[:,0]==itype])
                cur_xyz = cur_xyz[:,1:4]

                for ii_O in range(ict):
                    cur_O = cur_xyz[ii_O*3]
                    for jj_H in range(1,3):     # find unit vectors for each Hydrogen (2 for each water)
                        cur_H = cur_xyz[ii_O*3+jj_H]
                        O2H_uv = adjustL_conv2UV(cur_H-cur_O,boxL)
                        cosQ = np.sum(uv*O2H_uv)
                        z_ind = int(round((cur_O[2]+boxL[2]/2)*scaling))
                        Q_ind = int(round((cosQ+1)*Qscaling))
                        cosQ_Hist[z_ind,Q_ind] += 1

                # skip this many timestep in between
                for ii in range(TskipB):
                    cur_xyz, timestamp, Ntotal, boxL = get_wrapped_coord_Tqxyz(posfile)
                cur_t += TskipB + 1

                # check whether to output time info
                if cur_t*100/(Tstep-TskipI) >= cur_prcnt:
                    if cur_prcnt == 0:
                        print("*** Outputting time data. ***\nTstep,TskipI,TskipB,lammpstrj,itype,jtype,uv,scaling,Qscaling",
                              Tstep,TskipI,TskipB,lammpstrj,itype,jtype,uv,scaling,Qscaling,
                              "\n#######################################\nPercent[%] Time[s]")
                    print("{0:3.2f}       {1:3.2f}".format(cur_t*100/(Tstep-TskipI), time.time()-t0))
                    cur_prcnt += 10
                    t0 = time.time()

            ### print cosQ matrix to file
            outname = "Qmat_ij_"+str(itype)+"_"+str(jtype)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"TskipB"+str(TskipB)+"uv"+np.array2string(uv,separator='')[1:-1]+"scaling"+str(scaling)+"Qscaling"+str(Qscaling)+".out"
            tmpstr = "# boxL = "+str(boxL)+"\n"
            with open(outname, "w") as f:
                f.write(tmpstr)
            with open(outname, "ab") as f:
                np.savetxt(f,cosQ_Hist)

            outname1 = "cosQ_ij_"+str(itype)+"_"+str(jtype)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"TskipB"+str(TskipB)+"uv"+np.array2string(uv,separator='')[1:-1]+"scaling"+str(scaling)+"Qscaling"+str(Qscaling)+".out"
            with open(outname1,"w") as outfile1:
                for ii in range(int(round(boxL[2]*scaling+1))):
                    for jj in range(int(round(2*Qscaling+1))):
                        if np.sum(cosQ_Hist[ii,:]) != 0:
                            tmpstr = str(ii/scaling-boxL[2]/2)+" "+str(jj/Qscaling-1)+" "+str(cosQ_Hist[ii,jj]/np.sum(cosQ_Hist[ii,:]))+"\n"
                        else:
                            tmpstr = str(ii/scaling-boxL[2]/2)+" "+str(jj/Qscaling-1)+" "+str(cosQ_Hist[ii,jj])+"\n"
                        outfile1.write(tmpstr)
                    outfile1.write("\n")

            outname2 = "zdist_ij_"+str(itype)+"_"+str(jtype)+"T"+str(Tstep)+"TskipI"+str(TskipI)+"TskipB"+str(TskipB)+"uv"+np.array2string(uv,separator='')[1:-1]+"scaling"+str(scaling)+"Qscaling"+str(Qscaling)+".out"
            with open(outname2,"w") as outfile2:
                bulkCt = np.sum(cosQ_Hist[int(round(boxL[2]/2*scaling))])   # water
                bulkip2Ct = ions_Hist[int(round(boxL[2]/2*scaling)),0]
                bulkjp2Ct = ions_Hist[int(round(boxL[2]/2*scaling)),1]
                bulkiCt = ions_Hist[int(round(boxL[2]/2*scaling)),2]
                bulkjCt = ions_Hist[int(round(boxL[2]/2*scaling)),3]
                outfile2.write("# z_pos i_prob i+2_ion_prob j+2_ion_prob")  # i is O from H2O, i+2 is cation, j+2 is anion
                for ii in range(int(round(boxL[2]*scaling+1))):
                    tmpstr = str(ii/scaling-boxL[2]/2)+" "+str(np.sum(cosQ_Hist[ii])/bulkCt)+" "+str(ions_Hist[ii,0]/bulkip2Ct)+" "+str(ions_Hist[ii,1]/bulkjp2Ct)+" "+str(ions_Hist[ii,2]/bulkiCt)+" "+str(ions_Hist[ii,3]/bulkjCt)+"\n"
                    outfile2.write(tmpstr)



