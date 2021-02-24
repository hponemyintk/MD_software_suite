import numpy as np
import pandas as pd
from lib.math_lib import unwrapBond
from lib.read_coord import get_wrapped_coord,get_Coul,get_Coul_coord, r_sqrd, manual_parsing

def getdE(lammpstrj,CoulDump,RefFile,skipline):
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
    file_name_tail = "_TskipI" +str(skipline)+"RefFile"+RefFile+".out"

    # read in atom IDs and sort them by timestamp
    IDs = manual_parsing(RefFile," ",skipline,dtype=np.float64)
    # get rid of distance/angles infos
    IDs = np.hstack( ( IDs[:,0:4],np.array( [[i+1] for i in range(np.shape(IDs)[0])] ) ) )
    IDs = IDs[np.argsort(IDs[:,0])]
    print(IDs,np.shape(IDs))

    timecounter = 0
    prevStamp = 0
    outfile = open("config2dE"+file_name_tail,"w")
    outfile.write( "# ConfigNum timestamp catID FerriID FerroID dis cur_dE_A cur_dE_B\n")
    for cur_IDs in IDs:
        timestamp,catID,FerriID,FerroID,configN = (int(ii) for ii in cur_IDs)
        #-1 here to account for the fact that we already read in prevStamp. So, it is excluding intial and final points
        # also if timestamp we want is ever 0, range(-1) will return nothing!
        skipN = int(timestamp - prevStamp) if prevStamp == 0 else int(timestamp - prevStamp-1)
        for ii in range(skipN):
            get_wrapped_coord(posfile)
            get_Coul(coulfile)
            # get_Coul_coord(coulfile,posfile)
            timecounter += 1

        cur_Coul = get_Coul(coulfile)
        cur_xyz, timestep, Ntotal, boxL = get_wrapped_coord(posfile)
        timecounter += 1
        # print(int(timestamp),prevStamp,skipN,timecounter)

        ### calc dE
        cur_xyz_ferri = cur_xyz[ (cur_xyz[:,0]==itype) ]
        cur_xyz_ferro = cur_xyz[ (cur_xyz[:,0]==jtype) ]
        cur_Coul_ferri = cur_Coul[ (cur_Coul[:,0]==itype) ]
        cur_Coul_ferro = cur_Coul[ (cur_Coul[:,0]==jtype) ]

        dis = np.sqrt(r_sqrd(cur_xyz_ferri[FerriID][1:4],cur_xyz_ferro[FerroID][1:4], boxL))
        ferri_ML = cur_Coul_ferri[FerriID][1]/icharge*2.
        ferro_ML = cur_Coul_ferro[FerroID][1]/jcharge*2.

        # calculate and store dE_A and dE_B
        cur_dE_A = dq *(E_conv_const * dq/dis + ferri_ML - ferro_ML)
        cur_dE_B = dq *(-E_conv_const * dq/dis + ferro_ML - ferri_ML)
        # print(configN,timestamp,catID,FerriID,FerroID,dis,cur_dE_A,cur_dE_B)

        # cat ferri dis
        cur_xyz_cat = cur_xyz[ (cur_xyz[:,0]==3) ][catID]
        disFerri = np.sqrt(r_sqrd(cur_xyz_ferri[FerriID][1:4],cur_xyz_cat[1:4], boxL))
        disFerro = np.sqrt(r_sqrd(cur_xyz_ferro[FerroID][1:4],cur_xyz_cat[1:4], boxL))
        print(timestamp,catID,FerriID,FerroID,dis,disFerri,disFerro,cur_dE_A,cur_dE_B)
        outfile.write(str(configN)+" "+str(timestamp)+" "+str(catID)+" "+str(FerriID)+" "+str(FerroID)+" "+str(dis)+" "+str(disFerri)+" "+str(disFerro)+" "+str(cur_dE_A)+" "+str(cur_dE_B)+"\n")



        # all done with current timestep so mark it as an old time stamp
        prevStamp = int(timestamp)
    coulfile.close()
    posfile.close()
    outfile.close()

