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
    from lib.read_coord import get_wrapped_coord, r_sqrd
    from lib.math_lib import nearest_neigh_rsqrd, adjustL_conv2UV
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def get3nrstBAngHist(nrstIDs,FullMol_pos,cur_s,boxL,nrst3BAngHist,QQcol):
    for ii in nrstIDs:
        Fe_pos = FullMol_pos[(ii*7)]
        tmp_pos =  FullMol_pos[(ii*7+1):(ii*7+7)]
        Nsdis = np.sqrt(nearest_neigh_rsqrd(tmp_pos-cur_s,boxL))
        NsIDs = np.argsort(Nsdis)[0:3]      # select nearest FeN distances
        for ii_closet in range(3):
            # print("here!!!",Fe_pos,tmp_pos[NsIDs[ii_closet]],cur_s,ii_closet)
            Ns_uv = adjustL_conv2UV(tmp_pos[NsIDs[ii_closet]]-Fe_pos,boxL)
            Fes_uv = adjustL_conv2UV(cur_s-Fe_pos,boxL)

            Ns_dot_Fes = np.sum(Ns_uv*Fes_uv)
            # print("Here!!!",Ns_dot_Fes,)
            # print(Ns_dot_Fes)
            if Ns_dot_Fes < -1:
                Ns_dot_Fes = -1.0
            elif Ns_dot_Fes>1:
                Ns_dot_Fes = 1.0
            QQ = np.arccos(Ns_dot_Fes) * 180/np.pi
            # print(QQ)
            Q_ind = int(round(QQ))
            nrst3BAngHist[Q_ind,QQcol] += 1

# for t* use 2ps and tau = 100ps should be more than enough
def findCationBridge(lammpstrj,Tstep,stype,itype,jtype,K_lim,scaling,subroutine,threshold,coord_num):
    xyz_in = open(lammpstrj,'r')
    Hist = np.zeros((30*scaling+1,3))           # 1 for i-s distance, 2 for j-s distance
    isj_angHist = np.zeros(180+1)           # 0 is i and 1 for j statistics
    nrst3BAngHist = np.zeros((180+1,2))         # 1 for angle beteween the 3 closest FeCN bonds to 

    for cur_t0 in range(Tstep):
        cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(xyz_in)

        s_pos = cur_pos[cur_pos[:,0]==stype][:,1:4]
        i_pos = cur_pos[cur_pos[:,0]==itype][:,1:4]
        j_pos = cur_pos[cur_pos[:,0]==jtype][:,1:4]
        iFullMol_pos = cur_pos[(cur_pos[:,0]==itype)| (cur_pos[:,0]==itype+2)][:,1:4]
        jFullMol_pos = cur_pos[(cur_pos[:,0]==jtype)| (cur_pos[:,0]==jtype+2)][:,1:4]
        # print(iFullMol_pos,jFullMol_pos)
        s_ct = len(s_pos)

        if cur_t0 == 0:
            print(s_ct,len(cur_pos[cur_pos[:,0]==stype+1]),len(cur_pos[cur_pos[:,0]==itype]),len(cur_pos[cur_pos[:,0]==jtype]))

        for ii_s in range(s_ct):
            cur_s = s_pos[ii_s]
            iIDs = np.nan          # keep track of IDs of atoms
            jIDs = np.nan        # keep track of type of each atom in ID list -- will be used to check if there is a bridge or not

            idis = np.sqrt(nearest_neigh_rsqrd(i_pos-cur_s,boxL))
            iIDs = np.where(idis<K_lim)[0]

            jdis = np.sqrt(nearest_neigh_rsqrd(j_pos-cur_s,boxL))
            jIDs = np.where(jdis<K_lim)[0]

            # print(iIDs,jIDs)
            if len(iIDs)>0 and len(jIDs)>0:
                # histrogram angle between the nearest Fe-N unit vector and Caton-Fe vector
                get3nrstBAngHist(iIDs,iFullMol_pos,cur_s,boxL,nrst3BAngHist,0)
                get3nrstBAngHist(jIDs,jFullMol_pos,cur_s,boxL,nrst3BAngHist,1)
                # print("s ID:",ii_s,"\ni ID Dist:",iIDs,idis[iIDs],"\nj ID Dist:",jIDs,jdis[jIDs])

                # add to histogram
                for ii in iIDs:
                    for jj in jIDs:
                        redox_dis = np.sqrt(r_sqrd(j_pos[jj],i_pos[ii], boxL))
                        ind = int(round(redox_dis*scaling))
                        Hist[ind,0] += 1

                        # histogram i-s, j-s distances
                        Hist[int(round(idis[ii]*scaling)),1] += 1
                        Hist[int(round(jdis[jj]*scaling)),2] += 1
                        # print(i_pos[ii],j_pos[jj],redox_dis)

                        # histogram ferrij-cation/ferro-cation angle here
                        is_uv = adjustL_conv2UV(i_pos[ii]-cur_s,boxL)
                        js_uv = adjustL_conv2UV(j_pos[jj]-cur_s,boxL)
                        QQ = np.arccos(np.sum(is_uv*js_uv)) * 180/np.pi
                        Q_ind = int(round(QQ))
                        isj_angHist[Q_ind] += 1

                        if subroutine == 1:
                            leeway = 0.3        # smallest difference bteween exp adjusted rxn distances
                            if (redox_dis > threshold-leeway) and (redox_dis < threshold+leeway):
                                # select water around cation
                                w_pos = cur_pos[cur_pos[:,0]==1][:,1:4]
                                wdis = np.sqrt(nearest_neigh_rsqrd(w_pos-cur_s,boxL))
                                wIDs = np.where(wdis<coord_num)[0]

                                print("*** subroutine 1 printing CurT, catID, FerriID, FerroID, WaterID:::", cur_t0, ii_s, ii, jj, redox_dis, wIDs)

    tmpstr = "CationBridgeHist"+"T"+str(Tstep)+"stype"+str(stype)+"itype"+str(itype)+"jtype"+str(jtype)+"scaling"+str(scaling)+"K_lim"+str(K_lim)+".out"
    np.savetxt(tmpstr,Hist)
    tmpstr = "CationBridgeisj_angHist"+"T"+str(Tstep)+"stype"+str(stype)+"itype"+str(itype)+"jtype"+str(jtype)+"scaling"+str(scaling)+"K_lim"+str(K_lim)+".out"
    np.savetxt(tmpstr,isj_angHist)
    tmpstr = "CationBridge3nrst_BangHist"+"T"+str(Tstep)+"stype"+str(stype)+"itype"+str(itype)+"jtype"+str(jtype)+"scaling"+str(scaling)+"K_lim"+str(K_lim)+".out"
    np.savetxt(tmpstr,nrst3BAngHist)


