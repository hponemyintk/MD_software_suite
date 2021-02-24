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

def get3nrstBAngHist(nrstIDs,FullMol_pos,cur_s,boxL,nrst1BAngHist,nrst3BAngHist,QQcol):
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
            if ii_closet == 0:
                nrst1BAngHist[Q_ind,QQcol] += 1

### return only IDs that staisfies given angle constraints
def nrst_constraint(IDs,FullMol_pos,cur_s,boxL,minl,maxl):
    QQ = 0.
    # create boolean arrays for masking
    IDs_mask = np.full(len(IDs),0,dtype=bool)
    for ind,ii in enumerate(IDs):
        Fe_pos = FullMol_pos[(ii*7)]
        tmp_pos =  FullMol_pos[(ii*7+1):(ii*7+7)]
        Nsdis = np.sqrt(nearest_neigh_rsqrd(tmp_pos-cur_s,boxL))
        NsIDs = np.argsort(Nsdis)[0:3]      # select nearest FeN distances
        for ii_closet in range(1):
            Ns_uv = adjustL_conv2UV(tmp_pos[NsIDs[ii_closet]]-Fe_pos,boxL)
            Fes_uv = adjustL_conv2UV(cur_s-Fe_pos,boxL)

            Ns_dot_Fes = np.sum(Ns_uv*Fes_uv)
            if Ns_dot_Fes < -1:
                Ns_dot_Fes = -1.0
            elif Ns_dot_Fes>1:
                Ns_dot_Fes = 1.0
            QQ = np.arccos(Ns_dot_Fes) * 180/np.pi

            if minl<QQ and QQ<maxl:
                IDs_mask[ind] = True
    return IDs[IDs_mask], QQ


# for t* use 2ps and tau = 100ps should be more than enough
def findCationBridge(lammpstrj,Tstep,stype,itype,jtype,K_lim,scaling,subroutine,iminl,imaxl,jminl,jmaxl,isminl,ismaxl,jsminl,jsmaxl,Tdisminl,Tdismaxl):
    xyz_in = open(lammpstrj,'r')
    Hist = np.zeros((30*scaling+1,3))           # 1 for i-s distance, 2 for j-s distance
    isj_angHist = np.zeros(180+1)           # 0 is i and 1 for j statistics
    nrst1BAngHist = np.zeros((180+1,2))         # 0 for ferri/1 for ferro for angle beteween the 1 closest FeCN bonds to 
    nrst3BAngHist = np.zeros((180+1,2))         # 0 for ferri/1 for ferro for angle beteween the 3 closest FeCN bonds to 

    print("lammpstrj,Tstep,stype,itype,jtype,K_lim,scaling,subroutine,iminl,imaxl,jminl,jmaxl,isminl,ismaxl,jsminl,jsmaxl,Tdisminl,Tdismaxl:\n",lammpstrj,Tstep,stype,itype,jtype,K_lim,scaling,subroutine,iminl,imaxl,jminl,jmaxl,isminl,ismaxl,jsminl,jsmaxl,Tdisminl,Tdismaxl)
    if (iminl != np.inf and imaxl != np.inf) & (jminl == np.inf and jmaxl == np.inf):
        print("here doing 1 nrst constraint")
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
                    ### select isj state that are specified by 1 nearest ferri angle
                    if subroutine == 3:
                        iIDs, i1nrstang = nrst_constraint(iIDs,iFullMol_pos,cur_s,boxL,iminl,imaxl)
                    if subroutine == 4:
                        jIDs, o1nrstang = nrst_constraint(jIDs,jFullMol_pos,cur_s,boxL,iminl,imaxl)

                    # print(ii_s,cur_t0,iIDs,len(iIDs),jIDs,len(jIDs))
                    if len(iIDs)>0 and len(jIDs)>0:             # to make sure it only do the analysis using the specified population
                        # histrogram angle between the nearest Fe-N unit vector and Caton-Fe vector
                        get3nrstBAngHist(iIDs,iFullMol_pos,cur_s,boxL,nrst1BAngHist,nrst3BAngHist,0)
                        get3nrstBAngHist(jIDs,jFullMol_pos,cur_s,boxL,nrst1BAngHist,nrst3BAngHist,1)
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


        tmpstr = "nrstnghbr_cont_CationBridgeHist"+"T"+str(Tstep)+"stype"+str(stype)+"itype"+str(itype)+"jtype"+str(jtype)+"scaling"+str(scaling)+"K_lim"+str(K_lim)+"subR"+str(subroutine)+"iminl"+str(iminl)+"imaxl"+str(imaxl)+".out"
        np.savetxt(tmpstr,Hist)
        tmpstr = "nrstnghbr_cont_CationBridgeisj_angHist"+"T"+str(Tstep)+"stype"+str(stype)+"itype"+str(itype)+"jtype"+str(jtype)+"scaling"+str(scaling)+"K_lim"+str(K_lim)+"subR"+str(subroutine)+"iminl"+str(iminl)+"imaxl"+str(imaxl)+".out"
        np.savetxt(tmpstr,isj_angHist)
        tmpstr = "nrstnghbr_cont_CationBridge3nrst_BangHist"+"T"+str(Tstep)+"stype"+str(stype)+"itype"+str(itype)+"jtype"+str(jtype)+"scaling"+str(scaling)+"K_lim"+str(K_lim)+"subR"+str(subroutine)+"iminl"+str(iminl)+"imaxl"+str(imaxl)+".out"
        np.savetxt(tmpstr,nrst3BAngHist)
        tmpstr = "nrstnghbr_cont_CationBridge1nrst_BangHist"+"T"+str(Tstep)+"stype"+str(stype)+"itype"+str(itype)+"jtype"+str(jtype)+"scaling"+str(scaling)+"K_lim"+str(K_lim)+"subR"+str(subroutine)+"iminl"+str(iminl)+"imaxl"+str(imaxl)+".out"
        np.savetxt(tmpstr,nrst1BAngHist)


    if (iminl != np.inf and imaxl != np.inf) & (jminl != np.inf and jmaxl != np.inf):
        print("here doing 3 constraints")
        postfix = "T"+str(Tstep)+"stype"+str(stype)+"itype"+str(itype)+"jtype"+str(jtype)+"scaling"+str(scaling)+"K_lim"+str(K_lim)+"iMinMaxl"+str(iminl)+"_"+str(imaxl)+"jMinMaxl"+str(jminl)+"_"+str(jmaxl)+"isMinMaxl"+str(isminl)+"_"+str(ismaxl)+"jsMinMaxl"+str(jsminl)+"_"+str(jsmaxl)+"TdisMinMaxl"+str(Tdisminl)+"_"+str(Tdismaxl)+".out"
        IDout = open("IDs"+postfix,"w")
        printit = True
        for cur_t0 in range(Tstep):
            cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(xyz_in)

            # comment this out if you don't watnt to print water IDs
            water_pos = cur_pos[cur_pos[:,0]==1][:,1:4]

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

                # do Ferri-cation, Ferro-cation distance constraints here
                if isminl!=np.inf and ismaxl!=np.inf and jsminl!=np.inf and jsmaxl!=np.inf:
                    iIDs = np.where( (idis>isminl) & (idis<ismaxl) & (idis<K_lim) )[0]
                    jIDs = np.where( (jdis>jsminl) & (jdis<jsmaxl) & (jdis<K_lim) )[0]
                    # if len(iIDs)>0 and len(jIDs)>0:
                    #     print(idis[iIDs],jdis[jIDs])

                # select ij total distance within the constrained range
                og_iIDs = np.copy(iIDs)
                og_jIDs = np.copy(jIDs)
                for cur_i in og_iIDs:
                    for cur_j in og_jIDs:
                        cur_dis = r_sqrd(i_pos[cur_i],j_pos[cur_j],boxL)
                        if cur_dis > Tdisminl**2 and cur_dis < Tdismaxl**2:
                            iIDs = np.array([cur_i])
                            jIDs = np.array([cur_j])
                            # print("here***",i_pos[cur_i],j_pos[cur_j],boxL,np.sqrt(cur_dis))

                            ### select isj state that are specified by 1 nearest ferri angle
                            iIDs, i1nrstang = nrst_constraint(iIDs,iFullMol_pos,cur_s,boxL,iminl,imaxl)
                            jIDs, j1nrstang = nrst_constraint(jIDs,jFullMol_pos,cur_s,boxL,jminl,jmaxl)

                            if len(iIDs)>0 and len(jIDs)>0:             # to make sure it only do the analysis using the specified population
                                # histrogram angle between the nearest Fe-N unit vector and Caton-Fe vector
                                get3nrstBAngHist(iIDs,iFullMol_pos,cur_s,boxL,nrst1BAngHist,nrst3BAngHist,0)
                                get3nrstBAngHist(jIDs,jFullMol_pos,cur_s,boxL,nrst1BAngHist,nrst3BAngHist,1)
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
                                        # print(i_pos[ii],j_pos[jj],redox_dis,ind)

                                        # histogram ferrij-cation/ferro-cation angle here
                                        is_uv = adjustL_conv2UV(i_pos[ii]-cur_s,boxL)
                                        js_uv = adjustL_conv2UV(j_pos[jj]-cur_s,boxL)
                                        QQ = np.arccos(np.sum(is_uv*js_uv)) * 180/np.pi
                                        Q_ind = int(round(QQ))
                                        isj_angHist[Q_ind] += 1

                                        if printit == True:
                                            printit = False
                                            IDout.write("### CurT, catID, FerriID, FerroID, Tdis, Tangle, Ferri1nrstAngle, Ferro1nrstAngle, FerriCatDis, FerroCatDis ###\n")
                                        # IDout.write(str(cur_t0)+" "+str(ii_s)+" "+str(ii)+" "+str(jj)+" "+str(redox_dis)+" "+str(QQ)+" "+str(i1nrstang)+" "+str(j1nrstang)+" "+str(idis[ii])+" "+str(jdis[jj])+"\n")

                                        ############################################################################
                                        ### comment this block of code out if you don't want to print water IDs ###
                                        ############################################################################
                                        Rcut = 4 # Angstrom
                                        waterdis = np.sqrt(nearest_neigh_rsqrd(water_pos-cur_s,boxL))
                                        waterIDs = np.where(waterdis<Rcut)[0]
                                        Rcut += 2 # Angstrom
                                        waterdis = np.sqrt(nearest_neigh_rsqrd(water_pos-i_pos[ii],boxL))
                                        waterIDs = np.hstack((waterIDs,np.where(waterdis<Rcut)[0]))
                                        waterdis = np.sqrt(nearest_neigh_rsqrd(water_pos-j_pos[jj],boxL))
                                        waterIDs = np.hstack((waterIDs,np.where(waterdis<Rcut)[0]))
                                        # to make sure there is no duplicate waters
                                        #***** double check this code before submitting!!!
                                        rmDuWater = []
                                        for ii in waterIDs:
                                            if ii not in rmDuWater:
                                                rmDuWater.append(ii)
                                        waterIDs = np.array(rmDuWater)
                                        waterIDs = np.array2string(waterIDs,max_line_width=10000)[1:-1]
                                        IDout.write(str(cur_t0)+" "+str(ii_s)+" "+str(ii)+" "+str(jj)+" "+str(redox_dis)+" "+str(QQ)+" "+str(i1nrstang)+" "+str(j1nrstang)+" "+str(idis[ii])+" "+str(jdis[jj])+" "+waterIDs+"\n")

        tmpstr = "nrstnghbr_cont_CationBridgeHist"+postfix
        np.savetxt(tmpstr,Hist)
        tmpstr = "nrstnghbr_cont_CationBridgeisj_angHist"+postfix
        np.savetxt(tmpstr,isj_angHist)
        tmpstr = "nrstnghbr_cont_CationBridge3nrst_BangHist"+postfix
        np.savetxt(tmpstr,nrst3BAngHist)
        tmpstr = "nrstnghbr_cont_CationBridge1nrst_BangHist"+postfix
        np.savetxt(tmpstr,nrst1BAngHist)


