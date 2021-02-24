import numpy as np
##############################################################################################
### return square of the distance between a reference atom (i) and selected group of atoms ###
##############################################################################################
# input wrapped position and this will return the R**2 distance between the 2 group of atoms
def r_iN(s_pos,selected_pos,boxL):
    diff = np.abs(s_pos - selected_pos)
    for ii in range(3):
        diff[diff[:,ii]>boxL[ii]/2.0,ii] -= boxL[ii]
    diff = np.sum(diff*diff, axis=1)
    return diff

#################################################
### return a unit vector of an inputed vector ###
#################################################
# will also adjust the unit vector's direction depending on whether it is longer (change sign) or shorter (same sign) than boxL
def adjustL_conv2UV(og,boxL):
    a = abs(og)
    for ii in range(3):
        if a[ii]>boxL[ii]/2.0:
            a[ii] -= boxL[ii]
            a[ii] = -abs(a[ii]) if og[ii]>=0 else abs(a[ii])
        else:
            a[ii] = abs(a[ii]) if og[ii]>=0 else -abs(a[ii])
    return (a)/np.sqrt(np.sum(a*a))



#################################################################
### return an array of squared distances to nearest neighbors ###
#################################################################
def nearest_neigh_rsqrd(dist_diff,boxL):
    dist_diff = np.abs(dist_diff)
    for ii in range(3):
        dist_diff[dist_diff[:,ii]>boxL[ii]/2.0,ii] -= boxL[ii]
    return np.sum(dist_diff*dist_diff,axis=1)


##############################
### xyz, rotation matrices ###
##############################
# tansformation matrices
def x_trans(alpha):
    return np.array([[1,0,0],[0,np.cos(alpha),-np.sin(alpha)],[0,np.sin(alpha),np.cos(alpha)]])
def y_trans(alpha):
    return np.array([[np.cos(alpha),0,np.sin(alpha)],[0,1,0],[-np.sin(alpha),0,np.cos(alpha)]])
def z_trans(alpha):
    return np.array([[np.cos(alpha),-np.sin(alpha),0],[np.sin(alpha),np.cos(alpha),0],[0,0,1]])


#####################
### unwrap bounds ###
#####################
# given a reference atoms (refpos), and a list of atom bonded to it (curpos), make sure they are not bigger than half the box size
def unwrapBond(curpos,refpos,boxL):                # input arrays containing xyz coord of current and ref positions--no other info!
    BondCheck = curpos-refpos
    for ijk in range(3):
        curpos[ (abs(BondCheck[:,ijk]) > boxL[ijk]/2) & (BondCheck[:,ijk]<0), ijk ] += boxL[ijk]
        curpos[ (abs(BondCheck[:,ijk]) > boxL[ijk]/2) & (BondCheck[:,ijk]>0), ijk ] -= boxL[ijk]
    return curpos,refpos


###########################################################
### get the target IDs within cutoff distance of RefLst ###
###########################################################
# here RefLst is the 2D array with rows as each separate reference atoms coords
def get_solvation_shell(target,RefLst,K_lim,boxL):
    for ilst in range(len(RefLst)):
        dis = np.sqrt(nearest_neigh_rsqrd(target-RefLst[ilst],boxL))
        targetIDs = np.where(dis<K_lim)[0] if ilst==0 else np.hstack((targetIDs,np.where(dis<K_lim)[0]))
    # to make sure there is no duplicate IDs
    rmDuIDs = []
    for ii in targetIDs:
        if ii not in rmDuIDs:
            rmDuIDs.append(ii)
    targetIDs = np.array(rmDuIDs)
    print("...",np.sort(targetIDs))
    return targetIDs


#########################################################################################
### Given a ID of the first atom in the molecule, get the data for the whole molecule ###
#########################################################################################
# the data for molecule can be anything (as long as the data array inputted only contained the atoms from the molecule of the same type)
def ID2MolData(ID,data,NatomInMol):
    return data[ (ID*NatomInMol) : (ID*NatomInMol+NatomInMol) ,:]



