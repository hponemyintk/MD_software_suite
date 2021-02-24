import numpy as np
#############################################################
### Truncate float to float with specified decimal points ###
#############################################################
def truncate(n, decimals=0):
    multiplier = 10 ** int(decimals)
    return round(n * multiplier) / float(multiplier)

##################################
### Calculate Wrapped Distance ###
##################################
def r_sqrd(cur_coord,init_coord, boxL):
        dx = abs(cur_coord[0]-init_coord[0])
        if dx>boxL[0]/2.0:
                dx = boxL[0]-dx
        dy = abs(cur_coord[1]-init_coord[1])
        if dy>boxL[1]/2.0:
                dy = boxL[1]-dy
        dz = abs(cur_coord[2]-init_coord[2])
        if dz>boxL[2]/2.0:
                dz = boxL[2]-dz
        r_sq = dx**2+dy**2+dz**2
        return r_sq

#####################
### get Coul data ###
#####################
# dump file format::: id mol type q x y z ix iy iz
def get_Coul(coulfile):
    ### READ IN COUL VALUES ###
    # read in header
    for i in range(0,3):
        coulfile.readline()
    Ntotal = int(coulfile.readline())
    for i in range(0,5):
        coulfile.readline()

    cur_Coul = np.zeros((Ntotal,2))
    # read in Coul values
    for i_coul in range(0,Ntotal):
        tmp = coulfile.readline().split()
        for j_coul in range(0,2):
            cur_Coul[i_coul][j_coul] = float(tmp[j_coul+1])
    return cur_Coul


##################################
### get Coul and position data ###
##################################
# dump file format::: id mol type q x y z ix iy iz
def get_Coul_coord(coulfile,posfile):
    ### READ IN LAMMPSTRJ FILE ###
    # read in header first
    posfile.readline()
    tcur_Tstep = int(posfile.readline())
    posfile.readline()
    Ntotal = int(posfile.readline())
    posfile.readline()
    boxL =  np.zeros(3)             # to store x,y,z box lengths
    for i in range(0,3):
        boxDim = posfile.readline().rsplit()
        boxL[i] = (float(boxDim[1]) - float(boxDim[0]))
    posfile.readline()

    # now read in coord
    cur_xyz = np.zeros((Ntotal,4))
    for ii in range(0,Ntotal):
        tmp = posfile.readline().split()
        cur_xyz[ii][0] = tmp[2]
        for jj in range(4,7):
            cur_xyz[ii][jj-3] = float(tmp[jj]) #+(float(tmp[jj+3])*boxL[jj-4])

    ### READ IN COUL VALUES ###
    # read in header
    for i in range(0,9):
        coulfile.readline()

    cur_Coul = np.zeros((Ntotal,2))
    # read in Coul values
    for i_coul in range(0,Ntotal):
        tmp = coulfile.readline().split()
        for j_coul in range(0,2):
            cur_Coul[i_coul][j_coul] = float(tmp[j_coul+1])

    return cur_xyz, cur_Coul, tcur_Tstep, Ntotal, boxL



###################################
### get unwrapped position data ###
###################################
# dump file format::: id mol type q x y z ix iy iz
def get_unwrapped_coord(posfile):
    ### READ IN LAMMPSTRJ FILE ###
    # read in header first
    posfile.readline()
    tcur_Tstep = int(posfile.readline())
    posfile.readline()
    Ntotal = int(posfile.readline())
    posfile.readline()
    boxL =  np.zeros(3)             # to store x,y,z box lengths
    for i in range(0,3):
        boxDim = posfile.readline().rsplit()
        boxL[i] = (float(boxDim[1]) - float(boxDim[0]))
    posfile.readline()

    # now read in coord
    cur_xyz = np.zeros((Ntotal,4))
    for ii in range(0,Ntotal):
        tmp = posfile.readline().split()
        cur_xyz[ii][0] = tmp[2]
        for jj in range(4,7):
            cur_xyz[ii][jj-3] = float(tmp[jj])+(float(tmp[jj+3])*boxL[jj-4])
    return cur_xyz, tcur_Tstep, Ntotal, boxL

##################################
### get wrapped position data ###
##################################
# dump file format::: id mol type q x y z ix iy iz
# return::: Txyz (wrapped), tcur_Tstep, Ntotal, boxL
def get_wrapped_coord(posfile):
    ### READ IN LAMMPSTRJ FILE ###
    # read in header first
    posfile.readline()
    tcur_Tstep = int(posfile.readline())
    posfile.readline()
    Ntotal = int(posfile.readline())
    posfile.readline()
    boxL =  np.zeros(3)             # to store x,y,z box lengths
    for i in range(0,3):
        boxDim = posfile.readline().rsplit()
        boxL[i] = (float(boxDim[1]) - float(boxDim[0]))
    posfile.readline()

    # now read in coord
    cur_xyz = np.zeros((Ntotal,4))
    for ii in range(0,Ntotal):
        tmp = posfile.readline().split()
        cur_xyz[ii][0] = tmp[2]
        for jj in range(4,7):
            cur_xyz[ii][jj-3] = float(tmp[jj]) #+(float(tmp[jj+3])*boxL[jj-4])
            # print(cur_xyz[ii])
    return cur_xyz, tcur_Tstep, Ntotal, boxL

######################################################################
### get wrapped position and Coulomb data from lammpstrj+Madel dump ###
######################################################################
# dump file format::: id mol type q x y z ix iy iz c_peratom (PerAtomCoulData)
def get_wrapped_coord_Madel(posfile):
    ### READ IN LAMMPSTRJ FILE ###
    # read in header first
    posfile.readline()
    tcur_Tstep = int(posfile.readline())
    posfile.readline()
    Ntotal = int(posfile.readline())
    posfile.readline()
    boxL =  np.zeros(3)             # to store x,y,z box lengths
    for i in range(0,3):
        boxDim = posfile.readline().rsplit()
        boxL[i] = (float(boxDim[1]) - float(boxDim[0]))
    posfile.readline()

    # now read in coord
    xyz_MadeL = np.zeros((Ntotal,5))     # 0 atom type, 1-3 xyz pos, 4 Coul potential value
    for ii in range(0,Ntotal):
        tmp = posfile.readline().split()
        xyz_MadeL[ii][0] = tmp[2]
        xyz_MadeL[ii][4] = float(tmp[10])/float(tmp[3])*2         # convert Coul potential to Madelung potential
        for jj in range(4,7):
            xyz_MadeL[ii][jj-3] = float(tmp[jj]) #+(float(tmp[jj+3])*boxL[jj-4])
            # print(xyz_MadeL[ii])
    return xyz_MadeL, tcur_Tstep, Ntotal, boxL

######################################################################
### get wrapped atom type & position ###
######################################################################
# dump file format::: type q x y z
def get_wrapped_coord_Tqxyz(posfile):
    # read in header first
    posfile.readline()
    tcur_Tstep = int(posfile.readline())
    posfile.readline()
    Ntotal = int(posfile.readline())
    posfile.readline()
    boxL =  np.zeros(3)             # to store x,y,z box lengths
    for i in range(0,3):
        boxDim = posfile.readline().rsplit()
        boxL[i] = (float(boxDim[1]) - float(boxDim[0]))
    posfile.readline()

    # now read in coord
    TqxyzCoul = np.zeros((Ntotal,4))     # 0 atom type, 1-3 xyz pos
    for ii in range(0,Ntotal):
        tmp = posfile.readline().split()
        TqxyzCoul[ii][0] = tmp[0]
        for jj in range(2,5):
            TqxyzCoul[ii][jj-1] = float(tmp[jj])
    return TqxyzCoul, tcur_Tstep, Ntotal, boxL

#######################
### get wrapped xyz ###
#######################
# dump file format::: type x y z
# return::: Txyz, tcur_Tstep, Ntotal, boxL
def get_wrapped_xyz(posfile,boxL):
    ### READ IN LAMMPSTRJ FILE ###
    # read in header first
    Ntotal = int(posfile.readline())
    tcur_Tstep = int(posfile.readline().split()[-1])

    # now read in coord
    cur_xyz = np.zeros((Ntotal,4))
    for ii in range(0,Ntotal):
        tmp = posfile.readline().split()
        for jj in range(4):
            cur_xyz[ii][jj] = float(tmp[jj])
    return cur_xyz, tcur_Tstep, Ntotal, boxL

##########################
### read qxyzCoul data ###
##########################
# dump file format::: type q x y z c_peratom (PerAtomCoulData)
# return::: type q x y z c_peratom (PerAtomCoulData)
def get_TqxyzCoul(posfile):
    # read in header first
    posfile.readline()
    tcur_Tstep = int(posfile.readline())
    posfile.readline()
    Ntotal = int(posfile.readline())
    posfile.readline()
    boxL =  np.zeros(3)             # to store x,y,z box lengths
    for i in range(0,3):
        boxDim = posfile.readline().rsplit()
        boxL[i] = (float(boxDim[1]) - float(boxDim[0]))
    posfile.readline()

    # now read in coord
    TqxyzCoul = np.zeros((Ntotal,6))     # 0 atom type, 1 charge, 2-4 xyz pos, 5 CoulPE
    for ii in range(0,Ntotal):
        tmp = posfile.readline().split()
        for jj in range(len(tmp)):
            TqxyzCoul[ii][jj] = float(tmp[jj])
    return TqxyzCoul, tcur_Tstep, Ntotal, boxL


###################################################
### load text file with uneven columns as array ###
###################################################
# this will pad the missing data field at the end of the columns as np.nan
def manual_parsing(filename,delim,skipline,dtype):
    out = list()
    lengths = list()
    with open(filename,'r') as ins:
        for ct,line in enumerate(ins):
            if ct <skipline:
                pass
            else:
                l = line.split()
                # l = line.split(delim)
                out.append(l)
                lengths.append(len(l))
    lim = np.max(lengths)
    for l in out:
        while len(l)<lim:
            l.append("nan")
    # print(out)
    return np.array(out,dtype=dtype)
