# python3 main.py --lammpstrj equ_nve_woWater.lammpstrj --Tstep 100 --T0ct 3 --corr_t 10 --tau 70
import numpy as np
import sys

# move these to main eventually
lammpstrj = "equ_nve_woCN.lammpstrj"
CoulDump = "Coul_NoBond_pe_atom_spce_woCN.dump"
Tstep = 10 # 2 for Codecheck # 10 for short debug # 10000 for 1ns run with snapshot every 50timestep # 100000 for 10ns with every 50timestep dump
TskipI = 0 # 500 for skipping 50ps of data
TskipB = 0 # skip between shanpshots
T0ct = 500
corr_t = 50
tau = 250  # usually people look at 10ps window [250*50=12.5ps]
x_1sideL = 200
K_lim = 0
stype = False
itype = False
jtype = False
uv = False
iminl,imaxl,jminl,jmaxl,isminl,ismaxl,jsminl,jsmaxl,Tdisminl,Tdismaxl = np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,0.,np.inf
for i in range(len(sys.argv)):
    if sys.argv[i] == "--routine":              # 1 MSD, 2 VACF, 10 Dis Corr
        routine = int(sys.argv[i+1])
    if sys.argv[i] == "--subroutine":           # to choose what kind of analysis we want in each routine
        subroutine = int(sys.argv[i+1])
    if sys.argv[i] == "--lammpstrj":            # lammps traj file contianing xyz positions and image flags
        lammpstrj = (sys.argv[i+1])
    if sys.argv[i] == "--CoulDump":             # dump file with Coulomb potential values from lammps
        CoulDump = (sys.argv[i+1])
    if sys.argv[i] == "--RefFile":              # reference file (for lambda phob analysis this is bulk madelung dist)
        RefFile = (sys.argv[i+1])
    if sys.argv[i] == "--Tstep":                # total of snapshot to process (this will be the total number of snapshots outputted from simulation)
        Tstep = int(sys.argv[i+1])
    if sys.argv[i] == "--TskipI":               # number of snapshots to skip initially
        TskipI = int(sys.argv[i+1])
    if sys.argv[i] == "--TskipB":               # number of snapshots to skip between each analysis (so, the code will be doing analysis every TskipB+1 step)
        TskipB = int(sys.argv[i+1])
    if sys.argv[i] == "--T0ct":                 # total number of t0 for time window averaging
        T0ct = int(sys.argv[i+1])
    if sys.argv[i] == "--corr_t":               # total number of t0 for time window averaging
        corr_t = int(sys.argv[i+1])
    if sys.argv[i] == "--tau":                  # the duration of time window [in term of snapshots] to average over
        tau = int(sys.argv[i+1])
    if sys.argv[i] == "--x_1sideL":             # the number of bins to process for each side in Marcus energy gap calculation
        x_1sideL = int(sys.argv[i+1])
    if sys.argv[i] == "--K_lim":                # The intermolecular distance between the redox species that we are considering
        K_lim = float(sys.argv[i+1])
    if sys.argv[i] == "--scaling":              # scaling factor for rdf binning
        scaling = int(sys.argv[i+1])
    if sys.argv[i] == "--stype":                # Do anaylysis on only 1 type of atom
        stype = int(sys.argv[i+1])
    if sys.argv[i] == "--itype":                # Do pair analysis
        itype = int(sys.argv[i+1])
    if sys.argv[i] == "--jtype":                # Do pair analysis
        jtype = int(sys.argv[i+1])
    if sys.argv[i] == "--charge":                  # the duration of time window [in term of snapshots] to average over
        charge = float(sys.argv[i+1])
    if sys.argv[i] == "--threshold":                # Do pair analysis
        threshold = float(sys.argv[i+1])
    if sys.argv[i] == "--coord_num":                # coordination number of s i j
        coord_num = float(sys.argv[i+1])
    if sys.argv[i] == "--uv":                   # input the unit vector to get Cos(theta)
        uv = str(sys.argv[i+1])
        uv = uv.split(',')
        uv = np.array([float(x) for x in uv])
    if sys.argv[i] == "--iminl":                  # the duration of time window [in term of snapshots] to average over
        iminl = float(sys.argv[i+1])
    if sys.argv[i] == "--imaxl":                  # the duration of time window [in term of snapshots] to average over
        imaxl = float(sys.argv[i+1])
    if sys.argv[i] == "--jminl":                  # the duration of time window [in term of snapshots] to average over
        jminl = float(sys.argv[i+1])
    if sys.argv[i] == "--jmaxl":                  # the duration of time window [in term of snapshots] to average over
        jmaxl = float(sys.argv[i+1])
    if sys.argv[i] == "--isminl":                  # the duration of time window [in term of snapshots] to average over
        isminl = float(sys.argv[i+1])
    if sys.argv[i] == "--ismaxl":                  # the duration of time window [in term of snapshots] to average over
        ismaxl = float(sys.argv[i+1])
    if sys.argv[i] == "--jsminl":                  # the duration of time window [in term of snapshots] to average over
        jsminl = float(sys.argv[i+1])
    if sys.argv[i] == "--jsmaxl":                  # the duration of time window [in term of snapshots] to average over
        jsmaxl = float(sys.argv[i+1])
    if sys.argv[i] == "--Tdisminl":                  # the duration of time window [in term of snapshots] to average over
        Tdisminl = float(sys.argv[i+1])
    if sys.argv[i] == "--Tdismaxl":                  # the duration of time window [in term of snapshots] to average over
        Tdismaxl = float(sys.argv[i+1])


if routine == 1:                                # make sure the lammpstrj you link is the velocity (not the xyz file)
    from lib import vacf
    vacf.doVACF(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau,stype)

if routine == 2:
    from lib import msd
    msd.doMSD(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau,stype)

if routine == 3:
    from lib import rdf
    rdf.doRDF(lammpstrj,uv,Tstep,TskipB,itype,jtype,K_lim,scaling)

if routine == 4:
    # python3 main.py --routine 4 --lammpstrj equ_nve_woWater.lammpstrj --Tstep 101 --T0ct 4 --corr_t 10 --tau 50 --stype 3 --itype 5 --jtype 8 --K_lim 5 --threshold .1
    # RefFile::: 0 for ferriRoLammpstrj (IMTqxyzixiyiz) ; 1 for TqxyzCoul
    from lib import msd_pop as msd_pop
    msd_pop.doMSD_pop(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau,stype,itype,jtype,K_lim,threshold)

if routine == 5:
    # python3 main.py --routine 5 --lammpstrj equ_nve_woWater.lammpstrj --Tstep 10001 --T0ct 250 --corr_t 30 --tau 1000 --stype 3 --itype 8 --jtype 8 --K_lim 5.5 --threshold 20 --coord_num 5  --RefFile 0
    coord_num = int(round(coord_num))
    print("lammpstrj,Tstep,T0ct,corr_t,tau,stype,itype,jtype,K_lim,threshold,coord_num,RefFile\n",lammpstrj,Tstep,T0ct,corr_t,tau,stype,itype,jtype,K_lim,threshold,coord_num,RefFile)
    from lib import resid_time as resid_time
    resid_time.doResid_time(lammpstrj,Tstep,T0ct,corr_t,tau,stype,itype,jtype,K_lim,threshold,coord_num,RefFile)

if routine == 6:
    # python3 main.py --routine 6 --lammpstrj equ_nve_woWater.lammpstrj --Tstep 10001 --T0ct 250 --corr_t 30 --tau 1000 --stype 3 --itype 8 --jtype 8 --K_lim 5.5 --threshold 20 --coord_num 5
    coord_num = int(round(coord_num))
    from lib import long_time_pair as long_time_pair
    long_time_pair.find_long_time_pair(lammpstrj,Tstep,T0ct,corr_t,tau,stype,itype,K_lim,threshold,coord_num,scaling)

if routine == 7:
    # python3 main.py --routine 7 --lammpstrj equ_nve_woWater.lammpstrj --Tstep 10001 --stype 3 --itype 8 --K_lim 5.5 --threshold 20 --scaling 10
    from lib import topology_hist as topology_hist
    topology_hist.DoHist(lammpstrj,Tstep,stype,itype,K_lim,threshold,scaling)

if routine == 8:
    # python3 main.py --routine 8 --lammpstrj equ_nve_wC.lammpstrj --RefFile Ref_O_dist --Tstep 10001 --T0ct 300 --corr_t 30 --tau 200 --stype 4 --itype 1 --scaling 10
    from lib import lambda_phob as lambda_phob
    lambda_phob.DoLambda(lammpstrj,RefFile,Tstep,T0ct,corr_t,tau,stype,itype,scaling)

if routine == 9:
    # threshold here is the n in (x-x_m)^n ;so, threshold=8 will go from power of 1 to 8 AND first row in the output will be mean of the dist
    # python3 main.py --routine 9 --lammpstrj equ_nve_woCN.lammpstrj --CoulDump Coul_NoBond_pe_atom_spce_woCN.dump --Tstep 50001 --TskipI 0 --TskipB 0 --x_1sideL 500 --threshold 8
    from lib import MarcusMeanVar as MarcusMeanVar
    MarcusMeanVar.DoMeanVar(lammpstrj,CoulDump,Tstep,TskipI,TskipB,x_1sideL,threshold)

if routine == 10:
    from lib import dis_corr as discor
    discor.dis_corr(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau)

if routine == 11:
    # python3 main.py --routine 11 --CoulDump Coul_NoBond_pe_atom_spce_woCN.dump --Tstep 50000 --itype 5 --charge -3
    from lib import MadeL_mean_var
    MadeL_mean_var.doMadelungMeanVar(CoulDump,Tstep,itype,charge)

if routine == 12:
    # python3 main.py --routine 12 --Tstep 10001 --lammpstrj HLseeds_old/newseed1 --itype 5 --jtype 7 --scaling 10 --K_lim 5
    from lib import bond_angle_check
    bond_angle_check.doBondAnlgeHist(Tstep,lammpstrj,itype,jtype,scaling,K_lim)

if routine == 13:
    # this code is using K_lim as the scaling factor for cosQ
    # python3 main.py --routine 13 --Tstep 10001 --TskipI 50 --TskipB 0 --lammpstrj tmp.lammpstrj --itype 1 --jtype 2 --uv 0,0,-1 --scaling 10 --K_lim 10 --subroutine 2 --threshold -33
    # python3 main.py --routine 13 --Tstep 100 --TskipI 0 --TskipB 0 --lammpstrj tmp.lammpstrj --itype 1 --jtype 2 --uv 0,0,1 --scaling 10 --K_lim 10 --subroutine 2 --threshold -33
    # python3 main.py --routine 13 --Tstep 100 --TskipI 0 --TskipB 0 --lammpstrj tmp.lammpstrj --itype 1 --jtype 2 --uv 0,0,1 --scaling 10 --K_lim 10 --subroutine 3 --threshold -3
    from lib import HER_angZ_dist
    HER_angZ_dist.doAngZhist(Tstep,TskipI,TskipB,lammpstrj,itype,jtype,uv,scaling,K_lim,subroutine,threshold)

if routine == 14:
    # python3 main.py --routine 14 --Tstep 10001 --lammpstrj HLseeds_old/newseed1 --itype 3 --jtype 4 --scaling 10 --K_lim 3
    from lib import bond_angle_check
    bond_angle_check.doBondAnlgeHist(Tstep,lammpstrj,itype,jtype,scaling,K_lim)

if routine == 15:
    # check i in j's rdf or not. Then, check how many s are around i. Using threshold as is cut-off, and K_lim as ji cutoff.
    # python3 main.py --routine 15 --lammpstrj Coul_wTopo_qxyz.dump --Tstep 10000 --stype 4 --jtype 3 --x_1sideL 1000 --K_lim 3
    from lib import MadeLShellBulk
    MadeLShellBulk.doMadeLHist(lammpstrj,Tstep,itype,jtype,stype,x_1sideL,K_lim,threshold)

if routine == 16:
    # K_lim is used as r_cut between ferri-cation and ferro-cation
    # s,i,j are cation, ferri, ferro
    # scaling is used as scaling factor in historgramming ferri-ferro intermol distance
    # subroutine 1 will print out redox centers withint threshold+-leeway distance (along with ID of water mol specified in coord_num radius)
    # python3 main.py --routine 16 --lammpstrj tmp.lammpstrj --Tstep 1 --stype 3 --itype 5 --jtype 8 --K_lim 5 --scaling 10 --subroutine 1 --threshold 10 --coord_num 5
    from lib import CationBridge
    CationBridge.findCationBridge(lammpstrj,Tstep,stype,itype,jtype,K_lim,scaling,subroutine,threshold,coord_num)

if routine == 17:
    # python3 main.py --routine 17 --lammpstrj DFT_Ferri --RefFile MD_Ferri
    from lib import MDconfig2DFTconfig
    MDconfig2DFTconfig.ConvertConfig(lammpstrj,RefFile)

if routine == 18:
    # K_lim is used as r_cut between ferri-cation and ferro-cation
    # s,i,j are cation, ferri, ferro
    # scaling is used as scaling factor in historgramming ferri-ferro intermol distance
    # subroutine 3 for ferri & 4 for ferro
    # iminl and imaxl limits are angle in degree (the code will only constraint one angle--the type specified in subroutine && if jminl,jmaxl is given, it will constraint 2 peaks together).
    # python3 main.py --routine 18 --lammpstrj tmp.lammpstrj --Tstep 100 --stype 3 --itype 5 --jtype 8 --K_lim 6.5 --scaling 10 --subroutine 3 --iminl 10 --imaxl 25
    # python3 main.py --routine 18 --lammpstrj tmp.lammpstrj --Tstep 100 --stype 3 --itype 5 --jtype 8 --K_lim 6.5 --scaling 10 --subroutine 3 --iminl 10 --imaxl 25
    from lib import nrstnghbrCdist
    nrstnghbrCdist.findCationBridge(lammpstrj,Tstep,stype,itype,jtype,K_lim,scaling,subroutine,iminl,imaxl,jminl,jmaxl,isminl,ismaxl,jsminl,jsmaxl,Tdisminl,Tdismaxl)

if routine == 19:
    # python3 main.py --routine 19 --lammpstrj tmp --TskipI 1 --RefFile tmp
    # fname=MD2DFT1.xyz;lct=$(wc -l $fname|awk '{print $1}');sed "s/abcd1/$lct/g" header|cat - $fname >tmp.xyz;vmdFerri tmp.xyz
    # Tskip is used here to determine how many lines to skip in the beginning of the RefFile
    from lib import getConfig
    getConfig.makeConfig(lammpstrj,RefFile,TskipI)

if routine == 20:
    # python3 main.py --routine 20 --lammpstrj equ_nve_woWater.lammpstrj --CoulDump Coul_NoBond_pe_atom_spce_woCN.dump ---TskipI 0 --RefFile ID2dEref
    # fname=MD2DFT1.xyz;lct=$(wc -l $fname|awk '{print $1}');sed "s/abcd1/$lct/g" header|cat - $fname >tmp.xyz;vmdFerri tmp.xyz
    from lib import IDs2dE
    IDs2dE.getdE(lammpstrj,CoulDump,RefFile,TskipI)

if routine == 21:
    """
        This code finds the cation bridged state with desired dE ranged (as specified in uv array)
        uv here is the list of dE ranges to check in the following order: dE1_min, dE1_max, dE2_min, dE2_max etc.
    """
    # python3 main.py --routine 21 --lammpstrj equ_nve_woWater.lammpstrj --CoulDump Coul_NoBond_pe_atom_spce_woCN.dump --Tstep 1000 --stype 3 --itype 5 --jtype 8 --uv -53,-45,-40,-37,-24,0 --K_lim 6.6
    # fname=MD2DFT1.xyz;lct=$(wc -l $fname|awk '{print $1}');sed "s/abcd1/$lct/g" header|cat - $fname >tmp.xyz;vmdFerri tmp.xyz
    from lib import dE2config
    dE2config.getIDs(lammpstrj,CoulDump,Tstep,stype,itype,jtype,uv,K_lim)

if routine == 22:
    """
        Given timestep, catID, FerriID, FerroID, this code will find the solvated complex (with or without cations depending on what ref file you are using).
        It will then move the complex to origin while wrapping bonds accordingly.
        This is different from routine 19 in that this will find the immediate sovlation shell around Ferri and Ferro automatically (including the cations and anions around them).

        keysword definitions:
            lammpstrj = snapshot with water but not with CN on ferri/ferro
            CoulDump  = snapshot without water but include full redox centers           # I have to used two files as I deleted the og files for these
    """
    # python3 main.py --routine 22 --lammpstrj equ_nve_woCN_1.lammpstrj --CoulDump equ_nve_woWater.lammpstrj --TskipI 1 --RefFile Ref1 --K_lim 8.0
    # Tskip is used here to determine how many lines to skip in the beginning of the RefFile
    from lib import configwWaterfromIDs
    configwWaterfromIDs.makeConfig(lammpstrj,CoulDump,RefFile,TskipI,K_lim)

if routine == 23:
    # threshold here is the n in (x-x_m)^n ;so, threshold=8 will go from power of 1 to 8 AND first row in the output will be mean of the dist
    # python3 main.py --routine 23 --lammpstrj equ_nve.lammpstrj.10 --CoulDump Coul_NoBond_pe_atom_spce.dump.10 --Tstep 10 --TskipI 0 --TskipB 0 --x_1sideL 500 --threshold 8
    from lib import MarcusMeanVarFullAtom as MarcusMeanVarFullAtom
    MarcusMeanVarFullAtom.DoMeanVar(lammpstrj,CoulDump,Tstep,TskipI,TskipB,x_1sideL,threshold)

if routine == 24:
    # threshold here is the n in (x-x_m)^n ;so, threshold=8 will go from power of 1 to 8 AND first row in the output will be mean of the dist
    # python3 main.py --routine 24 --lammpstrj equ_nve.lammpstrj.10 --K_lim 14
    from lib import RealCoul as RealCoul
    RealCoul.DoMeanVar(lammpstrj,K_lim)

if routine == 25:
    # threshold here is the n in (x-x_m)^n ;so, threshold=8 will go from power of 1 to 8 AND first row in the output will be mean of the dist
    # python3 main.py --routine 25 --lammpstrj Ereorg.lammpstrj --CoulDump Ereorg.dump --Tstep 10 --TskipI 0 --TskipB 0 --x_1sideL 500 --threshold 8
    from lib import heterodE as heterodE
    heterodE.doHeterodE(lammpstrj,CoulDump,Tstep,TskipI,TskipB,x_1sideL,threshold)
