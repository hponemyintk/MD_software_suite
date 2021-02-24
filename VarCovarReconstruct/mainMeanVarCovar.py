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
for i in range(len(sys.argv)):
    if sys.argv[i] == "--routine":              # 1 MSD, 2 VACF, 10 Dis Corr
        routine = int(sys.argv[i+1])
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
        coord_num = int(sys.argv[i+1])
    if sys.argv[i] == "--uv":                   # input the unit vector to get Cos(theta)
        uv = str(sys.argv[i+1])
        uv = uv.split(',')
        uv = np.array([float(x) for x in uv])


if routine == 1:                                # make sure the lammpstrj you link is the velocity (not the xyz file)
    from libMeanVarCovar import vacf
    vacf.doVACF(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau,stype)

if routine == 2:
    from libMeanVarCovar import msd
    msd.doMSD(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau,stype)

if routine == 3:
    from libMeanVarCovar import rdf
    rdf.doRDF(lammpstrj,uv,Tstep,TskipB,itype,jtype,K_lim,scaling)

if routine == 4:
    # python3 main.py --routine 4 --lammpstrj equ_nve_woWater.lammpstrj --Tstep 101 --T0ct 4 --corr_t 10 --tau 50 --stype 3 --itype 5 --jtype 8 --K_lim 5 --threshold .1
    from libMeanVarCovar import msd_pop as msd_pop
    msd_pop.doMSD_pop(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau,stype,itype,jtype,K_lim,threshold)

if routine == 5:
    # python3 main.py --routine 5 --lammpstrj equ_nve_woWater.lammpstrj --Tstep 10001 --T0ct 250 --corr_t 30 --tau 1000 --stype 3 --itype 8 --jtype 8 --K_lim 5.5 --threshold 20 --coord_num 5
    from libMeanVarCovar import resid_time as resid_time
    resid_time.doResid_time(lammpstrj,Tstep,T0ct,corr_t,tau,stype,itype,jtype,K_lim,threshold,coord_num)

if routine == 6:
    # python3 main.py --routine 6 --lammpstrj equ_nve_woWater.lammpstrj --Tstep 10001 --T0ct 250 --corr_t 30 --tau 1000 --stype 3 --itype 8 --jtype 8 --K_lim 5.5 --threshold 20 --coord_num 5
    from libMeanVarCovar import long_time_pair as long_time_pair
    long_time_pair.find_long_time_pair(lammpstrj,Tstep,T0ct,corr_t,tau,stype,itype,K_lim,threshold,coord_num,scaling)

if routine == 7:
    # python3 main.py --routine 7 --lammpstrj equ_nve_woWater.lammpstrj --Tstep 10001 --stype 3 --itype 8 --K_lim 5.5 --threshold 20 --scaling 10
    from libMeanVarCovar import topology_hist as topology_hist
    topology_hist.DoHist(lammpstrj,Tstep,stype,itype,K_lim,threshold,scaling)

if routine == 8:
    # python3 main.py --routine 8 --lammpstrj equ_nve_wC.lammpstrj --RefFile Ref_O_dist --Tstep 10001 --T0ct 300 --corr_t 30 --tau 200 --stype 4 --itype 1 --scaling 10
    from libMeanVarCovar import lambda_phob as lambda_phob
    lambda_phob.DoLambda(lammpstrj,RefFile,Tstep,T0ct,corr_t,tau,stype,itype,scaling)

if routine == 9:
    # threshold here is the n in (x-x_m)^n ;so, threshold=8 will go from power of 1 to 8 AND first row in the output will be mean of the dist
    # python3 main.py --routine 9 --lammpstrj equ_nve_woCN.lammpstrj --CoulDump Coul_NoBond_pe_atom_spce_woCN.dump --Tstep 50001 --TskipI 0 --TskipB 0 --x_1sideL 500 --threshold 8
    from libMeanVarCovar import MarcusMeanVar as MarcusMeanVar
    MarcusMeanVar.DoMeanVar(lammpstrj,CoulDump,Tstep,TskipI,TskipB,x_1sideL,threshold)

if routine == 10:
    from libMeanVarCovar import dis_corr as discor
    discor.dis_corr(lammpstrj,Tstep,TskipI,TskipB,T0ct,corr_t,tau)

if routine == 11:
    # python3 main.py --routine 11 --CoulDump Coul_NoBond_pe_atom_spce_woCN.dump --Tstep 50000 --itype 5 --charge -3
    from libMeanVarCovar import MadeL_mean_var
    MadeL_mean_var.doMadelungMeanVar(CoulDump,Tstep,itype,charge)

if routine == 12:
    # python3 main.py --routine 12 --Tstep 10001 --lammpstrj HLseeds_old/newseed1 --itype 5 --jtype 7 --scaling 10 --K_lim 5
    from libMeanVarCovar import bond_angle_check
    bond_angle_check.doBondAnlgeHist(Tstep,lammpstrj,itype,jtype,scaling,K_lim)

if routine == 13:
    # this code is using K_lim as the scaling factor for cosQ
    # python3 main.py --routine 13 --Tstep 10001 --TskipI 50 --TskipB 0 --lammpstrj HLseeds_old/newseed1 --itype 1 --jtype 2 --uv 0,0,-1 --scaling 10 --K_lim 10
    from libMeanVarCovar import HER_angZ_dist
    HER_angZ_dist.doAngZhist(Tstep,TskipI,TskipB,lammpstrj,itype,jtype,uv,scaling,K_lim)
