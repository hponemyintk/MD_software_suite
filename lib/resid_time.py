# formula for claulating total = (T0ct-1)*corr_t+501
# python dis_corr.py --lammpstrj equ_nve_woWater.lammpstrj --Tstep 100 --T0ct 3 --corr_t 10 --tau 70
### doubled checked with excel and passed! Look at 20191001_MSD_FECN6.xlsx
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
    from lib.read_coord import get_wrapped_coord,get_wrapped_coord_Tqxyz
    from lib.math_lib import r_iN
except ImportError as import_error:
    print(import_error)
    print("cannot import the library")

def build_solvation_list(s_pos,redox_pos,boxL,R_lim,coord_num):
    # return 2D array containing solvation ion list for each redox species
    sol_list = np.full((len(redox_pos),coord_num*10),np.inf)
    pos_list = np.full((len(redox_pos),coord_num*10),np.inf)
    for cur_redox in range(len(redox_pos)):
        diff = r_iN(s_pos,redox_pos[cur_redox],boxL)
        ind_list = np.where(diff<R_lim**2)
        for ii in range(len(ind_list[0])):
            sol_list[cur_redox][ii]=ind_list[0][ii]
            pos_list[cur_redox][ii]=np.sqrt(diff[ind_list[0][ii]])
        # print("printing solvation data:::\n",cur_redox,len(ind_list[0]),ind_list)#,sol_list)
    return sol_list,pos_list

# for t* use 2ps and tau = 100ps should be more than enough
def doResid_time(lammpstrj,Tstep,T0ct,corr_t,tau,stype,itype,jtype,R_lim,threshold,coord_num,RefFile):
    if (T0ct-1)*corr_t+tau+1 > Tstep:
        print("The T0ct, corr_t and tau you have chose in larger than Ttep. Choose again!!!")
    else:
        hist = np.zeros(tau+1)                              # store sum of corr values for each t0 which will be averaged over by number of t0 in the end
        for cur_t0 in range(T0ct):
            xyz_in = open(lammpstrj,"r")
            print("***** doing t0#",cur_t0,xyz_in.tell())
            for ii in range(cur_t0*corr_t):                 # skip this many timestep before the new t0
                if RefFile == "0":
                    # print('getting xyz wrapped')
                    get_wrapped_coord(xyz_in)
                else:
                    # print('getting Tqxyz wrapped')
                    get_wrapped_coord_Tqxyz(xyz_in)

            for cur_t in range(tau+1):
                if RefFile == "0":
                    cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord(xyz_in)
                else:
                    cur_pos, tcur_Tstep, Ntotal, boxL = get_wrapped_coord_Tqxyz(xyz_in)

                s_pos=cur_pos[cur_pos[:,0]==stype][:,1:4]
                redox_pos=cur_pos[(cur_pos[:,0]==itype) | (cur_pos[:,0]==jtype)][:,1:4]
                cur_sol_list,pos_list = build_solvation_list(s_pos,redox_pos,boxL,R_lim,coord_num)

                if cur_t == 0:
                    # print("og sol list is ",cur_sol_list,"\npos_list is ",pos_list)
                    og_sol_list = np.copy(cur_sol_list)
                    out_list = np.copy(cur_sol_list)
                    out_list[out_list<np.inf] = 0
                prob_list = np.copy(og_sol_list)
                prob_list[prob_list<np.inf] = 1

                # here compare current sol_list to og and work out the probabilities. Remember the ions has to be in sol_list both at the beginning and end of any given time window to have probability of 1
                for ii in range(len(og_sol_list)):
                    for jj in range(len(og_sol_list[ii, og_sol_list[ii]<np.inf])):
                        if og_sol_list[ii][jj] in cur_sol_list[ii] and out_list[ii][jj]<threshold:
                            # print("here in the code1",ii,jj,og_sol_list[ii][jj], cur_sol_list[ii],out_list[ii][jj],og_sol_list[ii][jj] in cur_sol_list[ii],out_list[ii][jj]<threshold)
                            out_list[ii][jj] = 0
                        elif og_sol_list[ii][jj] not in cur_sol_list[ii]:
                            # print("here in the code2",ii,jj,og_sol_list[ii][jj],cur_sol_list[ii],og_sol_list[ii][jj] not in cur_sol_list[ii],len(og_sol_list[ii, og_sol_list[ii]<np.inf]))
                            prob_list[ii][jj] = 0
                            out_list[ii][jj] += 1
                        if prob_list[ii][jj] == 1 and out_list[ii][jj]>threshold:
                            # print("here in the code3",ii,jj,prob_list[ii][jj],out_list[ii][jj],prob_list[ii][jj] == 1,out_list[ii][jj]>threshold)
                            prob_list[ii][jj] = 0
                # if cur_t%1==0:
                #     print("printing lists")
                #     print(prob_list)
                #     print(out_list)
                #     print(cur_sol_list)
                    # print(pos_list)

                # now histogram the probabilities
                # print(np.sum(prob_list[prob_list==1]),prob_list[prob_list==1],len(prob_list))
                hist[cur_t] += np.sum(prob_list[prob_list==1])/len(prob_list)
        # print("final list",cur_sol_list,"\npos_list is ",pos_list)
        hist /= float(T0ct)

        # write to file
        outfile="resid_time_sijtype"+str(stype)+"_"+str(itype)+"_"+str(jtype)+"R_lim"+str(R_lim)+"Tthreshold"+str(threshold)+"T"+str(Tstep)+"T0ct"+str(T0ct)+"corr_t"+str(corr_t)+"tau"+str(tau)+".out"
        out_fname=open(outfile,"w")
        for i in range(tau+1):
            tmpstr=str(i)+" "+str(hist[i])+"\n"
            out_fname.write(tmpstr)
        out_fname.close()
