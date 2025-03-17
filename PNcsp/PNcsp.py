#!/usr/bin/env python

###############################################################
#                                                             #
#                          PNcsp                              #
#                                                             #
###############################################################

import numpy as np
import re
import os
import itertools
import time
import sys
import argparse
from db import OQMDoffline
from db import OQMDonline


def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def get_Symbol(PN):
    with open("./db/data/Z_PN_Elem_extended_MD.csv",'r') as file:
        Comp_list={}
        for line in file:
            if 'PN' in line:
                continue
            info=line.strip().replace('\ufeff','').split(',')
            Comp_list[info[0]]=info[1]
    # print(Comp_list)
    return(Comp_list[str(PN)])

def get_PN(Symb):
    with open("./db/data/Z_PN_Elem_extended_MD.csv",'r') as file:
        Comp_list={}
        for line in file:
            if 'PN' in line:
                continue
            info=line.strip().replace('\ufeff','').split(',')
            Comp_list[info[1]]=int(info[0])
    # print(Comp_list)
    return(Comp_list[Symb])

def separate_elems(formula):
  """Separates a chemical formula using regular expressions"""
  matches = re.findall(r'[A-Z][a-z]*|\d+', formula)
  elems = [match for match in matches if match.isalpha()]
  counts = [int(match) for match in matches if match.isdigit()]
  return elems, counts

def convert_formula(formula):
    elems, counts=separate_elems(formula)
    PNs=[]
    for Symb in elems:
        PNs.append(get_PN(Symb))
    print("Separated Formula:",elems)
    print("PNs:",PNs)
    print("Ratios:",counts)
    return elems, counts, PNs

def dist_classifier(PNs,res):
    dists=[]
    for i in range(len(res)):
        dist_local=[]
        for j in range(len(PNs)):
            # for k in range(len(res[i])):
            dist_local.append(abs(PNs[j]-res[i][j]))
        dist_max=max(dist_local)
        dists.append(dist_max)
    return dists

def get_Neig(formula, N_neig):
    elems, counts, PNs=convert_formula(formula)

    # In case constituent elemens are too close 
    for i in range(len(PNs)-1):
        for j in range(i+1,len(PNs)):
            if(abs(PNs[i]-PNs[j])<=N_neig):
                # print("\nWARNING: PN distance between some of constituent elemens are lower than N_neig. This may cause peculiar or incomplete outputs.\n")
                print("\nWARNING: PN distance between",elems[i],"and",elems[j] ,"are lower than N_neig. This may cause peculiar or incomplete outputs.")
    print("\n")

    PN_new_all=[]
    for i in range(len(PNs)):
        PN_new=[]
        for j in np.arange(-1*N_neig,N_neig+1):
            if((PNs[i]+j>0) and (PNs[i]+j<119)):
                PN_new.append(PNs[i]+j)

        PN_new_all.append(PN_new)
    print("PN_new_all: ",PN_new_all)
    
    # PN_new_all=[]
    # for i in range(len(PNs)):
    #     PN_new=[]
    #     for j in np.arange(-1*N_neig,N_neig+1):
    #         if((PNs[i]+j>0) and (PNs[i]+j<119)):
    #             ind_dup=0
    #             for k in range(len(PNs)):
    #                 if(k==i):
    #                     continue
    #                 if(PNs[i]+j == PNs[k]):
    #                     ind_dup=1
    #                     break
    #             if(ind_dup==0):
    #                 PN_new.append(PNs[i]+j)
    #     PN_new_all.append(PN_new)
    # print("Neigboring PNs: ",PN_new_all)
    
    # Exchange Look-up Table
    exchange_dict={}
    for i in range(len(PN_new_all)):
        for j in range(len(PN_new_all[i])):
            key=get_Symbol(PN_new_all[i][j])
            # print(key,":",elems[i])
            exchange_dict[key]=elems[i]
    print("exchange_dict: ",exchange_dict)

    res=np.array(list(itertools.product(*PN_new_all)))
    
    # # Drop original formula
    # for i in range(len(res)):
    #     if(np.array_equal(res[i], PNs)):
    #         res=np.concatenate((res[:i],res[i+1:]), axis=0)
    #         break

    drop_list=[]
    for i in range(len(res)):
        # Drop original formula
        if(np.array_equal(sorted(res[i]), sorted(PNs))):
            drop_list.append(i)
        # Drop elemental dublication
        if(len(set(res[i]))<len(res[i])):
            drop_list.append(i)

    res=np.delete(res, drop_list,axis=0)

    # Distance Detector
    dist_list=dist_classifier(PNs,res)

    # New order
    res_ordered=[]
    dist_list_ordered=[]
    for dst in range(1,N_neig+1):
        for i in range(len(dist_list)):
            if (dist_list[i]==dst):
                res_ordered.append(list(res[i]))
                dist_list_ordered.append(dist_list[i])

    # print("\nOrder_list:", dist_list_ordered)
    # print("\nPN  combinations:\n",res_ordered)

    Symbol_list=[]
    for i in range(len(res_ordered)):
        Symbols=[]
        for j in range(len(res_ordered[i])):
            Symbols.append(get_Symbol(res_ordered[i][j]))
        Neig_formula=""    
        for k in range(len(Symbols)):
            Neig_formula+=Symbols[k]+str(counts[k])
        Symbol_list.append(Neig_formula)

    for i in range(len(dist_list_ordered)):
        if(dist_list_ordered[i]==N_neig):
            res_ordered=res_ordered[i:]
            dist_list_ordered=dist_list_ordered[i:]
            Symbol_list=Symbol_list[i:]
            break

    # print("\nOrder_list:", dist_list_ordered)
    print("\nPN  combinations and Symbol list")
    for i in range(len(res_ordered)):
        print(res_ordered[i],Symbol_list[i])
    return Symbol_list,dist_list_ordered,exchange_dict

def categorize(N_neig,formula):
    import shutil  
    import os
    path="./output_"+formula+"/"+str(N_neig)+"_Neigh/"
    # cifs=os.listdir(path)
    cifs=[f for f in os.listdir(path) if ".cif" in f]
    print("Total Number of CIF files:",len(cifs))
    for cif in cifs:
        source_path=path+cif

        sym=cif.split("_")[2]
        dest_path=path+sym+"/"

        if not os.path.exists(dest_path):
            os.makedirs(dest_path)

        dest = shutil.move(source_path, dest_path)  

def show_config(formula,N_neig,E_filter,timer,online):
    print("\nProgram Configuration")
    print("---------------------")
    print("Query formula:\t",formula,"\nNeighbor order:\t",N_neig,"\nEnergy filter:\t",E_filter,"\nSleep timer:\t",timer,"\nOnline: \t",online)
    print("---------------------\n")

def main():
    parser = argparse.ArgumentParser(prog="PNcsp",description= "PNcsp: A PN similarty based initial structure generator.")
    parser.add_argument('formula')
    parser.add_argument('-n','--neighbor',default=1,help="Order of neighbors to be considered in the similarity search.")
    parser.add_argument('-f','--filter',default=0,help="Selected neighbors are limited to those below the energy filter value. (default: 0) unit: [eV/atom]. Use \"none\" for no filter.") 
    parser.add_argument('-t','--time_sleep',default="none",help="Set sleep time between queries. Excessive number of queries may cause the server to halt.(default: \"none\")")
    parser.add_argument('-o','--online',default=False,help="Sets online (True) or offline (False) search in OQMD (default: False). For offline seach, you should download and set up offline OQMD database. See https://oqmd.org/download/.")
    args = parser.parse_args()

    formula=args.formula
    N_neig=int(args.neighbor)
    online=bool(args.online)

    if(args.time_sleep =="none"):
        time_sleep=args.time_sleep
    else:
        time_sleep=int(args.time_sleep)

    if(args.filter =="none"):
        E_filter=args.filter
    else:
        E_filter=int(args.filter)

    show_config(formula=formula,N_neig=N_neig,E_filter=E_filter,timer=time_sleep,online=online)
    res,neigh_list,exchange_dict=get_Neig(formula=formula,N_neig=N_neig)

    if(online==True):
        All_list=OQMDonline.get_data_OQMD(res,neigh_list,Energy_filter=E_filter,timer=time_sleep)
        OQMDonline.create_prototype_OQMD(All_list,exchange_dict,formula=formula)
    else:
        All_list=OQMDoffline.get_data_OQMD(res,neigh_list,Energy_filter=E_filter,timer=time_sleep)
        OQMDoffline.create_prototype_OQMD(All_list,exchange_dict,formula=formula)
        
    print("TERMINATED SUCCESFULLY!")
    categorize(N_neig=N_neig,formula=formula)

if __name__=='__main__':
    main()
