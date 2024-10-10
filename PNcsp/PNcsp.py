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
import qmpy_rester as qr
import time
from ase.spacegroup import crystal
from ase.io import read, write
import sys
import argparse

def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def get_Symbol(PN):
    with open("./Z_PN_Elem_extended_MD.csv",'r') as file:
        Comp_list={}
        for line in file:
            if 'PN' in line:
                continue
            info=line.strip().replace('\ufeff','').split(',')
            Comp_list[info[0]]=info[1]
    # print(Comp_list)
    return(Comp_list[str(PN)])

def get_PN(Symb):
    with open("./Z_PN_Elem_extended_MD.csv",'r') as file:
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
    print("Separated Formula: ",elems, counts)
    print("PNs: ",PNs)
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
                print("\nWARNING: PN distance between some of constituent elemens are lower than N_neig. This may cause peculiar outputs or redundant operations.\n")

    PN_new_all=[]
    for i in range(len(PNs)):
        PN_new=[]

        for j in np.arange(-1*N_neig,N_neig+1):
            if((PNs[i]+j>0) and (PNs[i]+j<119)):
                PN_new.append(PNs[i]+j)

        PN_new_all.append(PN_new)
    print("PN_new_all: ",PN_new_all)
    
    # Exchange Look-up Table
    exchange_dict={}
    for i in range(len(PN_new_all)):
        for j in range(len(PN_new_all[i])):
            key=get_Symbol(PN_new_all[i][j])
            print(key,":",elems[i])
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

    print("\nOrder_list:", dist_list_ordered)
    print("\nPN  combinations:\n",res_ordered)

    Symbol_list=[]
    for i in range(len(res_ordered)):
        Symbols=[]
        for j in range(len(res_ordered[i])):
            Symbols.append(get_Symbol(res_ordered[i][j]))
        Neig_formula=""    
        for k in range(len(Symbols)):
            Neig_formula+=Symbols[k]+str(counts[k])
        Symbol_list.append(Neig_formula)

    print("\nSymbol list:\n",Symbol_list)

    for i in range(len(dist_list_ordered)):
        if(dist_list_ordered[i]==N_neig):
            res_ordered=res_ordered[i:]
            dist_list_ordered=dist_list_ordered[i:]
            Symbol_list=Symbol_list[i:]
            break

    # print("\nOrder_list:", dist_list_ordered)
    # print("\nPN  combinations:\n",res_ordered)
    # print("\nSymbol list:\n",Symbol_list)

    return Symbol_list,dist_list_ordered,exchange_dict

# OQMD
def get_data_OQMD(Comp_list,neigh_list,Energy_filter,timer):
    All_list=[]
    for i in range(len(Comp_list)):
        with qr.QMPYRester() as q:
            if(Energy_filter=="none"):
                kwargs = {
                    'composition': {Comp_list[i]},
                    'format': 'json',
                    }
            else:
                kwargs = {
                    'composition': {Comp_list[i]},
                    'format': 'json',
                    'delta_e': "<"+str(Energy_filter),
                    }
            list_of_data = q.get_oqmd_phases(False,**kwargs)


            if list_of_data is None:
                print("!!! time exceed !!!")
                print("!!! Wait for a while and use timer with higher value !!!")
                break
            
            if(list_of_data['data']==[]):
                print(Comp_list[i],"--> no structure")
                if(timer!="none"):
                    time.sleep(timer)
                continue
            
            for ind in range(len(list_of_data['data'])):
                list_of_data['data'][ind]['Neigh']=neigh_list[i]
                # list_of_data['data'][ind]['Original']=''.join([i for i in Comp_list[i] if not i.isdigit()])
                # list_of_data['data'][ind]['Original']=Comp_list[i].replace("1","")
                list_of_data['data'][ind]['Original']=Comp_list[i]


            All_list.append(list_of_data)
            print(Comp_list[i])
        if(timer!="none"):
            time.sleep(timer)
    if(All_list==[]):
        print("\nWARNING!")
        print("--------")
        print("** No compound could be found **")
    return All_list

def create_prototype_OQMD(All_list,exchange_dict,formula):
    path="./output_"+formula+"/"
    for num1 in range(len(All_list)):

        neigh=All_list[num1]['data'][0]['Neigh']
        dest_path=path+str(neigh)+"_Neigh/"

        compound=All_list[num1]
        for num2 in range(len(All_list[num1]['data'])):
            name=compound['data'][num2]['name']
            # elems=re.findall('[A-Z][^A-Z]*', name)
            spacegroup=compound['data'][num2]['spacegroup']
            unit_cell=compound['data'][num2]['unit_cell']

            sites=compound['data'][num2]['sites']
            site_list=[]
            elem_list=[]
            for i in range(len(sites)):
                site=tuple(sites[i].split(' @ ')[1].split(' '))
                site_list.append(site)
                elem=sites[i].split(' @ ')[0]
                elem_list.append(elem)

            # IF exchange_dict[elem_list[i]] values are the same, detect it and do realated operations.
            for i in range(len(elem_list)):
                elem_list[i]=exchange_dict[elem_list[i]]

            struct = crystal(elem_list, site_list, cell=unit_cell,size=(1,1,1))
            if not os.path.exists(dest_path):
                os.makedirs(dest_path)
            
            write(dest_path+formula+"_"+name+"_"+spacegroup.replace("/","")+"_"+str(num2)+'.cif',struct)

def categorize(N_neig,formula):
    import shutil  
    import os
    path="./output_"+formula+"/"+str(N_neig)+"_Neigh/"
    # cifs=os.listdir(path)
    cifs=[f for f in os.listdir(path) if ".cif" in f]
    for cif in cifs:
        source_path=path+cif

        sym=cif.split("_")[2]
        dest_path=path+sym+"/"

        if not os.path.exists(dest_path):
            os.makedirs(dest_path)

        dest = shutil.move(source_path, dest_path)  

def show_config(formula,N_neig,E_filter,timer):
    print("\nProgram Configuration")
    print("---------------------")
    print("Query formula:\t",formula,"\nNeighbor order:\t",N_neig,"\nEnergy filter:\t",E_filter,"\nSleep timer:\t",timer)
    print("---------------------\n")

def main():
    parser = argparse.ArgumentParser(prog="PNcsp",description= "PNcsp: A PN similarty based initial structure generator.")
    parser.add_argument('formula')
    parser.add_argument('-n','--neighbor',default=1,help="Order of neighbors to be considered in the similarity search.")
    parser.add_argument('-f','--filter',default=0,help="Selected neighbors are limited to those below the energy filter value. (default: 0) unit: [eV/atom]. Use \"none\" for no filter.") 
    parser.add_argument('-t','--time_sleep',default="none",help="Set sleep time between queries. Excessive number of queries may cause the server to halt.(default: \"none\")")
    args = parser.parse_args()

    formula=args.formula
    N_neig=int(args.neighbor)

    if(args.time_sleep =="none"):
        time_sleep=args.time_sleep
    else:
        time_sleep=int(args.time_sleep)

    if(args.filter =="none"):
        E_filter=args.filter
    else:
        E_filter=int(args.filter)

    show_config(formula=formula,N_neig=N_neig,E_filter=E_filter,timer=time_sleep)
    res,neigh_list,exchange_dict=get_Neig(formula=formula,N_neig=N_neig)
    All_list=get_data_OQMD(res,neigh_list,Energy_filter=E_filter,timer=time_sleep)
    create_prototype_OQMD(All_list,exchange_dict,formula=formula)
        
    print("TERMINATED SUCCESFULLY!")
    categorize(N_neig=N_neig,formula=formula)

if __name__=='__main__':
    main()
