import sys
import numpy as np
import re
import itertools
import qmpy_rester as qr
import time
import os
from ase.spacegroup import crystal
from ase.io import read, write

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
    print("****")
    print("Separated Formula: ",elems, counts)
    print("PNs: ",PNs)
    print("****")
    return elems, counts, PNs

def get_Neig(compound, N_neig=3):
    elems, counts, PNs=convert_formula(compound)
    PN_new_all=[]
    for i in range(len(PNs)):
        PN_new=[]
        for j in np.arange(-1*N_neig,N_neig+1):
            PN_new.append(PNs[i]+j)

        PN_new_all.append(PN_new)
    # print(PN_new_all)
    res=np.array(list(itertools.product(*PN_new_all)))
    for i in range(len(res)):
        if(np.array_equal(res[i], PNs)):
            res=np.concatenate((res[:i],res[i+1:]), axis=0)
            break
    Symbol_list=[]
    for i in range(len(res)):
        Symbols=[]
        for j in range(len(res[i])):
            Symbols.append(get_Symbol(res[i][j]))
        Neig_formula=""    
        for k in range(len(Symbols)):
            Neig_formula+=Symbols[k]+str(counts[k])
        Symbol_list.append(Neig_formula)
    return Symbol_list

def get_Structure(Comp_list,Energy_filter="<0"):
    All_list=[]
    print("**** Detected Neighbors ****")
    for i in range(4):
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
                    'delta_e': Energy_filter,
                    }
            list_of_data = q.get_oqmd_phases(False,**kwargs)
            if(list_of_data['data']==[]):
                print(Comp_list[i],"--> no structure")
                continue

            All_list.append(list_of_data)
            print(Comp_list[i])
    if(All_list==[]):
        print("\nWARNING!")
        print("--------")
        print("** No compound could be found **")
    return All_list

def create_Prototype(All_list,formula,dest_path):
    for num1 in range(len(All_list)):
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

            print(name,elem_list,site_list,spacegroup,unit_cell)

            elem_unique=list(set(elem_list))
            compound_elems, counts=separate_elems(formula)
            
            for i in range(len(elem_unique)):
                for j in range(len(elem_list)):
                    if(elem_list[j]==elem_unique[i]):
                        elem_list[j]=compound_elems[i]

            struct = crystal(elem_list, site_list, cell=unit_cell,size=(1,1,1))
            if not os.path.exists(dest_path):
                os.makedirs(dest_path)
            write(dest_path+name+"_sym"+spacegroup.replace("/","")+"_"+str(num2)+'.cif',struct)

def main():
    formula=sys.argv[1]
    print(sys.argv)
    
    if(len(sys.argv)==3):
        N_neig=int(sys.argv[2])
        res=get_Neig(formula,N_neig)
        All_list=get_Structure(res)
        create_Prototype(All_list,formula,"./output/")
    elif(len(sys.argv)==4):
        N_neig=int(sys.argv[2])
        Energy_filter=sys.argv[3]
        res=get_Neig(formula,N_neig)
        All_list=get_Structure(res,Energy_filter)
        create_Prototype(All_list,formula,"./output/")
    else:
        res=get_Neig(formula)
        All_list=get_Structure(res)
        create_Prototype(All_list,formula,"./output/")

if __name__ == "__main__":
    main()
    print("SUCCESS")