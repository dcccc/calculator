import numpy as np
import re,copy,os

def get_pos_restraint_line(para):
    pos_restraint_line = ["[ position_restraints ]",
                          "; atom  type  fx    fy   fz"]
    for i in para:
        pos_restraint_line.append("  ".join([str(j) for j in i])) 
    
    pos_restraint_str = "\n".join(pos_restraint_line)
    return(pos_restraint_str)


def itp_read(itp_file, ff_type="gaff"):
    
    itp_txt = open(itp_file, "r").read()
    
    part_list = ["atomtypes", "moleculetype", "defaults"]
    
    ff_dict = {"atomtypes":[], "mol_name":[],"other_str": [], "ff_type":ff_type}
    
    part_str_list = []
    for part in part_list:        
        tt = re.findall(r"(\[\s+{}\s+\].*?)\[".format(part), itp_txt, flags=re.DOTALL)
        if len(tt) == 1 :
            part_str_list.append(tt[0])
        else:
            part_str_list.append("")
   
    ifdef_num = len(re.findall(r"(#ifndef)", part_str_list[0]))
    endif_num = len(re.findall(r"(#endif)", part_str_list[0]))
    if ifdef_num == endif_num +1 :
        part_str_list[0] = re.findall(r"(\[\s+atomtypes\s+\].*?)#ifndef", part_str_list[0], flags=re.DOTALL)
        
    for n in [0,2]:
        if len(part_str_list[n]) > 0:
            idx = list(re.finditer("\[\s+{}\s+\]".format(part_list[n]), itp_txt))[0]
            idx = idx.start()
            itp_txt = itp_txt[:idx] + itp_txt[idx+len(part_str_list[n]):]
                
    ff_dict["other_str"] = itp_txt
    ff_dict["atomtypes"] = re.sub("\[\s+atomtypes\s+\]","", part_str_list[0])
    ff_dict["atomtypes"] = re.sub(';.*', '' , ff_dict["atomtypes"])
    ff_dict["atomtypes"] = re.sub('\n\n', '\n' , ff_dict["atomtypes"])
    
    mol_name_line = [ i for i in  part_str_list[1].split("\n") if len(i)>3 and i.strip()[0] != ";" ][1]            
    ff_dict["mol_name"] = mol_name_line.strip().split()[0]
    
    return(ff_dict)
    
def get_itp_line(ff_dict, pos_restraint = None): 
        
    if pos_restraint:
        pos_restraint_line = get_pos_restraint_line(pos_restraint)
    else:
        pos_restraint_line = ""

    
    itp_txt = [ff_dict["other_str"],
               pos_restraint_line ]
    
    
    
    itp_txt = "\n\n".join(itp_txt)
    
    return(itp_txt)
    

def write_top(structure, ff_dict_list, pos_restraint = None, write_itp = True):
    
    mol_types = structure.mol_types
    
    mol_list = [mol_types[0]]
    mol_num  = [1]
    for i in mol_types[1:]:
        if i == mol_list[-1]:
            mol_num[-1] = mol_num[-1] + 1
        else:
            mol_list.append(i)
            mol_num.append(1)
            
    
    ff_type = list(set([ff_dict.get("ff_type","gaff") for ff_dict in ff_dict_list]))
    assert len(ff_type) == 1, "Forcefield types should be same."
    ff_type = ff_type[0]

    if ff_type == "gaff" :
        default_line = ["[ defaults ]",
                    "  ;nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
                    "    1        2          yes        0.5      0.8333"]
    elif ff_type == "opls":
        default_line = ["[ defaults ]",
                    "  ;nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
                    "    1        3          yes        0.5      0.5"]
    default_line = "\n".join(default_line)


    atomtypes_line = "[ atomtypes  ]\n;   nr   type  resi   res  atom  cgnr     charge      mass"
    
    include_line = ["; Include itp topology file"]

    mol_restrant_list = list(pos_restraint.keys()) if pos_restraint else []
    for ff_dict in ff_dict_list:
        if pos_restraint and ff_dict["mol_name"] in mol_restrant_list:            
            restraint_para = pos_restraint[ff_dict["mol_name"]]
        else:
            restraint_para = None
        include_line.append('#include "%s.itp" '  %(ff_dict["mol_name"]))
        itp_file = open("%s.itp" %(ff_dict["mol_name"]), "w")
        itp_file.write(get_itp_line(ff_dict, pos_restraint = restraint_para))
        itp_file.close()
        
        atomtypes_line = atomtypes_line + ff_dict["atomtypes"]
    

    include_line = "\n".join(include_line)
    
    system_line = ["[ system ]",
                   "system"     ] 
    
    system_line = "\n".join(system_line)
    
    molecules_line = ["[ molecules ]",
                      "; Compound        nmols", ]
    
        
    for mol_type, mol_num in zip(mol_list, mol_num )  :
        molecules_line.append("  %s              %d" %(mol_type, mol_num ))
    
        
    molecules_line = "\n".join(molecules_line)
    

    top_line = [default_line,
                atomtypes_line,
                include_line,
                system_line,
                molecules_line,

               ]

    top_line = "\n\n".join(top_line)
    
    topology_file = open("topol.top", "w")
    
    topology_file.write(top_line)
    topology_file.close()