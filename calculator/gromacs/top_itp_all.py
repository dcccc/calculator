import numpy as np
import re, copy

# read itp file

def get_atom_para(ff_dict, ff_para):
    
    atom_types   = dict([[i[0],i[-2:]] for i in ff_dict["atomtypes"]])
    atom_type   = [i[1] for i in ff_dict["atoms"]]
    atom_name   = [i[4] for i in ff_dict["atoms"]]
    vdw_para    = np.array([[atom_types[i[1]] for i in ff_dict["atoms"]]], dtype=np.float32).tolist()
    atom_charge = np.array([i[-2] for i in ff_dict["atoms"]], dtype=np.float32).tolist()
    atom_mass   = np.array([i[-1] for i in ff_dict["atoms"]], dtype=np.float32).tolist()
    
    ff_para["vdw_para"]    = vdw_para[0]
    ff_para["atom_name"]   = atom_name
    ff_para["atom_type"]   = atom_type
    ff_para["mass"]        = atom_mass
    ff_para["charge"]      = atom_charge
    ff_para["mol_name"]    = ff_dict["moleculetype"][0][0]
    
    tmp_pairs = ff_dict.get("pairs", [])
    if ff_dict["pairs"] != []:
        tmp_pairs = np.array(ff_dict["pairs"], dtype=np.int)        
        ff_para["pairs"]["type"]   = tmp_pairs[:,2].tolist()
        ff_para["pairs"]["pair"]   = tmp_pairs[:,:2].tolist()
    
    
    exclusions = ff_dict.get("exclusions", [])
    if exclusions != []:
        exclusions = [[int(j)for j in i] for i in exclusions]    
        ff_para["exclusions"] = exclusions
    
    return(ff_para)

def get_bond_para(ff_dict, ff_para):
    
    bond_type = np.array([i[2]  for i in ff_dict["bonds"]], dtype=np.int).tolist()
    bond_pair = np.array([i[:2] for i in ff_dict["bonds"]], dtype=np.int).tolist()
    bond_para = np.array([i[3:] for i in ff_dict["bonds"]], dtype=np.float32).tolist()
    
    ff_para["bonds"]["type"] = bond_type
    ff_para["bonds"]["pair"] = bond_pair
    ff_para["bonds"]["para"] = bond_para
    
    
    return(ff_para)
    

def get_angle_para(ff_dict, ff_para):
    
    if ff_dict.get("angles", False) :
        angle_type = np.array([i[3]  for i in ff_dict["angles"]], dtype=np.int).tolist()
        angle_pair = np.array([i[:3] for i in ff_dict["angles"]], dtype=np.int).tolist()
        angle_para = np.array([i[4:] for i in ff_dict["angles"]], dtype=np.float32).tolist()


        ff_para["angles"]["type"] = angle_type
        ff_para["angles"]["pair"] = angle_pair
        ff_para["angles"]["para"] = angle_para
    
    return(ff_para)

def get_dihe_para(ff_dict, ff_para):
    
    if ff_dict.get("dihedrals", False) :

        tmp_order = ["1", "2", "3", "4", "6"]
        dihe_line = copy.deepcopy(ff_dict["dihedrals"])

        tmp_line = [i[:5] + [i[-1]]  for i in dihe_line if i[4]=="9"]

        for line in tmp_line:
            for i in tmp_order:
                if line[4] == "9" and not line[:5] +[i] in tmp_line :
                    dihe_line.append(line[:5] + ["0.0", "0.0"] +[i])
                    tmp_line.append(line[:5] + [i] )


        dihe_line = sorted(dihe_line, key=lambda x : " ".join(x[:5]+[x[-1]])) 


        dihe_type  = np.array([i[4]    for i in dihe_line if i[4]=="9" ], dtype=np.int).tolist()
        dihe_pair  = np.array([i[:4]   for i in dihe_line if i[4]=="9" ], dtype=np.int).tolist()
        dihe_para  = np.array([i[5:-1] for i in dihe_line if i[4]=="9" ], dtype=np.float32).tolist()
        dihe_order = np.array([i[-1]   for i in dihe_line if i[4]=="9" ], dtype=np.int).tolist()

        impo_type  = np.array([i[4]    for i in dihe_line if i[4]!="9" ], dtype=np.int).tolist()    
        impo_pair  = np.array([i[:4]   for i in dihe_line if i[4]!="9" ], dtype=np.int).tolist()
        impo_para  = np.array([i[5:-1] for i in dihe_line if i[4]!="9" ], dtype=np.float32).tolist()
        impo_order = np.array([i[-1]   for i in dihe_line if i[4]!="9" ], dtype=np.int).tolist() 

    
        ff_para["dihedrals"]["type"]  = dihe_type 
        ff_para["dihedrals"]["pair"]  = dihe_pair 
        ff_para["dihedrals"]["para"]  = dihe_para 
        ff_para["dihedrals"]["order"] = dihe_order

        ff_para["impropers"]["type"]  = impo_type 
        ff_para["impropers"]["pair"]  = impo_pair 
        ff_para["impropers"]["para"]  = impo_para 
        ff_para["impropers"]["order"] = impo_order
    

    return(ff_para)
 

def get_pos_restrain_para(ff_dict, ff_para):
    
    if ff_dict.get("position_restraints", False) :
        pos_restrain_line = copy.deepcopy(ff_dict["position_restraints"])

        atom_id  = [int(i[0]) for i in pos_restrain_line ]
        pos_type = [int(i[1]) for i in pos_restrain_line ]
        pos_para = [[flost(j) for j in i[2:]] for i in pos_restrain_line ]

        ff_dict["position_restraints"]["type"] = pos_type
        ff_dict["position_restraints"]["id"]   = atom_id 
        ff_dict["position_restraints"]["para"] = pos_para
    
    return(ff_dict)

def itp_read(itp_file, ff_type = "gaff"):
    
    itp_txt = open(itp_file, "r").read()
    itp_txt = re.sub(';.*', '' , itp_txt)
    itp_txt = re.sub('[\[\]]', '' , itp_txt)

    itp_line = itp_txt.split("\n")
    
    part_list = ["atomtypes", "moleculetype", "atoms","bonds","pairs","angles",
                 "dihedrals", "exclusions", "" ]
    
    ff_dict = {"atomtypes":[], "moleculetype":[], "atoms":[],     "bonds":[],
               "pairs":[],     "angles":[],       "dihedrals":[],  "exclusions":[], # "position_restraints":[]
              }
    
    addition_part = []
    for line in itp_line:
        line_split = line.split()
        
        if len(line_split) == 1 :
            part = line_split[0].lower()
        elif part not in part_list :
            addition_part.append(line)

        elif len(line_split) > 0 and part in part_list :
            ff_dict[part].append(line_split)
    

    ff_para = {"vdw_para":[], "bonds":{}, "angles":{}, "dihedrals":{},
               "impropers":{}, "pairs":{}, "mol_name":"", "ff_type":"gaff","charge":[],"mass":[],
               "atom_type":[],"atom_name":[], "exclusions":[],
               "position_restraints":[]}
    
    ff_para["addition_part"] = addition_part

    ff_para = get_atom_para(ff_dict,  ff_para)
    ff_para = get_bond_para(ff_dict,  ff_para)
    ff_para = get_angle_para(ff_dict, ff_para)
    ff_para = get_dihe_para(ff_dict,  ff_para)
    #ff_para = get_pos_restrain_para(ff_dict, ff_para)
    ff_para["ff_type"] = ff_type
    
    return(ff_para)
    


# write itp file

def get_moleculetype_line(ff_para):
    
    moleculetype = [ "[ moleculetype ]",
                     " ;name            nrexcl",
                     "   %s              3" %(ff_para["mol_name"])]
    
    moleculetype = "\n".join(moleculetype)
    
    return(moleculetype)

def get_atomtypes_line(ff_para):
    
    atom_type   = ff_para["atom_type"]
    vdw_para    = ff_para["vdw_para"]

    
    txt = ["[ atomtypes ]",
           ";name   bond_type     mass     charge   ptype   sigma         epsilon       Amb"]
    
    
    for n in range(len(atom_type)):   
        tmp_line = " %5s %5s          0.000      0.000    A   %10.7f  %10.7f " %(atom_type[n], atom_type[n],
                                                                                  vdw_para[n][0], vdw_para[n][1])
        if not tmp_line in txt:
            txt.append(tmp_line)
        
    return(txt)
    

def get_atom_line(ff_para):
    
    atom_type   = ff_para["atom_type"]
    atom_name   = ff_para["atom_name"]
    atom_mass   = ff_para["mass"]
    atom_charge = ff_para["charge"]
    mol_title   = ff_para["mol_name"]
    
    txt = ["[ atoms ]",
           ";   nr   type  resi   res  atom  cgnr     charge      mass"]
    
    
    for n in range(len(ff_para["mass"])):        
        

        
        tmp_line = " %5d %5s %5d  %5s%5s %5d   %10.7f  %10.7f " %(n+1,atom_type[n], 1, mol_title, atom_name[n],
                                                                 n, atom_charge[n], atom_mass[n])
        txt.append(tmp_line)
        
    return("\n".join(txt))
    
   
    

def get_bond_line(ff_para):
    
    txt = ["[ bonds ]",
           " ;   ai   aj   funct     r          k"]

    if ff_para.get("bonds", {}) != {}:
        atom_name  = ff_para["atom_name"]
        bond_type  = ff_para["bonds"]["type"]
        bond_pair  = ff_para["bonds"]["pair"]
        bond_para  = ff_para["bonds"]["para"]       
    
        
        for n in range(len(bond_pair)):
            
            tmp_pair   = bond_pair[n]
            bond_atom  = "%s - %s" %(atom_name[tmp_pair[0] - 1], atom_name[tmp_pair[1] - 1])
            tmp_line = "  %4d %4d   %4d    %8.6f   %e  ; %s" %(tmp_pair[0], tmp_pair[1],bond_type[n],
                                                               bond_para[n][0], bond_para[n][1],
                                                               bond_atom)
            txt.append(tmp_line)

    return("\n".join(txt))
  
    

def get_angle_line(ff_para):
    
    
    txt = ["[ angles ]",
           ";   ai   aj  ak    funct    theta          cth"]

    if ff_para.get("angles", {}) != {}:
        atom_name   = ff_para["atom_name"]
        angle_type  = ff_para["angles"]["type"]
        angle_pair  = ff_para["angles"]["pair"]
        angle_para  = ff_para["angles"]["para"]


        for n in range(len(angle_type)):

            tmp_pair   = angle_pair[n]
            angle_atom = " - ".join([atom_name[i - 1] for i in tmp_pair])

            tmp_line = "  %4d %4d %4d %4d    %e   %e  ; %s" %tuple(tmp_pair + [angle_type[n]] +
                                                               angle_para[n] + [angle_atom])
            txt.append(tmp_line)
    return("\n".join(txt))
 
    

def get_dihe_line(ff_para):

    txt = ["[ dihedrals ] ; propers",
           ";    i   j    k    l   func     phase          kd         pn"]
    if ff_para.get("dihedrals", {}) != {}:
        atom_name  =  ff_para["atom_name"]
        dihe_type  =  ff_para["dihedrals"]["type"]  
        dihe_pair  =  ff_para["dihedrals"]["pair"]  
        dihe_para  =  ff_para["dihedrals"]["para"]  
        dihe_order =  ff_para["dihedrals"]["order"] 

        for n in range(len(dihe_order)):

            tmp_pair   = dihe_pair[n]
            dihe_atom = " - ".join([atom_name[i - 1] for i in tmp_pair])

            tmp_line = "  %4d %4d %4d %4d %4d   %e   %e  %d  ; %s" %tuple(tmp_pair + [dihe_type[n]] +
                                                                          dihe_para[n] + [dihe_order[n], dihe_atom])
            txt.append(tmp_line)
    
    return("\n".join(txt))

def get_impo_line(ff_para):
    
    txt = ["[ dihedrals ] ; impropers",
           ";    i   j    k    l   func     phase          kd         pn"]


    if ff_para.get("impropers", {}) != {}:
        atom_name  =  ff_para["atom_name"]
        impo_type  =  ff_para["impropers"]["type"]   
        impo_pair  =  ff_para["impropers"]["pair"]   
        impo_para  =  ff_para["impropers"]["para"]   
        impo_order =  ff_para["impropers"]["order"]  
        
       
        for n in range(len(impo_order)):
            
            tmp_pair   = impo_pair[n]
            impo_atom = " - ".join([atom_name[i - 1] for i in tmp_pair])
            
            tmp_line = "  %4d %4d %4d %4d %4d   %e   %e  %d  ; %s" %tuple(tmp_pair + [impo_type[n]] +
                                                                          impo_para[n] + [impo_order[n], impo_atom])
            txt.append(tmp_line)
    
    return("\n".join(txt))

def get_pairs_line(ff_para):

    txt = ["[ pairs ]",
           ";   ai     aj    funct "]

    if ff_para.get("pairs", {}) != {}:
        atom_name = ff_para["atom_name"]
        pair_type = ff_para["pairs"]["type"]
        pair_pair = ff_para["pairs"]["pair"]
        
    
        
        for n in range(len(pair_type)):
            
            tmp_pair   = pair_pair[n]
            pair_atom  = "%s - %s" %(atom_name[tmp_pair[0] - 1], atom_name[tmp_pair[1] - 1])
            tmp_line = "  %4d  %4d   %4d    ; %s" %(tmp_pair[0], tmp_pair[1], pair_type[n], pair_atom)
            txt.append(tmp_line)

    return("\n".join(txt))

def get_exclusions_line(ff_para):

    txt = ["[ exclusions ]",
          "; i  j"]

    if ff_para.get("exclusions", {}) != {}:
        exclusions = ff_para["exclusions"]
        exclusions = [[str(j)for j in i]for i in exclusions]    



        for i in exclusions:
            tmp_line = "  ".join(i)
            txt.append(tmp_line)
    
    return("\n".join(txt))
    

def get_pos_restraint_line(para):
    pos_restraint_line = ["[ position_restraints ]",
                          "; atom  type  fx    fy   fz"]
    for i in para:
       pos_restraint_line.append("  %d  %d   %f    %f    %f" %tuple(i)) 
    
    pos_restraint_str = "\n".join(pos_restraint_line)
    return(pos_restraint_str)

def get_itp_line(ff_para, pos_restraint = None, write_atomtypes = False):
    
    moleculetype_line = get_moleculetype_line(ff_para)
    atomtypes_line    = get_atomtypes_line(ff_para)    
    atom_line         = get_atom_line(ff_para)
    bond_line         = get_bond_line(ff_para)
    angle_line        = get_angle_line(ff_para)
    dihe_line         = get_dihe_line(ff_para)
    impo_line         = get_impo_line(ff_para)
    pairs_line        = get_pairs_line(ff_para)
    exclusions_line   = get_exclusions_line(ff_para)
    
    

    if pos_restraint:
        pos_restraint_line = get_pos_restraint_line(pos_restraint)
    else:
        pos_restraint_line = ""
    
    ff_type = ff_para.get("ff_type","gaff")
    if ff_type == "gaff" :
        default_line = ["[ defaults ]",
                    "  ;nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
                    "    1        2          yes        0.5      0.8333"]
    elif ff_type == "opls":
        default_line = ["[ defaults ]",
                    "  ;nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
                    "    1        3          yes        0.5      0.5"]
    
    default_line = "\n".join(default_line)
    atomtypes_line = "\n".join(atomtypes_line)

    addition_line = "\n".join(ff_para["addition_aprt"])
                    

    if write_atomtypes:
        itp_txt = [default_line      ,
                   moleculetype_line ,
                   atomtypes_line    ,
                   atom_line         ,
                   bond_line         ,
                   angle_line        ,
                   dihe_line         ,
                   impo_line         ,
                   pairs_line        ,
                   exclusions_line   ,
                   addition_line,     ]
    else:
        itp_txt = [moleculetype_line ,
                   atom_line         ,
                   pos_restraint_line,
                   bond_line         ,
                   angle_line        ,
                   dihe_line         ,
                   impo_line         ,
                   pairs_line        ,
                   exclusions_line   ,
                   addition_line,     ]
        
        
    itp_txt = "\n\n".join(itp_txt)
    
    return(itp_txt)

# write top file

def write_top(structure, ff_para_list, pos_restraint = None, write_itp = True):
    
    mol_types = structure.mol_types
    
    mol_list = [mol_types[0]]
    mol_num  = [1]
    for i in mol_types[1:]:
        if i == mol_list[-1]:
            mol_num[-1] = mol_num[-1] + 1
        else:
            mol_list.append(i)
            mol_num.append(1)
            
    
    ff_type = list(set([ff_para.get("ff_type","gaff") for ff_para in ff_para_list]))
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


    atomtypes_line = ["[ atomtypes  ]",
                      ";   nr   type  resi   res  atom  cgnr     charge      mass"]
    
    include_line = ["; Include itp topology file"]

    mol_restrant_list = list(pos_restraint.keys()) if pos_restraint else []
    for ff_para in ff_para_list:
        if pos_restraint and ff_para["mol_name"] in mol_restrant_list:            
            restraint_para = pos_restraint[ff_para["mol_name"]]
        else:
            restraint_para = None
        include_line.append('#include "%s.itp" '  %(ff_para["mol_name"]))
        itp_file = open("%s.itp" %(ff_para["mol_name"]), "w")
        itp_file.write(get_itp_line(ff_para, pos_restraint = restraint_para, write_atomtypes = False))
        itp_file.close()
        
        atomtypes_line = atomtypes_line + get_atomtypes_line(ff_para)[2:]
    
    atomtypes_line = "\n".join(atomtypes_line)
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
