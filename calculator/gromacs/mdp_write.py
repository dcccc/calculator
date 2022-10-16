import copy


def get_mdp_line(gromacs_mdp, indoc):
    
    base_key_list = [i[0][:-1].strip() for i in gromacs_mdp if len(i) == 2 ]
    indoc_key = indoc.keys()
    for n, i in enumerate(gromacs_mdp):
        if len(i) == 2 and i[0][:-1].strip() in indoc_key:
            gromacs_mdp[n][1] = indoc[i[0][:-1].strip()]
    
    mdp_str = [[str(j) for j in i ] for i in gromacs_mdp]
    mdp_str = ["  ".join(i) for i in mdp_str]
    
    not_include = ["{:<23s}  =  {}".format(i, indoc[i]) for i in indoc_key 
                                                        if i not in base_key_list]
    
    mdp_str = "\n".join(mdp_str + ["\n"] + not_include)
    
    return(mdp_str)
    
    
    

def write_gromacs_mdp(indoc, mdp_name = "./grompp.mdp"):
    
    gromacs_mdp = copy.deepcopy(gromacs_mdp_base)  
    mdp_line = get_mdp_line(gromacs_mdp, indoc)
    
    mdp_file = open(mdp_name, "w")
    mdp_file.write(mdp_line)
    mdp_file.close()
    

def print_default_mdp(task_type=" ", mdp_para={}):
    if task_type in ["minimization", "npt", "nvt", "free_energy", "pull",\
                                    "anneal","normal_mode"] :
        mdp_para_defalut = copy.deepcopy(eval(task_type))
    elif task_type == "cmd":
        print("No defalut parameter for this task type.")
    elif task_type == "rerun" and mdp_para == {}:
        print("Defalut parameter for this task type is same with the step we do rerun on.\n\n ")
    elif task_type == "rerun" :
        print("Defalut parameter for this task type is same with minimization.\n\n ")
    elif task_type == "replica_exchange":
        print("Defalut parameter for this task type is same with nvt.\n\n")
        mdp_para_defalut = copy.deepcopy(nvt)
    else:
        assert False , print("Unsupported task type.")

    for i in mdp_para.keys():
        mdp_para_defalut[i] = mdp_para[i]
        
    mdp_line = get_mdp_line(gromacs_mdp_base, mdp_para_defalut)

    return(mdp_line)



gromacs_mdp_base =[

["; RUN CONTROL PARAMETERS "                           ],
["integrator               =" ,   "sd"                 ],  
["tinit                    =" ,   0                    ],
["dt                       =" ,   0.0000001            ],
["nsteps                   =" ,   0                    ],
["comm-mode                =" ,   "Linear"             ],
[                                                      ],
["; OUTPUT CONTROL OPTIONS"                            ],       
["nstlog                   =" ,   1000                 ],
["nstenergy                =" ,   1000                 ],
["nstcalcenergy            =" ,   1                    ],
["nstxout                  =" ,   1000                 ],
["nstxout-compressed       =" ,   0                    ],
["compressed-x-precision   =" ,   10000                ],
[                                                      ],
["; NEIGHBORSEARCHING PARAMETERS "                     ],        
["cutoff-scheme            =" ,   "Verlet"             ],
["nstlist                  =" ,   10                   ],
["ns-type                  =" ,   "grid"               ],
["pbc                      =" ,   "xyz"                ],
["rlist                    =" ,   1.2                  ],
[                                                      ],
["; OPTIONS FOR ELE AND VDW"                           ],         
["coulombtype              =" ,  "PME-Switch"          ],
["rcoulomb-switch          =" ,   1.18                 ],
["rcoulomb                 =" ,   1.2                  ],
["vdwtype                  =" ,  "PME"                 ],
["vdw-modifier             =" ,  "Potential-Shift"     ],
["rvdw                     =" ,   1.2                  ],
["DispCorr                 =" ,  "no"                  ],
["fourierspacing           =" ,   0.1                  ],
["pme-order                =" ,   4                    ],
["ewald-rtol               =" ,  "1e-05"               ],
["ewald-geometry           =" ,  "3d"                  ],
[                                                      ],
[                                                      ],
["; Temperature and pressure coupling"                 ],           
["nsttcouple               =" ,   1                    ],
["nstpcouple               =" ,   1                    ],
["tc-grps                  =" ,   "System"             ],
["tcoupl                   =" ,   "no"                 ],
["tau-t                    =" ,   1.0                  ],
["ref-t                    =" ,   300                  ],
["pcoupl                   =" ,   "no"                 ],
["pcoupltype               =" ,   "anisotropic"        ],
["tau-p                    =" ,   3.0                  ],
["ref-p                    =" ,   "1.0 1.0 1.0 1.0 1.0 1.0 "        ],
["compressibility          =" ,   "4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5"  ],
["refcoord-scaling         =" ,   "No"                 ],
[                                                      ],
[                                                      ],
["; SIMULATED ANNEALING  "                             ],
["annealing                =" ,  ""                    ], 
["annealing-npoints        =" ,  ""                    ], 
["annealing-time           =" ,  ""                    ], 
["annealing-temp           =" ,  ""                    ], 
[                                                      ],
["; velocity Generation"                               ],
["gen-vel                  =" ,   "no"                 ],
["gen-temp                 =" ,   300                  ],
["gen-seed                 =" ,   -1                   ],
[                                                      ],
["; OPTIONS FOR BONDS "                                ],
["constraints              =" ,   "none"               ],
["constraint-algorithm     =" ,   "lincs"              ],
["lincs-order              =" ,   12                   ],
["lincs-iter               =" ,   1                    ],
["continuation             =" ,   "no"                 ],                                            
[                                                      ],
["; Free energy variables "                            ],
["free-energy              =" ,   "no"                 ],
["couple-moltype           =" ,   ""                   ],
["couple-lambda0           =" ,   "vdw-q"              ],
["couple-lambda1           =" ,  "vdw-q"               ],
["init-lambda-state        =" ,   0                    ],
["delta-lambda             =" ,   0                    ],
["calc-lambda-neighbors    =" ,   -1                   ],
["nstdhdl                  =" ,   100                  ],
["dhdl-print-energy        =" ,   "total"              ],
]

minimization      = {
"integrator"         :   "steep",
"emtol"              :   2.0,
"emstep"             :   0.01,
"nsteps"             :   80000,
"nstcomm"            :   2,
"nstxout"            :   0,
"nstxout-compressed" :   1000,
"tcoupl"             :   "no",
"gen-vel"            :   "no",
"pcoupl"             :   "no",
"tau-p"              :   10.0,   
}


nvt                = {
"integrator"         :   "sd",
"dt"                 :   0.001   ,
"nsteps"             :   1000000,
"nstxout"            :   0,
"nstxout-compressed" :   10000,
    
"ref-t"              :   300   ,
"tcoupl"             :   "Berendsen",
"tau-t"              :   1.0,     
"gen-vel"            :   "yes",
"pcoupl"             :   "no",
"tau-p"              :   3.0,     
}

npt                = {
"integrator"         :   "sd",
"dt"                 :   0.001   ,
"nsteps"             :   1000000,
"nstxout"            :   0,
"nstxout-compressed" :   10000,    
"ref-t"              :   300   ,
"tcoupl"             :   "Berendsen",
"tau-t"              :   1.0,     
"gen-vel"            :   "yes",
"pcoupl"             :   "Berendsen",
"tau-p"              :   3.0,     
}


anneal             = {  
"integrator"         :   "sd",
"dt"                 :   0.001   ,
"nsteps"             :   2000000,
"nstxout"            :   0,
"nstxout-compressed" :   10000,    
"ref-t"              :   300   ,
"tcoupl"             :   "Berendsen",
"tau-t"              :   1.0,     
"gen-vel"            :   "yes",
"pcoupl"             :   "Berendsen",
"tau-p"              :   3.0,   
    
    
"annealing"          : "yes"    , 
"annealing-npoints"  : "2"   , 
"annealing-temp"     : "300  350"  , 
"annealing-time"     : "0  2000" , 
}

normal_mode = {"integrator" : "nm"}


pull               = {
"integrator"         :   "sd",
"dt"                 :   0.001   ,
"nsteps"             :   1000000,
"nstxout"            :   0,
"nstxout-compressed" :   10000,    
"ref-t"              :   300   ,
"tcoupl"             :   "no",
"tau-t"              :   1.0,     
"gen-vel"            :   "yes",
"pcoupl"             :   "no",
"tau-p"              :   3.0,    

"pull"                  : "yes" ,
"pull-ngroups"          : "2" ,
"pull-ncoords"          : "1" ,
"pull-group1-name"      : "group_A" ,
"pull-group2-name"      : "group_B" ,
"pull-coord1-type"      : "umbrella"  ,     
"pull-coord1-geometry"  : "distance" ,
"pull-coord1-groups"    : "1 2 "  ,
"pull-coord1-dim"       : "N N Y" ,
"pull-coord1-start"     : "yes" ,
"pull-coord1-rate"      : "0.001" ,         
"pull-coord1-k"         : "1000"  ,       
"pull-nstxout"          : "50" ,
"pull-nstfout"          : "50" }


free_energy         = {

"integrator"         :   "sd",
"dt"                 :   0.001   ,
"nsteps"             :   1000000,
"nstxout"            :   0,
"nstxout-compressed" :   10000,    
"ref-t"              :   300   ,
"tcoupl"             :   "no",
"tau-t"              :   1.0,     
"gen-vel"            :   "yes",
"pcoupl"             :   "no",
"tau-p"              :   3.0,    


"free-energy"              : "yes",
"init-lambda-state"        : "0",
"delta-lambda"             : "0",
"calc-lambda-neighbors"    : "1" , 
"vdw-lambdas"              : "0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00",
"coul-lambdas"             : "0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00",
"bonded-lambdas"           : "0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00",
"restraint-lambdas"        : "0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00",
"mass-lambdas"             : "0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00",
"temperature-lambdas"      : "0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00",
                           
"sc-alpha"                 : "0.5",
"sc-coul"                  : "no" ,      
"sc-power"                 : "1"  ,
"sc-sigma"                 : "0.3", 
"couple-moltype"           : "MOL" , 
"couple-lambda0"           : "vdw" ,     
"couple-lambda1"           : "none" ,    
"couple-intramol"          : "no",
"nstdhdl"                  : "10",
}
