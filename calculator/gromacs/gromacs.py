import os, copy, subprocess,glob
import numpy as np

from calculator.gromacs.mdp_write import *
from calculator.gromacs.pdb_gro import *
from calculator.gromacs.top_itp import *

def read_xvg(xvg_file, sel_list = []):
    
    sel_data_dict = {}
    xvg_data = np.loadtxt(xvg_file, comments=["@", "#"]).T.tolist()
    
    if "energy" in xvg_file:
    
        xvg_data = np.loadtxt(xvg_file, comments=["@", "#"]).T.tolist()
        item_list  = re.findall(r'@ s\d+ legend "(.*)"\n', open(xvg_file, "r").read())
        item_list  = ["Time"] + list(item_list)
        
        xvg_data_dict = dict(zip(item_list, xvg_data))
        
        tmp_list = [item for item in sel_list if item not in item_list]
        assert tmp_list == [], print(f"energy items {tmp_list} not in the energy! ")
        if sel_list != []:
            for i in sel_list:
                sel_data_dict[i] = xvg_data_dict[i]
        else:
            sel_data_dict = xvg_data_dict
    elif "rmsd" in xvg_file:
        item_list  = ["Time", "RMSD"]
        sel_data_dict = dict(zip(item_list, xvg_data))   
    else:
        # item_list  = ["Time", xvg_file[:-4]]
        sel_data_dict = xvg_data   
    
    
    return(sel_data_dict)

def write_index_file(st, index_doc, index_file = "index.ndx" ):
    
    all_atom_idx = np.arange(1, st.atom_num+1)
    index_str = ["[ system ]",  "\n".join(map(str, all_atom_idx)) ]
    for group in index_doc.keys():
        group_atom_idx_str = []
        for tp in index_doc[group].keys():
            atom_idx = []
            if tp in ["atom", "mol" ] :
                mol_idx1 = np.array([ i  for i in index_doc[group][tp] if isinstance(i, int)])
                mol_idx2 = [i.split("-") for i in index_doc[group][tp] if isinstance(i, str) and "-" in i]
                mol_idx2 = [np.arange(int(i[0]),int(i[1]) + 1 ) for i in mol_idx2]
                mol_idx  = np.hstack(mol_idx2 + [mol_idx1])
             
                
                if tp == "mol":
                    atom_idx = []
                    tmp_mask = all_atom_idx < 0
                    for i in mol_idx:
                        tmp_mask = np.logical_or(tmp_mask, st.mol_idx == i)
                    tmp_atom_idx =  all_atom_idx[tmp_mask]
                else:
                    tmp_atom_idx = mol_idx 
            elif tp in ["mol_name"]:
                tmp_mask = all_atom_idx < 0
                for mol_name in index_doc[group][tp]:
                    tmp_mask = np.logical_or(tmp_mask, st.mol_titles == i)
                    
                tmp_atom_idx =  all_atom_idx[tmp_mask]
            atom_idx.extend(tmp_atom_idx)

        atom_idx = "\n".join(map(str, sorted(atom_idx)))        
        index_str += [ f'[{group}]',atom_idx,"\n" ]
    
    index_str = "\n".join(index_str)
    index_file = open(index_file, "w")
    index_file.write(index_str)
    index_file.close()
    

def run_cmd(cmd，task_name=""):

    if task_name =="":
        task_name = "task"

    cmd = cmd + f" 1>{task_name}_run.out  2>{task_name}_run.err"
    status, stdout = subprocess.getstatusoutput(cmd)
    stderr_str = open(f"./{task_name}_run.err", "r").read()            
    assert status == 0, stderr_str


class Gromacs():
    def __init__(self, gmx_exe={}):
        self.gmx_exe = gmx_exe       
        self.gmx_exe_sp = self.gmx_exe.get("gmx", None)
        self.openmpi = self.gmx_exe.get("openmpi", None)

        # set the environment variables for openmpi
        if self.openmpi:
            openmpi_root = os.path.dirname(os.path.dirname(orca_exe["openmpi"]))
            os.environ["LD_LIBRARY_PATH"] = os.environ.get("LD_LIBRARY_PATH", "") + ":" + openmpi_root + "/lib"
            os.environ["OPAL_PREFIX"]     = openmpi_root
            os.environ["PATH"]            = os.environ.get("PATH", "") + ":" + openmpi_root + "/bin"
        

 
    def gen_tpr(self, mdp_para, structure, ff_list, replica_exchange = None, pos_restraint = None, 
                index_para = None, conf_name = "conf.gro", ref_st = None ):
        
        # get para_list
        para_list = [[i, mdp_para[i]] for i in mdp_para.keys() if isinstance(mdp_para[i], list) ]
        para_list_len = [len(i[1]) for i in para_list] 

        if replica_exchange and isinstance(structure, list):
            print("Error! Variable replica_exchange is a temperature list, and variable struture is a list of structures")
            exit()
            for n, i in enumerate(replica_exchange):
                os.mkdir("./remd{}".format(n))
                os.chidr("./remd{}".format(n))
                mdp_para["ref-t"] = i


                self.gen_tpr(mdp_para, structure[n], ff_list, pos_restraint = pos_restraint,
                             index_para = index_para, conf_name = "conf.gro" )
                os.chdir("../")   

        elif (replica_exchange == None and (isinstance(structure, list)) or para_list != []):
            
            if isinstance(structure, list) and para_list != []:
                assert len(set(para_list_len)) == 1, "Length of para_list should be the same"
                assert para_list_len[0] == len(structure), "Length of para_list and length of structures should be the same"
            elif isinstance(structure, list) and para_list == []:
                para_list = []
            elif not isinstance(structure, list) and para_list != []:
                structure = [structure] * para_list_len[0]

            if isinstance(index_para, dict):
                index_para = [index_para] * len(structure)
            elif isinstance(index_para, list):
                assert len(structure) == len(index_para), "Number of structure and length of index_para should be the same"
            else:
                index_para = [None]*len(structure)
            
            if not ref_st is None:
                if isinstance(ref_st, list) and isinstance(structure, list):
                    assert len(structure) == len(ref_st), "Length of ref_st and length of structures should be the same"
                else:
                    ref_st = [ref_st]*len(structure)
            else:
                ref_st = [None]*len(structure)

            
            for n in range(len(structure)):
                tmp_mdp_para = copy.deepcopy(mdp_para)
                os.mkdir("./md_{}".format(n))
                os.chdir("./md_{}".format(n))
                for para in para_list:
                    tmp_mdp_para[para[0]] = para[1][n]
                
                self.gen_tpr(tmp_mdp_para, structure[n], ff_list, pos_restraint = pos_restraint,
                             index_para = index_para[n], conf_name = conf_name, ref_st = ref_st[n])
                os.chdir("../")  
        else:
        
            # write index file
            index_file = "./index.ndx"
            if index_para:
                write_index_file(structure, index_para, index_file = "index.ndx" )
                index_file = os.path.abspath("./index.ndx")
            else:
                index_file = None

            if ref_st is None:
                ref_st = conf_name


            write_gromacs_mdp(mdp_para)
            structure.write(f"./{conf_name}")
            topol_str = write_top(structure, ff_list, pos_restraint = pos_restraint, 
                                  write_itp = True)
        
            cmd = f"{self.gmx_exe_sp} grompp -f grompp.mdp  -c {conf_name}  -p topol.top  -o md.tpr -maxwarn 10"
            
            if pos_restraint :
                cmd = cmd + f" -r {ref_st}"
            if index_file:
                cmd = cmd + f" -n {index_file}"
            
            run_cmd(cmd, task_name="gmx_grompp")

    def energy(self, edr_file = "md.edr", item = None, begin = 0 , end = None):

        if item:
            item_cmd  =  "  ".join(item)
        else:
            item_cmd = " ".join([str(i) for i in range(100)])
        edr_dir = os.path.dirname(os.path.abspath(edr_file))
        cmd = f" echo '{item_cmd}' |  {self.gmx_exe_sp}  energy -f  {edr_file} -o {edr_dir}/energy.xvg"
        
        if begin:
            cmd = cmd + f" -b {begin}"
        if end:
            cmd = cmd + f" -e {end}  "
        
        
        run_cmd(cmd, task_name="gmx_energy")

    def trjconv(self, traj_file_in = "md.xtc", tpr_file = "md.tpr", out_file ="md1.xtc",
                group = 0, begin = None, end = None, skip = None, pbc = None,
                ndec = None, dump = None, index_file = None):   
        
        cmd = f" echo {group} |  {self.gmx_exe_sp}  trjconv -f  {traj_file_in} -s {tpr_file} -o {out_file}"
        if index_file:
            index_cmd = f" -n {index_file} "
            cmd = cmd + index_cmd        
        if begin:
            cmd = cmd + f" -b    {begin}"
        if end:
            cmd = cmd + f" -e    {end}   "
        if skip:
            cmd = cmd + f" -skip {skip}  "
        if pbc:
            cmd = cmd + f" -pbc  {pbc}   "
        if ndec:
            cmd = cmd + f" -ndec {ndec}  "
        if dump:
            cmd = cmd + f" -dump {dump}  "

            
        run_cmd(cmd, task_name="gmx_trjconv")

    def rmsd(self, group = 0, index_file = None):
        
        # fix broken molecule by converting the xtc trajctory file
        cmd = f" echo 0 |  {self.gmx_exe_sp}  trjconv -f  md.xtc -s md.tpr -o rmsd.xtc -pbc whole"
        status, stdout = subprocess.getstatusoutput(cmd)
        assert status == 0, print(stdout)
        
        
        group_cmd = " {0} {0} ".format(group)

        cmd = f" echo {group_cmd} |  {self.gmx_exe_sp}  rms -f  rmsd.xtc -s md.tpr -o rmsd.xvg"
        if index_file:
            index_cmd = f" -n {index_file} "
            cmd = cmd + index_cmd
        
        run_cmd(cmd, task_name="gmx_rms")
        
        # delete rmsd.xtc
        os.remove("./rmsd.xtc")
        
        
        rmsd = np.loadtxt("rmsd.xvg", comments=["@", "#"]).T.tolist()        
        return(rmsd)

    def run(self, dir_path = "./", task_type="minimization",  run_type="ntmpi",
            replica_exchange = None, core_num=1, omp_num = 1, mpi_num=1, 
            precision = "single", continue_run = False, append = False, replex = 400, 
            gpu = None, gpu_id = None):
       
        pwd = os.getcwd()
        os.chdir(dir_path)      
       
        
        if precision == "single":
            gmx_exe = self.gmx_exe["gmx"]   
        else:
            gmx_exe = self.gmx_exe["gmx_d"] 
        
        if gpu and "openmpi" not in run_type:
            gmx_exe = self.gmx_exe["gmx_gpu"]
        
        
        assert not (gpu and precision == "double") , " gpu and double is not compatible"

        run_type_list = ["ntmpi", "ntomp", "ntmpi_ntomp", "openmpi", "openmpi_ntomp"]

        assert run_type in run_type_list , print("run_type is error")
        if   run_type == "ntmpi":
            os.environ["OMP_NUM_THREADS"] = str(1)
            cmd = f"{gmx_exe}  mdrun -ntmpi {core_num} -ntomp 1 "
        elif run_type == "ntomp":
            os.environ["OMP_NUM_THREADS"] = str(core_num)
            cmd = f"{gmx_exe}  mdrun -ntmpi  1 -ntomp {core_num}"    
        elif run_type == "ntmpi_ntomp" :
            os.environ["OMP_NUM_THREADS"] = str(omp_num)
            cmd = f"{gmx_exe}  mdrun -ntmpi  {mpi_num} -ntomp {omp_num}"
       
        elif "openmpi" in run_type : 
            mpirun = self.openmpi["mpirun"]
            
            if presion == "single":
                gmx_exe = self.gmx_exe["gmx_mpi"]
            else:
                gmx_exe = self.gmx_exe["gmx_mpi_d"]

            if gpu and "openmpi" in run_type:
                gmx_exe = self.gmx_exe["gmx_gpu_mpi"]
                
            if run_type == "openmpi" : 
                os.environ["OMP_NUM_THREADS"] = str(1)
                cmd = f"{mpirun} -np {core_num}  {gmx_exe} mdrun -ntomp 1"
            elif run_type == "openmpi_ntomp" : 
                os.environ["OMP_NUM_THREADS"] = str(omp_num)
                cmd = f"{mpirun} -np {mpi_num}  {gmx_exe} mdrun -ntomp {omp_num}"

 
       
        if task_type in ["minimization", "nvt", "npt", "anneal"]:
            cmd = cmd + " -deffnm md"
        elif task_type == "pull":
            cmd = cmd + " -deffnm md -px md_pullx.xvg -pf md_pullf.xvg"

        elif replica_exchange and task_type in ["remd", "replica_exchange"]:
            dir_list = ["remd{}".format(i) for i in range(len(replica_exchange))]
            dir_str = " ".join(dir_list)
            cmd = cmd + f" -multidir  {dir_str} -replex {replex} "
        elif task_type == "normal_mode":
            cmd = cmd + " -deffnm md -mtx nm.mtx"
        elif task_type == "free_energy":
            cmd = cmd +  " -deffnm md  -dhdl dhdl.xvg "
        elif task_type == "rerun" :
            cmd = cmd +  " -s md.tpr -rerun md.xtc -e md.edr "
        else:
            assert False, f"{task_type} is not supported! Exit."


        if gpu and gpu_id is None:
            cmd = cmd + " -gpu_id 0"
        elif gpu and gpu_id :
            cmd = cmd + f" -gpu_id {gpu_id}"
        
        if continue_run:
            if os.path.exists("./md.cpt"):
                cmd = cmd + " -cpi md.cpt "
                
        if append:
            cmd = cmd + " -append "
       
        cmd = cmd + " 1>stdout  2>stderr"
        status, stdout = subprocess.getstatusoutput(cmd)
        stderr_str = open("./stderr", "r").read()
        
        if status != 0 and run_type in ["openmpi", "ntmpi"]  and "There is no domain decomposition" in stderr_str:
            if task_type != "remd":
                os.environ["OMP_NUM_THREADS"] = str(core_num)                
                cmd = f"{gmx_exe}  mdrun -ntmpi  1 -ntomp {core_num}"
            else:
                omp_num = int(core_num / len(replica_exchange))
                os.environ["OMP_NUM_THREADS"] = str(omp_num)
                replica_exchange_num = len(replica_exchange)
                cmd = f"{mpirun} -np {replica_exchange_num} {gmx_exe} mdrun -ntomp {omp_num}"
                cmd = cmd + f" -multidir  {dir_str} -replex {replex} "
        
                run_cmd(cmd, task_name="gmx_md")
        
    
    def nmeig(self, run_type = "ntmpi", last = 10000):
        
        if run_type == "ntmpi":
            gmx_exe = self.gmx_exe["gmx_d"]
        else:
            gmx_exe = self.gmx_exe["gmx_mpi_d"] 
            
        cmd = f"{gmx_exe} nmeig -f nm.mtx -s md.tpr -last {last} "            
        
        run_cmd(cmd, task_name="gmx_nmeig")

        # eigenfreq = np.loadtxt("./eigenfreq.xvg", comments=["@", "#"]).T[1]
        # eigenval  = np.loadtxt("./eigenval.xvg",  comments=["@", "#"]).T[1]
        

        os.remove("./nm.mtx")
        os.remove("./eigenvec.trr")
        
        # return(eigenfreq, eigenval)
    
    def gmx_bar(self, xvg_file="md_*/dhdl.xvg", b = None, e = None, nbin = None, temp = None ) :
        
        xvg_file = "   ".join(glob.glob("md_*/dhdl.xvg"))

        cmd = f"{self.gmx_exe_sp} bar -f {xvg_file} -o -oi -oh" 

        for i in ["b" , "e" , "temp" , "nbin"  ]:
            cmd = cmd + eval(" f' -{i}  %s ' " ) %(eval(i)) if eval(f"{i}") else cmd

        run_cmd(cmd, task_name="gmx_bar")


    def gmx_wham(self, b = None, bins = None, temp = None, tol = None, max = None, min = None) :
        
        tpr_str = "\n".join(glob.glob("md_*/md.tpr"))
        pullf_str = "\n".join(glob.glob("md_*/md_pullf.xvg"))
        tpr_file = open("tpr.dat", "w")
        tpr_file.write(tpr_str)
        tpr_file.close()

        pullf_file = open("pullf.dat", "w")
        pullf_file.write(pullf_str)
        pullf_file.close()


        cmd = f"{self.gmx_exe_sp} wham -if pullf.dat -it tpr.dat  -o -hist  " 


        for i in ["b" , "bins" , "temp" , "tol" , "max" ,"min" ]:
            cmd = cmd + eval(" f' -{i}  %s ' " ) %(eval(i)) if eval(f"{i}") else cmd

        run_cmd(cmd, task_name="gmx_wham")
        
    
    def convert_tpr(self, tpr_file = "./md.tpr", extend = False, time = 0, index_file = None, group=None):
        
        
        cmd = f" {self.gmx_exe_sp}  convert-tpr -s {tpr_file} -o md.tpr "
        
        if extend:
            extend_cmd = f" -extend {time} "
            cmd = cmd + extend_cmd      
        
        if index_file:
            cmd = f"echo {group} | {cmd} -n {index_file} "
        
        run_cmd(cmd, task_name="gmx_convert_tpr")
        
    
    def extend_run(self, dir_path = "./", time = 100, task_type="md", 
                  replica_exchange = None, run_type="ntmpi", core_num=1, omp_num = 1, mpi_num=1, 
                  presion = "single", continue_run = False, append = False, replex=100,
                  gpu = None, gpu_id = None):

        os.chdir(dir_path)
        
        if replica_exchange:
            for n,i in enumerate(replica_exchange):
                os.chdir("{}/remd{}".format(pwd, n))
                self.convert_tpr(extend = True, time = time, index = index)
                
        else:
            self.convert_tpr(extend = True, time = time,index = index)

        
        self.run(task_type="md", replica_exchange = replica_exchange, 
                 run_type=run_type, core_num = core_num, omp_num = omp_num, mpi_num = mpi_num, 
                 presion = presion, continue_run = continue_run, append = append, replex = replex,
                 gpu = None, gpu_id = None)
         
    def rerun(self, traj_file = None, tpr_file = None, core_num = 1):
        
        assert traj_file != None, "Trajctory file is not setted"
        assert tpr_file  != None, "Tpr file is not setted"

        cmd = f"{self.gmx_exe_sp}  mdrun -nt {core_num} -s {tpr_file} -rerun {traj_file} -e md.edr"

        run_cmd(cmd, task_name="gmx_rerun")
                
