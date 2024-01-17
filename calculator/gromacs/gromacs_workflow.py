import os, copy, subprocess, glob,sys,json
import numpy as np

from calculator.gromacs.mdp_write import *
from calculator.gromacs.pdb_gro   import *
from calculator.gromacs.top_itp   import *
from calculator.gromacs.gromacs   import *
from calculator.util.logger       import *

def read_st(st_file):
    if "\n" not in st_file < 256 and os.path.exists(st_file):
        structure = Structure.read(st_file)
    elif isinstance(st_file, str)  :
        if "POSITION" in st_file:
            structure = Structure.read_g96_str(st_file)
        else:
            structure = Structure.read_gro_str(st_file)
    return(structure)



class Gromacs_workflow():
    def __init__(self, workflow_json = {}, gmx_exe = {}):
        self.workflow_json = workflow_json
        self.task_step_list  =  np.array(workflow_json.get("task_step", []))

        self.supported_task_type = ["minimization", "npt", "nvt", "replica_exchange", "free_energy", "pull",
                                    "anneal","normal_mode", "rerun","cmd"]

        if "openmpi_lib" in gmx_exe.keys():
            os.environ["LD_LIBRARY_PATH"] = os.environ.get("LD_LIBRARY_PATH", "") + ":" + gmx_exe["openmpi_lib"]
            os.environ["OPAL_PREFIX"] = os.path.dirmane(gmx_exe["openmpi_lib"]) 
        elif "openmpi" in orca_exe.keys():
            os.environ["LD_LIBRARY_PATH"] = os.environ.get("LD_LIBRARY_PATH", "") + ":" + orca_exe["openmpi"] + "/lib"
            os.environ["OPAL_PREFIX"]     = os.path.dirmane(orca_exe["openmpi"])  + "/lib"

        self.gmx = Gromacs(gmx_exe = gmx_exe)
        self.base_path = os.getcwd()
        self.logger = Logger()


    def task_run(self):

        if len(self.task_step_list) == 0:
            pass 
        else:
            # workflow steps task list and step task type list
            self.task_type_list = [self.workflow_json[task_step].get("task_type", None) for task_step in self.task_step_list]
            assert len(self.task_step_list) == len(set(self.task_step_list)), print("Task_step duplicated ! ")
            None_type_idx       = [n for n, i in enumerate(self.task_type_list) if i == None ]
            unknow_type_idx     = [n for n, i in enumerate(self.task_type_list) if i not in self.supported_task_type and \
                                                                                i != None ]

            # verify if all task_types of task_steps are supported  
            if  len(None_type_idx) != 0 or len(unknow_type_idx) != 0  :
                if len(None_type_idx) != 0 : 
                    None_type_task_step = self.task_step_list[None_type]
                    print(f"Task type of {None_type_task_step} is not seted !")
                if len(unknow_type_idx) != 0 : 
                    unknow_type_task_step = self.task_step_list[unknow_type_idx]
                    print(f"Task type of {unknow_type_task_step} is not recognized !")
                
                sys.exit()

            # start_step and start_mode
            start_step = self.workflow_json.get("start_step", 0)
            self.start_mode = self.workflow_json.get("start_mode", "from_scratch")

            # start step index
            start_step_idx = 0
            if isinstance(start_step , int) and start_step > len(self.task_type_list) :
                print("Start_step idx should be less then length of task_type_list ! Exit.")
                sys.exit()
            elif isinstance(start_step , int) and start_step < len(self.task_type_list):
                start_step_idx = start_step
            elif isinstance(start_step , str) and start_step in self.task_step_list:
                start_step_idx = self.task_step_list.index(start_step)
            self.start_step_idx = start_step_idx
            self.start_step     = self.task_step_list[start_step_idx]

            # run task steps
            for step_idx, task_step in enumerate(self.task_step_list):
                if step_idx < self.start_step_idx:
                    self.logger.logger.info(f"Task_step {task_step} is skipped")
                else:
                    self.logger.logger.info(f"Task_step {task_step} starts")

                self.step_idx = step_idx
                if step_idx >= self.start_step_idx and self.task_step_list[step_idx] != "cmd":

                    # extend run or continue run task steps
                    if self.start_mode in ["extend_run", "continue_run"] and step_idx == self.start_step_idx:
                        assert os.path.exists(os.path.abspath(task_step)), f"Step directory {task_step} doesn't exists"
                        os.chdir(os.path.join(self.base_path,task_step))
                        self.run_step(task_step)
                        self.result_pp(task_step)
                    else:
                        tmp = 1
                        step_path = os.path.join(self.base_path, task_step)
                        # handle old files
                        if os.path.exists(step_path):
                            while(os.path.exists(step_path + f"_bak_{tmp}")):
                                tmp += 1
                            os.rename(f"{step_path}", f"{step_path}_bak_{tmp}")

                        os.mkdir(step_path)
                        os.chdir(step_path)

                        # write itp, top, mdp, gro, files, and generate tpr file
                        self.gen_tpr(task_step)
                         # run gromacs
                        self.run_step(task_step)
                        # post precess 
                        self.result_pp(task_step)

                # cmd type canmmands
                elif step_idx >= self.start_step_idx and self.task_type_list[step_idx] == "cmd":
                     self.run_cmd(task_step)

                # back to base directory
                os.chdir(self.base_path)

                # save workflow result json file
                result_json_str = json.dumps(self.workflow_json)
                result_json = open("result.json","w")
                result_json.write(result_json_str)
                result_json.close()
                if step_idx >= self.start_step_idx:
                     self.logger.logger.info(f"Task_step {task_step} finish runing")

        self.logger.logger.info(f"All task steps in workflow finish runing")
                
    def gen_tpr(self, task_step ):

        step_para = self.workflow_json[task_step]        
        pos_restraint = step_para.get("pos_restraint", None)
        index_para    = step_para.get("index", None)
        replica_exchange =  step_para.get("replica_exchange", None)
        st_file = step_para.get("structure", None)
        ff_type = step_para.get("ff_type", "gaff")
        mdp_para = step_para.get("mdp_para", {})
        precision = step_para.get("precision", "single")
        ref_st    = step_para.get("ref_structure", None)


        # get default mdp parameter for task type 
        if step_para["task_type"] in ["minimization", "npt", "nvt", "free_energy", "pull",\
                                    "anneal","normal_mode"] :
            mdp_para_defalut = copy.deepcopy(eval(step_para["task_type"]))
        else: 
            mdp_para_defalut = copy.deepcopy(nvt)

        for i in mdp_para.keys():
            mdp_para_defalut[i] = mdp_para[i]

        mdp_para = mdp_para_defalut
        
        # get struture list
        if st_file == None and step_para["task_type"] not in ["rerun", "cmd"] :

            # use the input structure         
            if self.task_step_list[self.start_step_idx] != task_step:
                st_file = step_para.get("structure", None)

            # use the result structure of last step as inital structure input for this step
            if st_file == None and self.step_idx > self.start_step_idx and \
               self.workflow_json[self.task_step_list[self.step_idx - 1]]['result'].get("structure", False):
                st_file = self.workflow_json[self.task_step_list[self.step_idx - 1]]['result']["structure"].get("data", None)

            # use the md simulaiton out structure of step before as inital structure input for this step
            if st_file == None and self.step_idx > self.start_step_idx and task_type_list[self.step_idx - 1] in ["rerun", "cmd"]:
                st_file = self.workflow_json[self.task_step_list[self.step_idx - 2]]['result'].get("final_st", None)

        elif st_file == None and step_para["task_type"] == "rerun":
            # get input struture of rerun task step as input for rerun step
            task_step_path = os.path.join(self.base_path, step_para["task_step"])
            
            # use the md simulaiton out structure of step before as inital structure input for rerun step
            subtask_step_path = [os.path.abspath(i) for i in glob.glob(task_step_path + "/md_*") if os.path.isdir(i)] 
            if len(subtask_step_path) == 0:
                st_file = task_step_path + "/md.gro"
            else:
                st_file = [f"{i}/md.gro" for i in subtask_step_path]

            if mdp_para == {}:
                mdp_para = self.workflow_json[step_para["task_step"]].get("mdp_para", None)

        assert st_file != None, f" Input structure for step {task_step} is not setted!"


        self.logger.logger.info(f"The task type is {step_para['task_type']}, and mdp parameters is {mdp_para}")

        # use different structure formats for task step using different precision of gromacs
        # g96 format has more significant figures 
        if precision == "double":
            conf_name = "conf.g96"
        else:
            conf_name = "conf.gro"

        is_start_step = task_step == self.task_step_list[self.start_step_idx]
        if not is_start_step or self.start_mode in ["rerun", "from_scratch"] :

            if st_file == None:
                print(f"Input structure of step {task_step} is empty! Exit.")
                sys.exit()

            # if st_file is a list of structures
            structure = []
            if isinstance(st_file, list):
                for st in st_file:
                    tmp = read_st(st)
                    structure.append(tmp)
            else:
                structure = read_st(st_file)


            # read forcefield files
            ff_list = step_para.get("ff_list", []) 
            if step_para["task_type"] == "rerun" and  ff_list == []:
                ff_list = self.workflow_json[step_para["task_step"]].get("ff_list", [])
                ff_list = [itp_read(ff, ff_type = ff_type) for ff in ff_list]
            else:
                ff_list = [itp_read(ff, ff_type = ff_type) for ff in step_para.get("ff_list",[])]
            
            assert ff_list != [], " Forcefield is empty!"

            # generate tpr file
            self.gmx.gen_tpr(mdp_para, structure, ff_list, replica_exchange = replica_exchange, \
                             pos_restraint = pos_restraint, index_para = index_para, conf_name = conf_name,\
                             ref_st = ref_st)

            

        elif is_start_step and self.start_mode == "extend_run":
            
            extend_time = self.workflow_json.get("extend_time", None)
            if extend_time == None :
                print("Start_mode is extend_run, but extend_time item is not setted! Exit. ")
                sys.exit()
            # convert tpr file with extend time 
            self.gmx.convert_tpr(extend = True, time = extend_time)
           


    def run_step(self, task_step):
        
        step_para = self.workflow_json[task_step]
        task_type = step_para["task_type"]
        run_type  = step_para.get("run_type", "ntmpi")
        gpu_id    = step_para.get("gpu_id",   None )
        core_num  = step_para.get("core_num", 1)
        mpi_num   = step_para.get("mpi_num", 1)
        omp_num   = step_para.get("omp_num", 1)
        replica_exchange = step_para.get("replica_exchange", None)
        precision = step_para.get("precision", "single")
        assert precision in ["single", "double"],  "precision should be single or double!"
        mdp_para = step_para.get("mdp_para", {})
        replex = step_para.get("replex", 400)
        group = step_para.get("group", 0)

        # runtype
        if "gpu" in run_type:
            gpu = True
        else:
            gpu = None
        
        if "openmpi" in run_type :
            run_type0 = "openmpi"
        elif "ntmpi" in run_type :
            run_type0 = "ntmpi"
        else:
            run_type0 = "ntomp"
        
        if run_type0 != "ntomp" and  "ntomp" in run_type:
            run_type0 = run_type0 +"_ntomp"




        # get para_list
        para_list_len = len([i for i in glob.glob("./md_*") if os.path.isdir(i)])
        if isinstance(group, list) and para_list_len > 0:
            assert len(group) == para_list_len, "length of group and para_list should be same"
        elif para_list_len > 0:
            group = [group] * para_list_len

        # whetehr a extend run or continue run task step or not 
        is_start_step = task_step == self.task_step_list[self.start_step_idx]
        if is_start_step and self.start_mode in["extend_run", "continue_run"] :
            continue_run = True
            append = True
        else:
            continue_run = None
            append = None

        # whetehr a rerun task step or not 
        if task_type == "rerun":
            task_step_path = os.path.join(self.base_path, step_para["task_step"])
            subtask_step_path = [os.path.abspath(i) for i in glob.glob(task_step_path + "/md_*") if os.path.isdir(i)]

        # if thre are subtasks in this in task step
        if para_list_len > 0 :
            
            for i in range(para_list_len):
                self.logger.logger.info(f"The subtask md_{i} of {task_step} starts")

                # if task type is rerun
                os.chdir(f"./md_{i}")
                if task_type == "rerun":
                    traj_file = os.path.join(subtask_step_path[i] ,"md.xtc")
                    if not os.path.exists(traj_file):
                        traj_file = os.path.join(subtask_step_path[i] ,"md.trr")
                    assert os.path.exists(traj_file) , f"trajctory file {traj_file} does not exist"
                    out_file = "./md." + traj_file[-3:]
                    index_file = "./index.ndx" if os.path.exists("./index.ndx") else None

                    # convert trajctory files
                    if group[i] != 0:
                        self.gmx.trjconv(traj_file_in = traj_file, tpr_file = "./md.tpr", out_file =out_file,
                                        group = group[i], index_file = index_file)
                    else:
                        out_file = traj_file
                    
                    # convert tpr files
                    if index_file != None:
                        self.gmx.convert_tpr(tpr_file = "./md.tpr", group = group[i], index_file = index_file)
                    
                    self.gmx.rerun(traj_file = out_file, tpr_file = "./md.tpr", core_num = core_num,  )
                
                # mormal mode task type
                else:
                    self.gmx.run(dir_path = "./", task_type = task_type,  run_type = run_type,
                        replica_exchange = replica_exchange, core_num = core_num, omp_num = omp_num, mpi_num = mpi_num, 
                        precision = precision, continue_run = continue_run, append = append, replex = replex,
                        gpu = gpu, gpu_id = gpu_id )
                                        
                    if task_type == "normal_mode":
                        self.gmx.nmeig()
                self.logger.logger.info(f"The subtask md_{i} of {task_step} finishes")
                os.chdir("../")
        else:
            # if task type is normal md task
            self.gmx.run(dir_path = "./", task_type = task_type,  run_type = run_type,
                    replica_exchange = replica_exchange, core_num = core_num, omp_num = omp_num, mpi_num = mpi_num, 
                    precision = precision, continue_run = continue_run, append = append, replex = replex,
                    gpu = gpu, gpu_id = gpu_id )
            if task_type == "normal_mode":
                self.gmx.nmeig()
     
    def run_cmd(self, task_step):

        cmd_list = self.workflow_json[task_step].get("cmd", [])
        assert isinstance(cmd_list, str) or isinstance(cmd_list, list), "cmd should be a sring or list of string"
        if isinstance(cmd_list, str):
            cmd_list = [cmd_list]

        for cmd in cmd_list:
            run_cmd(cmd)


    def result_pp(self, task_step):
        
        # get para_list
        para_list_len = len([i for i in glob.glob("md_*") if os.path.isdir(i)])

        task_type = self.workflow_json[task_step]["task_type"]

        reuslt_para = self.workflow_json[task_step].get("result", {})
        self.workflow_json[task_step]["result"] = reuslt_para
        result_list = list(reuslt_para.keys())

        # no subtasks
        if (reuslt_para != {} and para_list_len == 0) or task_type == "normal_mode" :
            # task step is normal md task
            for result in result_list:
                self.workflow_json[task_step]['result'][result]['data'] = eval(f"self.result_{result}(task_step)")
            
            # normal mode task type 
            if task_type == "normal_mode":
                eigenfreq_data  = read_xvg("./eigenfreq.xvg")
                eigenvalue_data = read_xvg("./eigenval.xvg")

        # multi_subtasks in task step
        elif (reuslt_para != {} and para_list_len > 0) or task_type == "normal_mode":
            rmsd_data       = []
            energy_data     = []
            structure_data  = []
            eigenfreq_data  = []
            eigenvalue_data = []
            for i in range(para_list_len):
                os.chdir(f"md_{i}")
                for j in result_list:
                    if j not in ["bar", "wham"]:
                        eval(f"{j}_data.append(self.result_{j}(task_step))")
                if task_type == "normal_mode":
                    eigenfreq_data.append(read_xvg("./eigenfreq.xvg"))
                    eigenvalue_data.append(read_xvg("./eigenval.xvg"))
                os.chdir("../")

            for result in result_list:
                self.logger.logger.info(f"{task_step} result process {result}")
                self.workflow_json[task_step]['result'][result]['data'] =  eval(f"{result}_data")
            

        if task_type == "normal_mode":
            self.workflow_json[task_step]["result"]["eigenfreq"]    = eigenfreq_data     
            self.workflow_json[task_step]["result"]["eigenvalue"]   = eigenvalue_data 

        # read output structure of task step 
        if "structure" not in result_list:
            if para_list_len == 0 and os.path.exists("./md.gro"):
                self.workflow_json[task_step]['result']["final_st"] = open("./md.gro", "r").read()
            else:                
                final_st = []
                for i in range(para_list_len):
                    if os.path.exists(f"md_{i}./md.gro"):
                        final_st.append(open(f"md_{i}./md.gro", "r").read())
                self.workflow_json[task_step]['result']["final_st"] = final_st

    def result_structure(self, task_step):

        # output structure format
        precision = self.workflow_json[task_step].get("precision", "single")
        if precision == "double":
            structure_format = "g96"
        else: 
            structure_format = "gro"
        out_structure = f"out.{structure_format}"

        result_para = self.workflow_json[task_step].get("result", {})

        if result_para != {}:
            # get output structure
            out_structure_para = result_para.get("structure", None)
            if out_structure_para :
                dump_time = result_para["structure"].get("dump_time", None)
                st_index_para = result_para["structure"].get("index", {})
                group = result_para["structure"].get("group", 0)
                
                # write index file
                st_index_file = None
                if st_index_para != {}:
                    write_index_file(st, st_index_para, index_file = "index_st.ndx" )
                    st_index_file = "index_st.ndx"

                # get structure at dump time
                if isinstance(dump_time, float) or isinstance(dump_time, int):   
                    self.gmx.trjconv(out_file = out_structure, index_file = st_index_file, \
                                     dump = dump_time, group = group, ndec = 12 )
                    out_structure_str = open(out_structure, "r").read()
                elif isinstance(dump_time, list):
                    out_structure_str = []
                    for tmp_dump_time in dump_time:
                        self.gmx.trjconv(out_file = out_structure, index_file = st_index_file, \
                                         dump = tmp_dump_time, group = group )
                        out_structure_str.append(open(out_structure, "r").read())

        return(out_structure_str)
                    
    def result_rmsd(self, task_step):
        result_para = self.workflow_json[task_step].get("result", {})
        out_rmsd_para = result_para.get("rmsd", None)
        st_index_para = result_para["rmsd"].get("index", {})
        group = result_para["rmsd"].get("group", 0)
        st_index_file = None
        if st_index_para:
            index_file(st, st_index_para, index_file = "index_st.ndx" )
            st_index_file = "index_st.ndx"
        self.gmx.rmsd(group = group , index_file = st_index_file)
        rmsd_data = read_xvg("./rmsd.xvg")
        return(rmsd_data) 

    def result_energy(self, task_step):
        result_para = self.workflow_json[task_step].get("result", {})
        selected_item = result_para["energy"].get("item", [])
        self.gmx.energy(item = selected_item)
        energy_data = read_xvg("./energy.xvg")
        return(energy_data)

    def result_eigendata(self, task_step):
        eigenfreq_data  = read_xvg("./eigenfrequency.xvg")
        eigenvalue_data = read_xvg("./eigenvalue.xvg")
        data = {}
        data["eigenfrequency"] = eigenfreq_data
        data["eigenvalue"]     = eigenvalue_data
        return(data)


    def result_wham(self, task_step):
        wham_para = self.workflow_json[task_step]["result"].get("wham", {})
        if wham_para != {}:
            self.gmx.gmx_wham(**wham_para)
            self.workflow_json[task_step]["result"]["wham"]["profile"] = read_xvg("./profile.xvg")
            self.workflow_json[task_step]["result"]["wham"]["histo"]   = read_xvg("./histo.xvg")

    def result_bar(self, task_step):
        bar_para = self.workflow_json[task_step]["result"].get("bar", {})
        if wham_para != {}:
            self.gmx.gmx_wham(**bar_para)
            self.workflow_json[task_step]["result"]["bar"]["bar"] = read_xvg("./bar.xvg")
            self.workflow_json[task_step]["result"]["bar"]["barint"]   = read_xvg("./barint.xvg")
            self.workflow_json[task_step]["result"]["bar"]["histogram"]   = read_xvg("./histogram.xvg")
        
        



    