import os, copy, subprocess, glob, re, json
import numpy as np
from calculator.util.logger  import *
from calculator.orca.orca import Orca, run_cmd

class Orca_workflow():
    def __init__(self, orca_exe={}, workflow_json = {}):
        self.orca_exe      = orca_exe    
        self.workflow_json = workflow_json
        self.orca          = Orca(orca_exe=orca_exe)
        self.base_path     = os.getcwd()
        self.logger        = Logger()

    def run_step(self, task_step):
        
        os.chdir(self.base_path)
        
        # run one step of the workflow
        run_type = self.workflow_json[task_step].get("run_type", None)

        # if the task step is skipped
        if run_type == "skipped":
            self.logger.logger.info(f"Task step {task_step} is skipped!")
        # if the task step is a cmd step
        elif run_type == "cmd":
            if os.path.exists(task_step+".err"):
                nn = 1
                backup_file = f"{task_step}_bak{nn}"
                while(os.path.exists(backup_file+".err")):
                    backup_dir = f"{task_step}_bak{nn}"
                    nn = nn + 1
                
                os.system(f"mv {task_step}.err {backup_file}.err")
                os.system(f"mv {task_step}.out {backup_file}.out")
                self.logger.logger.info(f"Task step file '{task_step}'.out and '{task_step}'.err already exist, backup it to '{backup_file}.out' and '{backup_file}'.err")
            self.run_cmd(task_step)        
        # run task step
        else:
            # get para
            para = self.workflow_json[task_step].get("para", None)
            # check para 
            if para is None:
                self.logger.logger.info(f"Task parameter of task step '{task_step}' is missing !")
                assert False
            else:
                # check para keyword 
                if para.get("keyword", None) is None:
                    self.logger.logger.info(f"Task keyword of task step '{task_step}' is missing !")
                    assert False

            # get mol
            mol  = self.workflow_json[task_step].get("mol", None)
            if mol is None:
                self.logger.logger.info(f"Mol of task step {task_step} is missing! Mol from last step result will be used !")
                
                task_step_last = self.valid_task_step[self.valid_task_step.index(task_step)-1]
                # get the mol from last step result
                xyz_line = self.workflow_json[task_step_last]["result"].get("xyz", None)
                if xyz_line != None:
                    xyz_line = xyz_line.split("\n")
                    xyz_line = "\n".join(xyz_line[2:])
                # if the last step is a not a optimize step, get the mol from input 
                else:
                    xyz_line = self.workflow_json[task_step_last]["mol"]["coord"]
                
                mol = self.workflow_json[task_step_last]["mol"]
                mol["coord_type"]  = "xyz"                
                mol["coord"]  = xyz_line
                
                self.workflow_json[task_step]["mol"] = mol 
            
            # create folder for the task step 
            task_step_path = os.path.join(self.base_path, task_step)
            if not os.path.exists(task_step_path):
                os.mkdir(task_step_path)
            elif run_type == "continue_run":
                pass
            else:
                # backup existed directory
                nn = 1
                backup_dir = f"{task_step_path}_bak{nn}"
                while(os.path.exists(backup_dir)):
                    backup_dir = f"{task_step_path}_bak{nn}"
                    nn = nn + 1
                os.system(f"mv {task_step_path} {backup_dir}")
                os.mkdir(task_step_path)
                self.logger.logger.info(f"Task step directory '{task_step}' already exist, backup it to '{backup_dir}' ")
            os.chdir(task_step_path)

            # write input file
            self.orca.write_input(mol=mol, para=para)

            # run orca
            self.logger.logger.info(f"Task step '{task_step}' is running ! orca task keyword is '{para['keyword']}' ")
            mpirun_para = self.workflow_json[task_step].get("mpirun_para", "")
            self.orca.run(mpirun_para=mpirun_para)
            
            # get the result of the task step
            if run_type != "cmd":
                result = self.get_result(task_step)
                self.workflow_json[task_step]["result"] = result

            # delete files
            self.delete_files(file_list=self.workflow_json[task_step].get("delete_files", []))

            self.logger.logger.info(f"Task step '{task_step}' finished!")
        
        os.chdir(self.base_path)
        
    def get_result(self, task_step):
        # get the result of the task step
        result = self.workflow_json[task_step].get("result", None)

        xyz_result    = self.orca.get_xyz()
        energy_result = self.orca.get_energy()
        runtime       = self.orca.get_runtime()
        if runtime:
            coretime  = runtime * self.orca.nprocs
        else:
            coretime = None
        default_result = {"xyz": xyz_result, "energy": energy_result,
                          "runtime":runtime, "coretime":coretime}
        if result is None or result == [] or result == {}:
            result = default_result
        else:
            result_list = list(result.keys()) if isinstance(result, dict) else result
            assert isinstance(result_list, list), print("result should be a dict or a list!")

            supported_result_list = ["xyz", "energy", "trajectory", "hessian", "ir", "nmr", "raman", 
                                     "uv",  "ecd", "dipole", "force", "soc"]
            
            result = default_result
            if "energy" in result.keys():
                result["energy"] = self.orca.get_energy(energy_list=result["energy"])

            # get other result items
            for result_item in result_list:
                if result_item not in ["xyz", "energy", "runtime", "coretime"]:
                    assert result_item in supported_result_list, print(f"result '{result_item}' is not supported!")

                    result[result_item] = eval(f'self.orca.get_{result_item}()')
                    if isinstance(result[result_item], np.ndarray):
                        result[result_item] = result[result_item].tolist()
                    elif isinstance(result[result_item], list):
                        result[result_item] = [i.tolist() for i in result[result_item] if isinstance(i, np.ndarray)]

        return(result)




    def task_run(self):
        # run the workflow

        # get the task step list
        self.task_step_all_list = self.workflow_json["task_step"]

        # check whether task steps are in the list or not 
        self.valid_task_step = [i for i in self.task_step_all_list if self.workflow_json.get(i, None) != None]
        self.unvalid_task_step = [i for i in self.task_step_all_list if self.workflow_json.get(i, None) == None]
        assert len(self.unvalid_task_step) == 0, print(f"Task steps '{self.unvalid_task_step}' is not valid!")

        # get the task step list that will be run
        self.task_step_run_list = [i for i in self.task_step_all_list 
                                           if self.workflow_json[i].get("run_type", None) != "skipped"]
        

        # if no task step in the workflow 
        if len(self.task_step_run_list) == 0:
            self.logger.logger.info("No task step is found in the task workflow !")
        # run task step
        else:
            for task_step in self.task_step_all_list:
                self.run_step(task_step)
                json.dump(self.workflow_json, open("workflow.json", "w"), indent=4)

    def delete_files(self, file_list=[] ):
        # delete files in the task step directory
        if file_list != []:
            file_list = " ".join(file_list)
            run_cmd(f"rm -rf {file_list}", task_name="delete_files")

    def run_cmd(self, task_step):

        cmd_list = self.workflow_json[task_step].get("cmd", [])
        assert isinstance(cmd_list, str) or isinstance(cmd_list, list), "cmd should be a sring or list of string"
        if isinstance(cmd_list, str):
            cmd_list = [cmd_list]

        for cmd in cmd_list:
            assert isinstance(cmd, str)  , "cmd should be a string !"
            if cmd != "":
                self.logger.logger.info(f"Task step '{task_step}' cmd '{cmd}' is running ! ")
                run_cmd(cmd, task_name=task_step)
        self.logger.logger.info(f"Task step '{task_step}' finished ! ")