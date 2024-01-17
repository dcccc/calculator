import os, copy, subprocess,glob,re
import numpy as np


def run_cmd(cmd, task_name=""):

    if task_name =="":
        task_name = "task"

    cmd = cmd + f" 1>{task_name}.out  2>{task_name}.err"
    status, stdout = subprocess.getstatusoutput(cmd)
    stderr_str = open(f"./{task_name}.err", "r").read()            
    assert status == 0, stderr_str


class Orca():
    def __init__(self, orca_exe={}):        
              
        self.orca_exe = orca_exe       
        self.orca     = self.orca_exe["orca"]
        self.openmpi  = self.orca_exe.get("openmpi", None)
        
        # set the environment variables for openmpi
        if self.openmpi:
            openmpi_root = os.path.dirname(os.path.dirname(orca_exe["openmpi"]))
            os.environ["LD_LIBRARY_PATH"] = os.environ.get("LD_LIBRARY_PATH", "") + ":" + openmpi_root + "/lib"
            os.environ["OPAL_PREFIX"]     = openmpi_root
            os.environ["PATH"]            = os.environ.get("PATH", "") + ":" + openmpi_root + "/bin"
        
        # set the environment variables for orca
        orca_root = os.path.dirname(orca_exe["orca"])
        os.environ["LD_LIBRARY_PATH"]     = os.environ.get("LD_LIBRARY_PATH", "") + ":" + orca_root
        os.environ["PATH"]                = os.environ.get("PATH", "") + ":" +  orca_root + "/bin"

        
#     @classmethod
    def write_input(self, mol, para, print_input=False, write_input=True):
        orca_input = ""

        # the size of memery used for every core
        maxcore = para.get("maxcore", None)
        if maxcore:
            orca_input += "%MaxCore  {:d} \n".format(maxcore)

        # the number of cores used in calculation
        self.nprocs  = para.get("nprocs", 1)
        if self.nprocs > 1:
            orca_input += "%PAL NPROCS {:d} END \n".format(self.nprocs)
        
        # whether set the number of cores using the keyword in orca input file
        if self.nprocs == 1:
            nprocs = re.findall("\spal(\d+)\s", para["keyword"].lower())
            if len(nprocs) > 0:
                self.nprocs = int(nprocs[0])


        # keyword line
        keyword_line = para.get("keyword", "")
        if keyword_line != "":
            orca_input += "! {:s}\n".format(keyword_line)
        else:
            assert False, print("keyword is necessary !")

        # detail and additioanl parameters 
        detail_para = para.get("detail_para", None)

        def json_para_to_str_line(item):
            # convert json para to str line
            str_line = ""
            for key in item.keys():
                str_line += "    {: <10s} {}  \n".format(key, item[key])
            return(str_line)


        if detail_para:
            # detail_para is a dict
            if isinstance(detail_para, dict):
                for item in detail_para.keys():
                    # ignore the keyword "pal"
                    if item.lower() != "pal":
                        para_item_line = json_para_to_str_line(detail_para[item])
                        orca_input += "%{:s}\n{:s}END\n".format(item, para_item_line)
            else:
                assert False, print("detail_para should be a dict !")

        para_after_coord = para.get("para_after_coord", None)
        para_after_coord_txt = ""
        if para_after_coord:
            # detail_para is a dict
            if isinstance(para_after_coord, dict):
                for item in para_after_coord.keys():
                    para_item_line = json_para_to_str_line(para_after_coord[item])
                    para_after_coord_txt += "%{:s}\n{:s}END\n".format(item, para_item_line)
            else:
                assert False, print("para_after_coord should be a dict !")


        # mol atom symbols and coordinates
        mol_lines = "* {} {: >2d}  {: >2d} \n {}*"
        coord_type = mol.get("coord_type", "xyz")
        assert coord_type in ["xyz", "int"], print("coord_type should be 'xyz' or 'int'")
        orca_input += mol_lines.format(coord_type, mol["charge"], mol["spin"], mol["coord"] )
        orca_input += "\n" + para_after_coord_txt

        # wirte orca input file
        if write_input:
            orca_input_file = open("./orca.in", "w")
            orca_input_file.write(orca_input)
            orca_input_file.close()

        if print_input:
            print(orca_input)


#     @classmethod
    def run(self, mpirun_para=""):
        
        # openmpi parallel running
        if self.nprocs :
            cmd = '{} orca.in  " {}  --oversubscribe "  > orca.out'.format(self.orca, mpirun_para)  
        else:
            cmd = '{}         orca.in > orca.out'.format(self.orca)  
        
        # run orca
        run_cmd(cmd, task_name="orca")
        

    @staticmethod
    def get_energy(energy_list=[]):
        # read orca property file, get energy items from orca property file 
        # according to the keywords in energy_list, and return a dict
        # data format: {"energy": final_energy, energy_item_name: energy_item, ...}
        # unit: Eh
        # scf energy is the "FINAL SINGLE POINT ENERGY" in orca.out file, which includes
        # with zero point energy, DFT-D correction, solvation effect, etc.

        energy_dict = {}
        
        orca_out = open("./orca.out","r").read()
        # get final scf energy
        final_energy = re.findall("FINAL SINGLE POINT ENERGY\s+-?(\d+\.*\d+)", orca_out)
        if len(final_energy) >0 :
            final_energy = np.float64(final_energy[-1])
            energy_dict["energy"] = final_energy
        
        orca_property = open("./orca_property.txt","r").read()
        # get other energy items
        if energy_list != []:
            for i in energy_list:
                pattern = re.sub(r"\(", "\\(", i)
                pattern = re.sub(r"\)", "\\)", pattern)
                energy_result = re.findall(f"{pattern}\s*:?\s+(-?\d+\.*\d+)", orca_property)
                if len(energy_result) > 0:
                    energy_dict[i] = np.float64(energy_result[-1])
                else:
                    energy_dict[i] = None

        return(energy_dict)


    @staticmethod
    def get_xyz():
        # read final structure xyz file
        xyz_txt = ""
        if os.path.exists("orca.xyz"):
            xyz_txt = open("orca.xyz", "r").read()
        return(xyz_txt)

    @staticmethod
    def get_trajectory():
        # read trajectory
        traj = ""
        if os.path.exists("./orca_trj.xyz"):
            traj = open("./orca_trj.xyz", "r").read()
        return(traj)
    
    @staticmethod
    def get_hessian():
        # read hessian matrix and vibrational frequencies
        # data format: [hessian matrix, vibrational frequencies]
        # hessian matrix data format : [[dE/dx1dx1, dE/dx1dy1, dE/dx1dz1, dE/dx1dx2, dE/dx1dy2, dE/dx1dz2, ...],
        #                               [dE/dx2dx1, dE/dx2dy1, dE/dx2dz1, dE/dx2dx2, dE/dx2dy2, dE/dx2dz2, ...]]
        # unit : Eh/bohr2

        # vibrational frequencies data format: [freq1, freq2, freq3, ...] unit : cm-1

        hessian = "None"
        freq    = "None"
        if os.path.exists("./orca.hess"):
            hess = open("./orca.hess", "r").read()

            # read hessian matrix
            hessian = re.findall("\$hessian(.*?)\$", hess, flags=re.DOTALL)[0]
            if len(hessian) != 0:

                hessian = hessian.split("\n")[1:]
                dim = int(hessian[0])
                hessian = [i.strip().split()[1:] for i in hessian[2:] if "E" in i]
                hessian = [np.hstack(hessian[i::dim]) for i in range(dim)]
                hessian = np.array(hessian, dtype=np.float64)

                # read vibrational frequencies
                freq = re.findall("\$vibrational_frequencies(.*?)\$", hess, flags=re.DOTALL)[0]
                freq = freq.split("\n")[2:]
                freq = [i.strip().split()[1] for i in freq if len(i)>0]
                freq = np.array(freq, dtype=np.float64)

        return([hessian, freq])

    @staticmethod
    def get_ir():
        # read IR spectrum
        # data format: [[wavenumber(cm-1), intensity(L/(mol*cm))]...]
        ir_spectrum = "None"
        if os.path.exists("./orca.hess"):
            hess = open("./orca.hess", "r").read()

            if "$ir_spectrum" in hess:
                ir_spectrum = re.findall("\$ir_spectrum(.*?)\$", hess, flags=re.DOTALL)[0]
                ir_spectrum = ir_spectrum.split("\n")[2:]
                ir_spectrum = [i.strip().split()[:2] for i in ir_spectrum if len(i)>0]
                ir_spectrum = np.array(ir_spectrum, dtype=np.float64)[:,:2]
                ir_spectrum = ir_spectrum[ir_spectrum[:,0] > 0]

        return(ir_spectrum)

    @staticmethod
    def get_nmr():
        # read NMR data
        # data format: [[atom index A, atom index B, spin-spin coupling constant(Hz)]...]
        orca_property = open("./orca_property.txt","r").read()

        nmr_result = "None"
        if "$ EPRNMR_SSCoupling" in orca_property:
            nmr_txt = re.findall("\$ EPRNMR_SSCoupling(.*?)#", orca_property, flags=re.DOTALL)[0]
            if len(nmr_txt) != 0:
                idxA  = re.findall("Index A:\s+(\d+)", nmr_txt)
                idxB  = re.findall("Index B:\s+(\d+)", nmr_txt)
                nmr = re.findall("Total Spin-Spin Coupling ISO:\s+ ([-]+\d+\.\d+)", nmr_txt)

            nmr_result =  [[idxA[i], idxB[i], nmr[i]] for i in range(len(nmr))]
            nmr_result = np.array(nmr_result, dtype=np.float64)

        return(nmr_result)

    @staticmethod
    def get_raman():
        # read Raman spectrum
        # data format: [[wavenumber(cm-1), activity(A^4/AMU))]...]
        raman_spectrum = "None"
        if os.path.exists("./orca.hess"):
            hess = open("./orca.hess", "r").read()

            if "$raman_spectrum" in hess:
                raman_spectrum = re.findall("\$raman_spectrum(.*?)\$", hess, flags=re.DOTALL)[0]
                raman_spectrum = raman_spectrum.split("\n")[2:]
                raman_spectrum = [i.strip().split()[:3] for i in raman_spectrum if len(i)>0]
                raman_spectrum = np.array(raman_spectrum, dtype=np.float64)[:,[0,2]]
                raman_spectrum = raman_spectrum[raman_spectrum[:,0] > 0]

        return(raman_spectrum)


    @staticmethod
    def get_uv():
        # read uv spectrum
        # data format: [[wavelength(ev), fosc]...]

        uv_spectrum = "None"
        orca_property = open("./orca_property.txt","r").read()
        if "The CIS Absorption Spectrum" in orca_property:
            uv_spectrum = re.findall("The CIS Absorption Spectrum(.*?)The CIS Absorption Spectrum", orca_property, flags=re.DOTALL)[0]
            uv_spectrum = uv_spectrum.split("\n")[2:]
            uv_spectrum = [i.strip().split()[:3] for i in uv_spectrum if len(i.strip())>2]
            uv_spectrum = np.array(uv_spectrum, dtype=np.float64)
            uv_spectrum = uv_spectrum[:,1:3 ]

        return(uv_spectrum)

    @staticmethod
    def get_ecd():
        # read ECD spectrum
        # data format: [[wavelength(ev), R(1e40*cgs)]...]

        ECD_spectrum = "None"
        uv_spectrum = Orca.get_uv()
        orca_property = open("./orca_property.txt","r").read()
        if "CIS_CD" in orca_property:
            ECD_spectrum = re.findall("\$ CIS_CD(.*?)#", orca_property, flags=re.DOTALL)[0]
            ECD_spectrum = ECD_spectrum.split("\n")[7:]
            ECD_spectrum = [i.strip().split()[:3] for i in ECD_spectrum if len(i.strip())>2]
            ECD_spectrum = np.array(ECD_spectrum, dtype=np.float64)
            ECD_spectrum = ECD_spectrum[:,1]

            ECD_spectrum = np.vstack((uv_spectrum[:,0], ECD_spectrum)).T

        return(ECD_spectrum)


    @staticmethod
    def get_dipole():
        # read dipole moment
        # data format: [dipole_x, dipole_y, dipole_z] unit: Debye
        dipole = "None"
        orca_property = open("./orca_property.txt","r").read()
        if "Total Dipole moment" in orca_property:
            dipole = re.findall("Total Dipole moment:(.*?)#", orca_property, flags=re.DOTALL)[0]
            dipole = dipole.split("\n")[2:]
            dipole = [i.strip().split()[1] for i in dipole if len(i)>0]
            dipole = np.array(dipole, dtype=np.float64)

        return(dipole)

    @staticmethod
    def get_force():
        # read atom force
        # data format: [[force_x, force_y, force_z]...] unit: Eh/bohr
        force = "None"
        if os.path.exists("./orca.engrad"):
            engrad = open("./orca.engrad","r").read()
            force = re.findall("Eh/bohr\n#(.*?)#", engrad, flags=re.DOTALL)[0]
            force = force.split("\n")
            force = [i.strip().split()[0] for i in force if len(i)>0]
            force = -1*np.array(force, dtype=np.float64).reshape((-1,3))

        return(force)

    @staticmethod
    def get_soc():
        '''
        read soc form orca output file
        result data is soc result of between triplet substates and singlets
        data format: [[T root, S root, real part(ms= 0), imaginary  part(ms= 0), 
                                       real part(ms=-1), imaginary  part(ms=-1), 
                                       real part(ms= 1), imaginary  part(ms= 1) ]...]
                                       unit : cm-1
        '''
        soc = "None"
        orca_out = open("./orca.out","r").read()
        if "CALCULATED SOCME BETWEEN TRIPLETS AND SINGLETS" in orca_out:
            soc = re.findall("T      S           MS= 0                  \-1                    \+1(.*?)CALCULATED REDUCED SOCME BETWEEN TRIPLETS AND SINGLETS", orca_out, flags=re.DOTALL)[0]
            soc = re.sub(r"[\(\)\,]?", "", soc)
            soc = soc.split("\n")[2:-3]
            soc = [i.strip().split() for i in soc]
            soc = np.array(soc, dtype=np.float64)

        return(soc)    

    @staticmethod
    def get_mo_energy():
        # read molecular orbital energy
        # data format: [occupation(e), mo_energy(Eh)]
        mo_energy = "None"
        if os.path.exists("./orca.out"):
            mo_energy = open("./orca.out","r").read()
            mo_energy = re.findall("ORBITAL ENERGIES(.*?)\*", mo_energy, flags=re.DOTALL)[-1]
            mo_energy = mo_energy.split("\n")[4:-2]
            mo_energy = [i.strip().split() for i in mo_energy if len(i)>0]
            mo_energy = np.array(mo_energy, dtype=np.float64)[:,1:3]

        return(mo_energy)
    
    @staticmethod
    def get_runtime():
        # get total run time, unit : s
        runtime_line = open("./orca.out","r").readlines()[-1]
        if "TOTAL RUN TIME:" in runtime_line:
            runtime_line = runtime_line.strip().split()
            runtime      = [runtime_line[n] for n in [3, 5, 7, 9, 11]]
            runtime      = np.array(runtime, dtype=np.float32)
            runtime      = np.sum(runtime * np.array([24*3600, 3600, 60, 1, 0.001]))
        else:
            runtime      = None

        return(runtime)