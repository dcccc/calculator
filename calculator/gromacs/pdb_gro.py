import numpy as np



def latt6_to_latt9(latt6):
    cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma = latt6

    alpha_r = np.radians(cell_alpha)
    beta_r = np.radians(cell_beta)
    gamma_r = np.radians(cell_gamma)
    v_ = 1 - (np.cos(alpha_r))**2 - (np.cos(beta_r))**2 - (np.cos(gamma_r))**2 + \
         2 * np.cos(alpha_r) * np.cos(beta_r) * np.cos(gamma_r)

    assert v_ > 0, "lattice error: alpha + beta < gamma"

    vector_a = [float(cell_a), 0., 0.]
    vector_b = [cell_b * np.cos(gamma_r), cell_b * np.sin(gamma_r), 0.]
    vector_c = [cell_c * np.cos(beta_r),
                cell_c * (np.cos(alpha_r) - np.cos(beta_r) * np.cos(gamma_r)) / np.sin(gamma_r),
                cell_c / np.sin(gamma_r) * np.sqrt(v_)]

    return np.array([vector_a, vector_b, vector_c])

def latt9_to_latt6(latt9):
    
    a, b, c = np.linalg.norm(latt9, axis = 1)
    
    alpha = np.rad2deg(np.arccos(np.dot(latt9[1], latt9[2]) / b / c))
    beta  = np.rad2deg(np.arccos(np.dot(latt9[0], latt9[2]) / a / c))
    gamma = np.rad2deg(np.arccos(np.dot(latt9[0], latt9[1]) / a / b))
    
    return(np.array([a,b,c,alpha,beta,gamma]))

def write_xyz(xyz_name, atom_symbols, coords):

    xyz_line = ["%d" %(len(atom_symbols)) ,""]
    for i in range(len(atom_titles)):
        line = "%-4s  %9.4f  %9.4f  %9.4f    " \
                %(atom_symbols[i], coord[i][0], coord[i][1], coord[i][2])
        xyz_line.append(line)

    
    xyz_file = open(xyz_name, "w")
    xyz_file.write("\n".join(pdb_line))
    xyz_file.close()
    
    

class Structure():
    def __init__(self, latt, coord, atom_titles, mol_titles, mol_idx):
        
        if len(latt) == 6:
            self.latt6   = latt
            self.latt9   = latt6_to_latt9(self.latt6)
        elif isinstance(latt, np.ndarray) and latt.shape == (3,3):
            self.latt9   = latt
            self.latt6   = latt9_to_latt6(self.latt9)
        

        self.coord       = np.array(coord)
        self.atom_titles = atom_titles
        self.mol_titles  = mol_titles
        self.mol_idx     = mol_idx
        self.atom_num    = len(self.coord)
        self.frac_coord  = np.dot(self.coord, np.linalg.inv(self.latt9))
        atom_types       = self.get_mol_types()
        

    def get_mol_types(self): 
        
        tmp_list         = [i+str(j) for i,j in zip(self.mol_titles, self.mol_idx)]
        self.mol_idx     = np.array(self.mol_idx, dtype=np.int) 
        self.mol_types   = [ self.mol_titles[0] ]
        for n in range(0, self.atom_num - 1):
            if tmp_list[n] != tmp_list[n+1]:
                self.mol_types.append(self.mol_titles[n+1])
    
    
    @classmethod
    def read_gro_str(cls, gro_str):
        
        gro_line = [i for i in gro_str.split("\n") if len(i) > 1]
        
        atom_num = int(gro_line[1].strip())
        temp_latt9    = gro_line[-1].split()
        if len(temp_latt9) == 9:
            temp_latt9    = [[temp_latt9[0], temp_latt9[3], temp_latt9[4]],
                             [temp_latt9[5], temp_latt9[1], temp_latt9[6]],
                             [temp_latt9[7], temp_latt9[8], temp_latt9[2]] ]
        elif len(temp_latt9) == 3:
            temp_latt9    = [[temp_latt9[0],           0.0,           0.0],
                             [          0.0, temp_latt9[1],           0.0],
                             [          0.0,           0.0, temp_latt9[2]] ]
        else:
            temp_latt9 = np.zeros((3,3))
            
        
        latt9    = np.array(temp_latt9, dtype=np.float64) * 10.0
        
        mol_idx = []
        mol_titles =[]
        atom_titles = []
        coord = []
        for  line in gro_line[2:-1]:
            if len(line)>2 :
                mol_titles.append(line[5:10].strip())
                atom_titles.append(line[10:15].strip())
                coord.append(line[20:].split()[:3])
                mol_idx.append(line[:5].strip())

        coord       = np.array(coord, dtype=np.float64)*10.0
        atom_titles = np.array(atom_titles)
        mol_titles  = np.array(mol_titles)
        atom_num    = len(coord)
        
        mol_idx     = np.array(mol_idx, dtype=np.int32) 
       
        st = cls(latt9, coord, atom_titles, mol_titles, mol_idx)
        
        return(st)
    
    @classmethod
    def read_gro(cls, gro_file):        
        gro_str = open(gro_file, "r").read()  
        return(cls.read_gro_str(gro_str))

    @classmethod
    def read_g96_str(cls, g96_str):
        
        g96_line = [i for i in g96_str.split("\n") if len(i) > 1]
     
        temp_latt9    = np.array([0.0] * 9, dtype=np.float64) * 10.0
        
        mol_idx = []
        mol_titles =[]
        atom_titles = []
        coord = []
        tmp_tag = ""
        for  line in g96_line:
            line_strip = line.strip()
            if line_strip == "POSITION":
                tmp_tag = "POSITION"
            if len(line_strip.split())>2  and tmp_tag == "POSITION":
                mol_titles.append(line[6:12].strip())
                atom_titles.append(line[12:18].strip())
                coord.append(line[25:].split()[:3])
                mol_idx.append(line[:6].strip())

            if line_strip == "BOX":
                tmp_tag = "BOX"
            if len(line_strip.split())>2  and tmp_tag == "BOX":
                temp_latt9 = np.array(line.split(), dtype=np.float64)
                if len(temp_latt9) == 9:
                    temp_latt9    = [[temp_latt9[0], temp_latt9[3], temp_latt9[4]],
                                     [temp_latt9[5], temp_latt9[1], temp_latt9[6]],
                                     [temp_latt9[7], temp_latt9[8], temp_latt9[2]] ]
                elif len(temp_latt9) == 3:
                    temp_latt9    = [[temp_latt9[0],           0.0,           0.0],
                                     [          0.0, temp_latt9[1],           0.0],
                                     [          0.0,           0.0, temp_latt9[2]] ]


            latt9 = np.array(temp_latt9) * 10.0


            if line_strip == "END":
                tmp_tag = ""

        coord       = np.array(coord, dtype=np.float64)*10.0
        atom_titles = np.array(atom_titles)
        mol_titles  = np.array(mol_titles)
        atom_num    = len(coord)
        
        mol_idx     = np.array(mol_idx, dtype=np.int32) 
        mol_num     = len(set(mol_idx))
       
        st = cls(latt9, coord, atom_titles, mol_titles, mol_idx)
        
        return(st)
    
    @classmethod
    def read_g96(cls, g96_file):        
        g96_str = open(g96_file, "r").read()  
        return(cls.read_g96_str(g96_str))



    @classmethod
    def read_pdb_str(cls, pdb_str): 
        pdb_line = pdb_str.split("\n")
        
        for line in pdb_line:
            line = line.split()
            if line[0] == "CRYST1":
                latt6 = np.array(list(map(float,  line[1:7])))
                break
                
        mol_idx = []
        for line in pdb_line:
            if len(line)>2 and  line[20:22].strip() in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" :
                line = line[:20]+"  " +line[22:]
            line = line.split()
                        
            if "ATOM" in line[0]  or "HETATM" in line[0]:
                tmp_coord = line[5:8]
                coord.append(tmp_coord)
                atom_titles.append(line[2])
                mol_titles.append(line[3])
                mol_idx.append(line[4])
            
          
        mol_idx     = np.array(mol_idx, dtype=np.int32)
        coord       = np.array(coord, dtype=np.float64)
        atom_titles = np.array(atom_titles)
        mol_titles  = np.array(mol_titles)
        latt9       = latt6_to_latt9(latt6)       
        atom_num    = len(coord)
    
        st = cls(latt9, coord, atom_titles, mol_titles, mol_idx)
        
        return(st)
        
    @classmethod
    def read_pdb(cls, pdb_file): 
        pdb_str = open(pdb_file, "r").read()
        return(cls.read_pdb_str(pdb_str))

    @classmethod
    def read(cls, st_file):
      
        file_format  = st_file[-3:]        
        assert file_format in ["g96", "gro"], print(f"{file_format} is not sopported")
        
        if file_format == "g96":
            st = cls.read_96(st_file)            
        elif file_format == "gro":
            st = cls.read_gro(st_file)            
        return(st)

    def write_pdb(self, pdb_name):
        
        
        if isinstance(self.latt6, np.ndarray) :
            pdb_line = ["CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1"  %tuple(self.latt6)]
        else:
            pdb_line = ["CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1"  %tuple([0., 0., 0. ,90, 90, 90])]
        for i in range(len(self.atom_titles)):
            line = "ATOM  %5d %-4s %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00            " \
                    %((i+1)%9999, self.atom_titles[i], self.mol_titles[i], self.mol_idx[i], 
                      self.coord[i][0], self.coord[i][1],   self.coord[i][2])
            pdb_line.append(line)
            
        pdb_line.append("END")
        
        pdb_file = open(pdb_name, "w")
        pdb_file.write("\n".join(pdb_line))
        pdb_file.close()
    
    def write_gro(self, gro_name):
      
        if isinstance(self.latt9, np.ndarray) :
            latt9 = self.latt9 / 10.
            latt_line = ["  %10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f" 
                     %(latt9[0][0],latt9[1][1],latt9[2][2],latt9[0][1],latt9[0][2],
                       latt9[1][0],latt9[1][2],latt9[2][0],latt9[2][1])]
        else:
            latt9 = np.zeros((3,3))
            latt_line = ["  %10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f" 
                     %(latt9[0][0],latt9[1][1],latt9[2][2],latt9[0][1],latt9[0][2],
                       latt9[1][0],latt9[1][2],latt9[2][0],latt9[2][1])]
            
            
        coord = self.coord.T / 10.
        gro_line = ["write by python code", str(len(coord[0]))]
        
    
        coord_line = ["% 5d%-5s%5s% 5s% 8.3f% 8.3f% 8.3f"  %(idx, mol_title, atom_title, (n+1)%9999, x,y,z) 
                      for mol_title, atom_title, idx,n, x,y,z  in zip(self.mol_titles, self.atom_titles, self.mol_idx,
                                                                    np.arange(self.atom_num), coord[0], coord[1], coord[2])]
        
    
        
        gro_line = gro_line + coord_line + latt_line 
    
        gro_file = open(gro_name, "w")
        gro_file.write("\n".join(gro_line))
        gro_file.close()

    def write_g96(self, g96_name):

        latt_line = []
        if isinstance(self.latt9, np.ndarray) :
            latt9 = self.latt9 / 10.
            latt_line = "%15.9f"*9  %(latt9[0][0],latt9[1][1],latt9[2][2],latt9[0][1],latt9[0][2],
                                      latt9[1][0],latt9[1][2],latt9[2][0],latt9[2][1])
            latt_line = ["BOX", latt_line, "END"]


        title_line  = ["TITLE", "write by python code", "END"]
        coord = self.coord.T / 10.
        coord_line = ["% 5d %-6s%-6s% 6s% 15.9f% 15.9f% 15.9f"  %(idx, mol_title, atom_title, (n+1)%9999, x,y,z) 
                      for mol_title, atom_title, idx,n, x,y,z  in zip(self.mol_titles, self.atom_titles, self.mol_idx,
                                                                    np.arange(self.atom_num), coord[0], coord[1], coord[2])]

        coord_line = ["POSITION"] + coord_line + ["END"]

        g96_line = title_line + coord_line + latt_line

        g96_file = open(g96_name, "w")
        g96_file.write("\n".join(g96_line))
        g96_file.close()

    

    def write(self, st_file):
      
        file_format  = st_file[-3:]        
        assert file_format in ["pdb", "gro", "g96"], print(f"{file_format} is not sopported, only pdb/gro supported!")
        
        if file_format == "pdb":
            st = self.write_pdb(st_file)            
        elif file_format == "gro":
            st = self.write_gro(st_file)
        elif file_format == "g96":
            st = self.write_96(st_file)
