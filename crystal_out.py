# Crystal output file parser: Structure and Raman intensities 

import glob
import sys
from pymatgen.core import Lattice, Structure

# future: atomic species oxidation states 

def readUntil(f, s):
    """
    f: filename
    s: string to read until
    """
    while True:
        line = f.readline()
        if len(line) == 0:
            return ""
        if s in line:
            return line.strip()

spacegroups = {
    "P 1": 1, "P -1": 2, "P 2": 3, "P 21": 4, "C 2": 5, "P M": 6, "P C": 7,
    "C M": 8, "C C": 9, "P 2/M": 10, "P 21/M": 11, "C 2/M": 12, "P 2/C": 13,
    "P 21/C": 14, "C 2/C": 15, "P 2 2 2": 16, "P 2 2 21": 17,
    "P 21 21 2": 18, "P 21 21 21": 19, "C 2 2 21": 20, "C 2 2 2": 21,
    "F 2 2 2": 22, "I 2 2 2": 23, "I 21 21 21": 24, "P M M 2": 25,
    "P M C 21": 26, "P C C 2": 27, "P M A 2": 28, "P C A 21": 29,
    "P N C 2": 30, "P M N 21": 31, "P B A 2": 32, "P N A 21": 33,
    "P N N 2": 34, "C M M 2": 35, "C M C 21": 36, "C C C 2": 37,
    "A M M 2": 38, "A B M 2": 39, "A M A 2": 40, "A B A 2": 41,
    "F M M 2": 42, "F D D 2": 43, "I M M 2": 44, "I B A 2": 45,
    "I M A 2": 46, "P M M M": 47, "P N N N": 48, "P C C M": 49,
    "P B A N": 50, "P M M A": 51, "P N N A": 52, "P M N A": 53,
    "P C C A": 54, "P B A M": 55, "P C C N": 56, "P B C M": 57,
    "P N N M": 58, "P M M N": 59, "P B C N": 60, "P B C A": 61,
    "P N M A": 62, "C M C M": 63, "C M C A": 64, "C M M M": 65,
    "C C C M": 66, "C M M A": 67, "C C C A": 68, "F M M M": 69,
    "F D D D": 70, "I M M M": 71, "I B A M": 72, "I B C A": 73,
    "I M M A": 74, "P 4": 75, "P 41": 76, "P 42": 77, "P 43": 78, "I 4": 79,
    "I 41": 80, "P -4": 81, "I -4": 82, "P 4/M": 83, "P 42/M": 84,
    "P 4/N": 85, "P 42/N": 86, "I 4/M": 87, "I 41/A": 88, "P 4 2 2": 89,
    "P 4 21 2": 90, "P 41 2 2": 91, "P 41 21 2": 92, "P 42 2 2": 93,
    "P 42 21 2": 94, "P 43 2 2": 95, "P 43 21 2": 96, "I 4 2 2": 97,
    "I 41 2 2": 98, "P 4 M M": 99, "P 4 B M": 100, "P 42 C M": 101,
    "P 42 N M": 102, "P 4 C C": 103, "P 4 N C": 104, "P 42 M C": 105,
    "P 42 B C": 106, "I 4 M M": 107, "I 4 C M": 108, "I 41 M D": 109,
    "I 41 C D": 110, "P -4 2 M": 111, "P -4 2 C": 112, "P -4 21 M": 113,
    "P -4 21 C": 114, "P -4 M 2": 115, "P -4 C 2": 116, "P -4 B 2": 117,
    "P -4 N 2": 118, "I -4 M 2": 119, "I -4 C 2": 120, "I -4 2 M": 121,
    "I -4 2 D": 122, "P 4/M M M": 123, "P 4/M C C": 124, "P 4/N B M": 125,
    "P 4/N N C": 126, "P 4/M B M": 127, "P 4/M N C": 128, "P 4/N M M": 129,
    "P 4/N C C": 130, "P 42/M M C": 131, "P 42/M C M": 132,
    "P 42/N B C": 133, "P 42/N N M": 134, "P 42/M B C": 135,
    "P 42/M N M": 136, "P 42/N M C": 137, "P 42/N C M": 138,
    "I 4/M M M": 139, "I 4/M C M": 140, "I 41/A M D": 141, "I 41/A C D": 142,
    "P 3": 143, "P 31": 144, "P 32": 145, "R 3": 146, "P -3": 147,
    "R -3": 148, "P 3 1 2": 149, "P 3 2 1": 150, "P 31 1 2": 151,
    "P 31 2 1": 152, "P 32 1 2": 153, "P 32 2 1": 154, "R 3 2": 155,
    "P 3 M 1": 156, "P 3 1 M": 157, "P 3 C 1": 158, "P 3 1 C": 159,
    "R 3 M": 160, "R 3 C": 161, "P -3 1 M": 162, "P -3 1 C": 163,
    "P -3 M 1": 164, "P -3 C 1": 165, "R -3 M": 166, "R -3 C": 167,
    "P 6": 168, "P 61": 169, "P 65": 170, "P 62": 171, "P 64": 172,
    "P 63": 173, "P -6": 174, "P 6/M": 175, "P 63/M": 176, "P 6 2 2": 177,
    "P 61 2 2": 178, "P 65 2 2": 179, "P 62 2 2": 180, "P 64 2 2": 181,
    "P 63 2 2": 182, "P 6 M M": 183, "P 6 C C": 184, "P 63 C M": 185,
    "P 63 M C": 186, "P -6 M 2": 187, "P -6 C 2": 188, "P -6 2 M": 189,
    "P -6 2 C": 190, "P 6/M M M": 191, "P 6/M C C": 192, "P 63/M C M": 193,
    "P 63/M M C": 194, "P 2 3": 195, "F 2 3": 196, "I 2 3": 197,
    "P 21 3": 198, "I 21 3": 199, "P M 3": 200, "P N 3": 201, "F M 3": 202,
    "F D 3": 203, "I M 3": 204, "P A 3": 205, "I A 3": 206, "P 4 3 2": 207,
    "P 42 3 2": 208, "F 4 3 2": 209, "F 41 3 2": 210, "I 4 3 2": 211,
    "P 43 3 2": 212, "P 41 3 2": 213, "I 41 3 2": 214, "P -4 3 M": 215,
    "F -4 3 M": 216, "I -4 3 M": 217, "P -4 3 N": 218, "F -4 3 C": 219,
    "I -4 3 D": 220, "P M 3 M": 221, "P N 3 N": 222, "P M 3 N": 223,
    "P N 3 M": 224, "F M 3 M": 225, "F M 3 C": 226, "F D 3 M": 227,
    "F D 3 C": 228, "I M 3 M": 229, "I A 3 D": 230,
}

class crystalOut():
    """
    Crystal output file parser v0.1
    """

    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename, "r")
        self.parsed_data = {
            "space_group": "", "space_group_num" : 0,
            "a": 0,"b": 0,"c": 0,
            "alpha": 0,"beta" : 0,"gamma": 0,
            "volume": 0,"density": 0,
            "total_energy" : 0,
            "atom_lines": "",
            "intensities": {
            "polycrystalline_intensities": {}
            }
        }

        self.parsed_data["space_group"], self.parsed_data["space_group_num"]  = self.get_space_group(self.file)
        cell, volume_density, energy, atom_lines = self.get_cell_params(self.file)
        self.parsed_data["a"], self.parsed_data["b"], self.parsed_data["c"], self.parsed_data["alpha"], self.parsed_data["beta"], self.parsed_data["gamma"] = float(cell[0]), float(cell[1]), float(cell[2]), float(cell[3]), float(cell[4]), float(cell[5])
        self.parsed_data["volume"], self.parsed_data["density"] = float(volume_density[0]), float(volume_density[1]) 
        self.parsed_data["total_energy"] = float(energy)
        self.parsed_data["atom_lines"] = atom_lines
        self.parsed_data["intensities"]["polycrystalline_intensities"] = self.get_intensities(self.file)
        
        self.atoms, self.coords = self.get_atoms_coords()
        self.lattice = Lattice.from_parameters(a=self.parsed_data['a'], b=self.parsed_data['b'], c=self.parsed_data['c'], 
                                        alpha =self.parsed_data['alpha'], 
                                        beta=self.parsed_data['beta'], 
                                        gamma=self.parsed_data['gamma'])
        self.structure = Structure(self.lattice, self.atoms, self.coords)
        self.file.close()

    def get_space_group(self, file):
        sg = readUntil(self.file, "SPACE GROUP")
        if len(sg) > 0:
            sg_name = sg.split(":")[-1].strip()
            sg_num = spacegroups[sg_name]
        else:
            self.file.seek(0)
            sg = readUntil(self.file, "SPACE  GROUP")
            if len(sg) > 0:
                sg_num = int(sg.split(":")[-1].strip())
                for i, j in spacegroups.items():
                    if sg_num == j:
                        sg_name = i
            else:
                print("Space group could not be found. Assuming P1. Please check.")
                sg_name = "P 1"
                sg_num = 1
                self.file.seek(0)
        return sg_name, sg_num

    def get_cell_params(self, file):
        pos = self.file.tell()
        if len(readUntil(self.file, "TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL")) > 0:
            crystCell = True
        else:
            crystCell = False
        self.file.seek(pos)
        
        counter = 1
        while True:
            if crystCell:
                line = readUntil(self.file, "CRYSTALLOGRAPHIC CELL (VOLUME=")
                volume_density = 0, 0 # placeholder (not dealt with this type of output yet)
            else:
                line = readUntil(self.file, "PRIMITIVE CELL -")
                volume_density = line.split()[7], line.split()[10]
            if len(line) == 0:
                break
                
            self.file.readline()
            cell = self.file.readline().split()
            if len(cell) != 6:
                print("Weird number of cell parameters. Aborting.")
                sys.exit()
                
            readUntil(self.file, "*************************************************")
            if not crystCell:
                readUntil(self.file, "*************************************************")
            
            n = 0
            atom_lines = ""
            while True:
                line = self.file.readline().split()
                if len(line) != 7:
                    break
                if line[1] != "T":
                    continue
                n = n + 1
                atom = line[3][0] + line[3][1:].lower()
                atom_lines += atom + str(n) + " " + atom + " " + " ".join(line[4:]) + "\n"
            
            energy_line = readUntil(self.file, "TOTAL ENERGY(DFT)(AU)")
            if len(energy_line) > 48:
                energy = energy_line[26:48].strip()
            return cell, volume_density, energy, atom_lines

    def get_atoms_coords(self):
        atom_lines = self.parsed_data["atom_lines"].split("\n") 
        all_atoms = []
        all_coords = []
        for line in atom_lines:
            if line:
                s = line.split(" ")
                atom = s[1]
                coords = []
                for coord in s[2:5]:
                    coords.append(float(coord))
                all_atoms.append(str(atom))
                all_coords.append(coords)
        return all_atoms, all_coords

    def get_intensities(self, file):
        """
        only polycrystalline isotropic intensities (I_tot only), 
        TODO: add others and add single crystal 
        """
        readUntil(self.file, "POLYCRYSTALLINE ISOTROPIC INTENSITIES (ARBITRARY UNITS)")
        for i in range(4):
            line = self.file.readline()

        data = {}
        while not line.isspace():
            fields = line.split()
            data[float(fields[2])] = float(fields[5])
            line = self.file.readline()
        
        return data 
