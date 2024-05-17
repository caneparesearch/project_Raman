# Crystal output file parser: Structure and Raman intensities 

from ast import literal_eval
import json
import sys
import numpy as np
import pandas as pd
from pymatgen.core import Lattice, Structure
import itertools
import os
import warnings
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from scipy.special import voigt_profile


# future: atomic species oxidation states

def readUntil(f, s):
    """
    f: filename
    s: string to read until
    """
    while True:
        line = f.readline()
        if s in line:
            return line.strip()
        if not line:
            raise ValueError(f"NOT FOUND: {s}")


xyz_to_num = {"X": 0, "Y": 1, "Z": 2}

class crystalOut():
    """
    Crystal output file parser v0.1
    """

    def __init__(self, filename):
        self.filename = filename
        self.filedir = os.path.dirname(filename)
        self.file = open(filename, "r")
        self.parsed_data = {"cell": {"a": 0, "b": 0, "c": 0, "alpha": 0, "beta": 0, "gamma": 0},
                            "crystCell": {"a": 0, "b": 0, "c": 0, "alpha": 0, "beta": 0, "gamma": 0}, "volume": 0,
                            "cryst_volume": 0, "density": 0, "total_energy": 0,
                            "atom_lines": {"all_atom_lines": [], "asymm_atom_lines": [], "cryst_all_atom_lines": [],
                                           "cryst_asymm_atom_lines": []}}

        self.raman_temp, self.raman_wavelength = self.get_raman_exp()
        cell, volume_density, energy, all_atom_lines, asymm_atom_lines, crystCell, crystVolume, cryst_all_atom_lines, cryst_asymm_atom_lines = self.get_cell_params()

        self.parsed_data["cell"]["a"], self.parsed_data["cell"]["b"], self.parsed_data["cell"]["c"], \
        self.parsed_data["cell"]["alpha"], self.parsed_data["cell"]["beta"], self.parsed_data["cell"]["gamma"] = float(
            cell[0]), float(cell[1]), float(cell[2]), float(cell[3]), float(cell[4]), float(cell[5])
        self.parsed_data["volume"], self.parsed_data["density"] = float(volume_density[0]), float(volume_density[1])
        self.parsed_data["total_energy"] = float(energy)
        self.parsed_data["atom_lines"]["all_atom_lines"] = all_atom_lines
        self.parsed_data["atom_lines"]["asymm_atom_lines"] = asymm_atom_lines

        self.atoms, self.coords = self.get_atoms_coords(self.parsed_data["atom_lines"]["all_atom_lines"])
        self.asymm_atoms, self.asymm_coords = self.get_atoms_coords(self.parsed_data["atom_lines"]["asymm_atom_lines"])

        self.lattice = Lattice.from_parameters(a=self.parsed_data["cell"]['a'], b=self.parsed_data["cell"]['b'],
                                               c=self.parsed_data["cell"]['c'],
                                               alpha=self.parsed_data["cell"]['alpha'],
                                               beta=self.parsed_data["cell"]['beta'],
                                               gamma=self.parsed_data["cell"]['gamma'])
        self.structure = Structure(self.lattice, self.atoms, self.coords)
        self.space_group, self.space_group_no = self.get_space_group()

        if crystCell:
            self.parsed_data["crystCell"]["a"], self.parsed_data["crystCell"]["b"], self.parsed_data["crystCell"]["c"], \
            self.parsed_data["crystCell"]["alpha"], self.parsed_data["crystCell"]["beta"], \
            self.parsed_data["crystCell"]["gamma"] = float(crystCell[0]), float(crystCell[1]), float(
                crystCell[2]), float(crystCell[3]), float(crystCell[4]), float(crystCell[5])
            self.parsed_data["cryst_volume"] = float(crystVolume)
            self.parsed_data["atom_lines"]["cryst_all_atom_lines"] = cryst_all_atom_lines
            self.parsed_data["atom_lines"]["cryst_asymm_atom_lines"] = cryst_asymm_atom_lines
            self.cryst_atoms, self.cryst_coords = self.get_atoms_coords(
                self.parsed_data["atom_lines"]["cryst_all_atom_lines"])
            self.cryst_asymm_atoms, self.cryst_asymm_coords = self.get_atoms_coords(
                self.parsed_data["atom_lines"]["cryst_asymm_atom_lines"])

        self.atomic_masses = self.get_atomic_mass()
        self.raman_intensities = self.get_raman_IR_intensities()
        self.dielectric_tensor = self.get_dielectric_tensor()
        self.vib_contributions_dielectric_sum, self.vib_contributions_dielectric = self.get_vibrational_contributions()
        self.second_electric_susceptibility = self.get_second_electric_susceptibiliy()
        self.third_electric_susceptibiliy = self.get_third_electric_susceptibiliy()
        self.born_charge, self.born_charge_normal_mode = self.get_born_charge()
        self.thermodynamic_terms = self.get_thermodynamic_terms()
        self.file.close()

    def get_space_group(self):
        spg = SpacegroupAnalyzer(self.structure)
        return spg.get_space_group_symbol(), spg.get_space_group_number()

    def get_cell_params(self):
        pos = self.file.tell()
        try:
            readUntil(self.file, "FINAL OPTIMIZED GEOMETRY - DIMENSIONALITY OF THE SYSTEM      3")
            optimized = True
            print(">>>>>>> Optimized structure found.")
        except:
            optimized = False
        self.file.seek(pos)

        try:
            readUntil(self.file, "TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL")
            isCrystCell = True
        except ValueError:
            isCrystCell = False
        self.file.seek(pos)

        if optimized:
            readUntil(self.file, "FINAL OPTIMIZED GEOMETRY - DIMENSIONALITY OF THE SYSTEM      3")
            for i in range(3):
                self.file.readline()
            volume_density_line = self.file.readline().split()
        else:
            # read the first structure, if file does not contain optimization sequence
            volume_density_line = readUntil(self.file, "PRIMITIVE CELL -").split()
        volume_density = volume_density_line[7], volume_density_line[10]

        self.file.readline()
        cell = self.file.readline().split()
        if len(cell) != 6:
            print("Weird number of cell parameters. Aborting.")
            sys.exit()

        readUntil(self.file, "*************************************************")
        readUntil(self.file, "*************************************************")

        n = 0
        all_atom_lines = ""
        asymm_atom_lines = ""
        while True:
            line = self.file.readline().split()
            if len(line) != 7:
                break
            n += 1
            atom = line[3][0] + line[3][1:].lower()
            all_atom_lines += atom + str(n) + " " + atom + " " + " ".join(line[4:]) + "\n"
            if line[1] == "T":
                asymm_atom_lines += atom + str(n) + " " + atom + " " + " ".join(line[4:]) + "\n"

        if isCrystCell:
            #print(">>>>>>> CRYSTALLOGRAPHIC CELL found.")
            readUntil(self.file, "*************************************************")
            crystVolume = self.file.readline().split()[3].strip(")")
            self.file.readline()
            crystCell = self.file.readline().split()
            if len(crystCell) != 6:
                print("Weird number of cell parameters. Aborting.")
                sys.exit()
            readUntil(self.file, "*************************************************")

            n = 0
            cryst_all_atom_lines = ""
            cryst_asymm_atom_lines = ""
            while True:
                line = self.file.readline().split()
                if len(line) != 7:
                    break
                n += 1
                atom = line[3][0] + line[3][1:].lower()
                cryst_all_atom_lines += atom + str(n) + " " + atom + " " + " ".join(line[4:]) + "\n"
                if line[1] == "T":
                    cryst_asymm_atom_lines += atom + str(n) + " " + atom + " " + " ".join(line[4:]) + "\n"
        else:
            crystCell, crystVolume, cryst_all_atom_lines, cryst_asymm_atom_lines = None, None, None, None

        energy_line = readUntil(self.file, "TOTAL ENERGY(DFT)(AU)")
        if len(energy_line) > 48:
            energy = energy_line[26:48].strip()
        self.file.seek(0)
        return cell, volume_density, energy, all_atom_lines, asymm_atom_lines, crystCell, crystVolume, cryst_all_atom_lines, cryst_asymm_atom_lines

    def get_raman_exp(self):
        try:
            readUntil(self.file, "RAMANEXP")
            temp, wavelength = self.file.readline().split(", ")
            self.file.seek(0)
            return int(temp), int(wavelength)
        except ValueError:
            warnings.warn("RAMANEXP not found. Using default value of 298K and 514nm")
            self.file.seek(0)
            return 298, 514

    def get_atoms_coords(self, atom_lines):
        atom_lines_ = atom_lines.split("\n")
        all_atoms = []
        all_coords = []
        for line in atom_lines_:
            if line:
                s = line.split(" ")
                atom = s[1]
                coords = []
                for coord in s[2:5]:
                    coords.append(float(coord))
                all_atoms.append(str(atom))
                all_coords.append(coords)
        return all_atoms, all_coords

    def get_raman_IR_intensities(self):
        """
        returns the intensities as a dict with format frequency: intensity
        only polycrystalline isotropic intensities (I_tot only) and the corresponding IRREP
        """
        readUntil(self.file, " CONVERSION FACTORS FOR FREQUENCIES:")
        for _ in range(8):
            line = self.file.readline()
        data = []
        while not line.isspace():
            line = line.replace("("," ")
            line = line.replace(")"," ")
            fields = line.split()
            data.append(fields[1:])
            line = self.file.readline()
        
        df = pd.DataFrame(data, columns=["Mode","EIGV(Ha**2)","FREQ(CM**-1)","FREQ(THZ)","IRREP","IR","INTENS","RAMAN"])
        df.set_index("Mode", inplace=True)
        df["RAMAN INTENSITIES"] = np.zeros(len(df))

        readUntil(self.file, "POLYCRYSTALLINE ISOTROPIC INTENSITIES (ARBITRARY UNITS)")
        for i in range(4):
            line = self.file.readline()

        while not line.isspace():
            fields = line.split()
            tot_intensity = float(fields[5])
            mode = fields[1]
            df.loc[mode, "RAMAN INTENSITIES"] = tot_intensity
            line = self.file.readline()

        self.file.seek(0)
        return df

    def get_atomic_mass(self):
        readUntil(self.file, "ATOMS ISOTOPIC MASS (AMU) FOR FREQUENCY CALCULATION")
        line = self.file.readline()
        line = self.file.readline()
        atomic_masses = []
        while not line.isspace():
            fields = line.split()
            i = 0
            while i < len(fields):
                atomic_masses.append((fields[i + 1], float(fields[i + 2])))
                i += 3
            line = self.file.readline()
        self.file.seek(0)
        return atomic_masses

    def get_dielectric_tensor(self):
        dielectric_tensor = np.zeros([3, 3])
        readUntil(self.file, "POLARIZABILITY (ALPHA), DIELECTRIC (EPSILON) AND FIRST-ORDER ELECTRIC")
        for i in range(3):
            line = self.file.readline()
        for i in range(6):
            line = self.file.readline().split()
            if len(line) != 7:
                if any(char.isdigit() for char in line[1]): # bunched up numbers in second field
                    line.insert(1, line[1][0:2])
                    line[2] = line[2][2:]
                else:
                    raise ValueError("Error in Dielectric Tensor")
            perm = list(itertools.permutations(line[1]))
            for p in perm:
                index1 = xyz_to_num[p[0]]
                index2 = xyz_to_num[p[1]]
                dielectric_tensor[index1, index2] = float(line[4])
        self.file.seek(0)
        return dielectric_tensor

    def get_vibrational_contributions(self):
        dielectric_sum = np.zeros([3, 3])
        vibrational_contributions = {}
        readUntil(self.file, "VIBRATIONAL CONTRIBUTIONS TO THE STATIC DIELECTRIC TENSOR (OSCILLATOR")
        for i in range(18):
            line = self.file.readline()
        while "SUM TENSOR OF THE VIBRATIONAL CONTRIBUTIONS" not in line:
            tensor = np.zeros([3, 3])
            line = line.split()
            mode = int(line[0])
            tensor[0, :] = [float(i) for i in line[1:]]
            for i in range(1, 3):
                line = self.file.readline().split()
                if len(line) != 3: # some numbers bunched together
                    string = "".join(line)
                    line = []
                    start = 0
                    for ind, c in enumerate(string):
                        if c == ".":
                            end = ind+6 if ind+6 < len(string) else None                               
                            line.append(string[start : end])
                            start = end
                tensor[i, :] = [float(i) for i in line]
            vibrational_contributions[mode] = tensor
            line = self.file.readline()
            line = self.file.readline()
        line = self.file.readline()
        for i in range(3):
            line = self.file.readline().split()
            if len(line) != 3: # some numbers bunched together
                string = "".join(line)
                line = []
                start = 0
                for ind, c in enumerate(string):
                    if c == ".":
                        end = ind+6 if ind+6 < len(string) else None                               
                        line.append(string[start : end])
                        start = end
            dielectric_sum[i, :] = [float(i) for i in line]
        self.file.seek(0)
        return dielectric_sum, vibrational_contributions

    def get_second_electric_susceptibiliy(self):
        tensor = np.zeros([3, 3, 3])
        readUntil(self.file, "FIRST HYPERPOLARIZABILITY (BETA) AND SECOND ELECTRIC SUSCEPTIBILITY (CHI(2))")
        for i in range(5):
            self.file.readline()
        for i in range(10):
            line = self.file.readline().split()
            perm = list(itertools.permutations(line[0]))
            for p in perm:
                index1 = xyz_to_num[p[0]]
                index2 = xyz_to_num[p[1]]
                index3 = xyz_to_num[p[2]]
                tensor[index1, index2, index3] = float(line[4])
        self.file.seek(0)
        return tensor

    def get_third_electric_susceptibiliy(self):
        tensor = np.zeros([3, 3, 3, 3])
        try:
            readUntil(self.file, "SECOND HYPERPOLARIZABILITY (GAMMA) AND THIRD ELECTRIC SUSCEPTIBILITY (CHI(3))")
        except:
            self.file.seek(0)
            return None
        for i in range(5):
            self.file.readline()
        for i in range(15):
            line = self.file.readline().split()
            if len(line) == 4:  # line has 6 fields
                line.insert(0, line[0][0:4])  # splits the component and gamma fields
                line[1] = line[1][5:]
            elif len(line) == 3:
                line.insert(0, line[0][0:4])  # splits the component and gamma fields
                line[1] = line[1][5:]
                gamma_chi = line.pop(1).split(")")
                for i in range(len(gamma_chi)):
                    if "*" in gamma_chi[i]:
                        warnings.warn("Warning! SECOND HYPERPOLARIZABILITY computations incomplete. Please check.")
                        line.insert(i + 1, 0)
                    else:
                        line.insert(i + 1, gamma_chi[i])
            perm = list(itertools.permutations(line[0]))
            for p in perm:
                index1 = xyz_to_num[p[0]]
                index2 = xyz_to_num[p[1]]
                index3 = xyz_to_num[p[2]]
                index4 = xyz_to_num[p[3]]
                tensor[index1, index2, index3, index4] = float(line[4])
        self.file.seek(0)
        return tensor

    def get_born_charge(self):
        born_charge = []
        readUntil(self.file, "ATOMIC BORN CHARGE TENSOR (UNITS OF e, ELECTRON CHARGE).")
        for i in range(3):
            line = self.file.readline()
        while "DYNAMIC CHARGE" in line:
            line = line.split()
            atom_sym = line[2]
            dynamic_charge = float(line[5])
            tensor = np.zeros([3, 3])
            self.file.readline()
            self.file.readline()
            for i in range(3):
                line = self.file.readline().split()
                tensor[i, :] = [float(i) for i in line[1:]]
            atom_dict = {}
            atom_dict["Atomic symbol"] = atom_sym
            atom_dict["Dynamic Charge"] = dynamic_charge
            atom_dict["Born Charge"] = tensor
            born_charge.append(atom_dict)
            line = self.file.readline()
            line = self.file.readline()
        readUntil(self.file, "BORN CHARGE VECTOR IN THE BASIS OF NORMAL MODES")
        normal_basis = []
        for i in range(4):
            line = self.file.readline()
        while not line.isspace():
            line = line.split()
            normal_basis.append([float(i) for i in line[1:]])
            line = self.file.readline()
        normal_basis = np.array(normal_basis)
        self.file.seek(0)
        return born_charge, normal_basis

    def get_thermodynamic_terms(self):
        thermo = {}
        readUntil(self.file, "HARMONIC VIBRATIONAL CONTRIBUTIONS TO THERMODYNAMIC FUNCTIONS AT GIVEN")
        for i in range(9):
            self.file.readline()
        line = self.file.readline().split()
        thermo["EL"] = float(line[3])
        line = self.file.readline().split()
        thermo["ZPE"] = float(line[3])
        line = readUntil(self.file, "PV").split()
        thermo["PV"] = float(line[3])
        line = self.file.readline().split()
        thermo["TS"] = float(line[3])
        line = readUntil(self.file, "HEAT CAPACITY").split()
        thermo["HEAT CAPACITY"] = float(line[4])
        return thermo

    def get_IR_tensor(self):
        """
        returns the IR tensor of size 3n x 3 in a numpy array. The file TENS_IR.DAT must be in the same 
        directory as the CRYSTAL .out file
        """
        IR_file = self.filedir + "\\TENS_IR.DAT"
        if os.path.isfile(IR_file):
            IR_tensor = np.loadtxt(IR_file)
            return IR_tensor
        else:
            raise FileNotFoundError("Please place TENS_IR.DAT file in the same directory as the CRYSTAL .out file")

    def get_raman_tensor(self):
        """
        returns the Raman tensor of size 3n x 6 in a numpy array. The file TENS_RAMAN.DAT must be in the same 
        directory as the CRYSTAL .out file
        """
        Raman_file = self.filedir + "\\TENS_RAMAN.DAT"
        if os.path.isfile(Raman_file):
            Raman_tensor = np.loadtxt(Raman_file)
            return Raman_tensor
        else:
            raise FileNotFoundError("Please place TENS_RAMAN.DAT file in the same directory as the CRYSTAL .out file")

    def get_dynamical_matrix(self):
        """
        returns the dynamical matrix of size 3n x 3n in a numpy array. The file HESSFREQ.DAT must be in the same 
        directory as the CRYSTAL .out file
        """
        hessfreq_file = self.filedir + "\\HESSFREQ.DAT"
        no_atoms = len(self.atoms)
        if os.path.isfile(hessfreq_file):
            dyn_matrix = np.loadtxt(hessfreq_file)
            dyn_matrix = np.reshape(dyn_matrix.T, (3 * no_atoms, 3 * no_atoms))
            return dyn_matrix
        else:
            raise FileNotFoundError("Please place HESSFREQ.DAT file in the same directory as the CRYSTAL .out file")

    def get_convoluted_spectra(self, sigma, gamma, wavenumber_range=(0,1000), resolution=1000):
        """
        Returns: 
            Frequencies (np.array), Intensities (np.array) of convolution of Raman spectra intensities calculated.
            Convolution is based on Voigt profile.

        Parameters:
            sigma (float): standard deviation of Gaussian distribution. If sigma = 0, a Cauchy distribution is used.
            gamma (float): half-width at half-maximum of Cauchy distribution. If gamma = 0, a Gaussian distribution is used.
            padding (int): Frequency range added to the max and min of the frequencies to display the full spectra. Default=30.
            resolution (int): Number of points to compute between the max and min. Default=1000.
        """
        data = np.array(list(self.intRaman.items()))
        frequencies = np.linspace(wavenumber_range[0], wavenumber_range[1], resolution)
        convoluted_intensities = np.zeros((resolution))
        for i in range(data.shape[0]):
            freq = data[i, 0]
            intensity = data[i, 1]
            convolution_single_peak = []
            for x in frequencies:
                convoluted_amplitude = voigt_profile(x - freq, sigma, gamma) * intensity
                convolution_single_peak.append(convoluted_amplitude)
            convoluted_intensities += np.array(convolution_single_peak)

        return frequencies, convoluted_intensities


def return_structure_data(df2, idx):
    """
    returns csv data generated from crystalOutCSV back to dict 
    """
    structure_db_dict = {"structure": Structure.from_dict(json.loads(df2["structure"][idx])),
                         "spaceGroup": df2["spaceGroup"][idx],
                         "thermodynamicTerms": json.loads(df2["thermodynamicTerms"][idx]), "dielectricTensor": [],
                         "vibContributionsDielectric": [], "secondElectricSusceptibility": [],
                         "thirdElectricSusceptibility": [], "bornChargeArray": [], "bornChargeNormalModeBasis": [],
                         "intRaman": json.loads(df2["intRaman"][idx])}

    dielectricTensor = np.frombuffer(
        df2["dielectricTensor_flattened"][idx].encode().decode('unicode-escape').encode('ISO-8859-1')[2:-1],
        dtype=np.float64)
    dielectricTensor = dielectricTensor.reshape(literal_eval(df2["dielectricTensor_shape"][idx]))
    structure_db_dict["dielectricTensor"] = dielectricTensor

    vibContributionsDielectric = np.frombuffer(
        df2["vibContributionsDielectric_flattened"][idx].encode().decode('unicode-escape').encode('ISO-8859-1')[2:-1],
        dtype=np.float64)
    vibContributionsDielectric = vibContributionsDielectric.reshape(
        literal_eval(df2["vibContributionsDielectric_shape"][idx]))
    structure_db_dict["vibContributionsDielectric"] = vibContributionsDielectric

    secondElectricSusceptibility = np.frombuffer(
        df2["secondElectricSusceptibility_flattened"][idx].encode().decode('unicode-escape').encode('ISO-8859-1')[2:-1],
        dtype=np.float64)
    secondElectricSusceptibility = secondElectricSusceptibility.reshape(
        literal_eval(df2["secondElectricSusceptibility_shape"][idx]))
    structure_db_dict["secondElectricSusceptibility"] = secondElectricSusceptibility

    thirdElectricSusceptibility = np.frombuffer(
        df2["thirdElectricSusceptibility_flattened"][idx].encode().decode('unicode-escape').encode('ISO-8859-1')[2:-1],
        dtype=np.float64)
    thirdElectricSusceptibility = thirdElectricSusceptibility.reshape(
        literal_eval(df2["thirdElectricSusceptibility_shape"][idx]))
    structure_db_dict["thirdElectricSusceptibility"] = thirdElectricSusceptibility

    bornChargeArray = np.frombuffer(
        df2["bornChargeArray_flattened"][idx].encode().decode('unicode-escape').encode('ISO-8859-1')[2:-1],
        dtype=np.float64)
    bornChargeArray = bornChargeArray.reshape(literal_eval(df2["bornChargeArray_shape"][idx]))
    structure_db_dict["bornChargeArray"] = bornChargeArray

    bornChargeNormalModeBasis = np.frombuffer(
        df2["bornChargeNormalModeBasis_flattened"][idx].encode().decode('unicode-escape').encode('ISO-8859-1')[2:-1],
        dtype=np.float64)
    bornChargeNormalModeBasis = bornChargeNormalModeBasis.reshape(
        literal_eval(df2["bornChargeNormalModeBasis_shape"][idx]))
    structure_db_dict["bornChargeNormalModeBasis"] = bornChargeNormalModeBasis

    return structure_db_dict
