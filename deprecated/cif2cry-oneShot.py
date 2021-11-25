# required in working directory: this script, icsd_cif.txt, basis_sets/, icsd_cif/
# Usage example: python cif2cry.py icsd_cif-PX.txt -ubs --basis heyd tzvp --functional HSE06 PBE0-D3
import yaml
import argparse
import os
from pymatgen import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from itertools import combinations


def yaml_basis_set(path_bs):
    dict_bs = {}
    for file_bs in os.listdir(path_bs):
        with open(os.path.join(path_bs, file_bs)) as f:
            bs = f.read()
        ele_str = file_bs.split('_')[0]
        ele = Element(ele_str)
        ele_num = ele.Z
        if ele_num not in dict_bs.keys():
            dict_bs[ele_num] = {}
        bs_type_raw = file_bs.split('_', 1)[1]
        bs_type = ''
        if 'Heyd' in bs_type_raw:
            bs_type = 'heyd'
        elif 'DZVP' in bs_type_raw:
            bs_type = 'dzvp'
        elif 'TZVP' in bs_type_raw:
            bs_type = 'tzvp'
        dict_bs[ele_num][bs_type] = bs
        with open('basis_sets.yaml', 'w') as yaml_file:
            yaml.dump(dict_bs, yaml_file, default_flow_style=False)


def cif_to_crystal(cif_file, dict_bs, basis_type, functional_type):
    s = Structure.from_file("./icsd_cif/" + cif_file)    # added directory for the cif_file
    sga = SpacegroupAnalyzer(s)
    sgn = sga.get_space_group_number()
    cry = ""
    cry += cif_file + "\n"
    cry += "CRYSTAL\n"
    cry += "0 0 0\n"
    cry += str(sgn) + "\n"
    format_lattice = '{:> 3.12f} {:> 3.12f} {:> 3.12f} {:> 3.12f} {:> 3.12f} {:> 3.12f}\n'
    cry += format_lattice.format(s.lattice.abc[0], s.lattice.abc[1], s.lattice.abc[2],
                                 s.lattice.angles[0], s.lattice.angles[1], s.lattice.angles[2])
    cry += str(s.num_sites) + "\n"
    list_atom_num = []
    for i in range(s.num_sites):
        atom_num = s.atomic_numbers[i]
        if atom_num not in list_atom_num:
            list_atom_num.append(atom_num)
        format_coord_block = '{:3d}     {:>.12f} {:>.12f} {:>.12f}\n'
        cry += format_coord_block.format(atom_num, s.frac_coords[i][0], s.frac_coords[i][1], s.frac_coords[i][2])
    cry += "FREQCALC\n"
    cry += "PREOPTGEOM\n"
    cry += "END\n"
    cry += "END\n"
    cry += "END\n"
    for ele_num in list_atom_num:
        cry += dict_bs[ele_num][basis_type]    # check if needs another \n -> No.
    cry += "99 0\n"
    cry += "END\n"
    cry += "DFT\n"
    cry += functional_type + "\n"
    cry += "XXLGRID\nCHUNKS\n200\nEND\n"
    cry += "TOLINTEG\n7 7 7 9 30\nSHRINK\n8 8\nTOLDEE\n11\nEXCHSIZE\n12012489\nEND"
    cry_file = ''.join([s.composition.reduced_formula, "_", str(sgn), ".d12"])
    cry_dir = ''.join(["calc-", cry_file.split('.')[0], "_icsd", cif_file.split('.')[0].split('Code')[1],
                       "_", basis_type, "_", functional_type])
    if cry_dir not in os.listdir():
        os.mkdir(os.path.join(os.getcwd(), cry_dir))
    with open(os.path.join(os.getcwd(), cry_dir, cry_file), 'w') as f:
        f.write(cry)


def read2list(datafile):
    list_out = []
    with open(datafile) as f:
        while True:
            line = f.readline()
            if len(line) == 0:
                break
            list_out.append(line.split())
    return list_out


parser = argparse.ArgumentParser()
parser.add_argument("icsd_data_file", help="txt from collected ICSD-cif excel, e.g., icsd_cif-PX.txt")
parser.add_argument("-bsp", "--basis_sets_path", default="./basis_sets", help="path for basis-set files")
parser.add_argument("-ubs", "--update_basis_sets", help="a flag whether to update basis_sets.yaml", action="store_true")
parser.add_argument("--basis", nargs='*', help="multiple, e.g., heyd dzvp tzvp")
parser.add_argument("--functional", nargs='*', help="multiple, e.g., HSE06 PBE0 PBE0-D3")
args = parser.parse_args()

if args.update_basis_sets:
    yaml_basis_set(args.basis_sets_path)
dict_basis_sets = yaml.safe_load(open('basis_sets.yaml'))

list_icsd = read2list(args.icsd_data_file)
for item in list_icsd:
    icsd = item[0]
    cif_f = ''.join(["ICSD_CollCode", icsd, ".cif"])
    for bst in args.basis:
        for fun in args.functional:
            cif_to_crystal(cif_f, dict_basis_sets, bst, fun)
