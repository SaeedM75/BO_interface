import os
import json
import bayes_opt
import subprocess
import numpy as np
import pandas as pd
import pymatgen as mg
import shutil
import xml.etree.ElementTree as ET
import copy

from ase.io import *
from ase import Atoms, Atom
from ase.calculators.vasp.vasp2 import Vasp2
from ase.dft.kpoints import *

from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar, Poscar
from pymatgen import Lattice, Structure, Molecule
from pymatgen.io.vasp.outputs import BSVasprun, Vasprun
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.surface import Slab, SlabGenerator, ReconstructionGenerator
from pymatgen.analysis.substrate_analyzer import SubstrateAnalyzer, ZSLGenerator
from pymatgen.symmetry.analyzer import *

from bayes_opt import UtilityFunction
from bayes_opt import BayesianOptimization
from string import ascii_lowercase


def interface_shifter(Iface, sub_coords, film_coords, shift_array):
    """
        Generate shifted interface structures
        Parameters:
            Iface Pymatgen Interface structures
            sub_coords 1-D array of substrate atoms coordinates
            film_coords 1-D array of film atoms coordinates
            shift_array list of shift values in x,y,z planes
        Returns
            Pymatgen structure of shifted interfce structure

    """

    film_coords[:, 0] = film_coords[:, 0] + shift_array[0]
    film_coords[:, 1] = film_coords[:, 1] + shift_array[1]
    film_coords[:, 2] = film_coords[:, 2] + shift_array[2]
    int_coords = np.concatenate((film_coords, sub_coords), axis=0)
    new_int = Structure(Iface.lattice.matrix, Iface.species, int_coords, coords_are_cartesian=True)
    new_int_space = SpacegroupAnalyzer(new_int)
    new_int_pri = new_int_space.get_primitive_standard_structure()

    return new_int_pri


class vasp_init(object):
    def __init__(self, input_path, olddir):
        with open(input_path, 'r') as f:
            self.input_dict = json.load(f)
        self.struct_info = self.input_dict['structure_info']
        self.general_flags = self.input_dict['general_flags']
        self.atoms = None
        self.olddir = olddir
        self.substrate = Structure.from_file(self.olddir + "/POSCAR_sub")
        self.film = Structure.from_file(self.olddir + "/POSCAR_film")
        self.Iface = Structure.from_file(self.olddir + "/POSCAR_interface")
    def init_atoms(self):
        if os.path.exists(self.olddir + "/POSCAR_shifted"):
            self.atoms = read(self.olddir + "/POSCAR_shifted")
        else:
            self.atoms = read(self.olddir + "/POSCAR_interface")
        return self.atoms

    def modify_poscar(self, path='./'):
        with open(path + '/POSCAR', 'r') as f:
            poscar = f.readlines()
            ele = []
            for i in list(self.atoms.symbols):
                if i not in ele:
                    ele.append(i)
            poscar.insert(5, ' '.join(x for x in ele ) + '\n')
            poscar[7] = 'cartesian\n'
            f.close()

        with open(path + '/POSCAR','w') as d:
            d.writelines(poscar)
            d.close()

    def generate_input(self, directory, step, xc, import_kpath):
        if "magmom" in self.general_flags.keys():
            add_magmom = True
            magmom_val = self.general_flags['magmom']
            del self.general_flags['magmom']
        else:
            add_magmom = False

        flags = {}
        flags.update(self.general_flags)
        flags.update(self.input_dict[step])
        if step == 'scf':
            if xc == 'pbe':
                flags.update(self.input_dict[xc])
            calc = Vasp2(self.atoms,directory=directory,kpts=self.struct_info['kgrid_'+xc],gamma=True,**flags)
            calc.write_input(self.atoms)
            if add_magmom:
                incar_scf = Incar.from_file(directory + '/INCAR')
                incar_scf['MAGMOM'] = magmom_val
                incar_scf.write_file(directory + '/INCAR')
                
            self.modify_poscar(path=directory)


def readenergy(outcar, poscar):
    outcar_file = read(outcar)
    for l in outcar_file.readlines():
        if "energy without entropy" in l:
            energy = float(l.split()[7])

    return energy

class bayesOpt_DFTU(object):
    def __init__(self, path, kappa= 4):
        self.input = path + 'BO.txt'
        self.energy = readenergy(path + '/dftu/scf/OUTCAR', path + '/dftu/scf/POSCAR')
        self.kappa = kappa
        self.bounds = bounds

    def loss(self, energy):
        return -energy

    def bo(self, x_shift_range = (-bounds[0],bounds[0]), y_shift_range = (-bounds[1],bounds[1]), z_shift_range = (
            -bounds[2],bounds[2]) ):

        data = pd.read_csv(self.input, header=None, delimiter="\s", engine='python')
        num_rows, d = data.shape

        variables_string  = ['dx', 'dy', 'dz']
        pbounds = {}
        pbounds['dx'] = x_shift_range 
        pbounds['dy'] = y_shift_range 
        pbounds['dz'] = z_shift_range 

        utility = UtilityFunction(kind="ucb", kappa=self.kappa, xi=0.0)
        optimizer = BayesianOptimization(
            f=None,
            pbounds=pbounds,
            verbose=2,
        )

        for i in range(num_rows):
            values = list()
            for j in range(3):
                values.append(data.iloc[i][j])
            params = {}
            for (value, variable) in zip(values, variables_string):
                params[variable] = value
            target = self.loss(self.energy)
            optimizer.register(
                params=params,
                target=target,
            )

        next_point_to_probe = optimizer.suggest(utility)
        points = list(next_point_to_probe.values())
        points = [round(elem, 5) for elem in points]

        with open("BO.txt", "a") as res_file:
            print(points[0], points[1], points[2], end = " ", file = res_file)


        shifted_Iface =  interface_shifter(self.Iface, self.substrate.cart_coords, self.film.cart_coords, points)
        shifted_Iface_space = SpacegroupAnalyzer(shifted_Iface)
        shifted_Iface_pri = shifted_Iface_space.get_primitive_standard_structure()
        Poscar(shifted_Iface_pri).write_file("POSCAR_shifted", direct=False)

        return points


def calculate(command, outfilename, method, import_kpath):
    olddir = os.getcwd()
    calc = vasp_init(olddir+'/input.json', olddir)
    calc.init_atoms()
    calc.generate_input(olddir + '/%s/scf' % method, 'scf', 'pbe', import_kpath)

    try:
        os.chdir(olddir+'/%s/scf' %method)
        errorcode_scf = subprocess.call('%s > %s' %(command, outfilename), shell=True)
        if method == 'hse':
            calc.generate_input(olddir+'/%s/band' %method,'band','hse', import_kpath)
    finally:
        os.chdir(olddir)





