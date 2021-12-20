import yaml
import MDAnalysis as mda
import netCDF4 as netcdf
import os
import sys

from yank.experiment import *

class YankAnalyzer:
    def __init__(self, full_path_to_yaml):
        #Get to directory that contains simulation
        idx = [index for index, element in enumerate(full_path_to_yaml) if element == '/']
        self.yank_dir = full_path_to_yaml[:idx[-1]+1]  # Trailing Slash included
        this_dir = os.getcwd()
        os.chdir(self.yank_dir)

        #Obtain Contents of generate yaml file ($OUTPUT/experiments/experiments.yaml)
        f = open(full_path_to_yaml[idx[-1] + 1:], 'r')
        yaml_script_contents = [line[:-2] if line.endswith('\n') else line for line in f.readlines()]
        f.close()
        try:
            self.out_dir = get_from_yaml('output_dir:', yaml_script_contents)
        except:
            self.out_dir = 'output'
        f = open(self.yank_dir + '%s/experiments/experiments.yaml'%self.out_dir, 'r')
        self.yaml_contents = [line[:-2] if line.endswith('\n') else line for line in f.readlines()]
        f.close()

        #Retrieve other information as needed
        try:
            self.dsl_lig = get_from_yaml('ligand_dsl:', self.yaml_contents)
        except:
            self.dsl_lig = 'resname LIG'
        try:
            self.s_complex = list(get_from_yaml('phase1_path', self.yaml_contents))
        except:
            self.s_complex = ['complex.prmtop','complex.inpcrd']
        try:
            self.s_solvent = list(get_from_yaml('phase2_path', self.yaml_contents))
        except:
            self.s_solvent = ['solvent.prmtop', 'solvent.inpcrd']

        self.yaml = YankLoader(self.yank_dir + '%s/experiments/experiments.yaml'%self.out_dir)  # Probably
        self.complex_nc = netcdf.Dataset(self.yank_dir + '%s/experiments/complex.nc'%self.out_dir)  # Probably
        self.solvent_nc = netcdf.Dataset(self.yank_dir + '%s/experiments/solvent.nc' % self.out_dir)  # Probably

        #MDA Universes reference (input files), complex(complex.nc -> state0), alchemical(complex.nc -> state-1),
        # solvent(solvent.nc -> state0), vacuum(solvent.nc -> state-1)
        self.ref = mda.Universe(self.s_complex[0], self.s_complex[1])
        self.complex = mda.Universe(self.s_complex[0], self.s_complex[1])
        self.alchemical = mda.Universe(self.s_complex[0], self.s_complex[1])
        self.solvent = mda.Universe(self.s_solvent[0], self.s_solvent[1])
        self.vacuum = mda.Universe(self.s_solvent[0], self.s_solvent[1])


        self.ref = mda.Universe()
        os.chdir(this_dir)


    def graph_alchemical_paths(self):


    def graph_overlap_matrix(self):


    def graph_rmsds(self):


def get_from_yaml(name, yaml_contents):
    truth_list = [True if name in line else False for line in yaml_contents]
    if True in truth_list:
        idj = truth_list.index(True)
        idk = yaml_contents[idj].index(':')
        return yaml_contents[idj][idk + 2:]  # Colon space beginseq
    else:
        raise Exception('%s not found in yaml_contents'%name)