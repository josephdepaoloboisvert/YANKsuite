import yaml
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import align
import netCDF4 as netcdf
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
        # self.complex = mda.Universe(self.s_complex[0], self.s_complex[1])
        # self.alchemical = mda.Universe(self.s_complex[0], self.s_complex[1])
        # self.solvent = mda.Universe(self.s_solvent[0], self.s_solvent[1])
        # self.vacuum = mda.Universe(self.s_solvent[0], self.s_solvent[1])

        # Setting up directory that data will be stored in
        for dir in ['trajs', 'overlap_matrices', 'rsmds']:
            if not os.path.exists('./%s/%s'%(self.out_dir, dir)):
                os.mkdir('./%s/%s'%(self.out_dir, dir))


    def extract_trajs(self, state_string, state_index, sel_str_store, memory=False):
        # if memory, this will also return to an MDA universe of the trajectory
        if state_string == 'complex':
            nc = self.complex_nc
        elif state_string == 'solvent':
            nc = self.solvent_nc
        else:
            raise Exception('state_string must be either complex or solvent')

        niterations = nc.variables['positions'].shape[0]
        nstates = nc.variables['states'].shape[1]
        natoms = nc.variables['positions'].shape[2]
        replica_indices = [list(nc.variables['states'][iteration, :]).index(0) for iteration in range(0, niterations)]

        # Align and store snapshots
        temp_u = MDA.universe(self.s_complex[0], self.s_complex[1])
        sel_store = temp_u.select_atoms(sel_str_store)
        store_dcd = './%s/trajs/%s_%s_%s.dcd'%(self.out_dir, state_string, state_index, sel_str_store)
        writer = mda.Writer(store_dcd, sel_store.n_atoms)
        for frame in range(0, len(replica_indices)):
            coords = nc.variables['positions'][frame, replica_indices[frame], :, :] * 10.0
            temp_u.load_new(coords, format=MemoryReader)
            align.alignto(temp_u, self.ref, select='name CA')
            writer.write(sel_store)

        # Write a PDB file with the stored atoms, storing AMBER atom types as a remark
        remarks = []
        remarks.append('REMARK    <-- AMBER ATOM TYPES ')
        for n in range(0, len(sel_store.types), 20):
            remarks.append('REMARK    ' + ' '.join([f'{t:2s}' for t in sel_store.types[n:n + 20]]))
        remarks.append('REMARK    AMBER ATOM TYPES -->')
        remarks = '\n' + '\n'.join(remarks)
        sel_store.write(store_dcd.replace('.dcd','.pdb'), remarks=remarks)


    def graph_overlap_matrix(self, state):
        if state == 'complex':
            nc = self.complex_nc
        elif state == 'solvent':
            nc = self.solvent_nc
        else:
            raise Exception('State must be either complex or solvent')
        proposed = np.array(tuple([np.array(nc.variables['proposed'][i, :, :]) for i in range(nc.variables['positions'].shape[0])]))
        accepted = np.array(tuple([np.array(nc.variables['accepted'][i, :, :]) for i in range(nc.variables['positions'].shape[0])]))
        matrix = np.sum(accepted, axis=0) / np.sum(proposed, axis=0)
        plt.clf()
        plt.colorbar(plt.matshow(matrix))
        plt.xlabel('Replica i')
        plt.ylabel('Replica j')
        plt.savefig('./%s/%s/%s.png'%(self.out_dir, 'overlap_matrices', state))


    def graph_rmsds(self):


def get_from_yaml(name, yaml_contents):
    truth_list = [True if name in line else False for line in yaml_contents]
    if True in truth_list:
        idj = truth_list.index(True)
        idk = yaml_contents[idj].index(':')
        return yaml_contents[idj][idk + 2:]  # Colon space beginseq
    else:
        raise Exception('%s not found in yaml_contents'%name)