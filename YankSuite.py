import yaml
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import align
import netCDF4 as netcdf
import sys
from yank.experiment import *

class YankMultiAnalyzer():
    """A Class for comparing multiple YankAnalyzers"""
    def __init__(self, in_dictionary):
        """Format of in_dictionary is keys = string name; value = list of YankAnalyzer objects that are repeats of each other
        EX: Protein_in_dict = {'LIG1': [YankAnalyzerObjectFromRepeat1, YankAnalyzerObjectFromRepeat2, YankAnalyzerObjectFromRepeat3], etc... }
        Every value-list should have at least two entries and the overall dictionary should have at least one entry."""


    def convergence_of_rpts(self):
        """Analyze the convergence of repeats (see init for how repeats are defined)"""

class YankAnalyzer():
    def __init__(self, full_path_to_yaml):
        """
        Python class for managing the information exracted from YANK simulations

        PARAMETERS:
        full_path_to_yaml: full path to the yaml file which ran the simulation type(str)

        RETURNS:
            does not return
        """
        #Get to directory that contains simulation
        idx = [index for index, element in enumerate(full_path_to_yaml) if element == '/']
        self.yank_dir = full_path_to_yaml[:idx[-1]+1]  # Trailing Slash included
        f = open(full_path_to_yaml, 'r')
        yaml_script_contents = [line[:-2] if line.endswith('\n') else line for line in f.readlines()]
        f.close()

        # Obtain Contents of generate yaml file ($OUTPUT/experiments/experiments.yaml)
        try:
            self.out_dir = self.yank_dir + get_from_yaml('output_dir:', yaml_script_contents)
        except:
            self.out_dir = self.yank_dir + 'output'

        f = open(self.out_dir + '/experiments/experiments.yaml', 'r')
        self.yaml_contents = [line[:-2] if line.endswith('\n') else line for line in f.readlines()]
        f.close()

        #Retrieve other information as needed
        try:
            self.dsl_lig = get_from_yaml('ligand_dsl:', self.yaml_contents)
        except:
            self.dsl_lig = 'resname LIG'
        try:
            self.s_phase1 = [self.yank_dir + element for element in get_from_yaml('phase1_path', self.yaml_contents).strip('][').split(', ')]
        except:
            self.s_phase1 = [self.yank_dir + element for element in ['complex.prmtop','complex.inpcrd']]
        try:
            self.s_phase2 = [self.yank_dir + element for element in get_from_yaml('phase2_path', self.yaml_contents).strip('][').split(', ')]
        except:
            self.s_phase2 = [self.yank_dir + element for element in ['solvent.prmtop','solvent.inpcrd']]

        # Setting up directory that data will be stored in
        for dir in ['trajs', 'overlap_matrices', 'rsmds']:
            if not os.path.exists('%s/%s'%(self.out_dir, dir)):
                os.mkdir('%s/%s'%(self.out_dir, dir))

        self.yaml = YankLoader(self.out_dir + '/experiments/experiments.yaml')  # Probably
        self.phase1_nc = netcdf.Dataset(self.out_dir + '/experiments/complex.nc')  # Probably
        self.phase2_nc = netcdf.Dataset(self.out_dir + '/experiments/solvent.nc')  # Probably

        #MDA Universes reference (input files), complex(complex.nc -> state0), alchemical(complex.nc -> state-1), solvent(solvent.nc -> state0), vacuum(solvent.nc -> state-1)
        self.ref = mda.Universe(self.yank_dir+self.s_phase1[0], self.s_phase2[1])

        # extract trajectories
        self.extract_traj('complex', 0, 'all')
        self.extract_traj('complex', -1, 'protein')
        self.extract_traj('complex', -1, self.dsl_lig)
        self.extract_traj('solvent', 0, 'all')
        self.extract_traj('solvent', -1, self.dsl_lig)

    def extract_traj(self, state_string, state_index, sel_str_store, sel_str_align):
        """
        get a thermodynamic state of the yank simulation
        :param state_string: must be either the string 'complex' or 'solvent'
        :param state_index: replica/state index to extract (typically 0 or -1)
        :param sel_str_store: MDAnalysis selection string to store trajectory of
        :param memory: if memory, the extracted trajectory will be returned as an MDAnalysis Universe Object

        extract_trajs('complex', 0, 'all')
        extract_trajs('complex', -1, 'protein')
        extract_trajs('complex', -1, self.dsl_lig)
        extract_trajs('solvent', 0, 'all')
        extract_trajs('solvent', -1, self.dsl_lig)
        """
        if state_string == 'complex':
            nc = self.phase1_nc
            temp_u = mda.Universe(self.s_phase1[0], self.s_phase1[1])
        elif state_string == 'solvent':
            nc = self.phase2_nc
            temp_u = mda.Universe(self.s_phase2[0], self.s_phase2[1])
        else:
            raise Exception('state_string must be either complex or solvent')

        niterations = nc.variables['positions'].shape[0]
        nstates = nc.variables['states'].shape[1]
        natoms = nc.variables['positions'].shape[2]


        i_state = state_index
        if state_index == -1:
            i_state = nstates-1

        replica_indices = [list(nc.variables['states'][iteration, :]).index(i_state) for iteration in range(0, niterations)]

        # Align and store snapshots
        sel_store = temp_u.select_atoms(sel_str_store)
        store_dcd = '%s/trajs/%s_%s_%s.dcd'%(self.out_dir, state_string, i_state, sel_str_store)
        writer = mda.Writer(store_dcd, sel_store.n_atoms)
        for frame in range(0, len(replica_indices)):
            coords = nc.variables['positions'][frame, replica_indices[frame], :, :] * 10.0
            temp_u.load_new(coords, format=MemoryReader)
            align.alignto(temp_u, temp_u, select=sel_str_align)
            writer.write(sel_store)

        # Write a PDB file with the stored atoms, storing AMBER atom types as a remark
        remarks = []
        remarks.append('REMARK    <-- AMBER ATOM TYPES ')
        for n in range(0, len(sel_store.types), 20):
            remarks.append('REMARK    ' + ' '.join([f'{t:2s}' for t in sel_store.types[n:n + 20]]))
        remarks.append('REMARK    AMBER ATOM TYPES -->')
        remarks = '\n' + '\n'.join(remarks)
        sel_store.write(store_dcd.replace('.dcd','.pdb'), remarks=remarks)


    def graph_overlap_matrix(self, state_string):
        if state_string == 'complex':
            nc = self.phase1_nc
        elif state_string == 'solvent':
            nc = self.phase2_nc
        else:
            raise Exception('state_string must be either complex or solvent')
        proposed = np.array(tuple([np.array(nc.variables['proposed'][i, :, :]) for i in range(nc.variables['positions'].shape[0])]))
        accepted = np.array(tuple([np.array(nc.variables['accepted'][i, :, :]) for i in range(nc.variables['positions'].shape[0])]))
        matrix = np.sum(accepted, axis=0) / np.sum(proposed, axis=0)
        plt.clf()
        plt.colorbar(plt.matshow(matrix))
        plt.xlabel('Replica i')
        plt.ylabel('Replica j')
        plt.savefig('./%s/%s/%s.png'%(self.out_dir, 'overlap_matrices', state))


    def graph_rmsds(self):


    def get_GSBA_energies(self):


def get_from_yaml(name, yaml_contents):
    truth_list = [True if name in line else False for line in yaml_contents]
    if True in truth_list:
        idj = truth_list.index(True)
        idk = yaml_contents[idj].index(':')
        return yaml_contents[idj][idk + 2:]  # Colon space beginseq
    else:
        raise Exception('%s not found in yaml_contents'%name)