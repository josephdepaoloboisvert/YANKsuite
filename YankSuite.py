import yaml
import os
import matplotlib.pyplot as plt
import MDAnalysis as mda
import yank.analyze
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
import netCDF4 as netcdf

from yank.experiment import *
from yank.analyze import *
from openmm_yank_obc2_all import *
from slice_yank_nc import *

class YankMultiAnalyzer():
    """A Class for comparing multiple YankAnalyzers"""
    def __init__(self, in_dictionary):
        """Format of in_dictionary is keys = string name; value = list of YankAnalyzer objects that are repeats of each other
        EX: Protein_in_dict = {'LIG1': [YankAnalyzerObjectFromRepeat1, YankAnalyzerObjectFromRepeat2, YankAnalyzerObjectFromRepeat3], etc... }
        Every value-list should have at least two entries and the overall dictionary should have at least one entry."""
        self.data_dict = in_dictionary


    def convergence_of_rpts(self, data_dict_key, discard = 0.2, res = 50):
        """Analyze the convergence of repeats (see init for how repeats are defined)"""
        #copy_yank_nc is used to create many ncs to pass into the YankBFEAnalysis method
        sims = [sim for sim in self.data_dict[data_dict_key]]
        analyzers = [sim.YankBFEAnalysis(auto=False) for sim in sims]
        #Validate that sims have the same parameters (such as length)
        truth_table = [analyzers[i].get_general_simulation_data() == analyzers[i+1].get_general_simulation_data()
                       for i in range(len(analyzers) - 1)]
        if False in truth_table:
            raise Exception('Yank simulations were found to have different parameters')

        #Will create netcdf files to analyze every res iterations beginning discard% along the simulation
        # Ex. res=50 discard=0.2 gives BFE every 50 iterations for the last 80% of a simulation
        os.mkdir('temp_ncs')
        max_iterations = analyzers[0].get_general_simulation_data()['complex']['iterations'] # I think
        first_iter = int(max_iterations * discard)
        nc_iters = np.arange(first_iter, max_iterations, res)

        j = 0
        for sim in sims:
            sim_num = j
            j += 1
            os.mkdir('temp_ncs/sim%s'%sim_num)
            for iter_num in nc_iters:
                #_generate_nc_slice(self, phase, up2iter, out_nc_name):
                for i in [1, 2]:
                    sim._generate_nc_slice(i, iter_num, 'temp_ncs/sim%s/%s_%s.nc'%(sim_num, i, iter_num))


class YankAnalyzer():
    def __init__(self, full_path_to_yaml, do_trajs=False):
        """
        Python class for managing the information exracted from YANK simulations

        PARAMETERS:
        full_path_to_yaml: full path to the yaml file which ran the simulation type(str)
        do_trajs: Bool for whether to extract the physically meaningful trajs on init

        RETURNS:
            does not return
        """
        #Get to directory that contains simulation
        idx = [index for index, element in enumerate(full_path_to_yaml) if element == '/']
        self.yank_dir = full_path_to_yaml[:idx[-1]+1]  # Trailing / included
        print('Yank_dir is %s'%self.yank_dir)
        f = open(full_path_to_yaml, 'r')
        yaml_script_contents = [line[:-1] if line.endswith('\n') else line for line in f.readlines()]
        f.close()

        # Obtain Contents of generate yaml file ($OUTPUT/experiments/experiments.yaml)
        try:
            self.out_dir = self.yank_dir + get_from_yaml('output_dir:', yaml_script_contents) # No Trailing /
        except:
            self.out_dir = self.yank_dir + 'output'
        print('out_dir is %s'%self.out_dir)
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
        for dir in ['trajs', 'overlap_matrices', 'rmsds', 'CoM_distances']:
            if not os.path.exists('%s/%s'%(self.out_dir, dir)):
                os.mkdir('%s/%s'%(self.out_dir, dir))

        self.yaml = YankLoader(self.out_dir + '/experiments/experiments.yaml')  # Probably
        self.phase1_nc = netcdf.Dataset(self.out_dir + '/experiments/complex.nc')  # Probably
        self.phase2_nc = netcdf.Dataset(self.out_dir + '/experiments/solvent.nc')  # Probably

        #MDA Universes reference (input files), complex(complex.nc -> state0), alchemical(complex.nc -> state-1), solvent(solvent.nc -> state0), vacuum(solvent.nc -> state-1)
        self.ref = mda.Universe(self.s_phase1[0], self.s_phase2[1])

        #End Initialization with performing all actions that need the NC files in memory
        # extract trajectories
        if do_trajs:
            self.extract_meaningful_trajs()

        #Make Replica Overlap Matrices
        self.graph_overlap_matrix('complex')
        self.graph_overlap_matrix('solvent')


    def _generate_nc_slice(self, phase, up2iter, out_nc_name):
        """
        :param phase: int 1 or 2 for phase 1 or 2 ('complex' or 'solvent')
        :param up2iter: will copy the src nc file up to and including this iteration
        :param out_nc_name: path to write out_nc to
        :return: does not return
        """
        if phase == 1:
            src_nc = self.phase1_nc
        elif phase == 2:
            src_nc = self.phase2_nc
        else:
            raise Exception('phase should be 1 or 2')

        slice_yank_nc(src_nc, out_nc_name, up2iter)

    def extract_traj(self, state_string, state_index, sel_str_store, sel_str_align):
        """
        extract a thermodynamic state of the yank simulation,, storing the atoms in mdanalysis selection sel_str_store and aligning the trajectory on sel_str_align
        :param state_string: must be either the string 'complex' or 'solvent'
        :param state_index: replica/state index to extract (typically 0 or -1)
        :param sel_str_store: MDAnalysis selection string to store trajectory of
        :param memory: if memory, the extracted trajectory will be returned as an MDAnalysis Universe Object

        "Meaningful" trajs
        extract_trajs('complex', 0, 'all') #protein and ligand together
        extract_trajs('complex', -1, 'protein') # protein independent of ligand's presence
        extract_trajs('complex', -1, self.dsl_lig) # ligand independent of protein's presence (eff. vacuum)
        extract_trajs('solvent', 0, 'all') # ligand in solvent
        extract_trajs('solvent', -1, self.dsl_lig) # ligand indepedent of solvent (eff. vacuum)
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

    def extract_meaningful_trajs(self):
        # trajs are at '%s/trajs/%s_%s_%s.dcd' % (self.out_dir, state_string, i_state, sel_str_store)
        for keys in [['complex', 0, 'all'], ['complex', -1, 'protein'], ['complex', -1, self.dsl_lig],
                     ['solvent', 0, 'all'], ['solvent', -1, self.dsl_lig]]:

            if keys[0] == 'complex' and keys[2] != self.dsl_lig:
                sel_str_align = 'backbone'
            else:
                sel_str_align = self.dsl_lig

            if not os.path.isfile('%s/trajs/%s_%s_%s.dcd'%(self.out_dir, keys[0], keys[1], keys[2])):
                self.extract_traj(keys[0], keys[1], keys[2], sel_str_align)
            else:
                print('Traj already found for %s, state=%s, state=%s, store=%s'%(self.out_dir, keys[0], keys[1], keys[2]))

    def load_traj_as_MDAU(self, state_string, state_index, sel_str_store):
        return mda.Universe('%s/trajs/%s_%s_%s.dcd'%(self.out_dir, state_string, state_index, sel_str_store))


    def graph_overlap_matrix(self, state_string):
        """
        Plot the replica exchange matrix for either the phase1 or phase2 nc
        :param state_string: must be either the string 'phase1' or 'phase2'
        :return:
        """
        if state_string == 'phase1':
            nc = self.phase1_nc
        elif state_string == 'phase2':
            nc = self.phase2_nc
        else:
            raise Exception("state_string must be either the string 'phase1' or 'phase2'")

        proposed = np.array(tuple([np.array(nc.variables['proposed'][i, :, :]) for i in range(nc.variables['positions'].shape[0])]))
        accepted = np.array(tuple([np.array(nc.variables['accepted'][i, :, :]) for i in range(nc.variables['positions'].shape[0])]))
        matrix = np.sum(accepted, axis=0) / np.sum(proposed, axis=0)

        plt.clf()
        plt.colorbar(plt.matshow(matrix))
        plt.xlabel('Replica i')
        plt.ylabel('Replica j')
        plt.savefig('%s/overlap_matrices/%s.png'%(self.out_dir, state_string))


    def graph_rmsd(self, mda_universe, rmsd_select, save_name=None, align_string=None, ref_uni=None):
        """
        Saves a graph of the RMSD of rmsd_select over the simulation provided as an mda_universe
        :param mda_universe: simulation to analyze (see self.load_traj())
        :param rmsd_select: mdanalysis string of the group to obtain rmsd of
        :param save_name: save name for the rmsd graph as string ending in .png
        :param align_string: Perform alignment if desired (self.extract_traj() already aligns)
        :param ref_uni: Allows user to pass in a reference MDAnalysis universe (default is self.ref)
        :return: Does not return; saves png to self.out_dir/rmsds/SAVE_NAME
        """
        sim = mda_universe
        # This if block to allow this function to be more applicable outside YANK
        if ref_uni is not None:
            ref = ref_uni
        else:
            ref = self.ref

        # Align each trajectory to the reference
        if align_string is not None:
            aligner = align.AlignTraj(sim, ref, select=align_string, in_memory=True).run()

        ref_select = ref.select_atoms(rmsd_select)
        sim_select = sim.select_atoms(rmsd_select)
        data = np.array([[ts.time, rms.rmsd(sim_select.positions, ref_select.positions)] for ts in sim.trajectory])

        plt.clf()
        plt.title('RMSD of %s'%rmsd_select)
        plt.xlabel('Time (ns)')
        plt.ylabel('RMSD (Angstroms)')
        plt.scatter(data[:, 0]/1000, data[:, 1])
        if save_name is None:
            plt.savefig('%s/rmsds/%s.png'%(self.out_dir, save_name))
        else:
            plt.savefig(save_name)


    def graph_distance(self, mda_universe, selection_1, selection_2, save_name=None, align_string=None, ref_uni=None):
        """
        graphs the distance between the center of mass of two groups of atoms over the trajectory provided
        :param mda_universe: MDAnalysis Universe object with trajectory loaded
        :param selection_1: MDAnalysis string selection for first atom group
        :param selection_2: MDAnalysis string selection for second atom group
        :param save_name: if not None, will be saved as the argument instead of default
        :param align_string: Perform alignment if desired (self.extract_traj() already aligns)
        :return: does not return, saves png of graph to self.out_dir/CoM_distances/SAVE_NAME
        """
        sim = mda_universe
        # This if block to allow this function to be more applicable outside YANK
        if ref_uni is not None:
            ref = ref_uni
        else:
            ref = self.ref

        # Align each trajectory to the reference
        if align_string is not None:
            aligner = align.AlignTraj(sim, ref, select=align_string, in_memory=True).run()
        #Evaluate the distance between selection 1 and selection 2 over sim
        atoms_1, atoms_2 = sim.select_atoms(selection_1), sim.select_atoms(selection_2)
        data = np.array([[ts.time, distance_between_points(atoms_1.center_of_mass(), atoms_2.center_of_mass())]
                         for ts in sim.trajectory])

        plt.clf()
        plt.title('Distance between %s and %s'%(selection_1, selection_2))
        plt.xlabel('Time (ps?)')
        plt.ylabel('Distance (Angstroms)')
        plt.scatter(data[:, 0], data[:, 1])
        if save_name is None:
            save_name = 'CoM_%s_%s'%(selection_1, selection_2)
        plt.savefig('%s/CoM_distances/%s.png'%(self.out_dir, save_name))

    def YankBFEAnalysis(self, auto=True):
        """
        Perform the YANK binding free energy analysis
        :return:
        """
        analysis_dir = '%s/experiments/analysis.yaml'%self.out_dir
        analyzer = yank.analyze.ExperimentAnalyzer(analysis_dir)
        if auto:
            return analyzer.auto_analyze()
        else:
            return analyzer


    def get_GSBA_energies(self, prmtops, fname_dat, my_platform='CPU'):
        print('Not Fully Implemented GBSA ENERGY')
        complex_prmtop, ligand_prmtop, protein_prmtop = prmtops[0], prmtops[1], prmtops[2]
        fname_NC = self.out_dir + '/experiments/complex.nc'
        E_com = run_complex(complex_prmtop, fname_NC, my_platform)
        E_lig = run_ligand(ligand_prmtop, fname_NC, my_platform)
        E_prt = run_protein(protein_prmtop, fname_NC, my_platform)
        fout = open(fname_dat, 'w', 1)
        print("#  E(com)[kcal/mol]   E_lig    E_prt", file=fout)
        for it in range(len(E_com)):
            print(it, E_com[it], E_lig[it], E_prt[it], file=fout)
        fout.close()


#Utility Functions
def get_from_yaml(name, yaml_contents):
    truth_list = [True if name in line else False for line in yaml_contents]
    if True in truth_list:
        idj = truth_list.index(True)
        idk = yaml_contents[idj].index(':')
        return yaml_contents[idj][idk + 2:]  # Colon space beginseq
    else:
        raise Exception('%s not found in yaml_contents'%name)

def distance_between_points(x_y_z1, x_y_z2):
    """
    Gets the distance (vector difference magnitude) between two points
    :param x_y_z1: three floats as list or tuple
    :param x_y_z2: three floats as list or tuple
    :return: distance between the two points as float
    """
    x_y_z1, x_y_z2 = np.array(x_y_z1), np.array(x_y_z2)
    return np.linalg.norm(x_y_z1 - x_y_z2)