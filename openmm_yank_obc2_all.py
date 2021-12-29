from __future__ import print_function
import simtk
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import numpy as np
import netCDF4 as netcdf


def run_complex(fname_prmtop, fname_NC, my_platform):

    prmtop = app.AmberPrmtopFile(fname_prmtop)
    nc = netcdf.Dataset(fname_NC)

    niterations, nstates, natoms, ndim = nc.variables['positions'].shape

    # inpcrd = app.AmberInpcrdFile('input.inpcrd')
    # print('crds_A', crds_A)
    # """Remove gbsaModel"""
    system_obc2 = prmtop.createSystem(nonbondedMethod=app.NoCutoff,
                                      constraints=app.HBonds,
                                      implicitSolvent=app.OBC2)

    dt = 1.0
    temp = 300.0

    integrator = mm.LangevinIntegrator(
        temp*unit.kelvin, 1.0/unit.picoseconds,
        dt*unit.femtoseconds)

    integrator.setConstraintTolerance(0.00001)

    if my_platform == 'OpenCL':
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {'OpenCLPrecision': 'mixed'}
        sim_obc2 = app.Simulation(prmtop.topology, system_obc2,
                                  integrator, platform, properties)
    elif my_platform == 'CUDA':
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        sim_obc2 = app.Simulation(prmtop.topology, system_obc2,
                                  integrator, platform, properties)
    else:
        sim_obc2 = app.Simulation(prmtop.topology, system_obc2,
                                  integrator)

    E_list = []
    for iter in range(niterations):
        for istate in [0]:
            idx = list(nc.variables['states'][iter, :]).index(istate)
            pos = nc.variables['positions'][iter, idx].data  # nm
            sim_obc2.context.setPositions(pos)
            state = sim_obc2.context.getState(getEnergy=True)
            E_obc2 = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

            E_list.append(E_obc2)

    return np.array(E_list)




def run_ligand(fname_prmtop, fname_NC, my_platform):

    prmtop = app.AmberPrmtopFile(fname_prmtop)
    nc = netcdf.Dataset(fname_NC)
    numAtoms = prmtop._prmtop.getNumAtoms()

    niterations, nstates, natoms, ndim = nc.variables['positions'].shape

    # inpcrd = app.AmberInpcrdFile('input.inpcrd')
    # print('crds_A', crds_A)
    # """Remove gbsaModel"""
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff,
                                 constraints=app.HBonds,
                                 implicitSolvent=app.OBC2)

    dt = 1.0
    temp = 300.0

    integrator = mm.LangevinIntegrator(
        temp*unit.kelvin, 1.0/unit.picoseconds,
        dt*unit.femtoseconds)

    integrator.setConstraintTolerance(0.00001)

    if my_platform == 'OpenCL':
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {'OpenCLPrecision': 'mixed'}
        simulation = app.Simulation(prmtop.topology, system,
                                    integrator, platform, properties)
    elif my_platform == 'CUDA':
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        simulation = app.Simulation(prmtop.topology, system,
                                    integrator, platform, properties)
    else:
        simulation = app.Simulation(prmtop.topology, system, integrator)

    E_list = []
    for iter in range(niterations):
        for istate in [nstates-1]:
            idx = list(nc.variables['states'][iter, :]).index(istate)
            pos = nc.variables['positions'][iter, idx].data  # nm
            simulation.context.setPositions(pos[-numAtoms:])

            state = simulation.context.getState(getEnergy=True)

            E = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
            E_list.append(E)

    return np.array(E_list)


def run_protein (fname_prmtop, fname_NC, my_platform):

    prmtop = app.AmberPrmtopFile(fname_prmtop)
    nc = netcdf.Dataset(fname_NC)
    numAtoms = prmtop._prmtop.getNumAtoms()

    niterations, nstates, natoms, ndim = nc.variables['positions'].shape

    # inpcrd = app.AmberInpcrdFile('input.inpcrd')
    # print('crds_A', crds_A)
    # """Remove gbsaModel"""
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff,
                                 constraints=app.HBonds,
                                 implicitSolvent=app.OBC2)

    dt = 1.0
    temp = 300.0

    integrator = mm.LangevinIntegrator(
        temp*unit.kelvin, 1.0/unit.picoseconds,
        dt*unit.femtoseconds)

    integrator.setConstraintTolerance(0.00001)

    if my_platform == 'OpenCL':
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {'OpenCLPrecision': 'mixed'}
        simulation = app.Simulation(prmtop.topology, system,
                                    integrator, platform, properties)
    elif my_platform == 'CUDA':
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        simulation = app.Simulation(prmtop.topology, system,
                                    integrator, platform, properties)
    else:
        simulation = app.Simulation(prmtop.topology, system, integrator)

    E_list = []
    for iter in range(niterations):
        for istate in [nstates-1]:
            idx = list(nc.variables['states'][iter, :]).index(istate)
            pos = nc.variables['positions'][iter, idx].data  # nm
            simulation.context.setPositions(pos[:numAtoms])

            state = simulation.context.getState(getEnergy=True)

            E = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
            E_list.append(E)

    return np.array(E_list)



if __name__ == "__main__":
    import getopt
    import json
    import sys

    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, "hi:",
                               ["help=", "input="])

    if len(opts) == 0:
        print('python openmm_yank.py -i <input json file> ')
        sys.exit(2)

    fname_json = 'input.json'
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('python openmm_yank.py -i <input json file> ')
            sys.exit(1)
        if opt in ("-i", "--ifile"):
            fname_json = arg

    with open(fname_json) as f:
        data = json.load(f)
        fname_NC = data['NC']
        fname_dat = 'yank_complex.dat'
        if 'fname_dat' in data:
            fname_dat = data['fname_dat']
        my_platform = data['Platform']

        fname_prmtop = data['complex_prmtop']
        E_com = run_complex(fname_prmtop, fname_NC, my_platform)
        fname_prmtop = data['ligand_prmtop']
        E_lig = run_ligand(fname_prmtop, fname_NC, my_platform)
        fname_prmtop = data['protein_prmtop']
        E_prt = run_protein(fname_prmtop, fname_NC, my_platform)

        fout = open(fname_dat, 'w', 1)
        print("#  E(com)[kcal/mol]   E_lig    E_prt", file=fout)
        for it in range(len(E_com)):
            print(it, E_com[it], E_lig[it], E_prt[it], file=fout)

        fout.close()