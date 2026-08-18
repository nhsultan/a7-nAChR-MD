"""
openmm_run_Ca_restraints.py

Description:
This is a modified CHARMM-GUI OpenMM production script. 
It contains custom `CustomBondForce` restraints to maintain the coordination 
of Calcium (CAL) ions to the GLU 44 residues during the simulation (applicable 
to models like 7KOX).

Base Script: CHARMM-GUI (http://www.charmm-gui.org)
"""

from __future__ import print_function
import argparse
import sys
import os
import numpy as np

from omm_readinputs import *
from omm_readparams import *
from omm_vfswitch import *
from omm_barostat import *
from omm_restraints import *
from omm_rewrap import *

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

parser = argparse.ArgumentParser()
parser.add_argument('--platform', nargs=1, help='OpenMM platform (default: CUDA or OpenCL)')
parser.add_argument('-i', dest='inpfile', help='Input parameter file', required=True)
parser.add_argument('-p', dest='topfile', help='Input topology file', required=True)
parser.add_argument('-c', dest='crdfile', help='Input coordinate file', required=True)
parser.add_argument('-t', dest='toppar', help='Input CHARMM-GUI toppar stream file (optional)')
parser.add_argument('-b', dest='sysinfo', help='Input CHARMM-GUI sysinfo stream file (optional)')
parser.add_argument('-ff', dest='fftype', help='Input force field type (default: CHARMM)', default='CHARMM')
parser.add_argument('-icrst', metavar='RSTFILE', dest='icrst', help='Input CHARMM RST file (optional)')
parser.add_argument('-irst', metavar='RSTFILE', dest='irst', help='Input restart file (optional)')
parser.add_argument('-ichk', metavar='CHKFILE', dest='ichk', help='Input checkpoint file (optional)')
parser.add_argument('-opdb', metavar='PDBFILE', dest='opdb', help='Output PDB file (optional)')
parser.add_argument('-orst', metavar='RSTFILE', dest='orst', help='Output restart file (optional)')
parser.add_argument('-ochk', metavar='CHKFILE', dest='ochk', help='Output checkpoint file (optional)')
parser.add_argument('-odcd', metavar='DCDFILE', dest='odcd', help='Output trajectory file (optional)')
parser.add_argument('-rewrap', dest='rewrap', help='Re-wrap the coordinates in a molecular basis (optional)', action='store_true', default=False)
args = parser.parse_args()

# Load parameters
print("Loading parameters...")
inputs = read_inputs(args.inpfile)

top = read_top(args.topfile, args.fftype.upper())
crd = read_crd(args.crdfile, args.fftype.upper())
if args.fftype.upper() == 'CHARMM':
    params = read_params(args.toppar)
    top = read_box(top, args.sysinfo) if args.sysinfo else gen_box(top, crd)

# Build system
nboptions = dict(nonbondedMethod=inputs.coulomb,
                 nonbondedCutoff=inputs.r_off*nanometers,
                 constraints=inputs.cons,
                 ewaldErrorTolerance=inputs.ewald_Tol)
if inputs.vdw == 'Switch': nboptions['switchDistance'] = inputs.r_on*nanometers
if inputs.vdw == 'LJPME':  nboptions['nonbondedMethod'] = LJPME

if   args.fftype.upper() == 'CHARMM': system = top.createSystem(params, **nboptions)
elif args.fftype.upper() == 'AMBER':  system = top.createSystem(**nboptions)

if inputs.vdw == 'Force-switch': system = vfswitch(system, top, inputs)
if inputs.lj_lrc == 'yes':
    for force in system.getForces():
        if isinstance(force, NonbondedForce): force.setUseDispersionCorrection(True)
        if isinstance(force, CustomNonbondedForce) and force.getNumTabulatedFunctions() != 1:
            force.setUseLongRangeCorrection(True)
if inputs.e14scale != 1.0:
    for force in system.getForces():
        if isinstance(force, NonbondedForce): nonbonded = force; break
    for i in range(nonbonded.getNumExceptions()):
        atom1, atom2, chg, sig, eps = nonbonded.getExceptionParameters(i)
        nonbonded.setExceptionParameters(i, atom1, atom2, chg*inputs.e14scale, sig, eps)

if inputs.pcouple == 'yes':      system = barostat(system, inputs)
if inputs.rest == 'yes':         system = restraints(system, crd, inputs)
integrator = LangevinIntegrator(inputs.temp*kelvin, inputs.fric_coeff/picosecond, inputs.dt*picoseconds)

# =====================================================================
# CUSTOM RESTRAINTS: Calcium (CAL) ions to GLU 44
# =====================================================================
print("Applying custom flat-bottom distance restraints to Calcium ions...")

crd_in_nm = np.array(crd.positions/nanometer)

# atom1: Calcium ions (resname CAL)
atom1 = [31060, 31061, 31062, 31063, 31064] 
# atom2: GLU 44 (name CA)
atom2 = [727, 6939, 13151, 19363, 25575] 

restr_k = 5.0 # kcal/mol/A^2, force constant for the distance restraints
ref_dist = 0.55 # Target distance in nm

force_constant = restr_k * kilocalories_per_mole / angstrom**2
# Flat-bottom harmonic potential
restraint_force = CustomBondForce("0.5*step(r-r0)*k*(r-r0)^2")

restraint_force.addGlobalParameter("k", force_constant.value_in_unit(kilojoules_per_mole / nanometer**2))
restraint_force.addPerBondParameter("r0")

for i in range(len(atom1)):
    atom_a = atom1[i]
    atom_b = atom2[i]
    restraint_force.addBond(atom_a, atom_b, [ref_dist])
    print(f"  -> Bond {i+1}: Atom A ({atom_a}) attached to Atom B ({atom_b}) | Target = {ref_dist} nm")

system.addForce(restraint_force)
# =====================================================================

# Set platform
DEFAULT_PLATFORMS = 'CUDA', 'OpenCL'
enabled_platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
if args.platform:
    if not args.platform[0] in DEFAULT_PLATFORMS:
        print("Unable to find OpenMM platform '{}'; exiting".format(args.platform[0]), file=sys.stderr)
        sys.exit(1)

    platform = Platform.getPlatformByName(args.platform[0])
else:
    for platform in DEFAULT_PLATFORMS:
        if platform in enabled_platforms:
            platform = Platform.getPlatformByName(platform)
            break        

prop = dict(CudaPrecision='single') if platform.getName() == 'CUDA' else dict()

# Build simulation context
simulation = Simulation(top.topology, system, integrator, platform, prop)
simulation.context.setPositions(crd.positions)
if args.icrst:
    charmm_rst = read_charmm_rst(args.icrst)
    simulation.context.setPositions(charmm_rst.positions)
    simulation.context.setVelocities(charmm_rst.velocities)
    simulation.context.setPeriodicBoxVectors(charmm_rst.box[0], charmm_rst.box[1], charmm_rst.box[2])
if args.irst:
    with open(args.irst, 'r') as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))
if args.ichk:
    with open(args.ichk, 'rb') as f:
        simulation.context.loadCheckpoint(f.read())

# Re-wrap
if args.rewrap:
    simulation = rewrap(simulation)

# Calculate initial system energy
print("\nInitial system energy:")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Energy minimization
if inputs.mini_nstep > 0:
    print(f"\nEnergy minimization: {inputs.mini_nstep} steps")
    simulation.minimizeEnergy(tolerance=inputs.mini_Tol*kilojoule/mole, maxIterations=inputs.mini_nstep)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Generate initial velocities
if inputs.gen_vel == 'yes':
    print("\nGenerating initial velocities...")
    if inputs.gen_seed:
        simulation.context.setVelocitiesToTemperature(inputs.gen_temp, inputs.gen_seed)
    else:
        simulation.context.setVelocitiesToTemperature(inputs.gen_temp)

# Production
if inputs.nstep > 0:
    print(f"\nMD run: {inputs.nstep} steps")
    if inputs.nstdcd > 0:
        if not args.odcd: args.odcd = 'output.dcd'
        simulation.reporters.append(DCDReporter(args.odcd, inputs.nstdcd))
    simulation.reporters.append(
        StateDataReporter(sys.stdout, inputs.nstout, step=True, time=True, potentialEnergy=True, temperature=True, progress=True,
                          remainingTime=True, speed=True, totalSteps=inputs.nstep, separator='\t')
    )
    # Simulated annealing?
    if inputs.annealing == 'yes':
        interval = inputs.interval
        temp = inputs.temp_init
        for i in range(inputs.nstep):
            integrator.setTemperature(temp*kelvin)
            simulation.step(1)
            temp += interval
    else:
        simulation.step(inputs.nstep)

# Write restart file
if not (args.orst or args.ochk): args.orst = 'output.rst'
if args.orst:
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(args.orst, 'w') as f:
        f.write(XmlSerializer.serialize(state))
if args.ochk:
    with open(args.ochk, 'wb') as f:
        f.write(simulation.context.createCheckpoint())
if args.opdb:
    crd = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(top.topology, crd, open(args.opdb, 'w'))