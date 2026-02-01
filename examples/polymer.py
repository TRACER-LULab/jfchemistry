"""Example of using the PubChemCID node to get a molecule from PubChem."""

from jobflow.core.flow import Flow
from jobflow.managers.local import run_locally
from numpy.random import choice

from jfchemistry.calculators.torchsim.fairchem_calculator import FairChemCalculator
from jfchemistry.inputs import PolymerInput
from jfchemistry.molecular_dynamics.torchsim import (
    TorchSimMolecularDynamicsNPTNoseHoover,
    TorchSimMolecularDynamicsNVTNoseHoover,
)
from jfchemistry.packing.packmol import PackmolPacking
from jfchemistry.polymers import GenerateFinitePolymerChain

chain_length = 25
number_chains = 25
finite_chain_jobs = []
finite_chain_structures = []
polymer_job = PolymerInput().make(head="C[*:1]", monomer="[*:1]C[*:2]", tail="C[*:2]")

for _ in range(number_chains):
    rotation_angles = choice([180, 240, -240], size=chain_length, p=[0.75, 0.125, 0.125]).tolist()
    job = GenerateFinitePolymerChain(dihedral_angles=rotation_angles, monomer_dihedral=0.0).make(
        polymer_job.output.structure
    )
    finite_chain_jobs.append(job)
    finite_chain_structures.append(job.output.structure)

pack_job = PackmolPacking(
    packing_mode="box",
    num_molecules=[1] * number_chains,
    density=0.2,  # g/cm^3
).make(finite_chain_structures)

nvt_job = TorchSimMolecularDynamicsNVTNoseHoover(
    calculator=FairChemCalculator(device="cuda"),
    duration=20,
    timestep=1.0,
    temperature=300.0,
    log_temperature=True,
    log_potential_energy=True,
    log_interval=10.0,
    logfile="traj_nvt_1",
    tau=10.0,
).make(pack_job.output.structure)

npt_job = TorchSimMolecularDynamicsNPTNoseHoover(
    calculator=FairChemCalculator(device="cuda", compute_stress=True),
    duration=20,
    timestep=1.0,
    temperature=300.0,
    log_temperature=True,
    log_potential_energy=True,
    log_interval=10.0,
    log_volume=True,
    logfile="traj_npt_2",
    b_tau=1000.0,
    t_tau=10.0,
).make(pack_job.output.structure)

flow = Flow(
    [
        polymer_job,
        *finite_chain_jobs,
        pack_job,
        nvt_job,
        npt_job,
    ]
)

run_locally(flow)
