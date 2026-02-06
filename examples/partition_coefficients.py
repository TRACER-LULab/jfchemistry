"""Example workflow for calculating partition coefficients."""

# from fireworks import LaunchPad

from jobflow.core.flow import Flow
from jobflow.managers.local import run_locally

# from jobflow.managers.fireworks import flow_to_workflow
from jfchemistry.generation import RDKitGeneration
from jfchemistry.inputs import Smiles
from jfchemistry.workflows.partition_coefficient import PartitionCoefficientWorkflow

SMILES = "C(=O)(C(F)(F)F)O"
pubchem_cid = Smiles(remove_salts=False).make(input=SMILES)

generate_structure = RDKitGeneration(num_conformers=1).make(pubchem_cid.output.structure)

pc_calculation = PartitionCoefficientWorkflow(
    threads=16,
    alpha_phase="chloroform",
    crest_executable="/home/carson/Downloads/crest-gnu-12-ubuntu-latest/crest/crest",
).make(
    generate_structure.output.structure,
)

## ----- FLOW -------
flow = Flow(
    [pubchem_cid, generate_structure, pc_calculation],
    name="Partition Coefficient - TFA - Octanol/Water",
)

response = run_locally(flow)
