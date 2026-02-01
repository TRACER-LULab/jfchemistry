"""Polymer Input Ndoes."""

from dataclasses import dataclass

from jobflow.core.job import Response
from pydantic import ConfigDict
from rdkit.Chem import rdmolfiles, rdmolops

from jfchemistry.core.jfchem_job import jfchem_job
from jfchemistry.core.makers.core_maker import CoreMaker
from jfchemistry.core.outputs import Output
from jfchemistry.core.properties import Properties
from jfchemistry.core.structures import Polymer, RDMolMolecule


class PolymerInputOutput(Output):
    """Polymer Input Output."""

    model_config = ConfigDict(arbitrary_types_allowed=True)
    structure: Polymer
    files: None = None
    properties: None = None


@dataclass
class PolymerInput(CoreMaker):
    """Polymer Input."""

    name: str = "Polymer Input"
    _output_model: type[PolymerInputOutput] = PolymerInputOutput
    _properties_model: type[Properties] = Properties

    @jfchem_job()
    def make(
        self, monomer: str, head: str | None = None, tail: str | None = None
    ) -> Response[_output_model]:
        """Make a polymer."""
        monomer_mol = rdmolops.AddHs(rdmolfiles.MolFromSmiles(monomer))
        head_mol = rdmolops.AddHs(rdmolfiles.MolFromSmiles(head)) if head else None
        tail_mol = rdmolops.AddHs(rdmolfiles.MolFromSmiles(tail)) if tail else None

        polymer = Polymer(
            head=RDMolMolecule(head_mol) if head_mol else None,
            monomer=RDMolMolecule(monomer_mol),
            tail=RDMolMolecule(tail_mol) if tail_mol else None,
        )
        return Response(output=PolymerInputOutput(structure=polymer))
