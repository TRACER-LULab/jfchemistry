"""Base class for makers that process single molecules."""

from dataclasses import dataclass

from jobflow.core.job import Response
from rdkit.Chem import rdchem

from jfchemistry.core.input_types import (
    RecursiveMoleculeList,
    RecursiveRDKitMoleculeList,
    RecursiveStructureList,
)
from jfchemistry.core.makers.jfchem_maker import JFChemMaker
from jfchemistry.core.outputs import Output
from jfchemistry.core.properties import Properties
from jfchemistry.core.structures import RDMolMolecule


@dataclass
class RDKitMaker[
    InputType: RecursiveRDKitMoleculeList,
    OutputType: RecursiveRDKitMoleculeList | RecursiveMoleculeList | RecursiveStructureList,
](JFChemMaker[InputType, OutputType]):
    """Base class for operations on molecules without 3D geometry.

    This Maker processes RDMolMolecule objects that do not yet have assigned 3D
    coordinates. It is primarily used for structure generation tasks that convert
    molecular representations (SMILES, SMARTS) into 3D structures.

    The class handles automatic job distribution for lists of molecules and
    molecules with multiple conformers, enabling parallel processing of multiple
    structures.

    Attributes:
        name: Descriptive name for the job/operation being performed.

    """

    name: str = "Single RDMolMolecule Maker"
    _output_model: type[Output] = Output
    _properties_model: type[Properties] = Properties

    def _handle_conformers(self, structure: RDMolMolecule) -> Response[_output_model]:
        """Distribute workflow jobs for each conformer in a molecule.

        Creates separate jobs for each conformer in an RDKit molecule, allowing
        parallel processing of multiple conformations. This is useful when a
        molecule has multiple embedded conformers that need to be processed
        independently (e.g., optimized separately).

        Args:
            maker: A Maker instance that will process each conformer.
            structure: RDMolMolecule containing one or more conformers.

        Returns:
            Response containing:
                - structures: List of processed structures from each job
                - files: List of output files from each job
                - properties: List of computed properties from each job
                - detour: List of jobs to be executed

        """
        jobs: list[Response] = []
        for confId in range(structure.GetNumConformers()):
            s = RDMolMolecule(rdchem.Mol(structure, confId=confId))
            jobs.append(self.make(s))

        return Response(
            output=self._output_model(
                structure=[job.output.structure for job in jobs],
                files=[job.output.files for job in jobs],
                properties=[job.output.properties for job in jobs],
            ),
            detour=jobs,  # type: ignore
        )

    def _run_job(self, input: InputType | list[InputType], **kwargs) -> Response[_output_model]:
        """Run the job for a single molecule or a list of molecules."""
        if not isinstance(input, list) and input.GetNumConformers() > 1:
            return self._handle_conformers(input)
        return super()._run_job(input, **kwargs)
