"""Input types for JFChemistry."""

from pymatgen.core.structure import Molecule, Structure

from jfchemistry.core.structures import RDMolMolecule

type RecursiveStructureList = Structure | list["RecursiveStructureList"]
type RecursiveMoleculeList = Molecule | list[Molecule] | list["RecursiveMoleculeList"]
type RecursiveRDKitMoleculeList = RDMolMolecule | list["RecursiveRDKitMoleculeList"]
