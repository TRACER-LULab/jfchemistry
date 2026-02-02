"""Maker for single structure or molecule operations."""

from dataclasses import dataclass
from typing import Type

from jfchemistry.core.input_types import RecursiveMoleculeList, RecursiveStructureList
from jfchemistry.core.makers.jfchem_maker import JFChemMaker
from jfchemistry.core.outputs import Output


@dataclass
class PymatGenMaker[
    InputType: RecursiveStructureList | RecursiveMoleculeList,
    OutputType: RecursiveStructureList | RecursiveMoleculeList,
](JFChemMaker[InputType, OutputType]):
    """Base class for makers that process single structures or molecules."""

    _output_model: Type[Output] = Output
