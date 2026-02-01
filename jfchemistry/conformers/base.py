"""Base class for conformer generation.

This module provides the abstract base class for implementing conformer
generation methods in jfchemistry workflows.
"""


class ConformerGeneration:
    """Base class for conformer generation methods.

    This abstract class defines the interface for conformer generation
    implementations. Subclasses should implement the operation() method
    to generate multiple conformations of the input structure.

    Conformer generation explores the potential energy surface to find
    multiple low-energy conformations of a molecule, which is essential
    for accurately predicting molecular properties that depend on
    conformational flexibility.

    Attributes:
        name: Descriptive name for the conformer generation method.

    Examples:
        >>> # Subclass implementation
        >>> class MyConformerGenerator(ConformerGeneration): # doctest: +SKIP
        ...     def operation(self, structure): # doctest: +SKIP
        ...         # Generate conformers using some method
        ...         conformers = generate_conformers(structure) # doctest: +SKIP
        ...         properties = {"num_conformers": len(conformers)} # doctest: +SKIP
        ...         return conformers, properties # doctest: +SKIP
        >>>
        >>> gen = MyConformerGenerator() # doctest: +SKIP
        >>> job = gen.make(molecule) # doctest: +SKIP
        >>> conformers = job.output["structure"] # doctest: +SKIP
    """

    name: str = "Conformer Generation"
