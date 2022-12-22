import itertools
import typing
from dataclasses import dataclass

import rdkit.Chem.AllChem as rdkit
import stk


def rdkit_from_smiles(smiles: str) -> rdkit.Mol:
    """
    Get an embeded :mod:`rdkit` molecule from SMILES.

    Examples:

        *Create an rdkit molecule*

        .. doctest:: create-an-rdkit-molecule

            >>> import smores
            >>> molecule = smores.rdkit_from_smiles("Br")
            >>> import rdkit.Chem.AllChem as rdkit
            >>> rdkit.MolToSmiles(molecule)
            '[H]Br'

    Parameters:
        smiles: The SMILES of the molecule.
    Returns:
        The :mod:`rdkit` molecule.

    """

    molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
    _optimize(molecule)
    rdkit.Kekulize(molecule)
    return molecule


@dataclass
class Combination:
    """
    A combination of a core and substituent molecule.
    """

    #: The core molecule.
    core: rdkit.Mol
    #: The substituent molecule.
    substituent: rdkit.Mol
    #: The combination of the :attr:`core` and :attr:`substituent`.
    product: rdkit.Mol
    #: The index of the dummy atom.
    dummy_index: int
    #: The index of the attached atom of the substituent.
    attached_index: int
    #: The indices of atoms in :attr:`product`,
    #: which are from the :attr:`core`.
    core_indices: list[int]
    #: The indices of atoms in :attr:`product`, which are
    #: from the :attr:`substituent`.
    substituent_indices: list[int]


def combine(
    cores: typing.Iterable[rdkit.Mol],
    substituents: typing.Iterable[rdkit.Mol],
    join_atom: typing.Literal["F", "Cl", "Br"] = "Br",
    optimize: bool = True,
) -> typing.Iterator[Combination]:
    """
    Yield a set of molecules by combining `cores` with `substituents`.

    This function combines core and substituent molecules by creating a
    bond where the `join_atom` is specified. For example, if the join atom is
    Br, the operation is::

        CCBr + BrCCN -> CCCCN

    The function produces every possible pairwise combination of
    `cores` and `substituents`.

    Examples:

        *Compare the steric parameters of various substituents*

        .. testcode:: compare-steric-parameters-of-substituents

            import smores
            import rdkit.Chem.AllChem as rdkit

            cores=[
                smores.rdkit_from_smiles("NCBr"),
            ]
            substituents=[
                smores.rdkit_from_smiles("CCBr"),
                smores.rdkit_from_smiles("CCCCBr"),
            ]
            for combo in smores.combine(cores, substituents):
                molecule = smores.Molecule.from_combination(combo)
                params = molecule.get_steric_parameters()
                print(
                    f"Combination of \
{rdkit.MolToSmiles(rdkit.RemoveHs(combo.core))} and "
                    f"{rdkit.MolToSmiles(rdkit.RemoveHs(combo.substituent))} "
                    f"has SMORES parameters of {params}."
                )

        .. testoutput:: compare-steric-parameters-of-substituents

            Combination of NCBr and CCBr has SMORES parameters of \
StericParameters(L=4.768506137873876, B1=1.7742858563716175, \
B5=3.430948246451143).
            Combination of NCBr and CCCCBr has SMORES parameters of \
StericParameters(L=6.713047422661813, B1=1.7856585950046786, \
B5=3.5556262756974535).

    Parameters:

        cores (list[Mol]):
            The core molecules.

        substituents (list[Mol]):
            The substituent molecules.

        join_atom:
            The atom in `cores` and `substituents` which specifies
            the location where they are joined.

        optimize:
            If ``True``, the generated molecules will have an optimized
            structure generated with ETKDG_.

    Yield:
        A combination of a core and substituent molecule.

    .. _ETKDG: https://www.rdkit.org/docs/source/\
rdkit.Chem.rdDistGeom.html#rdkit.Chem.rdDistGeom.ETKDGv3


    """

    for core, substituent in itertools.product(cores, substituents):
        core_bb = _get_building_block(core, join_atom)
        subsituent_bb = _get_building_block(substituent, join_atom)
        construct = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(core_bb, subsituent_bb),
                repeating_unit="AB",
                num_repeating_units=1,
            ),
        )
        steric_atoms = _get_steric_atoms(construct, core_bb)

        product = construct.to_rdkit_mol()
        if optimize:
            _optimize(product)

        core_indices = []
        substituent_indices = []
        for atom_info in construct.get_atom_infos():
            if atom_info.get_building_block() is core_bb:
                core_indices.append(atom_info.get_atom().get_id())
            if atom_info.get_building_block() is subsituent_bb:
                substituent_indices.append(atom_info.get_atom().get_id())

        yield Combination(
            core=core,
            substituent=substituent,
            core_indices=core_indices,
            substituent_indices=substituent_indices,
            product=product,
            dummy_index=steric_atoms.core_atom,
            attached_index=steric_atoms.substituent_atom,
        )


def _get_building_block(
    molecule: rdkit.Mol,
    join_atom: str,
) -> stk.BuildingBlock:

    factories = {
        "F": stk.FluoroFactory(),
        "Br": stk.BromoFactory(),
        "I": stk.IodoFactory(),
    }
    rdkit.SanitizeMol(molecule)
    rdkit.Kekulize(molecule)
    return stk.BuildingBlock.init_from_rdkit_mol(
        molecule=molecule,
        functional_groups=[factories[join_atom]],
    )


@dataclass(frozen=True, slots=True)
class _StericAtoms:
    core_atom: int
    substituent_atom: int


def _get_steric_atoms(
    construct: stk.ConstructedMolecule,
    core: stk.BuildingBlock,
) -> _StericAtoms:
    added_bond = next(
        bond_info.get_bond()
        for bond_info in construct.get_bond_infos()
        if bond_info.get_building_block() is None
    )
    added_bond_atoms = (
        added_bond.get_atom1().get_id(),
        added_bond.get_atom2().get_id(),
    )
    core_atom = next(
        atom_info.get_atom()
        for atom_info in construct.get_atom_infos()
        if atom_info.get_building_block() is core
        and atom_info.get_atom().get_id() in added_bond_atoms
    )
    substituent_atom = (
        added_bond.get_atom1()
        if core_atom.get_id() == added_bond.get_atom2().get_id()
        else added_bond.get_atom2()
    )
    return _StericAtoms(
        core_atom=core_atom.get_id(),
        substituent_atom=substituent_atom.get_id(),
    )


def _optimize(molecule: rdkit.Mol) -> None:
    etkdg = rdkit.ETKDGv3()
    etkdg.randomSeed = 4
    rdkit.EmbedMolecule(molecule, etkdg)
