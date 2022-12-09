import itertools
import typing
from dataclasses import dataclass

import rdkit.Chem.AllChem as rdkit
import stk


def rdkit_from_smiles(smiles: str) -> rdkit.Mol:
    """
    Get an embeded :mod:`rdkit` molecule from SMILES.

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
    core: rdkit.Mol
    substituent: rdkit.Mol
    product: rdkit.Mol
    dummy_index: int
    attached_index: int
    core_indices: list[int]
    substituent_indices: list[int]


def combine(
    cores: typing.Iterable[rdkit.Mol],
    substituents: typing.Iterable[rdkit.Mol],
    join_atom: str = "Br",
    optimize: bool = True,
) -> typing.Iterator[Combination]:

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
