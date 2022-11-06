import itertools
import typing
from dataclasses import dataclass

import rdkit.Chem.AllChem as rdkit
import stk

from smores._internal.bond import Bond
from smores._internal.esp_molecule import EspMolecule
from smores._internal.molecule import Molecule


@dataclass
class Combination:
    core: Molecule | EspMolecule
    substituent: Molecule | EspMolecule
    product: Molecule
    dummy_index: int
    attached_index: int


def combine(
    cores: typing.Iterable[Molecule | EspMolecule],
    substituents: typing.Iterable[Molecule | EspMolecule],
    dummy_atom: str = "Br",
    optimize: bool = True,
) -> typing.Iterator[Combination]:

    for core, substituent in itertools.product(cores, substituents):
        core_bb = _get_building_block(core, dummy_atom)
        subsituent_bb = _get_building_block(substituent, dummy_atom)
        construct = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(core_bb, subsituent_bb),
                repeating_unit="AB",
                num_repeating_units=1,
            ),
        )
        steric_atoms = _get_steric_atoms(construct, core_bb)

        product = Molecule(
            atoms=(str(atom) for atom in construct.get_atoms()),
            bonds=(
                Bond(
                    atom1=bond.get_atom1().get_id(),
                    atom2=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                )
                for bond in construct.get_bonds()
            ),
            positions=construct.get_position_matrix(),
        )

        if optimize:
            product = _optimize(product)

        yield Combination(
            core=core,
            substituent=substituent,
            product=product,
            dummy_index=steric_atoms.core_atom,
            attached_index=steric_atoms.substituent_atom,
        )


def _get_building_block(
    molecule: Molecule | EspMolecule,
    dummy_atom: str,
) -> stk.BuildingBlock:

    factories = {
        "F": stk.FluoroFactory(),
        "Br": stk.BromoFactory(),
        "I": stk.IodoFactory(),
    }
    atoms = tuple(
        vars(stk)[atom](id) for id, atom in enumerate(molecule.atoms)
    )
    return stk.BuildingBlock.init(
        atoms=atoms,
        bonds=tuple(
            stk.Bond(
                atom1=atoms[bond.atom1],
                atom2=atoms[bond.atom2],
                order=bond.order,
            )
            for bond in molecule.bonds
        ),
        position_matrix=molecule.positions,
        functional_groups=[factories[dummy_atom]],
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


def _optimize(molecule: Molecule) -> Molecule:

    rdkit_molecule = molecule.to_rdkit()
    etkdg = rdkit.ETKDGv3()
    etkdg.randomSeed = 4
    conformer_id = rdkit.EmbedMolecule(rdkit_molecule, etkdg)
    return Molecule.from_rdkit(rdkit_molecule, conformer_id)
