import itertools
import typing
from dataclasses import dataclass

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

    factories = {
        "F": stk.FluoroFactory(),
        "Br": stk.BromoFactory(),
        "I": stk.IodoFactory(),
    }

    for core, substituent in itertools.product(cores, substituents):
        core_atoms = tuple(
            vars(stk)[atom](id) for id, atom in enumerate(core.atoms)
        )
        core_bb = stk.BuildingBlock.init(
            atoms=core_atoms,
            bonds=tuple(
                stk.Bond(
                    atom1=core_atoms[bond.atom1],
                    atom2=core_atoms[bond.atom2],
                    order=bond.order,
                )
                for bond in core.bonds
            ),
            position_matrix=core.positions,
            functional_groups=[factories[dummy_atom]],
        )
        substituent_atoms = tuple(
            vars(stk)[atom](id) for id, atom in enumerate(core.atoms)
        )
        subsituent_bb = stk.BuildingBlock.init(
            atoms=substituent_atoms,
            bonds=tuple(
                stk.Bond(
                    atom1=core_atoms[bond.atom1],
                    atom2=core_atoms[bond.atom2],
                    order=bond.order,
                )
                for bond in core.bonds
            ),
            position_matrix=core.positions,
            functional_groups=[factories[dummy_atom]],
        )
        construct = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(core_bb, subsituent_bb),
                repeating_unit="AB",
                num_repeating_units=1,
            ),
        )
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
            if atom_info.get_building_block() is core_bb
            and atom_info.get_atom().get_id() in added_bond_atoms
        )
        substituent_atom = (
            added_bond.get_atom1()
            if core_atom.get_id() == added_bond.get_atom2().get_id()
            else added_bond.get_atom2()
        )

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

        yield Combination(
            core=core,
            substituent=substituent,
            product=product,
            dummy_index=core_atom.get_id(),
            attached_index=substituent_atom.get_id(),
        )
