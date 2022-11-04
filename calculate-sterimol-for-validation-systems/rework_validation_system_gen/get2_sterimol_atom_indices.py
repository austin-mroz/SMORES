import argparse
import logging
import pathlib
import smores
import streusel
import morfeus
import dbstep.Dbstep as db
import seaborn as sns
import pandas as pd
import glob
from openbabel import openbabel
import rdkit.Chem.AllChem as rdkit
from rdkit.Chem import rdDepictor, Draw
from rdkit.Chem import rdRGroupDecomposition


def main() -> None:
    core = rdkit.MolFromSmiles('C[*:1]')
    Draw.MolToFile(core, 'core.png')

    subt = rdkit.MolFromSmiles('[*:1]CC')
    Draw.MolToFile(subt, 'subt.png')

    combo = rdkit.CombineMols(core, subt)
    Draw.DrawingOptions.includeAtomNumbers = True
    Draw.MolToFile(combo, 'combo.png')

    edcombo = rdkit.EditableMol(combo)

    for i,atom in enumerate(combo.GetAtoms()):
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    Draw.MolToFile(combo, 'combo2.png')


if __name__ == '__main__':
    main()
