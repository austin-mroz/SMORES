import rdkit.Chem.AllChem as rdkit

mol = rdkit.MolFromXYZFile('/home/he/work/SMORES/validation/catalytic_reactions/one_substituent/3_output/allylation/allylation_catalyst_Me/geom.xyz')

rdkit.SanitizeMol(mol)
mol = rdkit.AddHs(mol)

rdkit.EmbedMolecule(mol, rdkit.ETKDGv3())

rdkit.MolToXYZFile(mol, 'etkdg.xyz')
