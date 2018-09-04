# Scaffold Network Generator
## Usage
```bash
python output.py --np=50 --file_input=data/datasets/input.txt --scaffolds_output=data/datasets/scaffolds.smi --file_output=data/datasets/scaffolds.bin 
```
## Read the entire scaffolds message from a file
```python
from data import *
dic = DicSmScaffoldLs()
with open('data/datasets/scaffolds.bin','rb') as f:
    dic.ParseFromString(f.read())

print(len(dic.smiles_scaffold))

print(scaffold_smiles_idx(6))

print(len(dic.smiles_scaffold[scaffold_smiles_idx(6)].dic_mol_atoms))

print(dic.smiles_scaffold[scaffold_smiles_idx(6)].dic_mol_atoms[0].idx_mol)

print(dic.smiles_scaffold[scaffold_smiles_idx(6)].dic_mol_atoms[0].ls_atom)

```