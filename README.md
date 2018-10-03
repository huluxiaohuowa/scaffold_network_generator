
# Scaffold Network Generator
## Usage
```bash
python scaffolds_output.py --np=50 --file_input=data/datasets/input.smi --scaffolds_output=data/datasets/scaffolds.smi --file_output=data/datasets/scaffolds.bin 
```
## Read the entire scaffold message from a file


```python
from data import *
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
dic = DicSmScaffoldLs()
with open('data/datasets/scaffolds.bin','rb') as f:
    dic.ParseFromString(f.read())
    
print(len(dic.smiles_scaffold))
```

    402344


### Get the 1990th scaffold


```python
print(scaffold_smiles_idx(1990))
Chem.MolFromSmiles(scaffold_smiles_idx(1990))
```

    c1n[nH]nc1-c1nnc(-c2cnns2)o1





![png](images/output_3_1.png)




```python
dic.smiles_scaffold[scaffold_smiles_idx(1990)]
```




    dic_mol_atoms {
      idx_mol: 111921
      ls_atom {
        idx_atom: 1
        idx_atom: 2
        idx_atom: 3
        idx_atom: 10
        idx_atom: 11
        idx_atom: 12
        idx_atom: 13
        idx_atom: 14
        idx_atom: 15
        idx_atom: 16
        idx_atom: 17
        idx_atom: 18
        idx_atom: 19
        idx_atom: 20
        idx_atom: 22
      }
      ls_nh {
        idx_atom: 3
      }
      ls_np {
      }
    }



### The 111,921th molecule containing the 6th scaffold


```python
Chem.MolFromSmiles(smiles_from_line(111921))
```




![png](images/output_6_0.png)



### Get the scaffold network of The 111,921th molecule


```python
scaffold_list = get_mol_graph(smiles_from_line(111921)).ls_mol_from_sng_u()
```


```python
Draw.MolsToGridImage(scaffold_list)
```




![png](images/output_9_0.png)


