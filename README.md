

# Scaffold Network Generator


[![LICENSE](https://img.shields.io/badge/license-NPL%20(The%20996%20Prohibited%20License)-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)
<a href="https://996.icu"><img src="https://img.shields.io/badge/link-996.icu-red.svg" alt="996.icu"></a>

## Contributers
[![](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/images/0)](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/links/0)[![](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/images/1)](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/links/1)[![](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/images/2)](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/links/2)[![](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/images/3)](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/links/3)[![](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/images/4)](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/links/4)[![](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/images/5)](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/links/5)[![](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/images/6)](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/links/6)[![](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/images/7)](https://sourcerer.io/fame/jach4/jach4/scaffold_network_generator/links/7)
>>>>>>> gh/master
## Usage

```bash
git clone git@github.com:jach4/scaffold_network_generator.git sng
cd sng
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

{% assert_imag Scaffold-network-generator/output_3_1.png %}


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

{% assert_imag Scaffold-network-generator/output_6_0.png %}

### Get the scaffold network of The 111,921th molecule


```python
scaffold_list = get_mol_graph(smiles_from_line(111921)).ls_mol_from_sng_u()
```


```python
Draw.MolsToGridImage(scaffold_list)
```

{% assert_imag Scaffold-network-generator/output_9_0.png %}





