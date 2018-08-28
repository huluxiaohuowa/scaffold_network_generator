from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from copy import deepcopy

from data import graph
m = graph.get_mol_graph('CCC(C1=CC(C=NC2)=C2C(CC3(CC4)CCC4CC3)=C1CC5CCC6(CCCCC6)C5)=O')
print(m.graph_list_to_list)
