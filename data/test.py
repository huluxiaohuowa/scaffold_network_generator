from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from data import MolGraph

from copy import deepcopy
# m = Chem.MolFromSmiles('CCC(C1=CC(C=NC2)=C2C(CC3(CC4)CCC4CC3)=C1CC5CCC6(CCCCC6)C5)=O')
# m.__class__ = MolGraph
from . import graph
graph.get_mol_graph()


