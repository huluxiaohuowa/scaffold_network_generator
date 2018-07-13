import data
from rdkit import Chem

m = Chem.MolFromSmiles('CCC(C1=CC(C=NC2)=C2C(CC3(CC4)CCC4CC3)=C1CC5CCC6(CCCCC6)C5)=O')
m.__class__=data.MolGraph
t = m.sng_unique
print(len(t))