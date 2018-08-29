from rdkit import Chem
from os import path

__all__ = [
    'MoleculeSpec',
    'get_mol_spec',
]


class MoleculeSpec(object):
    def __init__(self,
                 file_name=path.join(path.dirname(__file__), 'datasets', 'atom_types.txt')
                 ):
        self.atom_types = []
        self.atom_symbols = []
        with open(file_name) as f:
            for line in f:
                atom_type_i = line.strip('\n').split(',')
                self.atom_types.append((atom_type_i[0], int(atom_type_i[1]), int(atom_type_i[2])))
                if atom_type_i[0] not in self.atom_symbols:
                    self.atom_symbols.append(atom_type_i[0])
        for atom_symbol_i in self.atom_symbols:
            if len([atom_type_i for atom_type_i in self.atom_types if atom_type_i[0]==atom_symbol_i]) != 1:
                # add unlabeled atom type
                self.atom_types.append((atom_symbol_i, None, None))
            else:
                continue

        self.bond_orders = [Chem.BondType.AROMATIC,
                            Chem.BondType.SINGLE,
                            Chem.BondType.DOUBLE,
                            Chem.BondType.TRIPLE]
        self.max_iter = 70

    def get_atom_type(self, atom, unlabeled=False):
        atom_symbol = atom.GetSymbol()
        if unlabeled and (atom_symbol, None, None) in self.atom_types:
            return self.atom_types.index((atom_symbol, None, None))
        atom_charge = atom.GetFormalCharge()
        atom_hs = atom.GetNumExplicitHs()
        return self.atom_types.index((atom_symbol, atom_charge, atom_hs))

    def get_bond_type(self, bond):
        return self.bond_orders.index(bond.GetBondType())

    def index_to_atom(self, idx):
        atom_symbol, atom_charge, atom_hs = self.atom_types[idx]
        assert atom_charge is not None and atom_hs is not None
        a = Chem.Atom(atom_symbol)
        a.SetFormalCharge(atom_charge)
        a.SetNumExplicitHs(atom_hs)
        return a

    def index_to_bond(self, mol, begin_id, end_id, idx):
        mol.AddBond(begin_id, end_id, self.bond_orders[idx])

    @property
    def num_atom_types(self):
        return len(self.atom_types)

    @property
    def num_bond_types(self):
        return len(self.bond_orders)


_mol_spec = None


def get_mol_spec():
    global _mol_spec
    if _mol_spec is None:
        _mol_spec = MoleculeSpec()
    return _mol_spec
