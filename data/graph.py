from .utils import *
from rdkit import Chem
import networkx as nx
from rdkit.Chem import rdmolops
from collections import Counter
from copy import deepcopy

__all__ = ['MolGraph',
           'get_mol_graph'
           ]


class MolGraph(object):
    def __init__(self, smiles):
        self.mol = Chem.MolFromSmiles(smiles)
        # self.aromatic_chained_nitro = []

    @property
    def sssr(self):
        list_sssr = []
        for ring in rdmolops.GetSymmSSSR(self.mol):
            list_sssr.append(list(ring))
            sssr_copy2 = deepcopy(list_sssr)
        if len(list_sssr) > self.mol.GetNumBonds() - self.mol.GetNumAtoms() + 1:
            sssr_copy, sssr_label = label_gen(list_sssr)
            counts = 0
            delta = len(list_sssr) - (self.mol.GetNumBonds() - self.mol.GetNumAtoms() + 1)
            for i in range(len(sssr_label)):
                if sum(sssr_label[i]) == 0:
                    sssr_copy2.remove(sssr_copy[i])
                    counts += 1
                    if counts >= delta:
                        break
        return sssr_copy2

    @property
    def sssr_list(self):
        list_sssr = []
        for ring in rdmolops.GetSymmSSSR(self.mol):
            list_sssr.append(list(ring))
        return list_sssr

    @property
    def graph(self):
        atom_types, bonds, bond_types = [], [], []
        for atom in self.mol.GetAtoms():
            atom_types.append(get_atom_type(atom))
        for bond in self.mol.GetBonds():
            idx_1, idx_2, bond_type = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), get_bond_type(bond)
            bonds.append([idx_1, idx_2])
            bond_types.append(bond_type)
        # build graph
        graph = nx.Graph()
        graph.add_nodes_from(range(self.mol.GetNumAtoms()))
        graph.add_edges_from(bonds)
        return graph

    @property
    def cracked_graph(self):
        graph = deepcopy(self.graph)
        for i in list(graph.edges):
            on = False
            for j in self.sssr_list:
                if i[0] in j and i[1] in j:
                    on = True
                    break
            if not on:
                graph.remove_edge(i[0], i[1])
        return graph

    @property
    def chains(self):
        sssr = self.sssr_list
        list_atom_idx = list(range(self.mol.GetNumAtoms()))
        sssr_tmp = []
        for i_sssr in sssr:
            sssr_tmp += i_sssr
        sssr_list = list(set(sssr_tmp))
        list_atom_removed_idx = list(set(list_atom_idx).difference(set(sssr_list)))
        return list_atom_removed_idx

    @property
    def ring_assemblies(self):
        graph = deepcopy(self.cracked_graph)
        graph.remove_nodes_from(self.chains)
        rings = list(nx.connected_component_subgraphs(graph))
        return rings

    def graph_list_to_list(self, graph_list=None):
        list_list_atom_idx_types = []
        list_list_bond_idx_types = []
        if graph_list is None:
            graph_list = self.sng_unique
        for graph in graph_list:
            list_atom_idx_types = []
            list_bond_idx_types = []
            i = 0
            for atom_idx in list(graph.nodes):
                if atom_idx not in self.hydro_nitro:
                    list_atom_idx_types.append((i,
                                                atom_idx,
                                                self.mol.GetAtomWithIdx(atom_idx).GetSymbol(),
                                                self.mol.GetAtomWithIdx(atom_idx).GetFormalCharge(),
                                                self.mol.GetAtomWithIdx(atom_idx).GetNumExplicitHs()
                                                ))
                else:
                    list_atom_idx_types.append((i,
                                                atom_idx,
                                                self.mol.GetAtomWithIdx(atom_idx).GetSymbol(),
                                                self.mol.GetAtomWithIdx(atom_idx).GetFormalCharge(),
                                                self.mol.GetAtomWithIdx(atom_idx).GetNumExplicitHs() + 1
                                                ))
                i += 1
            for edge in list(graph.edges):
                list_bond_idx_types.append((edge[0],
                                            edge[1],
                                            get_bond_type(self.mol.GetBondBetweenAtoms(edge[0], edge[1]))
                                            ))
            list_list_atom_idx_types.append(list_atom_idx_types)
            list_list_bond_idx_types.append(list_bond_idx_types)
        return list_list_atom_idx_types, list_list_bond_idx_types

    def ls_mol_from_sng_u(self,
                          ls_ls_atom=None,
                          ls_ls_bond=None,
                          sanitize=True
                          ):
        ls_mol = []
        if not any((ls_ls_atom, ls_ls_bond)):
            ls_ls_atom, ls_ls_bond = self.graph_list_to_list()
        for ls_atom, ls_bond in zip(ls_ls_atom, ls_ls_bond):
            i = 0
            dic = {}
            ls_atom_new = []
            ls_bond_new = []
            for atom in ls_atom:
                dic[atom[1]] = atom[0]
                ls_atom_new.append((atom[2], atom[3], atom[4]))
                i += 1
            for bond in ls_bond:
                ls_bond_new.append((dic[bond[0]], dic[bond[1]], bond[2]))
            mol = get_mol_from_graph(ls_atom_new, ls_bond_new, sanitize)
            ls_mol.append(mol)
        return ls_mol

    # remove side chains
    def get_murko_graph(self, graph=None):
        if graph is None:
            graph = self.graph
        murko = deepcopy(graph)
        if nx.is_connected(murko):
            while True:
                i = 0
                nodes = []
                for edge in list(murko.edges):
                    for item in edge:
                        nodes.append(item)
                counter = Counter(nodes)
                for node in counter:
                    if counter[node] == 1:
                        i += 1
                        murko.remove_node(node)
                if i == 0:
                    break
        else:
            raise ValueError
        return murko

    def get_next_level_ring_assemblies_graph(self, igraph=None, irings=None):
        graph_next_level = []
        rings_next_level = []
        if igraph is None:
            graph = deepcopy(self.graph)
        else:
            graph = deepcopy(igraph)
        if irings is None:
            rings = deepcopy(self.ring_assemblies)
        else:
            rings = deepcopy(irings)
        if nx.is_connected(graph):
            murko = self.get_murko_graph(graph)
            if len(rings) == 0:
                return [None], [None]
            else:
                for i in range(len(rings)):
                    rings_copy = deepcopy(rings)
                    murko_copy = deepcopy(murko)
                    murko_copy.remove_nodes_from(rings_copy[i])
                    if len(murko_copy) > 0 and nx.is_connected(murko_copy):
                        murko_copy = self.get_murko_graph(murko_copy)
                        graph_next_level.append(murko_copy)
                        rings_copy.remove(rings_copy[i])
                        rings_next_level.append(rings_copy)
        else:
            raise ValueError
        return graph_next_level, rings_next_level

    def get_next_level_ra_graph_from_list(self, list_graph, list_list_rings):
        list_graph_next_level = []
        list_rings_next_level = []
        if len(list_graph) > 0:
            for i in range(len(list_graph)):
                if list_graph[i] is not None:
                    l1, l2 = self.get_next_level_ring_assemblies_graph(list_graph[i],
                                                                       list_list_rings[i]
                                                                       )
                    list_graph_next_level += l1
                    list_rings_next_level += l2
            pass
        else:
            return [None], [None]
        return list_graph_next_level, list_rings_next_level

    def get_next_level_graph_from_list(self, list_graph, list_sssr):
        list_graph_next_level = []
        list_sssr_next_level = []
        if len(list_graph) > 0:
            for d_graph, d_sssr in zip(list_graph, list_sssr):
                if d_graph is not None:
                    l1, l2 = self.get_next_level_graph(d_graph, d_sssr)
                    list_graph_next_level += l1
                    list_sssr_next_level += l2
        else:
            return [None], [None]
        return list_graph_next_level, list_sssr_next_level

    # remove side rings
    def get_next_level_graph(self, igraph=None, isssr=None):
        graph_next_level = []
        sssr_next_level = []
        if igraph is None:
            graph = deepcopy(self.graph)
        else:
            graph = deepcopy(igraph)
        if isssr is None:
            sssr = deepcopy(self.sssr)
        else:
            sssr = deepcopy(isssr)
        d_sssr = d_sssr_single(sssr)
        if nx.is_connected(graph):
            murko = self.get_murko_graph(graph)
            if len(d_sssr) == 0:
                return [None], [None]
            else:
                for dd_sssr in d_sssr:
                    murko_copy = deepcopy(murko)
                    murko_copy.remove_nodes_from(dd_sssr)
                    if len(murko_copy) > 0 and nx.is_connected(murko_copy):
                        murko_copy = self.get_murko_graph(murko_copy)
                        graph_next_level.append(murko_copy)
                        sssr_next_level.append(next_sssr(sssr, dd_sssr))
        else:
            raise ValueError
        return graph_next_level, sssr_next_level

    @property
    def sng(self):
        graph_sng = []
        graph = self.graph
        sssr = self.sssr
        if not nx.is_connected(graph):
            raise ValueError
        else:
            graph_sng.append(self.get_murko_graph(graph))
            list_graph = [graph]
            list_sssr = [sssr]
            while True:
                list_graph, list_sssr = self.get_next_level_graph_from_list(list_graph, list_sssr)
                for i in list_graph:
                    if i is not None:
                        graph_sng.append(i)
                if not any(list_graph):
                    break
            return graph_sng

    @property
    def sng_unique(self):
        sng_u = []
        for i_graph in self.ra_sng:
            counts = 0
            for i_sng in sng_u:
                if graph_eq(i_graph, i_sng):
                    counts += 1
            if counts == 0:
                sng_u.append(i_graph)
        return sng_u

    @property
    def ra_sng(self):
        graph_sng = []
        graph = self.graph
        rings = self.ring_assemblies
        if not nx.is_connected(graph):
            raise ValueError
        else:
            graph_sng.append(self.get_murko_graph(graph))
            list_graph = [graph]
            list_rings = [rings]
            while True:
                list_graph, list_rings = self.get_next_level_ra_graph_from_list(list_graph,
                                                                                list_rings)
                for i in list_graph:
                    if i is not None:
                        graph_sng.append(i)
                if not any(list_graph):
                    break
            graph_sng += self.ring_assemblies
            return graph_sng

    @property
    def hydro_nitro(self):
        aromatic_chained_nitro = []
        aromatic_nitros= [a.GetIdx() for a in self.mol.GetAromaticAtoms() if a.GetSymbol() == 'N']
        for nitro_index in aromatic_nitros:
            bond_num = 0
            for atom_index in range(self.mol.GetNumAtoms()):
                if (nitro_index, atom_index) in self.graph.edges:
                    bond_num += 1
            if bond_num >= 3:
                aromatic_chained_nitro.append(nitro_index)
        return aromatic_chained_nitro


def get_mol_graph(smiles):
    return MolGraph(smiles)

