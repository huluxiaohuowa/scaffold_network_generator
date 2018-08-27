from .utils import *
from rdkit import Chem
import networkx as nx
from rdkit.Chem import rdmolops
from collections import Counter
from copy import deepcopy

__all__ = [
    'MolGraph',
    'get_mol_graph'
]


class MolGraph(object):
    def __init__(self, smiles):
        self.mol = Chem.MolFromSmiles(smiles)

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
            sssr_copy2 = deepcopy(list_sssr)
        return sssr_copy2

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

    @property
    def ring_assemblies_list(self):
        list_list_atom_idx_types = []
        list_list_bond_idx_types = []
        list_ring_assemblies = self.ring_assemblies
        for graph in list_ring_assemblies:
            list_atom_idx_types = []
            list_bond_idx_types = []
            for atom_idx in list(graph.nodes):
                list_atom_idx_types.append((
                    atom_idx,
                    self.mol.GetAtomWithIdx(atom_idx).GetSymbol(),
                    self.mol.GetAtomWithIdx(atom_idx).GetFormalCharge(),
                    self.mol.GetAtomWithIdx(atom_idx).GetNumExplicitHs()
                ))
                pass
            for edge in list(graph.edges):
                list_bond_idx_types.append((
                    edge[0],
                    edge[1],
                    get_bond_type(self.mol.GetBondBetweenAtoms(edge[0], edge[1]))
                ))
                pass
            list_list_atom_idx_types.append(list_atom_idx_types)
            list_list_bond_idx_types.append(list_bond_idx_types)
        return list_list_atom_idx_types, list_list_bond_idx_types

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
        for i_graph in self.sng:
            counts = 0
            for i_sng in sng_u:
                if graph_eq(i_graph, i_sng):
                    counts += 1
            if counts == 0:
                sng_u.append(i_graph)
        return sng_u


def get_mol_graph(smiles):
    return MolGraph(smiles)

# class MolGraph2(Chem.rdchem.Mol):
#     def __init__(self):
#         super().__init__()
#
#     @property
#     def sssr(self):
#         list_sssr = []
#         for ring in rdmolops.GetSymmSSSR(self):
#             list_sssr.append(list(ring))
#             sssr_copy2 = deepcopy(list_sssr)
#         if len(list_sssr) > self.GetNumBonds() - self.GetNumAtoms() + 1:
#             sssr_copy, sssr_label = label_gen(list_sssr)
#             counts = 0
#             delta = len(list_sssr) - (self.GetNumBonds() - self.GetNumAtoms() + 1)
#             for i in range(len(sssr_label)):
#                 if sum(sssr_label[i]) == 0:
#                     sssr_copy2.remove(sssr_copy[i])
#                     counts += 1
#                     if counts >= delta:
#                         break
#         return sssr_copy2
#
#     @property
#     def sssr_list(self):
#         list_sssr = []
#         for ring in rdmolops.GetSymmSSSR(self):
#             list_sssr.append(list(ring))
#             sssr_copy2 = deepcopy(list_sssr)
#         return sssr_copy2
#
#     @property
#     def graph(self):
#         atom_types, bonds, bond_types = [], [], []
#         for atom in self.GetAtoms():
#             atom_types.append(get_atom_type(atom))
#         for bond in self.GetBonds():
#             idx_1, idx_2, bond_type = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), get_bond_type(bond)
#             bonds.append([idx_1, idx_2])
#             bond_types.append(bond_type)
#         # build graph
#         graph = nx.Graph()
#         graph.add_nodes_from(range(self.GetNumAtoms()))
#         graph.add_edges_from(bonds)
#         return graph
#
#     @property
#     def cracked_graph(self):
#         graph = deepcopy(self.graph)
#         for i in list(graph.edges):
#             on = False
#             for j in self.sssr_list:
#                 if i[0] in j and i[1] in j:
#                     on = True
#                     break
#             if not on:
#                 graph.remove_edge(i[0], i[1])
#         return graph
#
#     @property
#     def chains(self):
#         sssr = self.sssr_list
#         list_atom_idx = list(range(self.GetNumAtoms()))
#         sssr_tmp = []
#         for i_sssr in sssr:
#             sssr_tmp += i_sssr
#         sssr_list = list(set(sssr_tmp))
#         list_atom_removed_idx = list(set(list_atom_idx).difference(set(sssr_list)))
#         return list_atom_removed_idx
#
#     @property
#     def ring_assemblies(self):
#         graph = deepcopy(self.cracked_graph)
#         graph.remove_nodes_from(self.chains)
#         rings = list(nx.connected_component_subgraphs(graph))
#         return rings
#
#     @property
#     def ring_assemblies_list(self):
#         list_list_atom_idx_types = []
#         list_list_bond_idx_types = []
#         list_ring_assemblies = self.ring_assemblies
#         for graph in list_ring_assemblies:
#             list_atom_idx_types = []
#             list_bond_idx_types = []
#             for atom_idx in list(graph.nodes):
#                 list_atom_idx_types.append((
#                     atom_idx,
#                     self.GetAtomWithIdx(atom_idx).GetSymbol(),
#                     self.GetAtomWithIdx(atom_idx).GetFormalCharge(),
#                     self.GetAtomWithIdx(atom_idx).GetNumExplicitHs()
#                 ))
#                 pass
#             for edge in list(graph.edges):
#                 list_bond_idx_types.append((
#                     edge[0],
#                     edge[1],
#                     get_bond_type(self.GetBondBetweenAtoms(edge[0], edge[1]))
#                 ))
#                 pass
#             list_list_atom_idx_types.append(list_atom_idx_types)
#             list_list_bond_idx_types.append(list_bond_idx_types)
#         return list_list_atom_idx_types, list_list_bond_idx_types
#
#     # remove side chains
#     def get_murko_graph(self, graph=None):
#         if graph is None:
#             graph = self.graph
#         murko = deepcopy(graph)
#         if nx.is_connected(murko):
#             while True:
#                 i = 0
#                 nodes = []
#                 for edge in list(murko.edges):
#                     for item in edge:
#                         nodes.append(item)
#                 counter = Counter(nodes)
#                 for node in counter:
#                     if counter[node] == 1:
#                         i += 1
#                         murko.remove_node(node)
#                 if i == 0:
#                     break
#         else:
#             raise ValueError
#         return murko
#
#     # remove side rings
#     def get_next_level_graph(self, igraph=None, isssr=None):
#         graph_next_level = []
#         sssr_next_level = []
#         if igraph is None:
#             graph = deepcopy(self.graph)
#         else:
#             graph = deepcopy(igraph)
#         if isssr is None:
#             sssr = deepcopy(self.sssr)
#         else:
#             sssr = deepcopy(isssr)
#         d_sssr = d_sssr_single(sssr)
#         if nx.is_connected(graph):
#             murko = self.get_murko_graph(graph)
#             if len(d_sssr) == 0:
#                 return [None], [None]
#             else:
#                 for dd_sssr in d_sssr:
#                     murko_copy = deepcopy(murko)
#                     murko_copy.remove_nodes_from(dd_sssr)
#                     if len(murko_copy) > 0 and nx.is_connected(murko_copy):
#                         murko_copy = self.get_murko_graph(murko_copy)
#                         graph_next_level.append(murko_copy)
#                         sssr_next_level.append(next_sssr(sssr, dd_sssr))
#         else:
#             raise ValueError
#         return graph_next_level, sssr_next_level
#
#     def get_next_level_graph_from_list(self, list_graph, list_sssr):
#         list_graph_next_level = []
#         list_sssr_next_level = []
#         if len(list_graph) > 0:
#             for d_graph, d_sssr in zip(list_graph, list_sssr):
#                 if d_graph is not None:
#                     l1, l2 = self.get_next_level_graph(d_graph, d_sssr)
#                     list_graph_next_level += l1
#                     list_sssr_next_level += l2
#         else:
#             return [None], [None]
#         return list_graph_next_level, list_sssr_next_level
#
#     @property
#     def sng(self):
#         graph_sng = []
#         graph = self.graph
#         sssr = self.sssr
#         if not nx.is_connected(graph):
#             raise ValueError
#         else:
#             graph_sng.append(self.get_murko_graph(graph))
#             list_graph = [graph]
#             list_sssr = [sssr]
#             while True:
#                 list_graph, list_sssr = self.get_next_level_graph_from_list(list_graph, list_sssr)
#                 for i in list_graph:
#                     if i is not None:
#                         graph_sng.append(i)
#                 if not any(list_graph):
#                     break
#             return graph_sng
#
#     @property
#     def sng_unique(self):
#         sng_u = []
#         for i_graph in self.sng:
#             counts = 0
#             for i_sng in sng_u:
#                 if graph_eq(i_graph, i_sng):
#                     counts += 1
#             if counts == 0:
#                 sng_u.append(i_graph)
#         return sng_u


