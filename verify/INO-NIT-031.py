# INO-NIT-031 検証：縮合リン酸の P のつながり（グラフ）を総当たりで数え，分子式・モル質量・滴定量を確かめる
import itertools, sys
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def graphs(n):
    nodes = range(n); pairs = list(itertools.combinations(nodes, 2)); reps = []
    for mask in range(1 << len(pairs)):
        es = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
        G = nx.Graph(); G.add_nodes_from(nodes); G.add_edges_from(es)
        if max(dict(G.degree()).values()) > 3 or not nx.is_connected(G): continue
        if any(nx.is_isomorphic(G, H) for H in reps if H.number_of_edges() == G.number_of_edges()): continue
        reps.append(G)
    return reps

def smiles_of(G):
    # P を原子，P–O–P を O で表した分子を RDKit で組み立てる
    m = Chem.RWMol(); idx = {}
    for v in G.nodes:
        p = m.AddAtom(Chem.Atom(15)); idx[v] = p
        o = m.AddAtom(Chem.Atom(8)); m.AddBond(p, o, Chem.BondType.DOUBLE)
    for u, v in G.edges:
        o = m.AddAtom(Chem.Atom(8)); m.AddBond(idx[u], o, Chem.BondType.SINGLE); m.AddBond(idx[v], o, Chem.BondType.SINGLE)
    for v in G.nodes:
        for _ in range(3 - G.degree(v)):
            o = m.AddAtom(Chem.Atom(8)); m.AddBond(idx[v], o, Chem.BondType.SINGLE)
    mol = m.GetMol(); Chem.SanitizeMol(mol); return mol

reps = graphs(6)
print('n=6 の連結なつながり方（次数3以下）:', len(reps))
def kind(G):
    deg = dict(G.degree()); E = sum(1 for d in deg.values() if d == 1); M = sum(1 for d in deg.values() if d == 2); B = sum(1 for d in deg.values() if d == 3)
    r = G.number_of_edges() - G.number_of_nodes() + 1
    return E, M, B, r
X = [G for G in reps if kind(G) == (3, 2, 1, 0)]
Y = [G for G in reps if kind(G) == (1, 4, 1, 1)]
print('X 候補', len(X), [sorted(G.edges) for G in X])
print('Y 候補', len(Y), [sorted(G.edges) for G in Y])
for G in X + Y:
    mol = smiles_of(G); E, M, B, r = kind(G)
    # 原子量 H1.0 O16 P31 で
    f = CalcMolFormula(mol); nH = sum(a.GetTotalNumHs() for a in mol.GetAtoms())
    Mw = nH*1.0 + 6*31 + sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'O')*16
    print(f, Mw, '強酸H', E + M, '弱酸H', E, Chem.MolToSmiles(mol))
# 一般式 M = 80n + 18 - 18r
for n in range(2, 12):
    for r in range(0, 4):
        M = 80*n + 18 - 18*r
        if M in (498, 480): print('M', M, 'n', n, 'r', r)
# 滴定量
for name, M, strong, weak in (('X', 498, 5, 3), ('Y', 480, 5, 1)):
    mmol = (0.249 if name == 'X' else 0.240)/M*1000
    v1 = mmol*strong/0.100; v2 = v1 + mmol*weak/0.100
    print(name, 'mmol', mmol, 'V1', v1, 'V2', v2)
# 末端を 1 個ずつ外す酵素：Y の各候補で H3PO4 がいくつ出て止まるか
for G in Y:
    H = G.copy(); k = 0
    while True:
        leaves = [v for v in H if H.degree(v) == 1]
        if not leaves: break
        H.remove_node(leaves[0]); k += 1
    print('Y 候補 環の大きさ', H.number_of_nodes(), '放出 H3PO4', k)
for e in ['H6P6O18 + 2H2O -> H4P4O12 + 2H3PO4', '6H3PO4 -> H8P6O19 + 5H2O', '6H3PO4 -> H6P6O18 + 6H2O']:
    print(e, check_equation(e) or 'OK')
