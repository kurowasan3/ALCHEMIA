# 検証スクリプト（東工大 g12）  python3 -I この.py
import sys, json, re, itertools
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
ITEM = '/home/claude/chemx/tokodai/out/items/POL-PRT-068/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
# 小問数と答えの数
assert it['problem'].count('\\item') == len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])) == 5
# アミノ酸（原子数 C,H,N,O）と等電点（教科書値）
AA = {'Gly':(2,5,1,2,6.0),'Ala':(3,7,1,2,6.0),'Ser':(3,7,1,3,5.7),'Lys':(6,14,2,2,9.7),'Asp':(4,7,1,4,2.8),'Glu':(5,9,1,4,3.2)}
pH = 6.0
cath = [a for a,v in AA.items() if v[4] > pH + 0.5]
anod = [a for a,v in AA.items() if v[4] < pH - 0.5]
neut = [a for a in AA if a not in cath + anod]
print('陰極', cath, '陽極', anod, 'ほぼ移動しない', neut)
assert cath == ['Lys'] and set(neut) == {'Gly','Ala','Ser'}
def tri(names):
    c=h=n=o=0
    for x in names:
        v=AA[x]; c+=v[0]; h+=v[1]; n+=v[2]; o+=v[3]
    return (c, h-4, n, o-2)
def mw(f): return 12*f[0]+1*f[1]+14*f[2]+16*f[3]
cands=[]
for x in neut:
    for combo in (['Lys','Lys',x], ['Lys',x,x]):
        f=tri(combo); cands.append((combo,f,mw(f)))
        print(combo, 'C%dH%dN%dO%d'%f, mw(f), 'C:H=1:2' if f[1]==2*f[0] else '')
ok=[c for c in cands if c[1][1]==2*c[1][0] and c[2]>300]
assert len(ok)==1 and sorted(ok[0][0])==['Lys','Ser','Ser'] and ok[0][1]==(12,24,4,6)
# 質量比 22:9 の確認
f=ok[0][1]; r=(44*f[0])/(18*f[1]/2); assert abs(r-22/9)<1e-9
# 解説中の候補の分子式と分子量が計算と一致するか
for combo,f,m in cands:
    s='C%dH%dN%dO%d'%f
    assert s in it['solution'] and str(m) in it['solution'], (s,m)
# RDKit で Lys-Ser-Ser とメチルエステル
P = Chem.MolFromSmiles('NCCCCC(N)C(=O)NC(CO)C(=O)NC(CO)C(=O)O')
Q = Chem.MolFromSmiles('NCCCCC(N)C(=O)NC(CO)C(=O)NC(CO)C(=O)OC')
print(rdMolDescriptors.CalcMolFormula(P), round(Descriptors.MolWt(P),1), rdMolDescriptors.CalcMolFormula(Q))
assert rdMolDescriptors.CalcMolFormula(P)=='C12H24N4O6' and rdMolDescriptors.CalcMolFormula(Q)=='C13H26N4O6'
# 配列の数
seqs=set(itertools.permutations(['Lys','Ser','Ser'])); print(len(seqs), seqs); assert len(seqs)==3
print('ALL OK')
