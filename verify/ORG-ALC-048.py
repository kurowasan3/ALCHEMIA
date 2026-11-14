# ORG-ALC-048 バイヤー・ビリガー酸化（置換基の多い側に O が入る）のモデルによる構造決定
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import product
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
def deg(atom):  # 炭素の級数（結合している炭素の数）
    return sum(1 for n in atom.GetNeighbors() if n.GetSymbol()=='C')
def BV(smi):
    """ケトン C(=O) の両側のうち，級数の大きい側の C と C=O の間に O を入れる"""
    m = Chem.MolFromSmiles(smi)
    patt = Chem.MolFromSmarts('[#6:1][CX3:2](=O)[#6:3]')
    res=set()
    for a,c,_o,b in m.GetSubstructMatches(patt):
        A, B = m.GetAtomWithIdx(a), m.GetAtomWithIdx(b)
        mig = a if deg(A) > deg(B) else b if deg(B) > deg(A) else None
        if mig is None: continue
        rw = Chem.RWMol(m); o = rw.AddAtom(Chem.Atom(8))
        rw.RemoveBond(c, mig); rw.AddBond(c, o, Chem.BondType.SINGLE); rw.AddBond(o, mig, Chem.BondType.SINGLE)
        p = rw.GetMol(); Chem.SanitizeMol(p); res.add(Chem.MolToSmiles(p))
    return res
def hydrolyze(ester):
    rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]')
    out=set()
    for ps in rxn.RunReactants((Chem.MolFromSmiles(ester),)):
        ss=[]
        for p in ps:
            Chem.SanitizeMol(p); ss.append(Chem.MolToSmiles(p))
        out.add('.'.join(sorted(ss)))
    return out
for k in ['CC(=O)C(C)C', 'O=C1CCCCC1', 'O=C1CCCCC1C', 'CC(=O)C(C)C(C)C']:
    e = BV(k); print(k, CalcMolFormula(Chem.MolFromSmiles(k)), '->', e, [hydrolyze(x) for x in e])
# 2-メチルシクロヘキサノンの逆の位置選択性（CH2 側）だった場合
print('other regio hydrolysis:', hydrolyze('O=C1OCCCCC1C'))
# 加水分解生成物のヨードホルム反応（CH3CH(OH)- をもつか）
iod = Chem.MolFromSmarts('[CH3][CH1]([OH])[#6]')
for s in ['CC(O)CCCCC(=O)O', 'CC(CCCCO)C(=O)O']:
    print(s, 'iodoform +' if Chem.MolFromSmiles(s).HasSubstructMatch(iod) else 'iodoform -')
# C5H12O のアルコールの網羅と不斉
alcs = {'1-ペンタノール':'CCCCCO','2-ペンタノール':'CCCC(C)O','3-ペンタノール':'CCC(O)CC','2-メチル-1-ブタノール':'CCC(C)CO',
        '3-メチル-1-ブタノール':'CC(C)CCO','2-メチル-2-ブタノール':'CCC(C)(C)O','3-メチル-2-ブタノール':'CC(C)C(C)O','2,2-ジメチル-1-プロパノール':'CC(C)(C)CO'}
for k,s in alcs.items():
    m = Chem.MolFromSmiles(s); cc = Chem.FindMolChiralCenters(m, includeUnassigned=True)
    print(k, CalcMolFormula(m), 'chiral' if cc else '')
# 光学活性な C5H12O アルコールを酢酸エステルとして与えるケトン CH3CO-R（R が移動）: R の級数 > メチル（1級以上）
for name, R in (('2-ペンタノール','C(C)CCC'), ('3-メチル-2-ブタノール','C(C)C(C)C'), ('2-メチル-1-ブタノール','CC(C)CC')):
    ket = 'CC(=O)'+R; m = Chem.MolFromSmiles(ket)
    print(name, '<- ketone', ket, CalcMolFormula(m), 'chiral' if Chem.FindMolChiralCenters(m, includeUnassigned=True) else 'achiral', BV(ket), hydrolyze(list(BV(ket))[0]))
# 脱水生成物（隣の炭素の H をとる）とシス-トランス異性
def dehyd(s):
    rxn = AllChem.ReactionFromSmarts('[C:1]([OH:3])[C;!H0:2]>>[C:1]=[C:2].[O:3]')
    out=set()
    for ps in rxn.RunReactants((Chem.MolFromSmiles(s),)):
        p=ps[0]; Chem.SanitizeMol(p); out.add(Chem.MolToSmiles(p))
    return out
for name in ('2-ペンタノール','3-メチル-2-ブタノール','2-メチル-1-ブタノール'):
    for a in dehyd(alcs[name]):
        m = Chem.MolFromSmiles(a)
        from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
        n = len(set(Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)))
        print(name, '->', a, 'stereoisomers', n)
print(check_equation('CH3COOCH(CH3)CH(CH3)2 + H2O -> CH3COOH + CH3CH(OH)CH(CH3)2'))
print(check_equation('C6H10O2 + H2O -> C6H12O3'))
# 環状エステル（ラクトン）の加水分解は手で開環して確認（RDKit テンプレートは環の開裂を正しく扱わないため）
for lac, acid in (('O=C1CCCCCO1','OCCCCCC(=O)O'), ('CC1CCCCC(=O)O1','CC(O)CCCCC(=O)O'), ('O=C1OCCCCC1C','OCCCCC(C)C(=O)O')):
    fl = CalcMolFormula(Chem.MolFromSmiles(lac)); fa = CalcMolFormula(Chem.MolFromSmiles(acid))
    print(lac, fl, '+ H2O ->', acid, fa)
print('cyclohexanone symmetric -> caprolactone', CalcMolFormula(Chem.MolFromSmiles('O=C1CCCCCO1')))
print(check_equation('C7H12O2 + H2O -> C7H14O3'))
