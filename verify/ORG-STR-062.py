ITEM = '/home/claude/chemx/tokodai/outh/items/ORG-STR-062/item.json'
# 検証スクリプト（東工大 g03）  python3 -I この.py
# 共通の検証ヘッダ（各検証スクリプトにそのまま埋め込む）
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j
eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
acids={'シュウ酸':'OC(=O)C(=O)O','マロン酸':'OC(=O)CC(=O)O','コハク酸':'OC(=O)CCC(=O)O','メチルマロン酸':'OC(=O)C(C)C(=O)O',
       'グルタル酸':'OC(=O)CCCC(=O)O','2-メチルコハク酸':'OC(=O)C(C)CC(=O)O','ジメチルマロン酸':'OC(=O)C(C)(C)C(=O)O','エチルマロン酸':'OC(=O)C(CC)C(=O)O'}
alk=['C','CC','CCC','C(C)C']
found=set()
for a in acids.values():
    # モノエステル：最初の O-H の O にアルキルをつける（2 通りの COOH を試す）
    for i in (0,1):
        for R in alk:
            s=a.replace('OC(=O)','OC(=O)',1)
            m=Chem.MolFromSmiles(a); rw=Chem.RWMol(m)
            ohs=[x.GetIdx() for x in m.GetAtoms() if x.GetSymbol()=='O' and x.GetTotalNumHs()==1]
            o=ohs[i]; r=Chem.MolFromSmiles(R)
            comb=Chem.RWMol(Chem.CombineMols(m,r)); comb.AddBond(o,m.GetNumAtoms(),Chem.BondType.SINGLE)
            mm=comb.GetMol(); Chem.SanitizeMol(mm)
            if CalcMolFormula(mm)=='C5H8O4': found.add(Chem.MolToSmiles(mm))
            for R2 in alk:
                c2=Chem.RWMol(Chem.CombineMols(mm,Chem.MolFromSmiles(R2)))
                o2=[x.GetIdx() for x in mm.GetAtoms() if x.GetSymbol()=='O' and x.GetTotalNumHs()==1]
                if not o2: continue
                c2.AddBond(o2[0],mm.GetNumAtoms(),Chem.BondType.SINGLE); m2=c2.GetMol(); Chem.SanitizeMol(m2)
                if CalcMolFormula(m2)=='C5H8O4': found.add(Chem.MolToSmiles(m2))
print(len(found),sorted(found)); assert len(found)==7
chiral=[s for s in found if Chem.FindMolChiralCenters(Chem.MolFromSmiles(s),includeUnassigned=True,useLegacyImplementation=False)]
print(chiral); assert len(chiral)==1 and Chem.MolToSmiles(Chem.MolFromSmiles('OC(=O)C(C)C(=O)OC'))==chiral[0]
for s in ('OC(=O)C(C)C(=O)O','OC(=O)CCC(=O)OC'):
    assert not Chem.FindMolChiralCenters(Chem.MolFromSmiles(s),includeUnassigned=True,useLegacyImplementation=False)
chk('DoU',(2*5+2-8)/2,2,1)
# P・Q・R の同定
props={'a':('di','シュウ酸',['メタノール','エタノール']),'b':('di','マロン酸',['メタノール']),'c':('mono','シュウ酸',['1-プロパノール']),
 'd':('mono','シュウ酸',['2-プロパノール']),'e':('mono','マロン酸',['エタノール']),'f':('mono','コハク酸',['メタノール']),'g':('mono','メチルマロン酸',['メタノール'])}
iodo={'エタノール','2-プロパノール'}; ketone={'2-プロパノール'}
P=[k for k,(t,a,al) in props.items() if t=='di' and any(x in iodo for x in al)]
Q=[k for k,(t,a,al) in props.items() if t=='mono' and a=='コハク酸']
R=[k for k,(t,a,al) in props.items() if t=='mono' and a=='シュウ酸' and any(x in ketone for x in al)]
print(P,Q,R); assert P==['a'] and Q==['f'] and R==['d']
print('ALL OK')
