# 検証スクリプト ORG-STR-024  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-STR-024/item.json'
import sys, json, re
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ORG-STR-024：MMFF でねじれ角を固定したエネルギーを比較 ----
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
def scan(smi, quad, angles):
    E = {}
    for ang in angles:
        m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=1); AllChem.MMFFOptimizeMolecule(m)
        mp = AllChem.MMFFGetMoleculeProperties(m)
        conf = m.GetConformer(); rdMolTransforms.SetDihedralDeg(conf, *quad, float(ang))
        ff = AllChem.MMFFGetMoleculeForceField(m, mp); ff.MMFFAddTorsionConstraint(*quad, False, float(ang), float(ang), 1e4)
        ff.Minimize(maxIts=2000); E[ang] = ff.CalcEnergy()
    return E
E = scan('CCCC', (0,1,2,3), [0, 60, 120, 180, 240, 300])
e0 = min(E.values()); print({k: round(v-e0,2) for k,v in E.items()}, 'kcal/mol')
assert E[0] > E[120] > E[60] > E[180]
assert abs(E[60]-E[300]) < 0.3
# 1,2-エタンジオール：O-C-C-O が 60° のほうが 180° より低い（分子内水素結合）
Eg = scan('OCCO', (0,1,2,3), [60, 180]); print('ethanediol', {k: round(v,2) for k,v in Eg.items()})
assert Eg[60] < Eg[180]
# エタン：極小 3 つが等しい
Ee = scan('CC', (2,0,1,5), [60,180,300]); print('ethane', {k: round(v,3) for k,v in Ee.items()})
assert max(Ee.values()) - min(Ee.values()) < 0.05
print('ALL OK')
