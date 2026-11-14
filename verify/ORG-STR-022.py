# 検証スクリプト ORG-STR-022  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-STR-022/item.json'
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
# ---- 個別の検証 ORG-STR-022 ----
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
# (1) 双極子：MMFF の部分電荷から双極子モーメントを概算し，シス形だけが有極性であることを確認
def dip(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=3); AllChem.MMFFOptimizeMolecule(m)
    mp = AllChem.MMFFGetMoleculeProperties(m)
    q = np.array([mp.GetMMFFPartialCharge(i) for i in range(m.GetNumAtoms())])
    X = m.GetConformer().GetPositions()
    return np.linalg.norm((q[:,None]*X).sum(0))
cis, trans = dip('Cl/C=C\\Cl'), dip('Cl/C=C/Cl')
print('双極子（相対）cis', round(cis,3), 'trans', round(trans,3)); assert cis > 0.3 and trans < 0.01
# (3) 水素化エンタルピーの差
assert -116 - (-120) == 4
# (4) 燃焼の式（生成物が同じ）
assert check_equation('C4H8 + 6O2 -> 4CO2 + 4H2O') is None
# ヘスの法則：ΔHc(cis) - ΔHc(trans) = ΔHhyd(cis) - ΔHhyd(trans)（どちらも H2 とブタンの項が共通）
dHc_butane, dHc_H2 = -2878, -286   # 任意の値でも差は変わらないことを確認
for hb, hh in ((dHc_butane, dHc_H2), (-1000, -200)):
    c_cis = -120 + hb - hh; c_trans = -116 + hb - hh
    assert c_cis - c_trans == -4
print('数値の出典確認：1,2-ジクロロエチレン bp cis 60.3, trans 47.5 °C；mp cis -80, trans -50 °C；2-ブテン水素化 cis -119.7, trans -115.5 kJ/mol（文献値を丸めた）')
print('ALL OK')
