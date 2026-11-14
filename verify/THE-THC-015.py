# 検証スクリプト（g03 作成）：python3 -I このファイル
import sys, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import *
import re
ID = __file__.replace('\\', '/').split('/')[-1][:-3]
ITEM = json.load(open(f'/home/claude/chem/out/items/{ID}/item.json', encoding='utf-8'))
ARR = re.compile(r'<=>|->')


def ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0:
            return out
        k, d = j + 4, 1
        while d:
            d += {'{': 1, '}': -1}.get(s[k], 0)
            k += 1
        out.append(s[j + 4:k - 1])
        i = k


def check_all_equations():
    """問題文・答え・解説の \\ce{…} のうち矢印1つの完全な式をすべて検査"""
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(ITEM[key]):
            if len(ARR.findall(b)) != 1 or '\\boxed' in b or re.search(r'(?<![A-Za-z])[Rn](?![a-z])', b):
                continue
            r = check_equation(b)
            assert r is None, (key, b, r)
            n += 1
    print(f'反応式 {n} 本：原子数・電荷のつり合い OK')


def near(a, b, rel=0.005):
    a, b = float(a), float(b)
    assert abs(a - b) <= rel * abs(b) + 1e-12, (a, b)


def rnd(x, n):
    """有効数字 n 桁に丸めた値"""
    x = float(x)
    return float(f'{x:.{n}g}')

# ---------------------------------------------------------------- 本体
check_all_equations()
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def bonds(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); c = {}
    for b in m.GetBonds():
        k = '-'.join(sorted([b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()])) + str(b.GetBondTypeAsDouble())
        c[k] = c.get(k, 0) + 1
    return CalcMolFormula(m), c
fb, bb = bonds('C=CCC'); fc, bc = bonds('C1CCC1'); print(fb, bb, fc, bc)
assert fb == fc == 'C4H8' and bb == {'C-C2.0':1, 'C-C1.0':2, 'C-H1.0':8} and bc == {'C-C1.0':4, 'C-H1.0':8}
assert check_equation('C4H8 + 6O2 -> 4CO2 + 4H2O') is None
P = 4*-394 + 4*-286; fB = P + 2717; fC = P + 2745; print(fB, fC, fB - fC); assert (fB, fC, fB-fC) == (-3, 25, -28)
at = 4*717 + 4*436; eB = at - (612 + 2*348 + 8*416); eC = at - (4*348 + 8*416)
print(eB, eC, fB - eB, fC - eC); assert (eB, eC) == (-24, -108)
