# 検証スクリプト（グループ g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def check_all_equations(it, extra=()):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|…', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]', b): continue
            if not (' + ' in b or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', b)): continue
            e = clean_ce(b)
            r = check_equation(e)
            assert r is None, (b, r)
            print('式OK:', b); n += 1
    for e in extra:
        r = check_equation(e); assert r is None, (e, r); print('式OK(追加):', e); n += 1
    return n

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def check_parts(it):
    k = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, k + 1)], (k, nums)
    print('小問数OK:', k)
ITEM = '/home/claude/chemk/out/items/POL-RES-030/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); check_all_equations(it)
from rdkit import Chem
from rdkit.Chem import rdchem
# X 側炭素に不対電子を置いたラジカル（R = CH3）の共鳴構造の数
rad = {'エチレン': 'CC[CH2]', '塩化ビニル': 'CC[CH]Cl', 'メタクリル酸メチル': 'CC[C](C)C(=O)OC',
       'イソプレン': 'CC[C](C)C=C', '酢酸ビニル': 'CC[CH]OC(C)=O', 'アクリロニトリル': 'CC[CH]C#N'}
res = {}
for k, s in rad.items():
    m = Chem.MolFromSmiles(s)
    sup = rdchem.ResonanceMolSupplier(m, rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS | rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS | rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS)
    n = len(set(Chem.MolToSmiles(x) for x in sup if x is not None))
    # ラジカルの非局在化を直接調べる：不対電子原子が π 系（共役）に入っているか
    a = [x for x in m.GetAtoms() if x.GetNumRadicalElectrons()][0]
    conj = any(b.GetIsConjugated() for b in a.GetBonds()) or any(nb.GetHybridization() == Chem.HybridizationType.SP2 or nb.GetHybridization() == Chem.HybridizationType.SP for nb in a.GetNeighbors() if nb.GetSymbol() == 'C')
    res[k] = conj; print(k, '共鳴構造数', n, '隣の炭素が sp2/sp:', conj)
assert res == {'エチレン': False, '塩化ビニル': False, 'メタクリル酸メチル': True, 'イソプレン': True, '酢酸ビニル': False, 'アクリロニトリル': True}
# Mayo–Lewis による組成のずれ（AN: r1=3.28, VC: r2=0.02）
r1, r2 = 3.28, 0.02
def F1(f1):
    f2 = 1 - f1; return (r1*f1*f1 + f1*f2)/(r1*f1*f1 + 2*f1*f2 + r2*f2*f2)
print('初期の共重合体の AN 分率', F1(0.5)); assert F1(0.5) > 0.5
N1, N2 = 1.0, 1.0; P1 = P2 = 0.0; hist = []
while N1 + N2 > 1e-3:
    f1 = N1/(N1 + N2); d = 1e-4; P1 += d*F1(f1); P2 += d*(1 - F1(f1)); N1 -= d*F1(f1); N2 -= d*(1 - F1(f1))
    N1 = max(N1, 0); hist.append(f1)
print('残り単量体の AN 分率: 初め', hist[0], '途中', hist[len(hist)//2], '終わり', hist[-1]); assert hist[-1] < hist[len(hist)//2] < hist[0]
print('全体の比 P1:P2 =', P1/P2); assert abs(P1/P2 - 1) < 0.01
print('全検証OK:', ITEM)
