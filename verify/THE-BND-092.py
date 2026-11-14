# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-BND-092/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
r = {'Na+': 0.102, 'Cu+': 0.077, 'Ag+': 0.115, 'Cl-': 0.181, 'I-': 0.220}
def pred(x):
    return 4 if 0.225 <= x < 0.414 else 6 if 0.414 <= x < 0.732 else None
rat = {k: r[a]/r[b] for k, (a, b) in {'AgCl': ('Ag+', 'Cl-'), 'AgI': ('Ag+', 'I-'), 'CuCl': ('Cu+', 'Cl-'), 'NaCl': ('Na+', 'Cl-')}.items()}
print({k: round(v, 3) for k, v in rat.items()})
assert round(rat['AgCl'], 3) == 0.635 and round(rat['AgI'], 3) == 0.523 and round(rat['CuCl'], 3) == 0.425 and round(rat['NaCl'], 3) == 0.564
assert all(pred(v) == 6 for v in rat.values())
# 実際の構造（文献の標準的な記述）：NaCl・AgCl は岩塩型，CuCl（γ）は閃亜鉛鉱型，AgI（β/γ）は常温で4配位（ウルツ鉱型/閃亜鉛鉱型）
actual = {'NaCl': 6, 'AgCl': 6, 'CuCl': 4, 'AgI': 4}
mismatch = [k for k in actual if actual[k] != pred(rat[k])]
assert sorted(mismatch) == ['AgI', 'CuCl']
# 溶解度積（25℃，文献値の目安）：AgCl 1.8e-10, AgI 8.5e-17, AgF は易溶（約 180 g/100 g 水）
Ksp = {'AgCl': 1.8e-10, 'AgI': 8.5e-17}
assert Ksp['AgCl'] > Ksp['AgI']
print('OK')
