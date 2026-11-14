# 検証スクリプト THE-ION-045  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-ION-045/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    return re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', s)
def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq
def should_check(eq):
    # 矢印が1つの完全な式はすべて（validate ツールが一般式とみなして飛ばす Zn・Sn なども含めて）検査する．
    # 一般式・同位体の式は矢印に空ラベル [\\text{}] を付けてあるので除く．
    if re.search(r'boxed|□|\?|？|…|\\square', eq): return False
    if re.search(r'(->|<=>|<-)\[\\text\{\}\]', eq): return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1: return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None
def check_all_equations():
    eqs = sorted({b for k in ('problem','answer','solution') for b in ce_blocks(it[k]) if should_check(b)})
    for e in eqs:
        r = check_equation(clean_ce(e)); print('式', e, '→', r or 'OK'); assert r is None, e
    return eqs
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    nums = [int(x) for x in re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, n_items+1)), nums
    print('小問数', n_items)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
EQS = check_all_equations()

# ---- 個別の検証 ----

from sympy import Rational as R
import math
Ksp = 1.1e-10
n = 0.233/233; c = n/0.0500; print(c); assert abs(c - 2.00e-2) < 1e-12
s = 1.05e-5; assert abs(math.sqrt(Ksp) - s) < 1e-7
loss = s*0.200*233; print('loss g', loss); assert sig(loss*1000,2) == 0.49
assert sig(loss/0.233*100,2) == 0.21
s2 = Ksp/1.0e-3; assert s2 < 1e-3/1000
loss2 = s2*0.200*233; print('loss2 mg', loss2*1000); assert sig(loss2*1000,2) == 5.1e-3
print('ratio', s/s2); assert round(s/s2) == 95
# 厳密解（SO4 の増加も含める）との比較
x = (-1e-3 + math.sqrt(1e-6 + 4*Ksp))/2; assert abs(x/s2 - 1) < 1e-3
one = 1/201; four = (1/51)**4; print(one, four); assert sig(one,2) == 5.0e-3 and sig(four,2) == 1.5e-7
ions = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json', encoding='utf-8'))['results']
assert ions['Ba^2+']['SO4']['kind'] == 'ppt' and ions['Ba^2+']['SO4']['color'] == '白'

print('全検証OK:', ITEM)
