# 検証スクリプト INO-PPT-032  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/INO-PPT-032/item.json'
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
assert 'Mg^2+ + 2OH- -> Mg(OH)2' in EQS
# 一般式（空ラベル付き）も原子数・電荷を確認
gen = [b for k in ('problem','solution') for b in ce_blocks(it[k]) if '[\\text{}]' in b]
for g in set(gen):
    r = check_equation(clean_ce(g)); print('一般式', g, r or 'OK'); assert r is None
c = R(1,100)
tot = c*R(12,1000); ca = c*R(8,1000); mg = tot - ca
Ca = ca/R(1,10); Mg = mg/R(1,10); print(float(Ca), float(Mg))
assert Ca == R(8,10000) and Mg == R(4,10000)
hard = (Ca+Mg)*100*1000; print('硬度', float(hard)); assert hard == 120
mgmg = Mg*R(243,10)*1000; print('Mg mg/L', float(mgmg)); assert sig(float(mgmg),3) == 9.72
v = R(1,1000)*R(10,1000)/c*1000; assert v == 1
ions = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json', encoding='utf-8'))['results']
assert ions['Mg^2+']['NaOH']['kind'] == 'ppt' and ions['Ca^2+']['NaOH']['kind'] == 'none'

print('全検証OK:', ITEM)
