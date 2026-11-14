# 検証スクリプト INO-PPT-033  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/INO-PPT-033/item.json'
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

import itertools
assert 'Al(OH)3 + OH- -> [Al(OH)4]-' in EQS and 'Zn(OH)2 + 4NH3 -> [Zn(NH3)4]^2+ + 2OH-' in EQS
C = ['HCl','NaOH','NH3','AgNO3','Zn(NO3)2','Al(NO3)3','Mg(NO3)2','BaCl2','Na2SO4']
# 規則（図鑑 ions.json と照合する）
cat = {'HCl':'Cl','BaCl2':'Cl','Na2SO4':'SO4','NaOH':'NaOH','NH3':'NH3'}
metal = {'AgNO3':'Ag+','Zn(NO3)2':'Zn^2+','Al(NO3)3':'Al^3+','Mg(NO3)2':'Mg^2+','BaCl2':'Ba^2+'}
ions = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json', encoding='utf-8'))['results']
def mix(x, y):
    """(観察, 過剰に加えると溶ける側の集合)．図鑑のデータから組み立てる"""
    for m, r in ((x, y), (y, x)):
        if m in metal and r in cat:
            key = cat[r]
            if metal[m] == 'Ba^2+' and key == 'Cl': continue
            d = ions[metal[m]][key]
            if d['kind'] == 'ppt':
                col = {'白':'白','褐色':'褐'}[d['color']]
                ex = set()
                if key in ('NaOH','NH3') and ions[metal[m]][key+'-ex']['kind'] == 'dissolve': ex.add(r)
                return (col, frozenset(ex))
    return ('−', frozenset())
# 手で確認した規則と一致すること
assert mix('AgNO3','NaOH') == ('褐', frozenset()) and mix('AgNO3','NH3') == ('褐', frozenset({'NH3'}))
assert mix('Zn(NO3)2','NH3') == ('白', frozenset({'NH3'})) and mix('Al(NO3)3','NH3') == ('白', frozenset())
assert mix('Al(NO3)3','NaOH') == ('白', frozenset({'NaOH'})) and mix('Mg(NO3)2','NaOH') == ('白', frozenset())
assert mix('AgNO3','HCl')[0] == '白' and mix('AgNO3','BaCl2')[0] == '白' and mix('BaCl2','Na2SO4')[0] == '白'
assert mix('AgNO3','Na2SO4')[0] == '−' and mix('BaCl2','NaOH')[0] == '−' and mix('NaOH','NH3')[0] == '−'
L = 'ABCDEF'; pairs = list(itertools.combinations(L, 2))
true = dict(zip(L, ['AgNO3','Mg(NO3)2','NaOH','Al(NO3)3','NH3','Zn(NO3)2']))
table = {('A','C'):'褐',('A','E'):'褐',('B','C'):'白',('B','E'):'白',('C','D'):'白',('C','F'):'白',('D','E'):'白',('E','F'):'白'}
for p in pairs: assert mix(true[p[0]], true[p[1]])[0] == table.get(p, '−'), p
exps = {'a':(('A','C'), ()), 'b':(('C','D'), ('C',)), 'c':(('E','F'), ('E',))}
def consistent(asg, use):
    for p in pairs:
        if mix(asg[p[0]], asg[p[1]])[0] != table.get(p, '−'): return False
    for k in use:
        (a, b), diss = exps[k]
        ex = mix(asg[a], asg[b])[1]
        if tuple(l for l in (a, b) if asg[l] in ex) != diss: return False
    return True
def count(use):
    return [p for p in itertools.permutations(C, 6) if consistent(dict(zip(L, p)), use)]
s0 = count(''); s_all = count('abc'); s_bc = count('bc')
print(len(s0), s_all, s_bc)
assert len(s0) == 12 and all(set(p) & {'HCl','BaCl2','Na2SO4'} == set() for p in s0)
assert s_all == [tuple(true[l] for l in L)]
assert sorted(s_bc) == sorted([tuple(true[l] for l in L), ('AgNO3','Mg(NO3)2','NH3','Zn(NO3)2','NaOH','Al(NO3)3')])

print('全検証OK:', ITEM)
