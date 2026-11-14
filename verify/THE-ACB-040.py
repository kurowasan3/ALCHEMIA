# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/out/items/THE-ACB-040/item.json'
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

# 酸化数の計算（各化学種の構成原子の酸化数を規則で与えて，前後で同じかを確認）
ox = {
 'Mg3N2': {'Mg': 2, 'N': -3}, 'H2O': {'H': 1, 'O': -2}, 'Mg(OH)2': {'Mg': 2, 'O': -2, 'H': 1}, 'NH3': {'N': -3, 'H': 1},
 'CaO': {'Ca': 2, 'O': -2}, 'SO2': {'S': 4, 'O': -2}, 'CaSO3': {'Ca': 2, 'S': 4, 'O': -2},
 'Zn(OH)2': {'Zn': 2, 'O': -2, 'H': 1}, 'OH-': {'O': -2, 'H': 1}, '[Zn(OH)4]2-': {'Zn': 2, 'O': -2, 'H': 1},
 'HSO4-': {'H': 1, 'S': 6, 'O': -2}, 'CH3COO-': {}, 'SO4 2-': {'S': 6, 'O': -2}, 'CH3COOH': {},
}
# 各化学種の電荷と酸化数の和が一致するか（与えた酸化数の整合性）
comp = {'Mg3N2': ({'Mg':3,'N':2},0), 'Mg(OH)2': ({'Mg':1,'O':2,'H':2},0), 'NH3': ({'N':1,'H':3},0), 'H2O': ({'H':2,'O':1},0),
        'CaO': ({'Ca':1,'O':1},0), 'SO2': ({'S':1,'O':2},0), 'CaSO3': ({'Ca':1,'S':1,'O':3},0),
        'Zn(OH)2': ({'Zn':1,'O':2,'H':2},0), 'OH-': ({'O':1,'H':1},-1), '[Zn(OH)4]2-': ({'Zn':1,'O':4,'H':4},-2),
        'HSO4-': ({'H':1,'S':1,'O':4},-1), 'SO4 2-': ({'S':1,'O':4},-2)}
for sp,(c,q) in comp.items():
    assert sum(ox[sp][e]*n for e,n in c.items()) == q, sp
reactions = {'a': (['Mg3N2','H2O'], ['Mg(OH)2','NH3']), 'b': (['CaO','SO2'], ['CaSO3']),
             'c': (['Zn(OH)2','OH-'], ['[Zn(OH)4]2-']), 'd': (['HSO4-'], ['SO4 2-'])}
for k,(L,R) in reactions.items():
    def states(side):
        s = {}
        for sp in side:
            for e,v in ox[sp].items(): s.setdefault(e,set()).add(v)
        return s
    l, r = states(L), states(R)
    for e in l:
        if e in r: assert l[e] == r[e], (k, e, l[e], r[e])
    print('酸化数不変:', k)
# H+ の移動の有無：(a)(d) は H をもつ化学種の H の数が変化，(b)(c) は H をもつ化学種間の H の受け渡しなし
assert '(b)，(c)' in it['answer'] and 'なし' in it['answer']

print('ALL OK')
