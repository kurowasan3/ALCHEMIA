# 検証スクリプト THE-GAS-032  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-GAS-032/item.json'
import sys, json, re, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = '原子の構造と電子配置 元素の周期表と周期性 化学結合と結晶 物質量と化学反応式 酸と塩基 酸化還元反応 電池と電気分解 気体の性質 溶液の性質 化学反応とエネルギー 反応速度と化学平衡 電離平衡と溶解平衡'.split()
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
def check_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and all(x.strip() for x in re.split(r'->|<=>|<-', b)):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK', e)
    return sorted(eqs)
def sig(x, n):
    return float(f'{x:.{n}g}')
def check_structure():
    assert it['field'] in FIELDS
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    for f in ('problem_fig', 'solution_fig'):
        if it[f]:
            import os
            p = os.path.join(os.path.dirname(ITEM), it[f]); assert os.path.exists(p), p
            t = open(p).read(); assert 'viewBox' in t and '<script' not in t and 'foreignObject' not in t
    print('小問数', n_items, '答え数', n_ans)
check_structure()
R = 8.31e3; T = 300
Z = {1e5: 1.00, 1e7: 0.85, 2e7: 0.82, 4e7: 1.00, 6e7: 1.20, 8e7: 1.40}
Vi = 1.00*R*T/2.0e7; Vr = Z[2e7]*Vi; print('(1)', Vi, Vr, (1 - Vr/Vi)*100)
assert sig(Vi, 2) == 0.12 and sig(Vr, 2) == 0.10 and round((1 - Z[2e7])*100) == 18
n3 = 8.0e7*10.0/(Z[8e7]*R*T); n3i = 8.0e7*10.0/(R*T); print('(3)', n3, n3i)
assert sig(n3, 2) == 230 and sig(n3i, 2) == 320 and n3 < n3i
n4 = 2.0e7*10.0/(Z[2e7]*R*T); print('(4)', n4, n3/n4); assert sig(n3/n4, 2) == 2.3
assert abs(n3/n4 - 4*0.82/1.40) < 1e-12
assert sig(n4, 3) == 97.8 and sig(n3, 3) == 229

check_equations()
print('全検証OK:', ITEM)
