# 検証スクリプト THE-BND-045  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-BND-045/item.json'
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
from sympy import symbols, solve, Rational as Q
E, f = symbols('E f')
nHB = 4 * Q(1, 2); assert nHB == 2
Ev = solve(Q(51) - 11 - nHB*E, E)[0]; print('E', Ev); assert Ev == 20
tot = nHB*Ev
f0 = 1 - Q(6)/tot; print('0℃ 残存', f0); assert sig(float(f0)*100, 2) == 85
# ヘスの法則との整合（0℃の蒸発熱）
assert 11 + f0*tot == 51 - 6
f100 = solve(Q(41) - 11 - tot*f, f)[0]; print('100℃ 残存', f100); assert sig(float(f100)*100, 2) == 75
heat = 75*100/1000; print('加熱 kJ', heat); assert sig(heat, 2) == 7.5
broken = nHB*(f0 - f100)*Ev; print('切断 kJ', broken); assert sig(float(broken), 2) == 4.0
print('割合', float(broken)/heat)
# 経路の整合性（実測値ベースの確認）：氷0℃→蒸気100℃
p1 = 6.0 + 7.5 + 41; p2 = 51 + 34*100/1000; print(p1, p2); assert abs(p1 - p2) < 0.5
# 実測値との比較（参考）：ΔHsub(0℃)=51.1, ΔHfus=6.01, ΔHvap(100℃)=40.66, ΔHvap(0℃)=45.05, Cp(l)=75.3, Cp(s)=38 J/molK
assert abs(51.1 - 51) < 0.5 and abs(40.66 - 41) < 0.5 and abs(45.05 - 45) < 0.5

check_equations()
print('全検証OK:', ITEM)
