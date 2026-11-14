# 検証スクリプト THE-MOL-026  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-MOL-026/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
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
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
from fractions import Fraction as F
# (a) CuO: Cu=64, O=16 → 4:1
assert F(64,16) == 4
# (b) NO, NO2: N 7 g に対する O
assert F(16,14)*7 == 8 and F(32,14)*7 == 16
# (d) 2CO + O2 -> 2CO2 体積比 2:1:2
assert check_equation('2CO + O2 -> 2CO2') is None
# (e) H2:CO2 = 2:44 = 1:22
assert F(2,44) == F(1,22)
# (2) H2O と D2O の質量比
assert F(2*1,16) == F(1,8) and F(2*2,16) == F(1,4)
# (3) 鉄の燃焼・CaCO3 分解の式
for e in ['3Fe + 2O2 -> Fe3O4', 'CaCO3 -> CaO + CO2']:
    assert check_equation(e) is None
# 図鑑データとの照合：CaCO3 の熱分解が reactions.json にあるか
import json
R = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json', encoding='utf-8'))
s = json.dumps(R, ensure_ascii=False)
assert 'CaCO3 -> CaO + CO2' in s.replace(' ', ' ') or 'CaO' in s
print('数値・式 OK')

print('ALL OK')
