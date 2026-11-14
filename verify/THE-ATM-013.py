# 検証スクリプト THE-ATM-013  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-ATM-013/item.json'
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
# 副殻モデル（構成原理＋Cr・Cu の例外）で殻ごとの電子数を計算し，陽イオンは n の大きい殻（同じ n では l の大きい順）から除く
ORDER = ['1s','2s','2p','3s','3p','4s','3d','4p']
CAP = {'s':2,'p':6,'d':10}
def config(Z):
    c = {}; left = Z
    for o in ORDER:
        k = min(left, CAP[o[1]]); 
        if k: c[o] = k
        left -= k
    if Z in (24, 29):   # Cr, Cu
        c['4s'] -= 1; c['3d'] += 1
    return c
def ionize(c, q):
    c = dict(c)
    for _ in range(q):
        o = max((x for x in c if c[x] > 0), key=lambda x: (int(x[0]), 'spd'.index(x[1])))
        c[o] -= 1
    return c
def shells(c):
    s = {}
    for o, k in c.items():
        if k: s['KLMN'[int(o[0])-1]] = s.get('KLMN'[int(o[0])-1], 0) + k
    return ' '.join(f'{k}{v}' for k, v in s.items())
assert shells(config(26)) == 'K2 L8 M14 N2' and shells(config(29)) == 'K2 L8 M18 N1'
assert shells(ionize(config(26),2)) == 'K2 L8 M14'
assert shells(ionize(config(26),3)) == 'K2 L8 M13'
parts = {'Ni':(28,0),'Cu+':(29,1),'Zn2+':(30,2),'Ga3+':(31,3),'Fe2+':(26,2)}
res = {k: shells(ionize(config(Z),q)) for k,(Z,q) in parts.items()}
print(res)
m18 = [k for k,v in res.items() if 'M18' in v]
assert m18 == ['Cu+','Zn2+','Ga3+'], m18
assert res['Ni'] == 'K2 L8 M16 N2'
assert all(Z-q == 28 for k,(Z,q) in parts.items() if k != 'Fe2+')
mn2 = shells(ionize(config(25),2)); print('Mn2+', mn2)
cands = [(Z,3) for Z in (25,26,28,29,30,31) if shells(ionize(config(Z),3)) == mn2]
assert cands == [(26,3)]
# 図鑑データ elements.json の殻ごとの電子数と照合
import json
E = {e['sym']: e for e in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/elements.json', encoding='utf-8'))['elements']}
for sym, Z in [('Mn',25),('Fe',26),('Ni',28),('Cu',29),('Zn',30),('Ga',31)]:
    assert shells(config(Z)) == ' '.join(f'{L}{n}' for L, n in zip('KLMN', E[sym]['shells'])), sym
print('elements.json と一致')

print('ALL OK')
