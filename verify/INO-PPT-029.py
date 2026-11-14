# 検証スクリプト（共通ヘッダ＋問題固有の検証）  python3 -I <this file>
import json, os, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation, parse_formula
ID = os.path.splitext(os.path.basename(__file__))[0]
ITEM = json.load(open(f'/home/claude/chem/out2/items/{ID}/item.json', encoding='utf-8'))
def _ce(s):
    out = []; i = 0
    while (i := s.find('\\ce{', i)) >= 0:
        j = i + 4; dep = 1
        while dep:
            dep += {'{': 1, '}': -1}.get(s[j], 0); j += 1
        out.append(s[i+4:j-1]); i = j
    return out
EQS = []
for _k in ('problem', 'answer', 'solution'):
    for _c in _ce(ITEM[_k]):
        if len(re.findall(r'<=>|->|<-', _c)) == 1 and not re.search(r'\\boxed|\\phantom|\\text|\bR\b|\bn\b', _c):
            EQS.append(_c)
for _c in EQS:
    _r = check_equation(_c)
    assert _r is None, (_c, _r)
print(f'{ID}: 反応式 {len(EQS)} 本はすべて原子数・電荷がつり合っている')
# 小問数と答えの番号の一致
_nq = ITEM['problem'].count('\\item')
_na = len(re.findall(r'(?:^|\s)\((\d+)\)', ITEM['answer']))
assert _nq == _na or (_nq == 0 and _na == 0), (_nq, _na)
print(f'{ID}: 小問 {_nq} 個と答え {_na} 個が対応')
M = dict(H=1.0, C=12.0, N=14.0, O=16.0, Na=23.0, Mg=24.0, Al=27.0, Si=28.0, P=31.0, S=32.0, Cl=35.5, K=39.0,
         Ca=40.0, Ti=48.0, Cr=52.0, Mn=55.0, Fe=56.0, Co=59.0, Cu=63.5, Zn=65.4, Br=80.0, Sr=88.0, Ag=108.0, I=127.0, Ba=137.0, Pt=195.0)
def mass(f, masses=None):
    m = masses or M
    c, _ = parse_formula(f)
    return sum(m[e] * float(n) for e, n in c.items())
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
# ---- 問題固有の検証 ----
# 平面正方形：位置 0,1,2,3（0-2, 1-3 がトランス）
TE = {'NO2': 3, 'Cl': 2, 'NH3': 1}
LEAVE = {'Cl': 0, 'NH3': 1, 'NO2': 2}   # 小さいほど脱離しやすい（同点のとき）
def step(c, L):
    cands = []
    for i in range(4):
        if c[i] == L: continue          # 変化しない置換は除く
        partner = c[(i+2) % 4]
        cands.append((LEAVE[c[i]], -TE[partner], i))   # R1: Cl を優先して置換, R2: トランス効果最大の向かい
    cands.sort()
    best = cands[0][:2]
    outs = set()
    for k in cands:
        if k[:2] == best:
            n = list(c); n[k[2]] = L; outs.add(canon(n))
    return outs
def canon(c):
    forms = []
    for r in range(4):
        x = c[r:] + c[:r]; forms.append(tuple(x)); forms.append(tuple(x[::-1]))
    return min(forms)
def run(start, seq):
    states = {canon(start)}
    for L in seq:
        new = set()
        for s in states: new |= step(list(s), L)
        states = new
    return states
def describe(c):
    c = list(c); out = []
    for i in range(2):
        out.append(f'{c[i]}-trans-{c[i+2]}')
    return ', '.join(out)
# 規則 1（Cl を優先して置換、変化のない置換は除く）・規則 2（向かい側のトランス効果最大）による全状態の列挙
PtCl4 = ['Cl']*4; PtA4 = ['NH3']*4
def res(st, seq):
    r = run(st, seq); assert len(r) == 1, r   # 生成物が一意に決まること
    return describe(next(iter(r)))
R = {k: res(*v) for k, v in {
    'A': (PtCl4, ['NH3','NH3']), 'B': (PtA4, ['Cl','Cl']),
    '2a': (PtCl4, ['NH3','NO2']), '2b': (PtCl4, ['NO2','NH3']),
    '3a': (PtA4, ['NO2','Cl']), '3b': (PtA4, ['Cl','NO2']),
    '4i': (PtCl4, ['NH3','NH3','NO2']), '4ii': (PtCl4, ['NH3','NO2','NH3']), '4iii': (PtCl4, ['NO2','NH3','NH3'])}.items()}
for k, v in R.items(): print(k, v)
assert R['A'] == 'Cl-trans-NH3, Cl-trans-NH3'          # シス形
assert R['B'] == 'Cl-trans-Cl, NH3-trans-NH3'          # トランス形
assert R['2a'] == 'Cl-trans-NH3, Cl-trans-NO2'         # NH3 と NO2 がシス
assert R['2b'] == 'Cl-trans-Cl, NH3-trans-NO2'         # NH3 と NO2 がトランス
assert R['3a'] == 'Cl-trans-NO2, NH3-trans-NH3'
assert R['3b'] == 'NH3-trans-NH3, NH3-trans-NO2'       # [Pt(NH3)3(NO2)]+
assert R['4i'] == R['4ii'] == R['4iii'] == 'Cl-trans-NH3, NH3-trans-NO2'   # NH3 どうしシス
# [PtCl(NH3)2(NO2)] の幾何異性体が 2 種類であること（平面正方形の全配置を正規化して数える）
from itertools import permutations
isos = {canon(list(p)) for p in permutations(['Cl','NH3','NH3','NO2'])}
print('isomers', isos); assert len(isos) == 2
assert canon(['NO2','NH3','Cl','NH3']) != canon(['NO2','NH3','NH3','Cl'])
print('ALL OK')
