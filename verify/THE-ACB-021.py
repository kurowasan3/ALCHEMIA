NAME = 'THE-ACB-021'
# ---- 共通チェック（単独で動く） ----
from sympy import *
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM_DIR = '/home/claude/chem/out2/items/' + NAME
it = json.load(open(ITEM_DIR + '/item.json', encoding='utf-8'))
def _ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k, d = j + 4, 1
        while d:
            d += {'{': 1, '}': -1}.get(s[k], 0); k += 1
        out.append(s[j + 4:k - 1]); i = k
def check_structure_and_equations():
    for key in ['src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig']:
        assert key in it, key
    assert it['difficulty'] in (4, 5) and 20 <= it['time_min'] <= 40
    fm = set('qa equation compare condition pathway calc structure'.split())
    assert all(f in fm for f in it['format'].split(',')) and 1 <= len(it['format'].split(',')) <= 2
    nitem = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i + 1) for i in range(nitem)], (nitem, nums)
    assert '\\hfill$\\square$' in it['solution']
    n, bad = 0, []
    for key in ('problem', 'answer', 'solution'):
        for b in _ce_blocks(it[key]):
            arr = re.findall(r'->|<=>', b)
            if len(arr) != 1 or '\\boxed' in b: continue
            l, r = re.split(r'->|<=>', b)
            if not l.strip() or not r.strip(): continue
            res = check_equation(b)
            if res and '読めません' in res and re.search(r'（(A|M|X|R|In)\s*の所）', res):
                continue   # 一般式（HA, M, X など）はつり合いの機械チェックの対象外
            n += 1
            if res: bad.append((b, res))
    assert not bad, bad
    for f in (it['problem_fig'], it['solution_fig']):
        if f:
            svg = open(ITEM_DIR + '/' + f, encoding='utf-8').read()
            assert 'viewBox' in svg and '<script' not in svg and 'foreignObject' not in svg and '<image' not in svg
            assert not re.search(r'\son[a-z]+=', svg) and 'href="http' not in svg
    print(f'{NAME}: 構造OK, 反応式 {n} 本つり合いOK')
def sig(v, n):
    return float(f'{float(v):.{n}g}')
# ---- ここまで共通 ----

K1,K2,K3,Kw=Rational(7,10**3),Rational(63,10**9),Rational(45,10**14),Rational(1,10**14)
a,b=symbols('a b'); s=solve([a+b-Rational(18,10), b-Rational(8,10)],[a,b]); print(s); assert s[a]/20==Rational(5,100) and s[b]/20==Rational(4,100)
x=symbols('x',positive=True); C=Rational(8,10)/30
xs=solve(x**2-K1*(C-x),x)[0]; print('frac',N(xs/C)); assert sig(xs/C*100,2)==40
Kh=Kw/K3; print('Kh',N(Kh)); Cp=Rational(8,10)/54; print('Cp',N(Cp))
xs=solve(x**2-Kh*(Cp-x),x)[0]; print('alpha',N(xs/Cp), 'pH',N(14+log(xs,10))); assert sig(xs/Cp*100,2)==69 and round(float(14+log(xs,10)),0)==12
# 厳密な滴定曲線で跳躍の有無
p=symbols('p'); h=Float(10)**(-p)
def pH_at(V):
    Vt=20+V; Na=Rational(1,10)*V/Vt; Cl=Rational(1,1000)/Vt*1000; P=Rational(8,10)/Vt
    D=h**3+K1*h**2+K1*K2*h+K1*K2*K3
    neg=P*(K1*h**2+2*K1*K2*h+3*K1*K2*K3)/D
    f=lambdify(p, Na+h-Cl-neg-Kw/h)
    lo,hi=0.0,14.0
    for _ in range(80):
        m=(lo+hi)/2
        if f(m)>0: lo=m
        else: hi=m
    return (lo+hi)/2
def slope(V): return (pH_at(V+0.5)-pH_at(V-0.5))
for V in [10,18,26,34]: print(V, pH_at(V), 'dpH/mL', slope(V))
assert slope(18)>2*slope(10) and slope(26)>3*slope(34)
assert Rational(8,10)/Rational(1,10)==8

check_structure_and_equations()
print('ALL OK')
