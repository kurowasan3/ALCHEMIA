NAME = 'THE-ION-032'
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

h,K1,K2=symbols('h K1 K2',positive=True)
f=K1*h/(h**2+K1*h+K1*K2)
crit=solve(diff(f,h),h); print(crit); assert simplify(crit[0]-sqrt(K1*K2))==0
fmax=simplify(f.subs(h,sqrt(K1*K2))); assert simplify(fmax-1/(1+2*sqrt(K2/K1)))==0
k1,k2=Float(10)**-3,Float(10)**-4.4
fm=fmax.subs({K1:k1,K2:k2}); print('fmax',fm); assert abs(fm-0.714)<0.001
r=(1/Float(0.714)-1)/2; ratio=1/r**2; print('K1/K2',ratio); assert round(float(ratio))==25
dpk=log(25,10); print(N(dpk)); pk1=3.70-dpk/2; pk2=3.70+dpk/2; print(N(pk1),N(pk2)); assert abs(pk1-3.00)<0.005 and abs(pk2-4.40)<0.005
D=1+1+Rational(1,25); print(N(1/D), N(Rational(1,25)/D)); assert sig(1/D,2)==0.49 and sig(Rational(1,25)/D,2)==0.020
p,hh=symbols('p hh',positive=True)
fm2=maximum(2*p*(1-p), p, Interval(0,1)); print('model max',fm2); assert fm2==Rational(1,2)
k1=hh*2*p*(1-p)/(1-p)**2; k2=hh*p**2/(2*p*(1-p)); print(simplify(k1/k2)); assert simplify(k1/k2)==4
r=solve(1/(1+2*sqrt(1/symbols('R',positive=True)))-Rational(1,2)); print(r); assert r==[4]

check_structure_and_equations()
print('ALL OK')
