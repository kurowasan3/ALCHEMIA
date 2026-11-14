NAME = 'THE-ION-031'
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

Ksp=Rational(2,10**10); beta=Rational(125,10)*10**6; Ka=Rational(5,10**10)
K=Ksp*beta; assert K==Rational(25,10**4)
assert Rational(1,10**6)/Rational(2,10)**2==Rational(25,10**6)
def onset(nAg, NT):
    Ag=Ksp/nAg; NH3=sqrt(nAg/(beta*Ag)); NH4=NT-2*nAg-NH3; return NH3, NH4, -log(Ka,10)+log(NH3/NH4,10)
NH3,NH4,pH=onset(Rational(1,1000), Rational(202,1000)); print(NH3,NH4,N(pH)); assert NH3==Rational(2,100) and abs(pH-8.35)<0.006 and NH4==Rational(18,100)
NH3b,NH4b,pHb=onset(Rational(1,10**4), Rational(2002,10000)); print(N(pHb)); assert abs(pHb-7.30)<0.006
# (3) 厳密: pH=7.30 固定, 未知 y (溶存銀=Cl-), NH3
y,n=symbols('y n',positive=True); H=Float(10)**(-7.30)
NT=Rational(202,1000)
eqs=[ y*y/(1+beta*n**2) - Ksp,                         # [Ag+][Cl-]=Ksp, [Ag+]=y/(1+beta n^2)
      n*(1+H/Ka) + 2*y*beta*n**2/(1+beta*n**2) - NT ]  # N 収支
sol=nsolve(eqs,[y,n],[1e-4,2e-3]); print(sol.T, 'pct', sol[0]/Rational(1,1000)*100); assert sig(sol[0]/Rational(1,1000)*100,2)==10

check_structure_and_equations()
print('ALL OK')
