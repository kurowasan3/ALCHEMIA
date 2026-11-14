NAME = 'THE-ION-028'
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

Kw=Rational(1,10**14); C=Rational(1,10**5)
Ka=Rational(1,10**7)*Rational(16,100)/(Rational(80,100)-Rational(16,100)); print('Ka',Ka); assert Ka==Rational(25,10**9)
pKa=-log(Ka,10); print(N(pKa)); assert abs(pKa-7.60)<0.003
pH2=pKa+log(Rational(56,24),10); print('pH2',N(pH2)); assert abs(pH2-7.97)<0.005
p=symbols('p'); h=Float(10)**(-p); In=C*Ka/(Ka+h)
ps=nsolve(h-In-Kw/h, p, 6.3); print('pH3 exact',ps); assert round(float(ps),1)==6.3
pa=-log(sqrt(Ka*C+Kw),10); print('approx',N(pa)); assert round(float(pa),1)==6.3
hh=Float(10)**(-ps); A=0.8*Ka/(Ka+hh); print('A',A)
Cmax=Rational(26,100)*Kw/Ka; print('Cmax',N(Cmax)); assert sig(Cmax,2)==1.0e-7
# 厳密にも pH 低下 <=0.05 か
ps2=nsolve(h-Cmax*Ka/(Ka+h)-Kw/h, p, 7); print('drop',7-ps2); assert 7-ps2<=0.0505

check_structure_and_equations()
print('ALL OK')
