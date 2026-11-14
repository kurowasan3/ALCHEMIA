NAME = 'THE-ION-033'
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

Ksp=Rational(1,10**20); c=Rational(1,10**22)
x=symbols('x',positive=True)
for c0,frac,pH in [(Rational(5,100),20,1.70),(Rational(5,1000),50,2.30)]:
    xs=solve((c0-x)*c/(2*x)**2-Ksp, x)[0]; print(xs, xs/c0*100, N(-log(2*xs,10)))
    assert xs/c0*100==frac and abs(-log(2*xs,10)-pH)<0.002
    # H2S 自身の H+ (sqrt(K1*0.1)=1e-4) を足しても大きくは変わらない
    xs2=[r for r in solve((c0-x)*c/Ksp-(2*x+Rational(1,10**4))**2, x) if r>0][0]; print(' with H2S H+', N(xs2/c0*100)); assert abs(xs2/c0*100-frac)<1.5
M=Ksp/(c/Rational(1,10**5)**2); print('M left',M); assert M==Rational(1,10**8)
a=symbols('a'); asol=solve(Rational(1,10)/(a-Rational(1,10))-Rational(1,2), a)[0]; print('NaOAc',asol); assert asol==Rational(3,10)
assert check_equation('Zn^2+ + H2S -> ZnS + 2H+') is None and check_equation('CH3COO- + H+ -> CH3COOH') is None

check_structure_and_equations()
print('ALL OK')
