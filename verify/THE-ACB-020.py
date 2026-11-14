NAME = 'THE-ACB-020'
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

Kw=Rational(1,10**14)
def E_exact(pH, c0, Ka):
    # c0: 酸・塩基の濃度, 体積 V0=20mL. 加えた塩基 Vb に対し pH が与えられた値になる Vb を厳密に解く
    H=Float(10)**(-pH); OH=Kw/H; Vb=symbols('Vb',positive=True); V0=20
    C=c0*V0/(V0+Vb); cB=c0*Vb/(V0+Vb)
    sol=nsolve(cB+H-(C*Ka/(Ka+H))-OH, Vb, 20)
    return (sol-V0)/V0
Ka=Rational(2,10**5)
for pH,c0,ref in [(9.70,Rational(1,10),0.099),(7.00,Rational(1,10),-0.50),(9.00,Rational(1,10),0.015),(9.00,Rational(1,100),0.20)]:
    e=E_exact(pH,c0,Ka)*100; print(pH,c0,N(e)); assert abs(sig(e,2)-ref)<1e-9 or abs(e-ref)<0.0006
h=symbols('h',positive=True)
for c,ref in [(Rational(1,20),8.70),(Rational(1,200),8.20)]:
    p=symbols('p'); hh=Float(10)**(-p)
    ps=nsolve(c*Ka/(Ka+hh)+Kw/hh-hh-c, p, 8.5)   # 中和点: Na+=C
    print('eq pH', ps); assert abs(ps-ref)<0.006
Kmin=Float(1e-9)*(1/Float(1.2e-3)-1); print('Kmin',Kmin); assert sig(Kmin,2)==8.3e-7
# Kmin で実際に E=-0.10% になること
e=E_exact(9.00,Rational(1,10),Kmin)*100; print('E at Kmin',N(e)); assert abs(e+0.10)<0.002
e=E_exact(9.00,Rational(1,10),Rational(1,10**10))*100; print('phenol',N(e)); assert abs(e+91)<1

check_structure_and_equations()
print('ALL OK')
