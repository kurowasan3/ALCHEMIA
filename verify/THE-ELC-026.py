NAME = 'THE-ELC-026'
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

k=Rational(59,1000)
AgL=Rational(1,10)/10**(Rational(472,1000)/k); print(N(AgL)); assert abs(AgL-Rational(1,10**9))<1e-12*1e-9*1e3
Ksp=AgL*Rational(1,10); print('Ksp',N(Ksp)); assert sig(Ksp,2)==1.0e-10
Ag3=Rational(1,10)/10**(Rational(531,1000)/k); beta=Rational(1,1000)/Ag3; print('beta',N(beta)); assert sig(beta,2)==1.0e7
E4=k*log(Rational(1,10)/(Rational(5,10**13)/Rational(1,10)),10); print('E AgBr',N(E4)); assert round(float(E4),2)==0.61

check_structure_and_equations()
print('ALL OK')
