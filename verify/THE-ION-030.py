NAME = 'THE-ION-030'
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

Ksp=Rational(1,10**10); K2=Rational(1,100)
s0=sqrt(Ksp); assert s0==Rational(1,10**5)
K=Ksp/K2; assert K==Rational(1,10**8)
s1=sqrt(Ksp*(1+1/K2)); print('s 1M',N(s1)); assert sig(s1,2)==1.0e-4
# [H+]一定の仮定の妥当性: 厳密に s を解く
s=symbols('s',positive=True); H=1-s*K2/(K2+1)  # 近似: HSO4- 生成で H+ がわずかに減る
ss=nsolve(Ksp-s*s/(1+H/K2), s, 1e-4); print('s exact',ss); assert abs(ss/s1-1)<0.001
loss=Ksp/Rational(1,100)*(1+Rational(1,10)/K2)*Rational(2,10)*233*1000; print('loss mg',N(loss)); assert sig(loss,2)==0.0051
s3=sqrt(Ksp*(1+Rational(1,10)/K2)); loss2=s3*Rational(2,10)*233*1000; print('loss eq mg',N(loss2), N(loss2/loss)); assert sig(loss2,2)==1.5
KB=Rational(5,10**9)/(Rational(4,10**7)*Rational(5,10**11)); print(KB); assert KB==250000000

check_structure_and_equations()
print('ALL OK')
