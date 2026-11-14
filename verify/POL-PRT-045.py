ITEM = '/home/claude/chem5/outh/items/POL-PRT-045/item.json'
# 検証スクリプト（グループ g09）  python3 -I この.py
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
assert [int(m.group(2)) for m in re.finditer(r'(^|\s)\((\d+)\)', it['answer'])] == list(range(1, np_+1))

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from itertools import product
from collections import Counter
def hydrolyze(seq, ribo):
    # i 番目の糖がリボースなら i と i+1 の間が切れる（最後の残基の後ろには結合がない）
    frags, cur, water = [], [], 0
    for i,(b,r) in enumerate(zip(seq,ribo)):
        cur.append(b)
        if i < len(seq)-1 and r:
            frags.append((tuple(sorted(cur)), True)); cur = []; water += 1
    frags.append((tuple(sorted(cur)), False))
    return Counter(frags), water
# (2) 全リボースの ACGU
prod, wat = hydrolyze('ACGU', [1,1,1,1])
print(prod, wat)
assert prod == Counter({(('A',),True):1, (('C',),True):1, (('G',),True):1, (('U',),False):1}) and wat == 3
# (3)(4) キメラ
obs = Counter({(('G',),True):1, (('A','C'),True):1, (('A','C','G'),False):1})
sols = []
for ribo in product([0,1], repeat=6):
    p,w = hydrolyze('GCAGCA', ribo)
    if p == obs: sols.append((ribo,w))
print(sols)
assert len(sols)==2 and all(s[0][:5]==(1,0,1,0,0) for s in sols) and {s[0][5] for s in sols}=={0,1}
assert all(w==2 for _,w in sols)
assert 5-2 == 3
print('ALL OK')
