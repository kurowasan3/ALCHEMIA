# 検証スクリプト（x2 g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x2/outh/items/THE-EQL-145/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 1), (nprob, nans)
print('小問数', nprob, '答え', nans)

from math import log10
R=8.31e-3; E1=180.0; E2=60.0
def lk1(T): return -6 - E1/(2.3*R)*(1/T-1/500)
def lk2(T): return -3 - E2/(2.3*R)*(1/T-1/500)
def lk(T): return log10(10**lk1(T)+10**lk2(T))
c=(E1-E2)/(2.3*R)
Tx=1/(1/500-3/c); print('交差温度', Tx); assert 650<Tx<662 and '6.6\\times10^{2}' in it['answer']
r=10**(lk1(800)-lk2(800)); f=1/(1+r); print('800K 触媒経路', f*100); assert abs(f*100-1.9)<0.05 and '1.9\\,\\%' in it['answer']
assert abs(10**0.709-5.12)<0.01
Ea=2.3*R*(lk(680)-lk(640))/(1/640-1/680); print('見かけEa', Ea); assert abs(Ea-120)<abs(Ea-60) and abs(Ea-120)<abs(Ea-180)
# 下に凸：傾き（-dlogk/d(1/T)）は高温ほど大きい
def slope(T,h=1e-6):
    x=1/T; return (lk(1/(x+h))-lk(1/(x-h)))/(2*h)
s_hi=slope(900); s_lo=slope(450); print('傾き', s_hi, s_lo); assert s_hi < s_lo < 0
Tp=1/(1/500-2/c); print('被毒後', Tp); assert 590<Tp<600 and '5.9\\times10^{2}' in it['answer']

print('ALL OK')
