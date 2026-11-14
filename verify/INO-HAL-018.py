# 検証スクリプト（x1 g03）  python3 -I この.py
import sys, json, re
from math import log10, floor, sqrt, pi, exp, log
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x1/out/items/INO-HAL-018/item.json'
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

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)

EA={'F':328,'Cl':349,'Br':325}; BE={'F':159,'Cl':243,'Br':193}; HY={'F':-515,'Cl':-381,'Br':-347}
d1=BE['F']+2*EA['Cl']-2*EA['F']-BE['Cl']; assert d1==-42
d2=BE['Cl']+2*EA['Br']-2*EA['Cl']-BE['Br']; assert d2==2
d3=d2-2*HY['Br']+2*HY['Cl']; assert d3==-66
# 水和の値を ±20 動かしても (3) は負
for a in (-20,0,20):
  for b in (-10,0,10): assert d2-2*(HY['Br']+a)+2*(HY['Cl']+a+b) < 0
print('ALL OK')
