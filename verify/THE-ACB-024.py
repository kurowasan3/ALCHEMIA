# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10, sqrt
ITEM = '/home/claude/chem2/out/items/THE-ACB-024/item.json'
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
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and ' + ' in b:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

# (1)
H = 2*0.0500*25.0; OH = 0.100*15.0; nNH3 = H - OH      # mmol
N = nNH3*14.0; prot = N/0.160; pct = prot/1000*100
print('(1)', nNH3, N, prot, pct); assert round(pct, 2) == 8.75
# (2) 終点の pH とフェノールフタレインでの過剰滴定
Ka = 5.7e-10
C = nNH3/40.0
print('終点 pH', -log10(sqrt(Ka*C)))
assert 4.2 < -log10(sqrt(Ka*C)) < 6.2
frac = Ka/(Ka + 10**-8.5); print('pH8.5 で NH3 になる割合', frac); assert frac > 0.1
# (3)
nN = 1.00*0.30*0.160/14*1000; nH2SO4 = nN/2; V = nH2SO4/0.0500
print('(3)', nN, nH2SO4, V); assert round(V) == 34
# (4)
Nmel = 10.0*6*14/126; dp = Nmel/0.160/1000*100
print('(4)', Nmel, dp, pct + dp); assert round(dp, 1) == 4.2
print('全検証OK:', ITEM)
