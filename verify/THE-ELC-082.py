ITEM = '/home/claude/chemx/x1/outh/items/THE-ELC-082/item.json'
# 検証スクリプト（x1 g04）  python3 -I この.py
import sys, json, re, math
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), ('小問数', nq, na)
print('小問数OK', nq)

# ---- 個別の検証 ----

for e_ in ['2H+ + 2e- -> H2','2H2O -> O2 + 4H+ + 4e-']:
    assert check_equation(e_) is None
s=1.36e-3; V=0.500; I=0.500; F=9.65e4; nCu=2.00e-2
tc=2*nCu*F/I/60; chk('tc',tc,129,3)
def xO2(t,nAr):
    e=I*t*60/F; O2=e/4; H2=max(0.0,(e-2*nCu)/2); return O2/(O2+H2+nAr)
for nAr,a,b in [(1e-2,3.40e-4,2.27e-4),(3e-2,1.70e-4,2.27e-4)]:
    chk('exh',s*V*xO2(tc,nAr),a,3); chk('lim',s*V*xO2(1e9,nAr),b,3)
    # 溶解量／気相の量が小さいことの確認
    assert s*V*xO2(tc,nAr)/(nCu/2) < 0.04
assert s*V*xO2(tc+10,1e-2) < s*V*xO2(tc,1e-2)
assert s*V*xO2(tc+10,3e-2) > s*V*xO2(tc,3e-2)
for nAr in [0.005,0.015,0.02,0.025,0.04]:
    d=xO2(tc+50,nAr)-xO2(tc+20,nAr)
    sgn=(nAr>nCu)-(nAr<nCu)
    assert (abs(d)<1e-12 if sgn==0 else d*sgn>0), (nAr,d)
chk('300min',s*V*xO2(300,1e-2),2.64e-4,3)
print('場合分けOK')

print('ALL OK')
