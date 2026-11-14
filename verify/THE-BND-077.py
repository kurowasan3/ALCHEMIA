ITEM = '/home/claude/chemx/tokodai/out/items/THE-BND-077/item.json'
# 検証スクリプト（東工大 g03）  python3 -I この.py
# 共通の検証ヘッダ（各検証スクリプトにそのまま埋め込む）
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
nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')
tab={'Na':('b',23.0,0.971),'Al':('f',27.0,2.70),'Cr':('b',52.0,7.19),'Cu':('f',63.5,8.96)}
given={'Na':(23.7,2.87),'Al':(10.0,2.15),'Cr':(7.23,1.93),'Cu':(7.09,1.92)}
for k,(s,M,D) in tab.items():
    assert abs(round(M/D,3 if M/D<10 else 1)-given[k][0])<0.006,(k,M/D); assert abs(round((M/D)**(1/3),2)-given[k][1])<0.006,(k,(M/D)**(1/3))
kf=4**(1/3)/2**0.5; kb=3**0.5/2*2**(1/3); kfg=1.59/1.41; kbg=1.73/2*1.26
for (kf_,kb_,lab) in ((kf,kb,'exact'),(kfg,kbg,'given')):
    dd={k:(kf_ if s=='f' else kb_)*((M/D)**(1/3) if lab=='exact' else given[k][1]) for k,(s,M,D) in tab.items()}
    order=sorted(dd,key=dd.get); print(lab,dd,order); assert order==['Cr','Cu','Al','Na']
    chk('ratio '+lab,dd['Na']/dd['Cr'],1.49,3)
print('packing',3**0.5*3.14159/8, 2**0.5*3.14159/6)
print('ALL OK')
