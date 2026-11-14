ITEM = '/home/claude/chem4/outh/items/THE-EQL-073/item.json'
# 検証スクリプト共通部（グループ g04）  python3 -I この.py
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
# 小問数と答えの番号の一致
nprob = it['problem'].count('\\item')
nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
if nprob:
    assert [int(n) for n in nums] == list(range(1, nprob + 1)), (nprob, nums)
print('小問数OK:', nprob)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
import sympy as sp, math
x=sp.symbols('x',positive=True); a=sp.Rational(3,10)
f=(a-x)*sp.sqrt(x); cr=sp.solve(sp.diff(f,x),x); print(cr); assert cr==[sp.Rational(1,10)]
Ka=2.0e-5; h=math.sqrt(Ka*0.10); pH=-math.log10(h); print('pH',pH, 3-0.5*0.30); assert abs(pH-2.85)<0.006
r=float(f.subs(x,sp.Rational(2,10))/f.subs(x,sp.Rational(1,10))); print(r); assert sig(r,2)==0.71
g=(a-x)*x; assert sp.solve(sp.diff(g,x),x)==[a/2]
hh=(0.010+math.sqrt(0.010**2+4*Ka*0.10))/2; print('h',hh); assert sig(hh,2)==0.010
# 数値積分
def run(h0, x0=1e-3, dt=0.002, T=60.0, k=1000.0):
    xs=x0; vmax=0; xm=None; t=0; vs=[]
    def vel(xx):
        H=(h0+math.sqrt(h0*h0+4*Ka*xx))/2
        return k*(0.30-xx)*H
    while t<T:
        v=vel(xs)
        if v>vmax: vmax=v; xm=xs
        vs.append(v)
        k1=vel(xs); k2=vel(xs+dt/2*k1); k3=vel(xs+dt/2*k2); k4=vel(xs+dt*k3)
        xs+=dt/6*(k1+2*k2+2*k3+k4); t+=dt
    return xm, vs
xm,_=run(0.0); print('weak xm',xm); assert abs(xm-0.10)<0.003
xm2,vs2=run(0.010,x0=0.0); print('HCl xm',xm2); assert xm2<1e-6 and all(vs2[i]>=vs2[i+1] for i in range(len(vs2)-1))
# 電離度
alpha=math.sqrt(Ka/0.10); print('alpha',alpha); assert alpha<0.05
# SVG 規格
svg=open(ITEM.replace('item.json','solution.svg'),encoding='utf-8').read()
assert 'viewBox' in svg and '<script' not in svg and 'foreignObject' not in svg and 'href' not in svg
print('全検証OK:', ITEM)
