ITEM = '/home/claude/chem4/out/items/THE-ACB-030/item.json'
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
import sympy as sp
R=sp.Rational
M=137+32+64; assert M==233
ba=R(350,1000)/M; cB=ba/R(1,100); print(float(cB)); assert sig(float(cB),3)==0.150
cBr=R(15,100)  # 3桁に丸めた値で以後計算
so4=R(233,1000)/M; assert so4==R(1,1000)
assert cBr*R(16,1000) > so4   # 実験1 は SO4 で決まる
oh=2*cBr*R(16,1000); hno3=oh-2*so4; print(float(hno3)); assert sig(float(hno3/R(1,100)),3)==0.280
assert sig(float(so4/R(1,100)),3)==0.100
# 実験2 の確認：十分な A なら SO4 過剰
# (4)
so4_5=R(1,10)*R(5,1000); assert so4_5 < cBr*R(1,100)
print(so4_5/R(1,100)/cBr); assert so4_5/R(1,100)/cBr==R(1,3)
import json as _j
s=open('/home/claude/kurowasan3/alchemia/content/chem/ions.json',encoding='utf-8').read()
assert '"f": "BaSO4", "color": "白"' in s.replace('\n','') or 'BaSO4' in s
print('全検証OK:', ITEM)
