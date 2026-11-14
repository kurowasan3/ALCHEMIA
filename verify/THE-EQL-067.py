ITEM = '/home/claude/chem4/out/items/THE-EQL-067/item.json'
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
dfH={'NO2':33.2,'N2O4':9.2,'HI':26.5,'I2g':62.4,'H2':0}
dH1=dfH['N2O4']-2*dfH['NO2']; dH2=2*dfH['HI']-dfH['I2g']-dfH['H2']
print('2NO2->N2O4',dH1,' H2+I2(g)->2HI',dH2); assert dH1<0 and dH2<0
assert dfH['N2O4']>0 and dfH['HI']>0   # 生成エンタルピーは実際には正 → 実験から負とは言えない
# Δn
assert (1-2)!=0 and (2-2)==0
# 実験Ⅰ：加熱で解離 → n 増 → P/T 増
n_before=1.0; n_after=1.2; assert n_after/n_before>1
ans=re.findall(r'\((\d)\)\s*([①②③])', it['answer']); print(ans)
assert [a for _,a in ans]==['①','③','①','②','①','①','③','②']
import json as _j
r=_j.load(open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json',encoding='utf-8'))
s=_j.dumps(r,ensure_ascii=False); assert '発熱反応なので冷やすと色が薄くなる' in s
print('全検証OK:', ITEM)
