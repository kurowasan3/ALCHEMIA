ITEM = '/home/claude/chem6/out/items/INO-HAL-016/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if should_check(b):
            r = check_equation(clean_ce(b))
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
print('検査した式の数:', n_eq)

# 小問数と答えの数の一致
n_items = it['problem'].count('\\item')
n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))


# 溶解性の表（教科書の記述）から観察結果と矛盾のない割り当てが一意であることを確認
import itertools
ag={'NaF':None,'NaCl':'白','NaBr':'淡黄','NaI':'黄'}
ca={'NaF':'白','NaCl':None,'NaBr':None,'NaI':None}
obsI={'A':'白','B':'黄','C':'淡黄','D':None}; obsII={'A':None,'B':None,'C':None,'D':'白'}
sols=[]
for p in itertools.permutations(ag):
    m=dict(zip('ABCD',p))
    if all(ag[m[k]]==obsI[k] and ca[m[k]]==obsII[k] for k in 'ABCD'): sols.append(m)
print(sols); assert sols==[{'A':'NaCl','B':'NaI','C':'NaBr','D':'NaF'}]
# 操作Ⅰだけでも一意か（Ⅱは確認）
s1=[p for p in itertools.permutations(ag) if all(ag[x]==obsI[k] for k,x in zip('ABCD',p))]
print('I only',len(s1))

print('全検証OK:', ITEM)
