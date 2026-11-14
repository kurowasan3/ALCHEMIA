# 検証スクリプト（グループ g08）  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/THE-RDX-043/item.json'
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
print('反応式のつり合い OK：', n_eq, '式')
# 小問数と答えの数の一致
m_en = re.search(r'\\begin\{enumerate\}(.*)\\end\{enumerate\}', it['problem'], re.S)
nq = m_en.group(1).count('\\item') if m_en else 0
assert '\\begin{itemize}' not in it['problem']
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), (nq, na)
print('小問数', nq, '答えの数', na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))


E = {'Fe3+':0.77,'Cu2+':0.34,'Ni2+':-0.25,'Fe2+':-0.44,'Zn2+':-0.76}
order = sorted(E, key=lambda k: -E[k]); assert order == ['Fe3+','Cu2+','Ni2+','Fe2+','Zn2+']
def run(metal_E, n_metal):
    # 電子を還元されやすい順に配分（金属のイオンより還元されやすいものだけ）
    e = 2*n_metal; Fe3, Cu2, Fe2 = 0.05, 0.10, 0.0; Cu = Fe = 0.0
    take = min(e, Fe3); Fe3 -= take; Fe2 += take; e -= take
    if E['Cu2+'] > metal_E:
        take = min(e, 2*Cu2); Cu += take/2; Cu2 -= take/2; e -= take
    if E['Fe2+'] > metal_E:
        take = min(e, 2*Fe2); Fe += take/2; Fe2 -= take/2; e -= take
    return Cu, Fe, e
cu, fe, rest = run(E['Ni2+'], 0.100); print(cu); assert abs(cu-0.075) < 1e-12 and sig(cu*63.5,3) == 4.76
assert sig(0.125*58.7,3) == 7.34
cu, fe, rest = run(E['Ni2+'], 0.200); assert fe == 0
for x in (0.125, 0.13, 0.175, 0.2):
    cu, fe, rest = run(E['Zn2+'], x); print(x, cu, fe)
assert run(E['Zn2+'],0.125)[1] == 0 and run(E['Zn2+'],0.13)[1] > 0 and abs(run(E['Zn2+'],0.175)[1]-0.05) < 1e-12

print('ALL OK')
