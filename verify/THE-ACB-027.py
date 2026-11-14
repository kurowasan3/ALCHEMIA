# 検証スクリプト THE-ACB-027  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/THE-ACB-027/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert len(re.findall(r'\((\d+)\)', it['answer'])) == len(nums), '答え中に余分な (n)'
    if n_items: assert nums == [str(i+1) for i in range(n_items)], (n_items, nums)
    else: assert not nums
    print('小問数', n_items, '答え', nums)
check_structure()
check_all_equations()
def bis(f, lo=-16, hi=1):
    """f(h) が h について増加 → f=0 の h（log 空間の二分法）"""
    for _ in range(200):
        m = (lo + hi) / 2
        if f(10**m) > 0: hi = m
        else: lo = m
    return 10**((lo + hi) / 2)
# ---- 個別の検証 ----
from math import log10
Ka=2.0e-5; Kw=1.0e-14
# (1) 中和点：Na 5.0 mmol, 酸 5.0 mmol, 100 mL
C=0.050
h = bis(lambda h: -(C + h - C*Ka/(h+Ka) - Kw/h))  # Na + H - A - OH = 0 ; decreasing in h? check sign
f = lambda h: C + h - C*Ka/(h+Ka) - Kw/h
print('f increasing?', f(1e-9) < f(1e-8))
h = bis(f)
print('中和点 pH', -log10(h)); assert round(-log10(h),1) == 8.7
# (3) pH 7.00 となる NaOH 体積 v
def pH_at(v):
    Ct = 5.0/(50+v); Na = 0.10*v/(50+v)
    g = lambda h: Na + h - Ct*Ka/(h+Ka) - Kw/h
    return -log10(bis(g))
lo, hi = 40.0, 50.0
for _ in range(100):
    m = (lo+hi)/2
    if pH_at(m) < 7: lo = m
    else: hi = m
print('pH7 の体積', m, '差', 50-m); assert sig(50-m, 2) == 0.25
# (4) 弱塩基 B
Kb=2.0e-5
def pH_B(v):
    Ca = 5.0/(50+v); Cb = 0.10*v/(50+v)
    g = lambda h: Cb*(h/(h+Kw/Kb)) + h - Ca*Ka/(h+Ka) - Kw/h   # [BH+] = Cb*h/(h+Ka') , Ka'=Kw/Kb
    return -log10(bis(g))
lo, hi = 40.0, 60.0
for _ in range(100):
    m = (lo+hi)/2
    if pH_B(m) < 7: lo = m
    else: hi = m
print('B の体積', m); assert abs(m - 50.0) < 0.01

print('ALL OK')
