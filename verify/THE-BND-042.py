# 検証スクリプト THE-BND-042  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-BND-042/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
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
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
import itertools
from fractions import Fraction as Fr
def hc(seq):
    n = len(seq); return ''.join('h' if seq[(i-1)%n] == seq[(i+1)%n] else 'c' for i in range(n))
def minper(s):
    n = len(s)
    for p in range(1, n+1):
        if n % p == 0 and all(s[i] == s[(i+p)%n] for i in range(n)): return p
def canon(seq):
    n = len(seq); forms = []
    for perm in itertools.permutations('ABC'):
        m = dict(zip('ABC', perm)); s = ''.join(m[x] for x in seq)
        for t in (s, s[::-1]):
            for k in range(n): forms.append(t[k:]+t[:k])
    return min(forms)
def frac_h(s): return Fr(hc(s).count('h'), len(s))
# (1)
assert frac_h('AB') == 1 and frac_h('ABC') == 0 and frac_h('ABAC') == Fr(1,2)
# (2) 符号の変わり目 = h
step = {('A','B'):1,('B','C'):1,('C','A'):1,('A','C'):-1,('C','B'):-1,('B','A'):-1}
for N in range(2, 9):
    for seq in itertools.product('ABC', repeat=N):
        if any(seq[i] == seq[(i+1)%N] for i in range(N)): continue
        st = [step[(seq[i], seq[(i+1)%N])] for i in range(N)]
        assert sum(st) % 3 == 0                    # (3)
        for i in range(N):
            assert (hc(seq)[i] == 'h') == (st[i-1] != st[i])
# (4)(5) 6 層周期の全列挙
classes = {}
for seq in itertools.product('ABC', repeat=6):
    s = ''.join(seq)
    if any(s[i] == s[(i+1)%6] for i in range(6)) or minper(s) != 6: continue
    classes.setdefault(canon(s), set()).add(s)
print({k: frac_h(k) for k in classes})
assert len(classes) == 2
fr = sorted(frac_h(k) for k in classes); assert fr == [Fr(1,3), Fr(2,3)]
# A,B で始まり h の割合 1/3 の 6 層周期の列（数え始めのずれを除いて 1 通り）
cands = sorted(s for k in classes for s in classes[k] if s[:4] == 'ABCA' and frac_h(s) == Fr(1,3))
print(cands)
assert cands == ['ABCACB']
assert hc('ABCACB') == 'hcchcc'
# (5) の例 ABCBCB が ABABAC と同型
assert canon('ABCBCB') == canon('ABABAC') and frac_h('ABCBCB') == Fr(2,3)
# 図の確認（SVG：viewBox あり，script なし）
svg = open('/home/claude/chem2/outh/items/THE-BND-042/problem.svg', encoding='utf-8').read()
assert 'viewBox' in svg and '<script' not in svg and 'foreignObject' not in svg

print('ALL OK')
