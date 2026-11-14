# 検証スクリプト THE-ATM-014  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-ATM-014/item.json'
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
# 核反応式（\text{} 付きで chemlib の対象外）は質量数・原子番号（電荷）の保存を独自に確認する
def nuc(tok):
    tok = tok.strip()
    if tok == 'e-': return (0, -1)
    m = re.fullmatch(r'\^\{(\d+)\}_\{(\d+)\}([A-Za-z]+)', tok)
    assert m, tok
    return (int(m.group(1)), int(m.group(2)))
nuc_eqs = [b for k in ('problem','answer','solution') for b in ce_blocks(it[k]) if '->' in b]
print(nuc_eqs)
for e in nuc_eqs:
    lhs, rhs = re.split(r'->\[\\text\{\}\]', e)
    L = [nuc(t) for t in lhs.split(' + ')]; R = [nuc(t) for t in rhs.split(' + ')]
    assert sum(a for a,_ in L) == sum(a for a,_ in R) and sum(z for _,z in L) == sum(z for _,z in R), e
assert len(nuc_eqs) == 2
# (3) α・β の回数を全探索で一意に決める
sols = [(a,b) for a in range(20) for b in range(20) if 232-4*a == 208 and 90-2*a+b == 82]
assert sols == [(6,4)], sols
# (4) 中性子数
assert (232-90)-(208-82) == 16 == 2*6+1*4
print('崩壊の回数・中性子数 OK')

print('ALL OK')
