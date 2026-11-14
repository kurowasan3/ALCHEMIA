# 検証スクリプト ORG-DET-025  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-DET-025/item.json'
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
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b) and (' + ' in b or re.search(r'(^|\s)\d+\s*[A-Z(\[]', b)):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert len(nums) == n_items or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items, '答え数', len(nums))
check_structure()
check_all_equations()

# ---- 個別の検証 ----

# 官能基 → 塩をつくる試薬
groups = {'A': {'amine'}, 'B': {'amine', 'cooh'}, 'C': {'phenol', 'ester'}, 'D': {'amide'}, 'E': {'phenol'}}
reacts = {'NaHCO3': {'cooh'}, 'NaOH': {'cooh', 'phenol'}, 'HCl': {'amine'}}
def run(order, comps):
    ether = set(comps); out = {}
    for lab, reag in order:
        moved = {c for c in ether if groups[c] & reacts[reag]}
        for c in moved: out[c] = lab
        ether -= moved
    for c in ether: out[c] = 'IV'
    return out
r = run([('I', 'NaHCO3'), ('II', 'NaOH'), ('III', 'HCl')], 'ABCDE')
assert r == {'A': 'III', 'B': 'I', 'C': 'II', 'D': 'IV', 'E': 'II'}, r
# (2) 希塩酸を先に
aq = {c for c in 'ABCDE' if groups[c] & reacts['HCl']}
assert aq == {'A', 'B'}
# 塩基性にすると酸性基（cooh, phenol）をもつものはイオンのまま水層に残る
to_ether = {c for c in aq if not (groups[c] & reacts['NaOH'])}
assert to_ether == {'A'}
print('層の割り当て OK')

print('ALL OK')
