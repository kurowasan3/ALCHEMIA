# 検証スクリプト POL-FIB-012  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/POL-FIB-012/item.json'
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

fiber = {'羊毛': {'NH3+'}, 'ナイロン': {'NH3+'}, '木綿': {'OH'}, 'ポリエステル': set(), 'アクリル': {'SO3-'}}
dye = {'ア': {'SO3-'}, 'イ': {'OH', 'NH2'}, 'ウ': set(), '陽イオン染料': {'NR3+'}}
def bond(f, d):
    F, D = fiber[f], dye[d]
    if any(x.endswith('+') for x in F) and any(x.endswith('-') for x in D): return 'イオン結合'
    if any(x.endswith('-') for x in F) and any(x.endswith('+') for x in D): return 'イオン結合'
    if F & {'OH', 'NH2'} and D & {'OH', 'NH2'}: return '水素結合'
    return '分子間力のみ'
assert bond('羊毛', 'ア') == 'イオン結合'
assert bond('木綿', 'ア') == '分子間力のみ'
assert bond('木綿', 'イ') == '水素結合'
assert bond('ポリエステル', 'ア') == '分子間力のみ' and bond('ポリエステル', 'イ') == '分子間力のみ'
assert bond('アクリル', '陽イオン染料') == 'イオン結合'
dyed = [f for f in ('羊毛', 'ポリエステル') if bond(f, 'ア') == 'イオン結合']
assert dyed == ['羊毛']
print('染料と繊維の対応 OK')

print('ALL OK')
