# 検証スクリプト ORG-ARO-043  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-ARO-043/item.json'
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
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ORG-ARO-043：pKa（水溶液中の代表値）で強弱と反応の向きを確認 ----
pKa = {'HCl': -7, 'CH3COOH': 4.76, 'H2CO3': 6.35, 'C6H5OH': 9.99, 'H2O': 15.7, 'C2H5OH': 15.9}
order = sorted(pKa, key=lambda k: pKa[k]); print(order)
assert order == ['HCl', 'CH3COOH', 'H2CO3', 'C6H5OH', 'H2O', 'C2H5OH']
base = {'CH3COO-': 'CH3COOH', 'C2H5O-': 'C2H5OH', 'C6H5O-': 'C6H5OH', 'OH-': 'H2O', 'HCO3-': 'H2CO3'}
border = sorted(base, key=lambda b: -pKa[base[b]]); print(border)
assert border == ['C2H5O-', 'OH-', 'C6H5O-', 'HCO3-', 'CH3COO-']
# (2) 反応する酸 HA と生じる酸 HB：pKa(HA) < pKa(HB) なら進む
cases = {'a': ('H2O', 'C2H5OH'), 'b': ('H2O', 'C6H5OH'), 'c': ('C6H5OH', 'H2CO3'), 'd': ('C6H5OH', 'C2H5OH'), 'e': ('C2H5OH', 'H2O'), 'f': ('H2CO3', 'C6H5OH')}
go = {k for k, (ha, hb) in cases.items() if pKa[ha] < pKa[hb]}; print(go); assert go == {'a', 'd', 'f'}
print('注：水・エタノールの pKa は値の取り方で差が小さいが，エトキシドが水で加水分解されること（NaOEt + H2O → EtOH + NaOH）は教科書・文献の標準的な記述')
print('ALL OK')
