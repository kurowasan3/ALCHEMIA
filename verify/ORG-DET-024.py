# 検証スクリプト ORG-DET-024  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-DET-024/item.json'
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
# ---- 個別の検証 ORG-DET-024 ----
for e in ['CH3CHO + 3OH- -> CH3COO- + 2H2O + 2e-', '2Cu^2+ + 2OH- + 2e- -> Cu2O + H2O', '[Ag(NH3)2]+ + e- -> Ag + 2NH3',
          'CH3CHO + 2Cu^2+ + 5OH- -> CH3COO- + Cu2O + 3H2O', 'HCHO + 4[Ag(NH3)2]+ + 6OH- -> CO3^2- + 4Ag + 8NH3 + 4H2O',
          'HCHO + 3OH- -> HCOO- + 2H2O + 2e-', 'HCOO- + 3OH- -> CO3^2- + 2H2O + 2e-']:
    assert check_equation(e) is None, (e, check_equation(e))
M_Cu2O = 63.5*2 + 16; assert M_Cu2O == 143
n = 1.43/143; m_ald = n*44; assert sig(m_ald, 3) == 0.440
m_ag = 2*n*108; assert sig(m_ag, 3) == 2.16
nf = 0.300/30; m_ag2 = 4*nf*108; assert sig(m_ag2, 3) == 4.32
print(m_ald, m_ag, m_ag2)
print('ALL OK')
