# 検証スクリプト THE-PER-016  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-PER-016/item.json'
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
import json
E = {e['sym']: e for e in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/elements.json', encoding='utf-8'))['elements']}
assert E['Na']['ie1'] == 496
# 電子親和力の文献値（CRC Handbook：kJ/mol）。Ne は安定な陰イオンをつくらない（負）
EA = {'O':141.0,'F':328.2,'S':200.4,'Cl':348.6}
assert round(EA['F']) == 328 and round(EA['Cl']) == 349
assert EA['F'] > EA['O'] and EA['Cl'] > EA['S'] and EA['Cl'] > EA['F']
# 殻の電子配置（Ne は閉殻）
assert E['Ne']['shells'] == [2,8] and E['F']['shells'] == [2,7]
assert E['F']['en'] > E['Cl']['en']
dE = 496 - 349; assert dE == 147
# NaCl の格子エネルギー（文献値 787 kJ/mol）が 700 を超えることの確認
assert 787 > 700 and 787 > dE
print('OK')

print('ALL OK')
