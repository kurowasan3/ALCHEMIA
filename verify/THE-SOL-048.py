import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_all_equations(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts(it):
    np_ = it['problem'].count('\\item')
    na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), (np_, na)
    print('小問数OK:', np_)
import glob
it = json.load(open('/home/claude/chem3/out/items/THE-SOL-048/item.json', encoding='utf-8'))
check_all_equations(it); check_parts(it)
hits = [f for f in glob.glob('/home/claude/kurowasan3/alchemia/content/chem/*.json') if 'コロイド' in open(f, encoding='utf-8').read()]
print('図鑑中のコロイド項目:', hits)
# 判定表（教科書記述に基づく）：疎水=少量で凝析, 親水=多量で塩析, 保護=少量で凝析しない
table = {'X': ('生じる', '生じる'), 'Y': ('生じない', '生じる'), 'Z': ('生じない', '生じる')}
for k, (a, b) in table.items():
    assert f'{k}：(a) {a}，(b) {b}' in it['answer'], k
print('ALL OK')
