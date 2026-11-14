ITEM = '/home/claude/chem4/out/items/INO-CSI-011/item.json'
# 検証スクリプト（グループ g07）  python3 -I この.py
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

# 小問数と答えの数の一致
n_items = it['problem'].count('\\item')
n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
print('小問数', n_items, '答え', n_ans)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

# 容器の割り当ての一意性（総当たり）
from itertools import permutations
cont = {'A': dict(glass=True, gstop=True, brown=False, rubber=False),
        'B': dict(glass=True, gstop=False, brown=False, rubber=True),
        'C': dict(glass=True, gstop=True, brown=True, rubber=False),
        'D': dict(glass=False, gstop=False, brown=False, rubber=False)}
need = {
 'ア HF':   lambda c: not c['glass'],                       # ガラスを溶かす
 'イ NaOH': lambda c: not c['gstop'],                       # ガラス栓が固着
 'ウ 濃HNO3': lambda c: c['brown'] and not c['rubber'],      # 光分解・ゴムを酸化
 'エ 濃HCl': lambda c: True,
}
names = list(need)
sols = []
for p in permutations('ABCD'):
    if all(need[n](cont[k]) for n, k in zip(names, p)):
        sols.append(dict(zip(names, p)))
print(sols)
assert len(sols) == 1 and sols[0] == {'ア HF': 'D', 'イ NaOH': 'B', 'ウ 濃HNO3': 'C', 'エ 濃HCl': 'A'}
print('すべて OK')
