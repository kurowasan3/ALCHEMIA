# INO-PPT-042 検証：キレート錯体の組成と電荷
# --- 共通：item.json の \ce{} 反応式のつり合い検査（単独で動く） ---
import json, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
def _ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k, d = j + 4, 1
        while d:
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k
def check_item_equations(path):
    it = json.load(open(path, encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for ce in _ce_blocks(it[key]):
            arrows = re.findall(r'<=>|->|<-', ce)
            if len(arrows) != 1 or '\\boxed' in ce or '\\text' in ce or '\\mathrm' in ce:
                continue
            r = check_equation(ce)
            assert r is None, (key, ce, r)
            n += 1
    print('反応式', n, '本 OK')
    return it
check_item_equations('/home/claude/chemx/todai/out/items/INO-PPT-042/item.json')
metal = {'Al': (3, 6), 'Fe': (3, 6), 'Cu': (2, 4)}
lig = {'Q': -1, 'en': 0, 'gly': -1, 'ox': -2}
def cx(m, l):
    q, cn = metal[m]; n = cn // 2
    return n, q + n*lig[l]
assert cx('Al', 'Q') == (3, 0) and cx('Cu', 'Q') == (2, 0)
res = {k: cx(*v) for k, v in {'a': ('Fe', 'Q'), 'b': ('Cu', 'en'), 'c': ('Cu', 'gly'), 'd': ('Fe', 'ox')}.items()}
print(res)
assert [k for k, v in res.items() if v[1] == 0] == ['a', 'c']
assert res['b'][1] == 2 and res['d'][1] == -3
print('OK')
