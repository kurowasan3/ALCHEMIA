# INO-PPT-053 検証：試薬の順序の数え上げと ALCHEMIA 図鑑データとの照合，反応式のつり合い
import sys, json, itertools
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ions = ['Ag+', 'Cu^2+', 'Al^3+', 'Ba^2+']
d = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))
res = d['results']
code = {'A': 'Cl', 'B': 'SO4', 'C': 'H2S-a', 'D': 'NH3-ex'}
ppt = {}
for k, rid in code.items():
    s = set()
    for i in ions:
        r = res[i].get(rid, {})
        print(k, i, r.get('kind'), r.get('f'), r.get('note', '')[:40])
        if r.get('kind') in ('ppt', 'remain'):  # remain = 少量で生じた沈殿が過剰でも残る
            s.add(i)
    ppt[k] = s
# 過剰のアンモニア水：図鑑では「少量で沈殿→過剰で溶解」の kind を確認
print(ppt)
expected = {'A': {'Ag+'}, 'B': {'Ba^2+'}, 'C': {'Ag+', 'Cu^2+'}, 'D': {'Al^3+'}}
assert ppt == expected, ppt
cnt = 0
for order in itertools.permutations('ABCD'):
    rem = set(ions); ok = True
    for k in order:
        p = ppt[k] & rem
        if len(p) != 1: ok = False; break
        rem -= p
    if ok and not rem:
        cnt += 1
print('valid orders =', cnt)
assert cnt == 12
for e in ['Ag+ + Cl- -> AgCl', 'Ba^2+ + SO4^2- -> BaSO4', '2Ag+ + H2S -> Ag2S + 2H+',
          'Cu^2+ + H2S -> CuS + 2H+', 'Al^3+ + 3NH3 + 3H2O -> Al(OH)3 + 3NH4+',
          'Ag+ + 2NH3 -> [Ag(NH3)2]+', 'Cu^2+ + 4NH3 -> [Cu(NH3)4]^2+']:
    r = check_equation(e); print(e, r); assert r is None
print('OK')
