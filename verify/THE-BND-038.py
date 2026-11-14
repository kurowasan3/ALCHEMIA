# 検証スクリプト THE-BND-038  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-BND-038/item.json'
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
# 価電子数・非共有電子対数・電子対の方向数から形を判定（VSEPR）
val = {'H':1,'N':5,'O':6,'C':4}
def species(center, ligands, charge, bonds, fc_center):
    """bonds: 各配位原子への結合次数，fc_center: 中心原子の形式電荷。中心原子の非共有電子数 = 価電子 - 結合に使う電子 - 形式電荷"""
    total = val[center] + sum(val[l] for l in ligands) - charge
    lone_e = val[center] - sum(bonds) - fc_center
    # 配位原子のオクテット：各 O は 8 - 2*結合次数 個の非共有電子（H は 0）→ 全電子数の整合
    lig_lone = sum(0 if l == 'H' else 8 - 2*b for l, b in zip(ligands, bonds))
    assert total == 2*sum(bonds) + lone_e + lig_lone, (center, ligands, charge)
    return total, lone_e
def shape(n_lig, lone_pairs):
    dirs = n_lig + lone_pairs
    return {(4,0):'正四面体形',(3,1):'三角錐形',(2,2):'折れ線形',(2,0):'直線形',(2,1):'折れ線形'}[(n_lig,lone_pairs)]
r = {}
r['NH4+'] = species('N',['H']*4, +1, [1,1,1,1], +1)
r['NH3']  = species('N',['H']*3, 0, [1,1,1], 0)
r['NH2-'] = species('N',['H']*2, -1, [1,1], -1)
r['NO2+'] = species('N',['O']*2, +1, [2,2], +1)
r['NO2-'] = species('N',['O']*2, -1, [2,1], 0)
r['NO2']  = species('N',['O']*2, 0, [2,1], +1)
r['CO2']  = species('C',['O']*2, 0, [2,2], 0)
print(r)
assert r['NH4+'] == (8,0) and r['NH3'] == (8,2) and r['NH2-'] == (8,4)
assert r['NO2+'] == (16,0) and r['NO2-'] == (18,2) and r['NO2'] == (17,1) and r['CO2'] == (16,0)
assert shape(4,0) == '正四面体形' and shape(3,1) == '三角錐形' and shape(2,2) == '折れ線形'
assert shape(2,0) == '直線形' and shape(2,1) == '折れ線形'
# 各 O 原子のオクテット確認（NO2+：O=N=O で O は二重結合＋非共有2組＝8，NO2-：O=N-O^- で単結合側 O は 3 組）
# 実測の結合角（文献値）：NH4+ 109.5, NH3 106.7, NH2- 約104, NO2+ 180, NO2 134.1, NO2- 115.4, O3 116.8
ang = {'NH4+':109.5,'NH3':106.7,'NH2-':104.4,'NO2+':180,'NO2':134.1,'NO2-':115.4}
assert ang['NH4+'] > ang['NH3'] > ang['NH2-'] and ang['NO2+'] > ang['NO2'] > ang['NO2-']
# O3 の価電子 18 = NO2-
assert 6*3 == r['NO2-'][0]

print('ALL OK')
