# 検証スクリプト THE-SOL-056  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/THE-SOL-056/item.json'
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
    FIELDS = '原子の構造と電子配置／元素の周期表と周期性／化学結合と結晶／物質量と化学反応式／酸と塩基／酸化還元反応／電池と電気分解／気体の性質／溶液の性質／化学反応とエネルギー／反応速度と化学平衡／電離平衡と溶解平衡／ハロゲン・貴ガス・水素／酸素・硫黄／窒素・リン／炭素・ケイ素／アルカリ金属・アルカリ土類金属／アルミニウム・亜鉛・スズ・鉛／遷移元素（鉄・銅・銀・クロム・マンガン）／金属イオンの沈殿・錯イオン・系統分析／気体の製法と性質／有機化合物の構造と元素分析／脂肪族炭化水素／アルコール・アルデヒド・ケトン／カルボン酸・エステル・油脂／芳香族化合物／有機化合物の検出反応と分離／高分子の基礎と合成樹脂／合成繊維／ゴム・機能性高分子／糖類／アミノ酸・タンパク質・核酸'.split('／')
    assert it['field'] in FIELDS, it['field']
    assert all(f in ('qa','equation','compare','condition','pathway','calc','structure') for f in it['format'].split(','))
    import os
    for k in ('problem_fig','solution_fig'):
        if it[k]: assert os.path.exists(os.path.join(os.path.dirname(ITEM), it[k])), it[k]
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
import os, re
# 臨界ミセル濃度（文献値，mmol/L, 25 ℃）：アルキル硫酸ナトリウム
cmc = {10: 33.0, 12: 8.2, 14: 2.1}
assert cmc[10] > cmc[12] > cmc[14]
r1, r2 = cmc[12]/cmc[10], cmc[14]/cmc[12]
print('炭素数 +2 ごとの比', round(r1, 3), round(r2, 3))
assert 0.2 < r1 < 0.3 and 0.2 < r2 < 0.3
cstar = 8e-3
P, Q = 1e-3, 2e-2
assert P < cstar < Q
assert abs(cstar*1000 - cmc[12]) < 0.5
# 図の確認：SVG の規格（viewBox・禁止要素なし・currentColor）
svg = open(os.path.join(os.path.dirname(ITEM), 'problem.svg'), encoding='utf-8').read()
assert 'viewBox' in svg and 'script' not in svg and 'foreignObject' not in svg and not re.search(r'\son\w+=', svg)
assert 'currentColor' in svg
print('図 OK')

print('ALL OK')
