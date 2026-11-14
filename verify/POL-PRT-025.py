# 検証スクリプト POL-PRT-025  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/POL-PRT-025/item.json'
import json, re
it = json.load(open(ITEM, encoding='utf-8'))
n_items = len(re.findall(r'\\item\b', it['problem'])); n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans == 4, (n_items, n_ans)
eqs = [b for k in ('problem', 'answer', 'solution') for b in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[k])
       if len(re.findall(r'->|<=>|<-', b)) == 1]
assert eqs == []
data = {'ア': dict(A=30, G=20, C=20, T=30, U=0), 'イ': dict(A=25, G=24, C=18, T=33, U=0),
        'ウ': dict(A=28, G=22, C=22, T=0, U=28), 'エ': dict(A=30, G=19, C=26, T=0, U=25)}
def judge(d):
    assert sum(d.values()) == 100
    kind = 'DNA' if d['T'] > 0 else 'RNA'
    pair = d['T'] if kind == 'DNA' else d['U']
    return kind, ('二本鎖' if d['A'] == pair and d['G'] == d['C'] else '一本鎖')
res = {k: judge(v) for k, v in data.items()}; print(res)
assert res == {'ア': ('DNA', '二本鎖'), 'イ': ('DNA', '一本鎖'), 'ウ': ('RNA', '二本鎖'), 'エ': ('RNA', '一本鎖')}
s1 = dict(A=35, G=15, C=25, T=25)
s2 = dict(A=s1['T'], T=s1['A'], G=s1['C'], C=s1['G']); print('相補鎖', s2)
assert s2 == dict(A=25, G=25, C=15, T=35)
assert {b: (s1[b] + s2[b]) / 2 for b in 'AGCT'} == {b: data['ア'][b] for b in 'AGCT'}
# (3) の反例：s1 と s2 をつないだ 1 本鎖（長さ 2 倍）の組成
joined = {b: (s1[b] + s2[b]) / 2 for b in 'AGCT'}; assert joined['A'] == joined['T'] and joined['G'] == joined['C']
e = data['エ']; tmpl = dict(T=e['A'], A=e['U'], C=e['G'], G=e['C']); print('鋳型', tmpl)
assert tmpl == dict(A=25, G=26, C=19, T=30) and sum(tmpl.values()) == 100
print('全検証OK:', ITEM)
