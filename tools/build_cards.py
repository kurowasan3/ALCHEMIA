"""暗記カード（data/cards.json）を作る

カードは2か所から作る：
  1. 図鑑のデータ（content/chem/）・法則（laws.json）・誤概念（misconceptions.json）から自動で作るデッキ
     → 図鑑を直すとカードも同じ内容に直る（図鑑と食い違わない）
  2. content/cards/*.json に手で書いたデッキ
     {"schema": "alchemia.cards/1", "id": "constants", "name": "…", "field": "THE", "desc": "…",
      "cards": [{"id": "avogadro", "f": "表（$...$ や \\ce が使える）", "b": "裏", "n": "補足（任意）"}]}

カードの ID は「デッキID/カードID」。一度使った ID は変えない（覚えた記録がカード ID に結びつく）。
自動で作るカードの ID は、図鑑のデータの ID（イオン・試薬・反応・化合物…）からできているので、図鑑の ID を変えなければ変わらない。
"""
import json
import re
from pathlib import Path

from content_lib import load_chem, load_knowledge, read

KIND_TXT = {'remain': '沈殿は溶けない', 'none': '沈殿しない（変化なし）'}


def ce(f):
    return f'$\\ce{{{f}}}$'


def card(cid, f, b, n=None, c=None, link=None):
    d = {'id': cid, 'f': f, 'b': b}
    if n:
        d['n'] = n
    if c:
        d['c'] = c
    if link:
        d['link'] = link
    return d


def lhs(eq):
    m = re.search(r'\s*(<=>|->|<-)', eq)
    return eq[:m.start()] if m else eq


def build_cards(content):
    C = load_chem(content)
    K = load_knowledge(content)
    decks = []

    def deck(did, name, field, desc, cards):
        decks.append({'id': did, 'name': name, 'field': field, 'desc': desc, 'cards': [{**c, 'id': f'{did}/{c["id"]}'} for c in cards]})

    # ---- 無機：沈殿・錯イオン・呈色
    I = C['ions']
    rg = {r['id']: r for r in I['reagents']}
    cs = []
    for cat in I['cations']:
        for rid, v in I['results'][cat['id']].items():
            if v['kind'] in ('na',) or (v['kind'] == 'none' and not v.get('note')):
                continue
            if v['kind'] == 'none' and rid in ('NaOH', 'NH3'):
                continue
            front = f'{ce(cat["id"])}（{cat["name"]}）に<br>{rg[rid]["name"]}'
            if v['kind'] == 'ppt':
                back = f'<b>{v["color"]}色の沈殿</b>' + (f'　{ce(v["f"])}' if v.get('f') else '') + (f'<br>{v["name"]}' if v.get('name') else '')
            elif v['kind'] == 'dissolve':
                back = f'<b>沈殿が溶ける</b>　{ce(v["f"])}（{v["color"]}）' + (f'<br>{v["name"]}' if v.get('name') else '')
            elif v['kind'] == 'solution':
                back = f'<b>{v["color"]}色の溶液</b>' + (f'　{ce(v["f"])}' if v.get('f') else '') + '（沈殿しない）'
            else:
                back = f'<b>{KIND_TXT[v["kind"]]}</b>'
            if v.get('eq'):
                back += f'<br>{ce(v["eq"])}'
            color = I['colors'].get(v.get('color')) if v['kind'] in ('ppt', 'solution', 'dissolve') else None
            cs.append(card(f'{cat["id"]}|{rid}', front, back, v.get('note'), color, f'#zukan=ions/{cat["id"]}'))
    deck('ions', '沈殿・錯イオン・呈色', 'INO', '金属イオンに試薬を加えたときに起こること（図鑑の反応表から）', cs)
    deck('flame', '炎色反応', 'INO', '炎色反応の色', [card(c['id'], f'{ce(c["id"].replace("^2+", "").replace("+", ""))} の炎色反応', f'<b>{c["flame"]}</b>', c=I['flameColors'][c['flame']])
                                                for c in I['cations'] if c.get('flame')])
    deck('systematic', '系統分析', 'INO', '各段階で加える試薬と沈殿するもの', [
        card(f'g{st["step"]}', f'系統分析の{st["group"]}：加える試薬と沈殿',
             f'{st["reagent"]}<br>' + ('沈殿：' + '，'.join(ce(p) for p in st['ppt']) if st['ppt'] else 'ろ液に残る：' + '，'.join(ce(i) for i in st['ions'])), st['note'], link='#zukan=ions')
        for st in I['systematic']])

    # ---- 無機：気体・製法・反応
    R = C['reactions']
    G = R['gases']
    gas = [r for r in R['reactions'] if r['cat'] == 'gas']
    deck('gas-prep', '気体の製法', 'INO', '気体の製法の反応式と条件', [
        card(r['id'], f'{r["title"]}<br><small>{r["cond"]}</small>', ce(r['eq']), r.get('note'), link=f'#zukan=reactions/{r["id"]}') for r in gas])
    deck('gas-props', '気体の性質', 'INO', '色・におい・水溶液の性質・捕集法・乾燥剤', [
        card(f, f'{ce(f)}（{g["name"]}）の性質と捕集法',
             f'{g["color"]}・{g["smell"]}・{g["nature"]}，{g["sol"]}<br>捕集：{g["collect"] or "—"}　乾燥剤：{"・".join(g["dry"]) or "—"}', g.get('test') or None)
        for f, g in G.items()])
    deck('industry', '工業的製法', 'INO', '工業的製法の条件（触媒）と反応', [
        card(r['id'], r['title'], (f'{r["cond"]}<br>' if r['cond'] else '') + '<br>'.join(ce(s['eq']) for s in r['steps']) + (f'<br>全体：{ce(r["total"])}' if r.get('total') else ''),
             link=f'#zukan=reactions/{r["id"]}')
        for r in R['reactions'] if r['cat'] == 'industry'])
    deck('inorg-rx', '無機の反応式', 'INO', '元素ごとの反応：左辺から右辺を思い出す', [
        card(r['id'], f'{r["title"]}<br>{ce(lhs(r["eq"]) + " ->")} ？', ce(r['eq']), r.get('note'), link=f'#zukan=reactions/{r["id"]}')
        for r in R['reactions'] if r['cat'] == 'reaction'])

    # ---- 有機
    O = C['organic']
    comp = {c['id']: c for c in O['compounds']}
    nm = lambda i: comp[i]['name'] + (f'　{ce(comp[i]["cf"])}' if comp[i].get('cf') else '')  # noqa: E731
    deck('org-rx', '有機の反応', 'ORG', '反応物と試薬から生成物を思い出す（図鑑の反応経路から）', [
        card(e['id'], f'{comp[e["from"]]["name"]} → ？<br><small>{e["reagent"]}（{e["type"]}）</small>', nm(e['to']) + (f'<br>{ce(e["eq"])}' if e.get('eq') else ''),
             e.get('note'), link=f'#zukan=organic/{e["to"]}')
        for e in O['reactions']])
    deck('org-names', '有機化合物の示性式', 'ORG', '化合物名から示性式（分子式）', [
        card(c['id'], c['name'], ce(c['cf']), c.get('note'), link=f'#zukan=organic/{c["id"]}')
        for c in O['compounds'] if c.get('cf') and c['cls'] not in ('無機', '高分子')])
    deck('tests', '検出反応', 'ORG', '操作・観察と、反応する構造', [
        card(t['id'], f'{t["name"]}<br><small>{t["reagent"]}</small>', f'<b>{t["obs"]}</b><br>反応する構造：{t["target"]}<br>陽性の例：' + '，'.join(comp[k]['name'] for k in t['pos']),
             link=f'#zukan=tests/{t["id"]}')
        for t in C['tests']['tests']])

    # ---- 高分子
    deck('polymers', '高分子の単量体と重合', 'POL', '単量体・重合の種類・結合', [
        card(p['id'], f'{p["name"]}' + (f'（{p["abbr"]}）' if p.get('abbr') else '') + 'の単量体と重合',
             '単量体：' + '＋'.join(m['name'] for m in p['monomers']) + f'<br>{p["poly"]}' + (f'　{p["bond"]}' if p.get('bond') else '') + (f'　{p["heat"]}' if p.get('heat') else ''),
             p.get('props'), link=f'#zukan=polymers/{p["id"]}')
        for p in C['polymers']['polymers']])

    # ---- 理論：元素・法則・誤概念
    E = C['elements']['elements']
    deck('elements', '元素記号と電子配置', 'THE', '原子番号 1〜36 の元素名・電子配置', [
        card(e['sym'], f'{ce(e["sym"])}（原子番号 {e["z"]}）', f'<b>{e["ja"]}</b>　{e["period"]}周期 {e["group"]}族<br>' + ' '.join(f'{"KLMN"[i]}{n}' for i, n in enumerate(e['shells'])),
             link=f'#zukan=elements/{e["sym"]}')
        for e in E if e['z'] <= 36])
    ORIGIN = {'fundamental': '出発点となる法則', 'empirical': '実験から得た法則', 'derived': '導出できる法則', 'definition': '定義'}
    deck('laws', '定義・法則', 'THE', '名前から式と成り立つ条件を思い出す（タップで導出・根拠へ）', [
        card(l['id'], f'{l["name"]}<br><small>{ORIGIN.get(l["origin"], "")}</small>', f'${l["tex"]}$<br><small>{l["conditions"]}</small>', link=f'#law={l["id"]}')
        for l in K['laws']])
    deck('misc', 'よくある誤解（○×）', 'THE', '正しいか、まちがいか', [
        card(m['id'], f'○か×か<br>{m["wrong"]}', f'<b>×</b>　{m["right"]}') for m in K['misconceptions']])

    # ---- 手で書いたデッキ
    d = Path(content) / 'cards'
    for f in sorted(d.glob('*.json')) if d.exists() else []:
        x = read(f)
        deck(x['id'], x['name'], x.get('field', 'THE'), x.get('desc', ''), [card(c['id'], c['f'], c['b'], c.get('n'), link=c.get('link')) for c in x['cards']])

    ids = [c['id'] for dk in decks for c in dk['cards']]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        raise SystemExit(f'カードの ID が重複しています：{sorted(dup)[:10]}')
    order = {'THE': 0, 'INO': 1, 'ORG': 2, 'POL': 3}
    decks.sort(key=lambda dk: order.get(dk['field'], 9))
    return decks
