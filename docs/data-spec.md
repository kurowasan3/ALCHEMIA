# ALCHEMIA データ仕様（v0.1）

正規データ `content/`（マスター）→ `tools/build_data.py` → アプリ用データ `data/`。アプリが読むのは `data/` だけ。
JSON Schema は `tools/schema/`、Schema で見られない関係の検査は `tools/validate_content.py`。

## ID の規則

| ID | 形 | 規則 |
| --- | --- | --- |
| 問題ID | `INO-PPT-001` | `<分野3字>-<単元略号3字>-<番号>`。一度使ったら変えない・再利用しない |
| 分野ID | `THE` `INO` `ORG` `POL` | 理論化学・無機化学・有機化学・高分子化合物 |
| 単元ID | `INO-08` | `taxonomy.json` で固定。使わなくなったら `retired` に残し、再利用しない |
| 概念ID・法則ID | `normal-force` `newton-2nd` | 英小文字とハイフン |
| 誤概念ID | `MC-EQ-STOP` | `MC-` で始まる英大文字 |
| 図ID | `q-setup` `a-forces` | 問題の中で一意。英小文字とハイフン |

## 問題（`content/problems/<分野>/<問題ID>/problem.json`、`alchemia.problem/1`）

| 項目 | 必須 | 内容 |
| --- | --- | --- |
| `id` `schema` `field` `unit` | ○ | 問題ID、`"alchemia.problem/1"`、分野ID、単元ID |
| `topic` `course` `minutes` | | テーマ、`化学基礎`／`化学`、想定時間（分） |
| `level` | ○ | 0=基礎 1=標準 2=入試基礎 3=国公立2次 4=難関 5=最難関 |
| `format` | ○ | 出題形式：`qa` 一問一答・`equation` 反応式補完・`compare` 比較・`condition` 条件判断・`pathway` 反応経路・`calc` 計算・`structure` 構造決定・`written` 記述ほか（複数可）。ホームの「形式ごとに演習」はこれで選ぶ |
| `examStyle` | | `common-test` `national-2nd` `private` `original` |
| `statement.blocks` | ○ | 本文（下の「ブロック」） |
| `figures` | ○ | 図の一覧（下の「図」）。無ければ `[]` |
| `answer.parts` | ○ | 小問ごとの答え `{id, tex, grading?, wrong?}`。小問が無い問題は `id: "main"` 1つ |
| `solutions` | ○ | `{id: "main"|"alt1"|"note", kind, blocks}`。別解・補足を足せる |
| `idea` | | 核となる考え方 `{tex}` |
| `zukan` | | 解いたあとに開く図鑑のページ `{page, key?, label?}`。page は `ions` `reactions` `organic` `tests` `polymers` `elements` `laws`、key はイオン（`Cu^2+`）・反応ID・化合物ID・高分子ID・元素記号・法則ID |
| `choices` | | 選択問題の選択肢。誤答の選択肢に `misconception` を付けると、選んだときに誤概念として記録できる（Phase 5） |
| `concepts` `laws` `misconceptions` | ○ | 概念・法則・誤概念の ID（無ければ `[]`） |
| `tags` `provenance` `quality` | ○ | `quality.status` は `draft`／`review`／`published`／`hidden`。アプリに出るのは `published` だけ |

### ブロック

```json
{"type": "text", "tex": "文章（LaTeX。$...$ が数式）"}
{"type": "figure", "ref": "図ID", "placement": "block | float-right | side-by-side"}
{"type": "parts", "items": [{"id": "1", "tex": "…", "figures": ["図ID"]}]}
```

並びがそのまま表示順。図の位置も、ブロックの位置で決まる（`order` を手で管理しない）。

### 入力用の答え（`answer.parts[].grading`）

1行に1小問（`(1) …`）、小問の中の複数の欄は ` ; ` で区切り、`ラベル：` で欄の名前を付ける（Lemma・CELESTIA と同じ）。化学で足した書き方：

| 書き方 | 入力 | 採点 |
| --- | --- | --- |
| `語{アンモニア\|NH3}` | 文字（化学式は入力中に下付き・上付きで表示） | `\|` で区切ったどれかと一致。全角半角・ひらがなカタカナ・空白・`^ { } _`・中黒・下付き上付き文字の違いは無視。大文字小文字は区別（Co と CO） |
| `選ぶ{塩酸\|硫酸\|硝酸}=2` | ボタンで1つ選ぶ | `=` の後は何番目か |
| `すべて選ぶ{A\|B\|C\|D}=1,3` | ボタンで複数選ぶ | 過不足なく選んだときだけ正解 |
| `並べる{HF\|HCl\|HBr}=3>2>1` | タップした順に並ぶ | `=` の後に正しい順の番号を `>` でつなぐ |
| `約{2.24\|3} [L]` | 数値（`1.2e-3`・`1.2×10^-3` も可） | `\|` の後の有効数字の桁で、最後の桁の ±0.55 まで正解（省略すると2桁） |
| `2` | 数（係数など） | 整数だけの答えは小さな欄・数字キーボードになる |
| `記述` | なし | 自己採点 |

選択肢・並べかえの札には `$\ce{…}$` の数式を書ける。化学式だけの札（`HCl` など）は自動で化学式として表示する。

### 反応式

文章中の反応式は KaTeX の mhchem で `$\ce{2Na + 2H2O -> 2NaOH + H2}$` と書く。2価以上の電荷は `^`（`Cu^2+`、`SO4^2-`）、1価は `Ag+`・`OH-`。
`tools/validate_content.py` が、矢印が1つで `+` か係数のある完全な式について、原子の数と電荷のつり合いを調べる（空欄 `\boxed`・一般式の R・重合度 n を含む式は調べない）。

## 図

```json
{"id": "a-product", "role": "answer", "kind": "structure", "format": "svg", "src": "a-product.svg",
 "alt": "生成物 X の構造式：ベンゼン環にヒドロキシ基が1つ", "caption": "…",
 "smiles": "Oc1ccccc1"}
```

| 項目 | 内容 |
| --- | --- |
| `role` | `question`（問題用）／`answer`（解答用）／`both`。**解答用の図は問題文から参照できない**（検査でエラー）。アプリは答え合わせまで表示もファイルの読み込みもしない |
| `kind` | `diagram` `structure`（構造式） `scheme`（反応経路） `apparatus`（装置） `graph` `table` `photo` |
| `format` | `svg`（標準）、`png`／`webp`（写真） |
| `src` | 問題フォルダ内のファイル名。Base64 の埋め込みは不可 |
| `alt` | 代替テキスト（必須。図が見えない人が問題を解けるように） |

### 構造式の図

図に `"smiles": "Oc1ccccc1"` を書くと，`tools/structgen.py` が RDKit で構造式の SVG（`src` のファイル）を作る。`*` は高分子の鎖が続く印（波線）になる。

### 手描き・AI 生成の SVG の決まり

- `viewBox` を付ける。`script`・`foreignObject`・`on*` 属性・外部参照・埋め込み画像は使わない（検査でエラー）
- 色は `currentColor`（明暗の設定に合わせて変わる）。色を直接指定した図はダークモードでも白い紙の上に出る

## 概念・法則・誤概念

| ファイル | 1件の形 |
| --- | --- |
| `concepts.json` | `{id, name, field, requires: [前提の概念], laws: [法則]}`（前提関係の循環は検査でエラー） |
| `laws.json` | `{id, name, tex, conditions, concepts, origin, basis?, derivation?}`。`origin` は `fundamental`／`empirical`／`derived`／`definition`。`derived` は `derivation`（1段ずつ）が必須、`fundamental`・`empirical` は `basis` が必須 |
| `misconceptions.json` | `{id, wrong, right, concepts, remedy: [補強問題の ID]}` |

## 図鑑のデータ（`content/chem/`）

| ファイル | 中身 |
| --- | --- |
| `ions.json` | `cations[{id, name, aq, flame, group}]`、`reagents[{id, name, ion}]`、`results[イオン][試薬] = {kind, f?, color?, eq?, note?, name?}`（kind：`ppt` 沈殿／`dissolve` 過剰で溶ける／`remain` 溶けない／`solution` 溶液の色が変わる／`none` 変化なし／`na` 扱わない）、`systematic`（系統分析の6段階）、`colors`・`flameColors` |
| `reactions.json` | `reactions[{id, unit, cat: gas｜industry｜reaction, title, eq｜steps[{eq, note}], total?, cond, type, gas?, note?}]`、`gases{化学式: {name, color, smell, nature, sol, collect, dry, test}}`、`types` |
| `organic.json` | `compounds[{id, name, cf: 示性式, smiles, cls}]`、`reactions[{id, from, to, reagent, type, eq?, note?, map}]`、`maps{aliphatic｜aromatic: {pos{化合物: [列, 行]}}}` |
| `tests.json` | `tests[{id, name, reagent, obs, target, scope: [対象の化合物], pos{化合物: 補足}, tri{…}, neg_note{…}}]`（scope のうち pos・tri 以外は陰性） |
| `polymers.json` | `polymers[{id, unit, name, abbr?, cls, monomers[{name, smiles}], poly, repeat?: SMILES(*), rf?: 繰り返し単位の式, bond?, heat?, props, uses}]` |
| `elements.json` | `elements[{z, sym, ja, mass, exam?, group, period, block, cat, conf, orbitals, shells, en, ie1, radius, state}]` |

## アプリ用データ（`data/`、`alchemia.data/1`）

| ファイル | 中身 |
| --- | --- |
| `index.json` | `format` `dataVersion` `levels` `subjects`（分野の並び）`fields`（単元）`problems`（一覧：`id f lv t th tg fm c m zk in ok`） |
| `f-<単元ID>.json` | 本文：`q a s k`（HTML）、`figs`（図の情報）、`in`（入力用の答え）、`zk` `choices` `laws` `alt`（別解） |
| `knowledge.json` | 概念・法則・誤概念 |
| `zukan.json` | 図鑑のデータ（content/chem/ をまとめたもの）。`struct` は構造式 SVG の場所 |
| `search.json` `collections.json` | 検索の索引（概念名・法則名も入る）、問題セット |
| `assets/<問題ID>/` | 図のファイル。名前に中身のハッシュ（`q-setup-35755270.svg`）が付くので、図を直すとキャッシュも入れ替わる |

本文の HTML では、図は `<figure class="pfig" data-fig="図ID"></figure>` の空の枠で、`js/problem/figures.js` が中身を入れる。

## 版は3つを分ける（Lemma と同じ）

| 版 | 場所 | いつ上げるか |
| --- | --- | --- |
| コードの版 | 各 JS の `FILE_VERSION`、CSS の `--alchemia-version`、`sw.js` の `VERSION` | アプリを更新したとき（食い違うと赤い警告） |
| データ形式の版 | `index.json` の `format` | `data/` の項目の意味を変えたとき |
| 学習記録の版 | localStorage `alchemia.v1` の `schema`（今は 1） | 学習記録の形を変えたとき |
