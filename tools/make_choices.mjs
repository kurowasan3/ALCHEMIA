// 選択式の選択肢を作る：入力用の答え（grading）の語句・化学式・数値の欄に、まちがいの選択肢を3つ足す（Lemma の tools/make_choices.mjs の化学版）
//   node tools/make_choices.mjs [content フォルダ]        … content/choices.json を作り直す
//   node tools/make_choices.mjs [content] --check          … 書き出さずに、今の choices.json と同じか調べる
//
// 作り方（AI ではなく決まった手順なので、何度実行しても同じ結果）：
//   語{ }  … 同じ単元（足りなければ同じ分野・全体）のほかの問題の答えから、正解と同じ種類（化学式どうし・日本語の語句どうし）を選ぶ
//   約{ }  … 正解の数値を 2倍・1/2・10倍・1/10 などにして、同じ有効数字で書く
//   数（係数など）… 正解の数を ±1・2倍などに変える
//   どれも、アプリの採点（answer.js の checkField）で「不正解」になり、候補どうし・正解と重ならないことを確かめる。
//   選ぶ{ }・並べる{ }・すべて選ぶ{ } と自己採点の欄は、そのままの形で出すので選択肢を作らない（null）。
// choices.json：{ 問題ID: { spec: 入力用の答え, f: [ [正解, 誤答, …] | null, …（欄の順） ] } }
//   選択肢の書き方：化学式・語句はそのまま、数式は $…$（アプリの optHTML が描く）。
// 入力用の答えが変わった問題は spec が合わなくなるので、build_data.py は使わない（ここを実行し直す）。
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { parseSpec, checkField, normWord, roundSig, sciTeX, exprTeX, unitTeX } from '../js/answer.js';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2).filter(a => !a.startsWith('--'));
const content = path.resolve(args[0] || path.join(ROOT, 'content'));
const CHECK = process.argv.includes('--check');

// ---- 決まった乱数（問題ID と欄の番号から）
function seedOf(s) { let h = 2166136261; for (const c of s) { h ^= c.charCodeAt(0); h = Math.imul(h, 16777619); } return h >>> 0; }
function rng(seed) { let x = seed || 1; return () => { x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0; return x / 4294967296; }; }
function shuffled(list, r) { const a = list.slice(); for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(r() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; } return a; }

// ---- content/problems を読む（build_data.py の grading_spec と同じ形の spec を作る）
function gradingSpec(parts) {
  const lines = parts.filter(p => 'grading' in p).map(p => (p.id === 'main' ? p.grading : `(${p.id}) ${p.grading}`));
  return lines.length ? lines.join('\n') : null;
}
const probs = [];
for (const subj of fs.readdirSync(path.join(content, 'problems')).sort()) {
  const dir = path.join(content, 'problems', subj);
  if (!fs.statSync(dir).isDirectory()) continue;
  for (const id of fs.readdirSync(dir).sort()) {
    const f = path.join(dir, id, 'problem.json');
    if (!fs.existsSync(f)) continue;
    const p = JSON.parse(fs.readFileSync(f, 'utf8'));
    if (p.quality && p.quality.status !== 'published') continue;
    const spec = gradingSpec(p.answer.parts);
    if (!spec) continue;
    let parsed;
    try { parsed = parseSpec(spec); } catch { continue; }
    probs.push({ id: p.id, unit: p.unit, field: p.field || p.unit.split('-')[0], spec, fields: parsed.parts.flatMap(pt => pt.fields) });
  }
}

// ---- 語句・化学式の候補の集まり
const isFormula = (w) => /^[A-Z\[(][A-Za-z0-9()\[\]+\-^·.=#{}]*$/.test(String(w).normalize('NFKC').replace(/\s+/g, ''));
const kindOfWord = (w) => (isFormula(w) ? 'formula' : /[぀-ヿ一-鿿]/.test(w) ? 'ja' : 'other');
const words = [];   // { id, unit, field, w, kind }
for (const p of probs) for (const f of p.fields) if (f.kind === 'word') words.push({ id: p.id, unit: p.unit, field: p.field, w: f.accept[0], kind: kindOfWord(f.accept[0]) });

const wrong = (f, text) => { const r = checkField(f, { mode: 'plain', text }); return r && r.ok === false && !r.syntax; };

function wordChoices(p, f, k) {
  const r = rng(seedOf(`${p.id}#${k}`));
  const kind = kindOfWord(f.accept[0]);
  const len = [...f.accept[0]].length;
  const seen = new Set(f.norm);
  const out = [];
  // 近い順：同じ単元 → 同じ分野 → 全体。長さが近いものを先に（見ただけで外せる選択肢を避ける）
  // 日本語の語句は、分野が違うと見ただけで外せるので、同じ分野の中だけから選ぶ
  const tiers = kind === 'ja' ? [(x) => x.unit === p.unit, (x) => x.field === p.field] : [(x) => x.unit === p.unit, (x) => x.field === p.field, () => true];
  for (const inTier of tiers) {
    const cand = shuffled(words.filter(x => x.id !== p.id && x.kind === kind && inTier(x)), r)
      .sort((a, b) => Math.abs([...a.w].length - len) - Math.abs([...b.w].length - len));
    for (const c of cand) {
      const n = normWord(c.w);
      if (seen.has(n) || !wrong(f, c.w)) continue;
      if (kind === 'ja' && Math.abs([...c.w].length - len) > Math.max(4, len)) continue;
      seen.add(n); out.push(c.w);
      if (out.length >= 3) return [f.accept[0], ...out];
    }
  }
  return out.length >= 2 ? [f.accept[0], ...out] : null;
}

function numTeX(v, sig, unit) {
  return `$${sciTeX(v, sig)}${unitTeX(unit)}$`;
}
function approxChoices(p, f, k) {
  const r = rng(seedOf(`${p.id}#${k}`));
  const v = f.value;
  const muls = shuffled([2, 0.5, 10, 0.1, 1.5, 4, 0.25, 3, 1 / 3], r);
  const seen = new Set([sciTeX(v, f.sig)]);
  const out = [];
  for (const m of muls) {
    const w = roundSig(v * m, f.sig);
    if (!Number.isFinite(w) || w === 0) continue;
    const t = sciTeX(w, f.sig);
    if (seen.has(t) || !wrong(f, String(w))) continue;
    seen.add(t); out.push(numTeX(w, f.sig, f.unit));
    if (out.length >= 3) break;
  }
  return out.length >= 2 ? [numTeX(v, f.sig, f.unit), ...out] : null;
}

function exprChoices(p, f, k) {
  const src = String(f.value.src).trim();
  const unit = f.value.unit || '';
  const tex = (s) => { const t = exprTeX(s); return t == null ? null : `$${t}${unitTeX(unit)}$`; };
  const r = rng(seedOf(`${p.id}#${k}`));
  let alts = [];
  if (/^-?\d+$/.test(src)) {   // 係数・価数・個数など
    const v = parseInt(src, 10);
    alts = shuffled([v + 1, v - 1, v + 2, v * 2, v - 2, v + 3], r).filter(w => (v > 0 ? w > 0 : true) && w !== v).map(String);
  } else if (/^-?\d+(\.\d+)?$/.test(src)) {
    const v = Number(src), d = (src.split('.')[1] || '').length;
    alts = shuffled([v * 2, v / 2, v * 10, v / 10, v + 10 ** -d, v - 10 ** -d], r).filter(w => w > 0 && w !== v).map(w => String(+w.toFixed(d + 1)));
  } else {
    // 数字を1つずつ変える（2/3 → 1/3・4/3 など）
    for (const m of src.matchAll(/\d+/g)) {
      const v = Number(m[0]);
      for (const w of [v + 1, v - 1, v * 2]) {
        if (w <= 0 || w === v) continue;
        alts.push(src.slice(0, m.index) + w + src.slice(m.index + m[0].length));
      }
    }
    alts = shuffled(alts, r);
  }
  const right = tex(src);
  if (!right) return null;
  const seen = new Set([right]), out = [];
  for (const a of alts) {
    const t = tex(a);
    if (!t || seen.has(t) || !wrong(f, a)) continue;
    seen.add(t); out.push(t);
    if (out.length >= 3) break;
  }
  return out.length >= 2 ? [right, ...out] : null;
}

const result = {};
let nFields = 0, nMade = 0;
for (const p of probs) {
  const f = p.fields.map((fl, k) => {
    if (fl.kind === 'word') { nFields++; return wordChoices(p, fl, k); }
    if (fl.kind === 'approx') { nFields++; return approxChoices(p, fl, k); }
    if (fl.kind === 'expr') { nFields++; return exprChoices(p, fl, k); }
    return null;
  });
  nMade += f.filter(Boolean).length;
  if (f.some(Boolean)) result[p.id] = { spec: p.spec, f };
}
const sorted = Object.fromEntries(Object.keys(result).sort().map(k => [k, result[k]]));
const text = JSON.stringify(sorted, null, 1) + '\n';
const out = path.join(content, 'choices.json');
if (CHECK) {
  const cur = fs.existsSync(out) ? fs.readFileSync(out, 'utf8') : '';
  if (cur !== text) { console.log('content/choices.json が古くなっています。node tools/make_choices.mjs を実行してください'); process.exit(1); }
  console.log(`choices.json は最新です（${Object.keys(sorted).length}問）`);
} else {
  fs.writeFileSync(out, text);
  console.log(`選択肢を作りました：${Object.keys(sorted).length}問（語句・化学式・数値の欄 ${nFields}のうち ${nMade}欄）→ ${path.relative(ROOT, out)}`);
}
