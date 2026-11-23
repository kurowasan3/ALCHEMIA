// 学習記録の保存・点検・バックアップ・移行（js/store.js。Lemma v0.27 と同じ仕組み）の単体テスト
//   node tests/store.test.mjs
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const S = await import('../js/store.js');

let fail = 0, pass = 0;
const ok = (cond, msg) => { if (cond) pass++; else { fail++; console.log('FAIL', msg); } };
const clone = (x) => JSON.parse(JSON.stringify(x));
const base = JSON.parse(fs.readFileSync(path.join(root, 'tests/fixtures/state-v1.json'), 'utf8'));   // v0.2（schema 1）の実際の形
const NOW = Date.UTC(2026, 10, 23);

function memStorage({ quota = Infinity } = {}) {
  const m = new Map();
  let failNext = false;
  return {
    m, failNext: () => { failNext = true; },
    getItem: k => (m.has(k) ? m.get(k) : null),
    setItem(k, v) {
      if (failNext) { failNext = false; throw new Error('write failed'); }
      const size = [...m.entries()].reduce((a, [kk, vv]) => a + (kk === k ? 0 : kk.length + vv.length), 0) + k.length + String(v).length;
      if (size > quota) { const e = new Error('quota'); e.name = 'QuotaExceededError'; throw e; }
      m.set(k, String(v));
    },
    removeItem: k => m.delete(k),
    key: i => [...m.keys()][i] ?? null,
    get length() { return m.size; },
  };
}

// ---------------------------------------------------------------- 形・移行
{
  ok(!S.structural(S.freshState()).length, '新しい記録は正しい形');
  ok(S.KEY === 'alchemia.v1' && S.BACKUP_FORMAT === 'alchemia-backup', 'ALCHEMIA のキー名と形式名');
  const m = S.migrate(clone(base));
  ok(m.state.schema === S.SCHEMA && m.applied.includes('schema 1→2'), 'schema 1 の記録を移行する');
  ok(m.state.log.length === base.log.length && S.stable(m.state.pstat) === S.stable(base.pstat), '移行しても解答の記録・復習状態は変わらない');
  ok(m.state.bookmarks && typeof m.state.bookmarks === 'object' && m.state.cards, '足りない項目（ブックマーク）を補う');
  ok(base.schema === 1 && !('bookmarks' in base), '元の記録は変えない');
  const nw = S.migrate({ ...clone(base), schema: 9, future: 1 });
  ok(nw.newer && nw.state.future === 1, '新しい版の記録は形を変えない');
}

// ---------------------------------------------------------------- 点検
{
  const st = S.migrate(clone(base)).state;
  st.cards = { 'ions-ag': { box: 2, due: NOW, last: NOW - 1000, n: 3 }, bad: { S: 1 } };
  st.bookmarks = { 'INO-HAL-001': { on: true, t: NOW } };
  const rep = S.audit(st, { problemIds: new Set(base.log.map(l => l.p).concat(Object.keys(base.pstat), Object.keys(base.problems), ['INO-HAL-001'])), cardIds: new Set(['ions-ag']), now: NOW });
  ok(rep.errors.some(e => e.code === 'card_invalid' && e.ids.includes('bad')), '暗記カードの記録の形（box・due・last・n）を調べる');
  ok(rep.warnings.some(w => w.code === 'card_unknown' && w.ids.includes('bad')), '今のデータに無いカードIDは注意');
  ok(!rep.errors.some(e => e.code === 'card_invalid' && e.ids.includes('ions-ag')), '正しいカードの記録は通る');
}

// ---------------------------------------------------------------- バックアップ
{
  const st = S.migrate(clone(base)).state;
  st.cards = { 'ions-ag': { box: 2, due: NOW, last: NOW - 1000, n: 3 } };
  const b = S.makeBackup(st, '0.3.0', NOW);
  ok(b.format === 'alchemia-backup' && b.version === 1 && b.counts.cards === 1 && /^fnv1a32:/.test(b.hash), 'バックアップに形式・件数・確認用の値');
  const back = S.readBackup(JSON.parse(JSON.stringify(b)), { now: NOW });
  ok(back.state.log.length === st.log.length && back.state.cards['ions-ag'].box === 2, '書き出した記録をそのまま読める');
  const broken = clone(b); broken.state.log.pop();
  let msg = ''; try { S.readBackup(broken, { now: NOW }); } catch (e) { msg = e.code; }
  ok(msg === 'hash_mismatch', '書き換えられたバックアップは断る');
  msg = ''; try { S.readBackup({ format: 'lemma-backup', version: 1 }, { now: NOW }); } catch (e) { msg = e.code; }
  ok(msg === 'not_lemma', 'Lemma など別のアプリのバックアップは断る');
  const old = S.readBackup(clone(base), { now: NOW });   // v0.2 の書き出し（記録そのもの）
  ok(old.meta.format === 0 && old.state.schema === S.SCHEMA, 'v0.2 までの書き出し（記録そのもの）も移行して読む');
  const cur = S.migrate(clone(base)).state; cur.settings = { ...cur.settings, theme: 'dark' };
  const rep = S.restoreState(cur, back.state, 'replace', NOW);
  ok(rep.resetAt === NOW && rep.settings.theme === 'dark', '置き換え：消した印を付け、設定は今のまま');
  const mg = S.restoreState(cur, back.state, 'merge', NOW);
  ok(mg.cards['ions-ag'] && mg.log.length === cur.log.length, '合わせる：両方の記録を残す');
}

// ---------------------------------------------------------------- 保存先
{
  const ls = memStorage();
  const st = S.createStore(ls, { now: () => NOW, app: 't' });
  let L = st.load();
  ok(L.status === 'new', '記録が無いときは新しく始める');
  L.state.log.push({ p: 'INO-HAL-001', t: NOW, g: 3 });
  ok(st.save(L.state).ok && JSON.parse(ls.getItem(S.KEY)).log.length === 1, '保存できる');
  // 別のタブ：先に保存された記録と合わせる
  const other = JSON.parse(ls.getItem(S.KEY)); other.log.push({ p: 'INO-HAL-002', t: NOW + 5, g: 1 }); ls.setItem(S.KEY, JSON.stringify(other));
  L.state.log.push({ p: 'INO-HAL-003', t: NOW + 9, g: 4 });
  const r = st.save(L.state);
  ok(r.ok && r.merged && r.merged.log.length === 3, '別のタブの保存とは上書きせずに合わせる');
  // 読めない記録：退避してから新しく始める
  const ls2 = memStorage(); ls2.setItem(S.KEY, '{壊れた');
  const st2 = S.createStore(ls2, { now: () => NOW });
  const L2 = st2.load();
  ok(L2.status === 'broken' && L2.quarantined && st2.brokenKeys().length === 1, '読めない記録は消さずに退避する');
  // 移行の前に控えを残す
  const ls3 = memStorage(); ls3.setItem(S.KEY, JSON.stringify(base));
  const st3 = S.createStore(ls3, { now: () => NOW });
  const L3 = st3.load();
  ok(L3.status === 'migrated' && st3.snapshots().length === 1 && JSON.parse(st3.snapshots()[0].raw).schema === 1, '移行の前の記録を控えに残す');
  // 書き込みに失敗したら前の記録のまま
  const before = ls3.getItem(S.KEY);
  L3.state.log.push({ p: 'X', t: NOW, g: 3 }); ls3.failNext();
  const w = st3.save(L3.state);
  ok(!w.ok && w.reason === 'write' && ls3.getItem(S.KEY) === before, '保存に失敗しても前の記録は残る');
  // 新しい版の記録には書かない
  const ls4 = memStorage(); ls4.setItem(S.KEY, JSON.stringify({ ...base, schema: 99 }));
  const st4 = S.createStore(ls4, { now: () => NOW });
  const L4 = st4.load();
  ok(L4.status === 'newer' && !st4.save(L4.state).ok, '新しい版のアプリの記録は保存しない');
  // 置き換え（fence）：ほかのタブの古い記録を持ち込ませない
  const rr = st.replace(L.state, { ...S.freshState() }, null, { fence: true });
  ok(rr.ok && rr.state.localResetAt === NOW && !rr.state.log.length, 'ログアウトして消す：localResetAt を付けて置き換える');
  const stale = { ...S.freshState(), log: [{ p: 'OLD', t: NOW - 1000, g: 3 }] };
  const m = S.mergeTabs(rr.state, stale);
  ok(!m.log.length, '消したあとに、古いタブの記録はよみがえらない');
}

// ---------------------------------------------------------------- 合わせ方（化学版の項目）
{
  const A = S.migrate(clone(base)).state, B = S.migrate(clone(base)).state;
  A.cards = { c1: { box: 3, due: NOW + 9, last: NOW, n: 4 } };
  B.cards = { c1: { box: 1, due: NOW, last: NOW - 50, n: 2 }, c2: { box: 0, due: NOW, last: NOW - 10, n: 1 } };
  A.bookmarks = { p1: { on: false, t: NOW } }; B.bookmarks = { p1: { on: true, t: NOW - 100 } };
  A.misc = { m1: { n: 2, last: 10 } }; B.misc = { m1: { n: 5, last: 5, cleared: true } };
  const M = S.mergeState(A, B);
  ok(M.cards.c1.box === 3 && M.cards.c2, '暗記カードはカードごとに最後に答えた方');
  ok(M.bookmarks.p1.on === false, 'ブックマークは最後に付け外しした方（外した印も残す）');
  ok(M.misc.m1.n === 5 && M.misc.m1.last === 10, '誤概念の検出は回数の多い方・最後の時刻の新しい方');
  const R = { ...S.migrate(clone(base)).state, resetAt: NOW };
  const M2 = S.mergeState(R, B);
  ok(!M2.cards.c1 && M2.log.every(l => l.t >= NOW), '消した印より前の記録は、印の古い側から持ち込まない');
}

console.log(`store.test: ${pass} 件 OK${fail ? `、${fail} 件失敗` : ''}`);
process.exit(fail ? 1 : 0);
