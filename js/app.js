import { mountFigures, setFigureEsc } from './problem/figures.js';
// 版の確認用（古いファイルに FILE_VERSION が無くても読み込みが止まらないよう名前空間で取り込む）
import * as M_FIG from './problem/figures.js';
import * as M_SCHED from './scheduler.js';
import * as M_PARSER from './parser.js';
import * as M_CAT from './catalog.js';
import * as M_ANS from './answer.js';
import * as M_UI from './ui.js';
import * as M_STORE from './store.js';
import { freshState, createStore, audit, makeBackup, readBackup, restoreState, counts, mergeState, BROKEN_PREFIX } from './store.js';
import { parseSpec, checkField, fieldHint, fieldUnit, fieldTeX, condTeX, exprTeX, fieldAnswerHTML, ceFix, CHEM_KINDS } from './answer.js';
import { ICON, play, setSound, applyTheme, MASCOT_SVG, MASCOT_NAME, mascotMessage } from './ui.js';
import { review, retrievability, dayKey } from './scheduler.js';
import {
  LEVELS, SUBJECT_ORDER, FIELDS, ENTRIES, META, FORMATS, fieldByKey, entryById, loadCatalog, ensureBodies, getBody,
  VALUES, subjectsWithFields, pool, pickPractice, pickReview, dueEntries, nextReviewTime, pickToday, similarTo, loadKnowledge, loadZukan,
} from './catalog.js';

const VERSION = '0.3.0';
let sync = null, syncApplying = false, syncStatus = null;   // アカウントと同期（initSync で用意）
const $ = (sel, el = document) => el.querySelector(sel);
const $$ = (sel, el = document) => [...el.querySelectorAll(sel)];

// ---------- 保存（Lemma v0.23 と同じ：読み書き・移行・控えは js/store.js） ----------
//   読めない記録は上書きせずに退避、保存に失敗したら前の記録のまま、別のタブの保存とは合わせてから保存する。
// localStorage そのものを使えない環境（一部のプライベートモード・保存を禁止した設定）だけ、この画面の中だけで動く。
let LS, noStorage = false;
try { LS = window.localStorage; if (!LS) throw new Error('no localStorage'); } catch {
  noStorage = true;
  const m = new Map();
  LS = { getItem: k => (m.has(k) ? m.get(k) : null), setItem: (k, v) => m.set(k, String(v)), removeItem: k => m.delete(k), key: i => [...m.keys()][i] ?? null, get length() { return m.size; } };
}
const store = createStore(LS, { app: VERSION });
const boot = store.load();
let state = boot.state;
// 設定・選んだ分野などが変わった時刻（端末どうしで食い違ったら新しい方を採る）
let lwwSig = null;
let saveWarned = false;
let persistAsked = false;
const SAVE_MSG = {
  quota: 'この端末の保存領域がいっぱいで、記録を保存できませんでした（前に保存した記録は残っています）。設定から書き出しておいてください。',
  write: 'この端末に記録を書き込めませんでした（前に保存した記録は残っています）。設定から書き出しておいてください。',
  read: 'この端末に保存した記録を読めないため、上書きしないよう保存を止めています。しばらくしてから再読み込みしてください。',
  verify: '記録を保存できたか確かめられませんでした。設定から書き出しておいてください。',
  verify_read: '記録を保存できたか確かめられませんでした。設定から書き出しておいてください。',
  conflict: 'ほかの画面と同時に保存が重なりました。この画面の記録は残っています。もう一度操作すると保存します。',
  invalid: '記録の形が正しくないため、保存を止めました（前に保存した記録は残っています）。設定から書き出しておいてください。',
  broken: '読めなかった記録を退避できなかったため、この端末への保存を止めています。設定から書き出しておいてください。',
  migrate: '記録を新しい形に移せなかったため、元の記録を守るよう保存を止めています。アプリを更新してから再読み込みしてください。',
  newer: '新しい版の ALCHEMIA で保存された記録なので、この版では保存しません（元の記録はそのまま）。ページを再読み込みしてアプリを更新してください。',
  serialize: '記録を保存できませんでした。設定から書き出しておいてください。',
};
const BIG = 3e6;   // 保存している記録がこの文字数を超えたら1回知らせる（localStorage は多くのブラウザで約500万文字）
let bigWarned = false;
const lwwOf = (st) => JSON.stringify([st.settings, st.practice, st.zukanPrefs, st.cardPrefs]);
function save() {
  const sig = lwwOf(state);
  if (lwwSig !== null && sig !== lwwSig) state.settingsAt = Date.now();
  lwwSig = sig;
  const r = store.save(state);
  // 別のタブが先に保存していた：その記録と合わせた結果を、この画面の記録にも入れる（同じオブジェクトのまま）
  if (r.merged) {
    for (const k of Object.keys(state)) if (!(k in r.merged)) delete state[k];
    for (const k of Object.keys(r.merged)) state[k] = r.merged[k];
    lwwSig = lwwOf(state);   // ほかのタブで変えた設定を、この画面で変えたことにしない
  }
  if (!r.ok && !saveWarned) {   // 黙って失うと気づけないので、1回だけ知らせる
    saveWarned = true; setTimeout(() => toast(SAVE_MSG[r.reason] || SAVE_MSG.serialize), 0);
  }
  if (r.ok) saveWarned = false;
  if (r.ok && r.size > BIG && !bigWarned) { bigWarned = true; setTimeout(() => toast('保存している学習記録が大きくなり、この端末の保存領域の上限に近づいています。設定から書き出しておいてください。'), 0); }
  if (sync && !syncApplying) sync.schedule();
  return r.ok;
}

// ---------- アカウントと同期（js/sync.js を読めないときは、この機能なしで動く）
function initSync() {
  import('./sync.js').then((M) => {
    if (M.FILE_VERSION !== VERSION) showBanner(`<b>古いファイルが混ざっています（v${VERSION}）</b>js/sync.js が更新されていません。リポジトリにzipの中身をすべて上書きアップロードしてください。`);
    if (!M.enabled()) return;
    sync = M.createSync({
      getState: () => state,
      applyState: (merged) => {
        syncApplying = true;
        // 通信している間にこの端末で解いた分を落とさないよう、今の記録ともう一度合わせてから入れる
        merged = mergeState(state, merged);
        for (const k of Object.keys(merged)) state[k] = merged[k];
        save();
        syncApplying = false;
        // 解いている途中でなければ、画面を描き直して新しい記録を見せる
        const tab = (location.hash || '#today').slice(1).split('=')[0];
        if (!session && ['today', 'map', 'settings'].includes(tab) && !document.activeElement?.matches('input, textarea')) route();
      },
      onStatus: (s) => { syncStatus = s; const el = $('#acct'); if (el) drawAccount(); },
    });
    syncStatus = sync.status();
    if ((location.hash || '').startsWith('#settings')) renderSettings();
    if (sync.loggedIn()) sync.syncNow({ auto: true });
    document.addEventListener('visibilitychange', () => { if (document.visibilityState === 'hidden') sync.syncNow({ auto: true }); });
    window.addEventListener('online', () => sync.syncNow({ auto: true }));
  }).catch(() => { /* 同期なし */ });
}
// ログイン直後：端末とアカウントの両方に記録があるときだけ、合わせるかを聞く（共用の端末で他人の記録が混ざらないように）
async function firstSyncAfterLogin() {
  let remote = null;
  try { remote = await sync.peek(); } catch { /* 取れなければいつもの同期に任せる */ }
  if (remote && remote.n && hasRecords(state)) {
    const ln = state.log.length, last = remote.last ? new Date(remote.last).toLocaleDateString('ja-JP') : '';
    const k = await choose({ title: 'この端末の記録をアカウントに合わせますか？',
      body: `この端末に解答の記録が ${ln}件、アカウントに ${remote.n}件（最後は ${esc(last)}）あります。`,
      options: [{ key: 'merge', label: '合わせる（両方の記録を残す）', cls: 'primary' },
        { key: 'remote', label: 'アカウントの記録だけを使う（この端末の分は書き出してから消す）' }] });
    if (k === 'remote') {
      download(`alchemia-この端末の記録-${dayKey()}.json`, JSON.stringify(makeBackup(state, VERSION), null, 1), 'application/json');
      // 端末の記録を空にしてから同期すると、アカウントの記録だけになる（設定は端末のまま）
      //   fence：ほかのタブに残っている古い記録が、あとで持ち込まれないようにする
      const next = { ...freshState(), ...keepPrefs() };
      const r = store.replace(state, next, null, { fence: true });
      if (!r.ok) { toast(`${SAVE_MSG[r.reason] || SAVE_MSG.write}（切り替えていません）`); await sync.logout(); return; }
      state = r.state; lwwSig = null; session = null; clearSavedSession();
      await sync.syncNow(); toast('アカウントの記録に切り替えました'); route(); return;
    }
    if (!k) { await sync.logout(); toast('ログインを取りやめました'); return; }
  }
  await sync.syncNow(); toast('ログインしました'); route();
}
// 記録を消す・切り替えるときにも残すもの（設定・図鑑とカードの表示）
const keepPrefs = () => ({ settings: state.settings, practice: state.practice, settingsAt: state.settingsAt, zukanPrefs: state.zukanPrefs, cardPrefs: state.cardPrefs });
function drawAccount() {
  const box = $('#acct');
  if (!box) return;
  if (!sync) { box.innerHTML = '<p class="muted" style="margin:0">この版では同期を使えません。</p>'; return; }
  const s = syncStatus || sync.status();
  if (!sync.loggedIn()) {
    box.innerHTML = `
      <p style="margin:0 0 10px">ログインすると、学習記録を複数の端末で使え、サーバーにも保存されます。ログインしなくても、今まで通りこの端末だけで使えます。</p>
      <div id="acct-step1">
        <input class="field-in" type="email" id="acct-email" placeholder="メールアドレス" autocomplete="email" inputmode="email">
        <button class="btn primary block" id="acct-send" style="margin-top:8px">ログイン用のコードを送る</button>
      </div>
      <div id="acct-step2" hidden>
        <p class="muted" style="margin:0 0 6px"><span id="acct-to"></span> に届いたコード（数字）を入れてください（届くまで1分ほどかかることがあります）。</p>
        <input class="field-in code" type="text" id="acct-code" placeholder="123456" inputmode="numeric" autocomplete="one-time-code" maxlength="10">
        <button class="btn primary block" id="acct-verify" style="margin-top:8px">ログイン</button>
        <button class="linkbtn" id="acct-back" style="margin-top:6px">メールアドレスを入れ直す</button>
      </div>
      <p class="acct-err" id="acct-err">${s.error ? esc(s.error) : ''}</p>
      <p class="muted small" style="margin:6px 0 0">保存するのはメールアドレスと学習記録だけです。<a href="privacy.html">プライバシーポリシー</a>・<a href="terms.html">利用規約</a></p>`;
    const err = (m) => { $('#acct-err').textContent = m || ''; };
    $('#acct-send').addEventListener('click', async (ev) => {
      const email = $('#acct-email').value.trim();
      if (!/^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(email)) { err('メールアドレスを確かめてください。'); return; }
      ev.target.disabled = true; err('');
      try {
        await sync.sendCode(email);
        $('#acct-step1').hidden = true; $('#acct-step2').hidden = false; $('#acct-to').textContent = email;
        $('#acct-code').focus();
      } catch (e) { err(e.message); }
      ev.target.disabled = false;
    });
    $('#acct-email').addEventListener('keydown', (ev) => { if (ev.key === 'Enter') $('#acct-send').click(); });
    $('#acct-code').addEventListener('keydown', (ev) => { if (ev.key === 'Enter') $('#acct-verify').click(); });
    $('#acct-back').addEventListener('click', () => { $('#acct-step1').hidden = false; $('#acct-step2').hidden = true; err(''); });
    $('#acct-verify').addEventListener('click', async (ev) => {
      ev.target.disabled = true; err('');
      try {
        await sync.verify($('#acct-email').value.trim(), $('#acct-code').value, { autoSync: false });
        await firstSyncAfterLogin();
        drawAccount();
      } catch (e) { err(e.message); ev.target.disabled = false; }
    });
    return;
  }
  const when = s.lastSync ? new Date(s.lastSync).toLocaleString('ja-JP', { month: 'numeric', day: 'numeric', hour: '2-digit', minute: '2-digit' }) : 'まだ';
  const label = { syncing: '同期しています…', ok: `最後の同期：${when}`, idle: `最後の同期：${when}`, offline: `オフラインのため同期待ち（最後の同期：${when}）`, error: '同期できませんでした' }[s.status] || '';
  box.innerHTML = `
    <div class="setting" style="padding-top:0"><span><b>${esc(s.email || '')}</b><br><span class="muted small" id="acct-state">${esc(label)}</span></span>
      <button class="btn" id="acct-now" ${s.status === 'syncing' ? 'disabled' : ''}>今すぐ同期</button></div>
    ${s.error && s.status !== 'offline' ? `<p class="acct-err">${esc(s.error)}${s.halted ? '<br>同じ失敗を繰り返さないよう、自動の同期を止めています。対処したら「今すぐ同期」を押してください。' : ''}</p>` : ''}
    <p class="muted small" style="margin:8px 0">アカウントは Lemma と共通です（記録は別々に保存）。この端末の記録とサーバーの記録は、解答の記録を足し合わせ、復習の予定は新しい方を使って合わせます。</p>
    <div class="btn-row" style="margin-top:4px"><button class="btn ghost" id="acct-out">ログアウト</button><button class="btn ghost danger" id="acct-del">退会する</button></div>
    <p class="muted small" style="margin:6px 0 0"><a href="privacy.html">プライバシーポリシー</a>・<a href="terms.html">利用規約</a></p>`;
  $('#acct-now').addEventListener('click', () => sync.syncNow());
  $('#acct-out').addEventListener('click', async () => {
    const k = await choose({ title: 'ログアウト', body: 'サーバーの記録はそのまま残ります。この端末の記録をどうしますか？（家族や学校の共用の端末なら、消すのがおすすめです）',
      options: [{ key: 'keep', label: 'この端末の記録は残してログアウト', cls: 'primary' }, { key: 'wipe', label: 'この端末の記録も消してログアウト' }] });
    if (!k) return;
    const synced = await sync.syncNow();
    // 消す前に、サーバーに送れたことを確かめる（送れていない記録を消さない）
    if (k === 'wipe' && !synced) { toast('サーバーと同期できなかったため、この端末の記録は消していません（ログアウトもしていません）'); drawAccount(); return; }
    await sync.logout();
    if (k === 'wipe') {
      const r = store.replace(state, { ...freshState(), settings: state.settings }, null, { fence: true });
      if (!r.ok) { toast(`${SAVE_MSG[r.reason] || SAVE_MSG.write}（この端末の記録は消していません）`); route(); return; }
      state = r.state; lwwSig = null; session = null; clearSavedSession();
    }
    toast(k === 'wipe' ? 'ログアウトし、この端末の記録を消しました' : 'ログアウトしました'); route();
  });
  $('#acct-del').addEventListener('click', async () => {
    if (!confirm('アカウントは Lemma（受験数学）と共通です。退会すると、アカウントと、ALCHEMIA・Lemma の両方でサーバーに保存した学習記録がすべて消えます（この端末の記録は残ります）。元に戻せません。')) return;
    try { await sync.deleteAccount(); toast('退会しました'); } catch (e) { toast(e.message); }
    drawAccount();
  });
}

// ---------- 表示ユーティリティ ----------
const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const rich = (s) => esc(s).replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>');
function math(el) {
  if (!window.renderMathInElement) { setTimeout(() => math(el), 60); return; }
  window.renderMathInElement(el, {
    delimiters: [
      { left: '$$', right: '$$', display: true }, { left: '\\[', right: '\\]', display: true },
      { left: '\\(', right: '\\)', display: false }, { left: '$', right: '$', display: false },
    ],
    throwOnError: false,
  });
}
function tex(el, src, display = false) {
  if (!window.katex) { el.textContent = src; return; }
  window.katex.render(src, el, { throwOnError: false, displayMode: display });
}
function toast(msg) {
  let t = $('.toast');
  if (!t) { t = document.createElement('div'); t.className = 'toast'; document.body.appendChild(t); }
  t.textContent = msg; t.classList.add('show');
  clearTimeout(toast.tm); toast.tm = setTimeout(() => t.classList.remove('show'), 1800);
}
function relTime(t) {
  const d = (t - Date.now()) / 86400000;
  if (d <= 0) return '今';
  if (d < 1 / 24) return `${Math.max(1, Math.round(d * 1440))}分後`;
  if (d < 1) return `${Math.round(d * 24)}時間後`;
  return `${Math.round(d)}日後`;
}

const view = $('#view');
let knowledgeCache = null;   // data/knowledge.json（概念・法則・誤概念）
let activeGraphs = [];   // 画面を移るときに止めるもの（図鑑の操作など。destroy() を持つ）
function clearGraphs() { activeGraphs.forEach(g => g.destroy()); activeGraphs = []; }

// ---------- ルーティング ----------
let session = null; // { steps, i, results: [], extra }
// 中断した演習：端末の中だけに保存（同期はしない）。答え終えた問題の次から再開する
const SKEY = 'alchemia.session.v1';
function saveSession() {
  if (!session) return;
  const cur = session.steps[session.i];
  const answered = cur && session.results.length && session.results[session.results.length - 1].pid === cur.pid;
  const i = session.i + (answered ? 1 : 0);
  try {
    if (i >= session.steps.length) LS.removeItem(SKEY);
    else LS.setItem(SKEY, JSON.stringify({ steps: session.steps, i, results: session.results, cfg: session.cfg, t: Date.now() }));
  } catch { /* 保存できない環境では再開できないだけ */ }
}
function clearSavedSession() { try { LS.removeItem(SKEY); } catch { /* なし */ } }
function savedSession() {
  try {
    const d = JSON.parse(LS.getItem(SKEY) || 'null');
    if (!d || !Array.isArray(d.steps) || !Number.isInteger(d.i) || d.i >= d.steps.length) return null;
    if (Date.now() - d.t > 14 * 86400000) { clearSavedSession(); return null; }   // 2週間たったら捨てる
    const left = d.steps.slice(d.i).filter(st => st.kind !== 'problem' || entryById.has(st.pid));
    return left.some(st => st.kind === 'problem') ? d : null;
  } catch { return null; }
}

let ready = false;
let greeted = false;
let timer = null;
function stopTimer() { if (timer) { clearInterval(timer); timer = null; } }
function route() {
  stopTimer();
  if (!ready) { view.innerHTML = '<p class="muted" style="padding:40px 4px">問題データを読み込んでいます…</p>'; return; }
  let tab = (location.hash || '#today').slice(1);
  let arg = '';
  if (tab.includes('=')) { arg = decodeURIComponent(tab.slice(tab.indexOf('=') + 1)); tab = tab.slice(0, tab.indexOf('=')); }
  // 共有リンク：#set=<問題セットID>、#ids=<問題ID,…>（開いたら「今日」に戻してから出題）
  if (tab === 'set' || tab === 'ids') {
    history.replaceState(null, '', '#today');
    if (tab === 'set') openCollection(arg); else openIds(arg.split(',').map(x => x.trim()).filter(Boolean));
    return;
  }
  const tabOf = { practice: 'today', search: 'today', concept: 'cards', concepts: 'cards', law: 'cards', deck: 'cards', study: 'cards' };
  $$('.tabbar a').forEach(a => a.classList.toggle('on', a.dataset.tab === (tabOf[tab] || tab)));
  document.body.dataset.tab = tab;
  clearGraphs();
  document.body.classList.remove('immersive');
  window.scrollTo(0, 0);
  $('#mascot')?.classList.remove('talk'); // 画面を移ったら吹き出しは閉じる
  if (tab === 'zukan') return renderZukanTab(arg);
  if (tab === 'cards' || tab === 'deck' || tab === 'study') return renderCardsTab(tab, arg);
  if (tab === 'concepts' || tab === 'concept') { if (tab === 'concept') history.replaceState(null, '', '#concepts'); return renderConcepts(); }
  if (tab === 'law') return renderLaw(arg);
  if (tab === 'map') return renderMap();
  if (tab === 'settings') return renderSettings();
  if (tab === 'practice') return renderPractice();
  if (tab === 'search') return renderSearch(arg);
  // 結果画面を見たあとに「今日」へ戻ったら、ホームを出す
  if (session && session.shown && !session.steps[session.i]) session = null;
  if (session) return renderStep();
  return renderToday();
}
window.addEventListener('hashchange', route);

// ---------- 今日 ----------
function levelToolbar(pr, onChange) {
  return {
    html: `<div class="toolbar">
        <div class="toggles" id="lv">${Object.entries(LEVELS).map(([k, v]) => `<button class="tog lv${k} ${pr.levels.includes(+k) ? 'on' : ''}" data-l="${k}">${v}</button>`).join('')}</div>
        <select id="tsize" aria-label="問題数">${[4, 6, 8, 10, 15, 20].map(k => `<option value="${k}" ${k === pr.size ? 'selected' : ''}>${k}問</option>`).join('')}</select>
      </div>`,
    bind() {
      $$('#lv .tog').forEach(b => b.addEventListener('click', () => {
        const k = +b.dataset.l, set = new Set(pr.levels);
        if (set.has(k) && set.size === 1) { toast('難易度は1つ以上選んでください'); return; }
        set.has(k) ? set.delete(k) : set.add(k);
        pr.levels = [...set].sort((x, y) => x - y); save(); onChange();
      }));
      $('#tsize').addEventListener('change', (e) => { pr.size = +e.target.value; save(); onChange(); });
    },
  };
}

// 今日の演習：復習（scheduler。自信なしで正解した問題は早めに戻る）＋苦手な単元＋新しい問題
function todayPlan(pr, weak) {
  return pickToday(state, { levels: pr.levels.slice(), size: pr.size, weak: weak ? weak.f.key : null, weakName: weak ? weak.f.name : '' });
}
function planSummary(plan) {
  const c = plan.count;
  return [c.due ? `復習 ${c.due}問` : '', c.weak ? `苦手な単元 ${c.weak}問` : '', c.fresh ? `新しい問題 ${c.fresh}問` : ''].filter(Boolean).join('・')
    + (plan.short ? `（${plan.short.text}）` : '');
}
function renderToday() {
  const pr = state.practice;
  const solved = (list) => list.filter(e => state.problems[e.id]).length;
  const all = pool({ levels: pr.levels });
  const due = dueEntries(state).length;
  const marked = markedIds();
  const bar = levelToolbar(pr, renderToday);
  const today = dayKey();
  const doneN = state.log.filter(l => dayKey(l.t) === today).length;
  const goal = pr.size;
  const resume = savedSession();
  const weak = weakestField();
  // 今日の演習の中身をここで決め、ボタンの説明と実際の出題を一致させる
  const plan = todayPlan(pr, weak);
  const firstRun = !state.log.length && !resume;   // 初めて：説明と「まず8問」だけを大きく出す
  const planText = planSummary(plan);
  view.innerHTML = `
    <div class="fade-in">
      <div class="kicker">ALCHEMIA · 高校化学</div>
      <h1>今日の学習</h1>
      ${firstRun ? `<div class="panel first-run">
        <p style="margin:0;line-height:1.8"><b>ALCHEMIA は、高校化学の問題を解きながら、反応・性質・法則を「使える知識」にするアプリです。</b>解いた結果から、問題ごとに解き直す時期を記録し、次に解く問題を選びます。覚えることは「暗記」タブのカード、確かめたいことは「図鑑」で。</p>
        <button class="hero" id="firstgo" style="margin-top:12px"><b>まず8問解いてみる</b><span>基礎〜入試基礎の問題を、単元を混ぜて出します。</span></button>
        <p class="muted small" style="margin:8px 0 0">慣れてきたら、下の難易度ボタンで範囲を変えられます。</p>
      </div>` : ''}
      ${firstRun ? '' : `<div class="goal" aria-label="今日の進み具合">
        <span>今日 <b>${doneN}</b> / ${goal}問${doneN >= goal ? '　目標達成' : ''}${state.streak.count ? `　連続 ${state.streak.count}日` : ''}</span>
        <span class="bar"><div style="width:${Math.min(100, Math.round(doneN / Math.max(1, goal) * 100))}%"></div></span>
      </div>`}
      ${resume ? `<div class="resume-row">
        <button class="hero alt" id="resume"><b>前回の続き</b><span>全${resume.steps.length}問の ${resume.i + 1}問目から再開します。</span></button>
        <button class="btn ghost small" id="dropresume">続きを破棄</button>
      </div>` : ''}
      <button class="${firstRun ? 'hero alt' : 'hero'}" id="todaygo" ${plan.steps.length ? '' : 'disabled'}>
        <b>今日の演習を始める</b>
        <span>${planText || '今は出題できる問題がありません'}。単元を混ぜて出すので、どの知識を使うかを見抜く練習になります。</span>
      </button>
      ${weak ? `<p class="weak-note">苦手な単元：<b>${esc(weak.f.name)}</b>（得意度 ${Math.round(weak.score * 100)}・${weak.n}問の結果${weak.n < 10 ? '・暫定' : ''}）<button class="linkbtn" id="weakgo">この単元だけ演習 ›</button></p>` : ''}
      ${bar.html}
      <div class="btn-row" style="margin-top:0">
        <button class="btn" id="random">全分野ランダム</button>
        ${due ? `<button class="btn" id="review">復習だけ ${due}問</button>` : ''}
        ${marked.length ? `<button class="btn" id="bmset">★ ブックマーク ${marked.length}問</button>` : ''}
      </div>
      <p class="sub" style="margin:12px 2px 0">解いた問題 ${solved(all)} / ${all.length}（選んだ難易度）</p>
      <form class="search" id="sform" role="search">
        <input type="search" id="sq" placeholder="問題を探す（例：錯イオン、アンモニアソーダ法、ヨードホルム）" enterkeyhint="search" autocomplete="off">
        <button class="btn" type="submit">検索</button>
      </form>
      <div id="cardsT"></div>
      <div id="sets"></div>
      <h2>形式ごとに演習</h2>
      <div class="fmt-grid">${FORMATS.map(fm => {
        const list = pool({ formats: [fm.id], levels: pr.levels });
        return `<button class="fmt-card" data-fm="${fm.id}" ${list.length ? '' : 'disabled'}><b>${esc(fm.name)}</b><span>${esc(fm.desc)}</span><i>${solved(list)} / ${list.length}</i></button>`;
      }).join('')}</div>
      <h2>単元ごとに演習</h2>
      ${subjectsWithFields().map(sub => {
        const subAll = pool({ fields: sub.fields.map(f => f.key), levels: pr.levels });
        return `<div class="subj-block">
          <button class="subj-row" data-s="${esc(sub.name)}" ${subAll.length ? '' : 'disabled'}>
            <b>${esc(sub.name)}</b><span>${solved(subAll)} / ${subAll.length}</span><i>まとめて ›</i>
          </button>
          ${sub.fields.map(f => {
            const list = pool({ fields: [f.key], levels: pr.levels }), d = solved(list);
            return `<button class="unit-row" data-f="${esc(f.key)}" ${list.length ? '' : 'disabled'}>
              <span class="nm">${esc(f.name)}</span><span class="ct">${d} / ${list.length}</span>
              <span class="bar"><div style="width:${list.length ? Math.round(d / list.length * 100) : 0}%"></div></span>
            </button>`;
          }).join('')}
        </div>`;
      }).join('')}
      <button class="btn ghost block" id="pick" style="margin-top:14px">複数の分野を組み合わせて選ぶ</button>
      <details class="howto">
        <summary>このアプリの学び方</summary>
        <ol style="margin:0;padding-left:1.2em">
          <li><strong>反応を「見える形」で思い浮かべてから解く。</strong>何色の沈殿か、どの結合が切れてどこにつくかを先に言葉にします。</li>
          <li><strong>覚えたことは図鑑で確かめる。</strong>解いたあと「図鑑で確かめる」で、反応表・経路の地図・周期表に戻れます。</li>
          <li><strong>ヒントは答えを出さない。</strong>詰まったら「考え方」だけを見られます。</li>
          <li><strong>答え合わせの前に自信度を選ぶ。</strong>「わかったつもり」とのズレを定着度ページで見られます。</li>
          <li><strong>忘れかけた頃に復習。</strong>解けなかった問題・自信なく正解した問題ほど早く、解けた問題は間をあけて復習に出します。</li>
        </ol>
      </details>
    </div>`;
  bar.bind();
  if (!greeted) { greeted = true; setTimeout(() => mascotSay(mascotMessage('greet'), { sound: false }), 700); }
  const go = (fields, formats = null) => {
    const cfg = { fields, formats, levels: pr.levels.slice(), size: pr.size };
    const st = pickPractice(state, cfg);
    if (st.length) startSession(st, cfg); else toast('この条件の問題はありません');
  };
  $('#random').addEventListener('click', () => go(null));
  $('#todaygo').addEventListener('click', () => {
    if (plan.steps.length) startSession(plan.steps, { today: true, levels: pr.levels.slice(), size: pr.size }); else toast('今は出題できる問題がありません');
  });
  $('#resume')?.addEventListener('click', () => resumeSession(resume));
  $('#firstgo')?.addEventListener('click', () => {
    const steps = pickPractice(state, { levels: [0, 1, 2], size: 8 }).map(s => ({ ...s, why: 'はじめての演習（基礎〜入試基礎）' }));
    if (steps.length) startSession(steps, { levels: [0, 1, 2], size: 8 });
  });
  $('#dropresume')?.addEventListener('click', () => { if (confirm('中断した演習の続きを破棄しますか？（解いた記録は残ります）')) { clearSavedSession(); renderToday(); } });
  $('#weakgo')?.addEventListener('click', () => go([weak.f.key]));
  $$('.fmt-card').forEach(b => b.addEventListener('click', () => go(null, [b.dataset.fm])));
  $$('.unit-row').forEach(b => b.addEventListener('click', () => go([b.dataset.f])));
  $$('.subj-row').forEach(b => b.addEventListener('click', () => go(FIELDS.filter(f => f.subject === b.dataset.s).map(f => f.key))));
  $('#pick').addEventListener('click', () => { location.hash = '#practice'; });
  $('#sform').addEventListener('submit', (ev) => { ev.preventDefault(); const q = $('#sq').value.trim(); if (q) location.hash = '#search=' + encodeURIComponent(q); });
  drawSets();
  cardsTeaser($('#cardsT'));
  $('#review')?.addEventListener('click', () => {
    const st = pickReview(state, { size: pr.size });
    if (st.length) startSession(st, { review: true, size: pr.size });
  });
  // ブックマーク：混ぜて最大30問（交互練習になるように、付けた順ではなく混ぜる）
  $('#bmset')?.addEventListener('click', () => openIds(shuffle(marked).slice(0, 30)));
}

// ---------- 問題セット・共有リンク・検索 ----------
const shareURL = (frag) => location.href.split('#')[0] + '#' + frag;
async function copyLink(url) {
  try { await navigator.clipboard.writeText(url); toast('リンクをコピーしました'); }
  catch { toast(url); }
}
async function drawSets() {
  const box = $('#sets');
  if (!box || !M_CAT.loadCollections) return;
  const cols = await M_CAT.loadCollections();
  if (!cols.length || !document.body.contains(box)) return;
  box.innerHTML = `<h2>問題セット</h2>${cols.map(c => {
    const done = M_CAT.collectionProgress(state, c), n = c.items.length;
    return `<div class="set-card">
      <button class="set-main" data-set="${esc(c.id)}">
        <span class="set-head"><b>${esc(c.title)}</b><span class="ct">${done} / ${n}問</span></span>
        ${c.description ? `<span class="muted set-desc">${esc(c.description)}</span>` : ''}
        <span class="bar"><div style="width:${Math.round(done / n * 100)}%"></div></span>
      </button>
      <div class="set-foot"><span class="set-go">${done >= n ? 'もう一周 ›' : done ? '続きから ›' : '始める ›'}</span><button class="linkbtn" data-share="${esc(c.id)}">リンクをコピー</button></div>
    </div>`;
  }).join('')}`;
  $$('[data-set]', box).forEach(b => b.addEventListener('click', () => openCollection(b.dataset.set)));
  $$('.set-go', box).forEach(x => x.addEventListener('click', () => openCollection(x.closest('.set-card').querySelector('[data-set]').dataset.set)));
  $$('[data-share]', box).forEach(b => b.addEventListener('click', () => copyLink(shareURL('set=' + b.dataset.share))));
}
async function openCollection(id) {
  const col = (await M_CAT.loadCollections()).find(c => c.id === id);
  if (!col) { toast('この問題セットは見つかりません'); route(); return; }
  const steps = M_CAT.pickCollection(state, col, state.practice.size);
  startSession(steps, { collection: id, size: state.practice.size });
}
function openIds(ids) {
  const ok = ids.filter(id => entryById.has(id)).slice(0, 50);
  if (!ok.length) { toast('リンクの問題が見つかりません'); route(); return; }
  startSession(ok.map(pid => ({ kind: 'problem', pid })), { ids: ok, size: ok.length });
}
let searchShown = 30;
async function renderSearch(q) {
  const pr = state.practice;
  const bar = levelToolbar(pr, () => renderSearch(q));
  view.innerHTML = `<div class="fade-in">
    <div class="kicker">問題を探す</div>
    <h1>検索</h1>
    <form class="search" id="sform" role="search">
      <input type="search" id="sq" value="${esc(q)}" placeholder="例：漸化式 確率、回転体、チェビシェフ" enterkeyhint="search" autocomplete="off">
      <button class="btn" type="submit">検索</button>
    </form>
    ${bar.html}
    <div id="sres"><p class="muted">探しています…</p></div>
  </div>`;
  bar.bind();
  $('#sform').addEventListener('submit', (ev) => { ev.preventDefault(); const v = $('#sq').value.trim(); searchShown = 30; location.hash = '#search=' + encodeURIComponent(v); });
  const box = $('#sres');
  if (!q) { box.innerHTML = '<p class="muted">テーマ・タグ・問題文の言葉で探せます。空白で区切ると、すべてを含む問題を探します。</p>'; return; }
  let res;
  try { res = await M_CAT.searchProblems(q, { levels: pr.levels }); }
  catch (e) { box.innerHTML = `<p class="muted">検索の索引を読めませんでした（${esc(e.message)}）。</p>`; return; }
  dataWarn();
  if (!res.length) { box.innerHTML = '<p class="muted">見つかりませんでした。言葉を短くするか、難易度の選択を広げてみてください。</p>'; return; }
  const fresh = res.filter(e => !state.problems[e.id]);
  const pickIds = [...fresh, ...res.filter(e => state.problems[e.id])].slice(0, pr.size).map(e => e.id);
  box.innerHTML = `<p class="sub" style="margin:4px 0 10px">${res.length}問 見つかりました（未着手 ${fresh.length}問）</p>
    <div class="btn-row" style="margin-top:0"><button class="btn primary" id="sgo">${pickIds.length}問 演習する</button><button class="btn ghost" id="sshare">リンクをコピー</button></div>
    <div class="results">${res.slice(0, searchShown).map(e => {
      const f = fieldByKey[e.f] || {};
      return `<button class="res-row" data-pid="${esc(e.id)}">
        <span class="nm">${e.th || esc(e.id)}</span>
        <span class="muted small">${esc(f.subject || '')}・${esc(f.name || '')}　<span class="chip lv${e.lv}">${LEVELS[e.lv]}</span>${state.problems[e.id] ? '　✓ 解いた' : ''}${(e.tg || []).length ? '　' + e.tg.slice(0, 3).map(esc).join('・') : ''}</span>
      </button>`;
    }).join('')}</div>
    ${res.length > searchShown ? `<button class="btn ghost block" id="smore">さらに表示（残り ${res.length - searchShown}問）</button>` : ''}`;
  math(box);
  $('#sgo').addEventListener('click', () => startSession(pickIds.map(pid => ({ kind: 'problem', pid })), { ids: pickIds, size: pickIds.length }));
  $('#sshare').addEventListener('click', () => copyLink(shareURL('search=' + encodeURIComponent(q))));
  $('#smore')?.addEventListener('click', () => { searchShown += 60; renderSearch(q); });
  $$('.res-row', box).forEach(b => b.addEventListener('click', () => openIds([b.dataset.pid])));
}

async function startSession(steps, cfg = null) {
  view.innerHTML = '<p class="muted" style="padding:40px 4px">問題を読み込んでいます…</p>';
  // 本文と、概念・法則の登録簿（knowledge.json）をそろえてから出題する
  try { await Promise.all([ensureBodies(steps.map(s => s.pid).filter(Boolean)), loadKnowledge().then(k => { knowledgeCache = k; }).catch(() => {})]); dataWarn(); }
  catch (e) { showBanner(`<b>問題を読み込めませんでした</b>${esc(e.message)}`); route(); return; }
  session = { steps, i: 0, results: [], cfg };
  saveSession();
  if (!studyStart || Date.now() - lastActivity > 15 * 60000) { studyStart = Date.now(); breakNotified = 0; }
  lastActivity = Date.now();
  const hr = new Date().getHours();
  if (hr >= 0 && hr < 4) setTimeout(() => mascotSay(mascotMessage('late')), 800);
  if (location.hash !== '#today') location.hash = '#today'; else route();
}

async function resumeSession(d) {
  const steps = d.steps.slice(d.i).filter(st => st.kind !== 'problem' || entryById.has(st.pid));
  view.innerHTML = '<p class="muted" style="padding:40px 4px">問題を読み込んでいます…</p>';
  try { await Promise.all([ensureBodies(steps.map(s => s.pid).filter(Boolean)), loadKnowledge().then(k => { knowledgeCache = k; }).catch(() => {})]); dataWarn(); }
  catch (e) { showBanner(`<b>問題を読み込めませんでした</b>${esc(e.message)}`); route(); return; }
  // 進み具合の表示を保つため、解き終えた分も steps に残して i から始める
  session = { steps: [...d.steps.slice(0, d.i), ...steps], i: d.i, results: Array.isArray(d.results) ? d.results : [], cfg: d.cfg || null };
  saveSession();
  studyStart = Date.now(); breakNotified = 0; lastActivity = Date.now();
  if (location.hash !== '#today') location.hash = '#today'; else route();
}

// ---------- セッション ----------
function renderStep() {
  document.body.classList.add('immersive');
  const step = session.steps[session.i];
  if (!step) { session.shown = true; clearSavedSession(); return renderSummary(); }
  return renderProblem(step);
}

function topBar() {
  const done = session.i, total = session.steps.length;
  return `<div class="player-top">
    <button class="icon-btn" id="quit" aria-label="中断">×</button>
    <div class="progress"><div style="width:${Math.round(done / total * 100)}%"></div></div>
  </div>`;
}
function bindQuit() {
  $('#quit').addEventListener('click', () => {
    if (session.results.length && !confirm('中断しますか？ここまでの結果は保存されています。ホームの「前回の続き」から再開できます。')) return;
    saveSession(); session = null; route();
  });
}
function nextStep() { lastActivity = Date.now(); stopTimer(); session.i++; saveSession(); clearGraphs(); window.scrollTo(0, 0); renderStep(); }

const CONF = [{ id: 'low', label: '自信なし' }, { id: 'mid', label: 'たぶん' }, { id: 'high', label: '確実' }];
const FEEL = [
  { k: 1, label: '手が出ない' }, { k: 2, label: '苦戦' }, { k: 3, label: 'ふつう' }, { k: 4, label: 'ほぼ解けた' }, { k: 5, label: '余裕' },
];

// 解説を読みやすく：文の終わりで改行し、∫・Σ・lim を大きく表示する
function bigOps(html) {
  return html.replace(/\$([^$]+)\$/g, (m, body) => /\\(int|iint|sum|prod|lim)\b/.test(body) && !/\\displaystyle/.test(body) ? `$\\displaystyle ${body}$` : m);
}
function formatSolution(html) {
  if (!html) return '';
  const keep = [];
  // 数式とタグを退避してから、文末の「．」「。」の後で改行
  let t = html.replace(/\$[^$]*\$|\\\[[\s\S]*?\\\]|<[^>]+>/g, (m) => { keep.push(m); return `\u0000${keep.length - 1}\u0000`; });
  t = t.replace(/([．。])(?=\s*[^\s<\u0000)）」』．。,，、]|\s*\u0000\d+\u0000\s*[^<\s])/g, '$1<br>');
  t = t.replace(/\u0000(\d+)\u0000/g, (_, k) => keep[+k]);
  // 数式の直後の句読点だけが次の行に送られないようにする
  t = t.replace(/\$([^$]+)\$([．。，、])/g, (m, body, p) => `$${body}\\text{${p}}$`);
  t = t.replace(/(<br>\s*){2,}/g, '<br>').replace(/<br>\s*(<\/?(div|ol|ul|li|table)\b)/g, '$1').replace(/<br>\s*$/, '');
  return bigOps(t);
}

// 大学で習う用語が出てきたら、高校の言葉で説明を添える
// ---------- マスコット
let mascotTimer = null;
function mascotSay(msg, { sound = true, ms = 9000 } = {}) {
  let el = $('#mascot');
  if (!el) return;
  const bub = $('#mascot .bubble');
  bub.innerHTML = `<b>${MASCOT_NAME}</b>${esc(msg).replace(/\n/g, '<br>')}`;
  el.classList.add('talk');
  if (sound && state.settings.sound !== false) play('mascot');
  clearTimeout(mascotTimer);
  mascotTimer = setTimeout(() => el.classList.remove('talk'), ms);
}
function mountMascot() {
  if ($('#mascot') || state.settings.mascot === false) return;
  const el = document.createElement('div');
  el.id = 'mascot';
  el.innerHTML = `<div class="bubble" role="status"></div><button class="mascot-btn" aria-label="${MASCOT_NAME}に話しかける">${MASCOT_SVG}</button>`;
  document.body.appendChild(el);
  $('.mascot-btn', el).addEventListener('click', () => {
    if (el.classList.contains('talk')) { el.classList.remove('talk'); return; }
    mascotSay(mascotMessage('random'));
  });
  $('.bubble', el).addEventListener('click', () => el.classList.remove('talk'));
}
let studyStart = null, breakNotified = 0, lastActivity = 0;
function mascotCheckBreak() {
  if (!studyStart) return;
  const min = Math.floor((Date.now() - studyStart) / 60000);
  const step = min >= 50 ? 2 : min >= 25 ? 1 : 0;
  if (step > breakNotified) { breakNotified = step; mascotSay(mascotMessage('break', { min }), { ms: 12000 }); }
}

// 共通の問題ビュー（catalog.js が古いときは一覧の項目から作る：一部だけ更新された場合でも動くように）
function problemView(e) {
  if (M_CAT.getProblem) return M_CAT.getProblem(e.id);
  const f = fieldByKey[e.f] || {};
  return { id: e.id, kind: e.kind === 'bank' ? 'written' : 'interactive', subject: f.subject, fieldName: f.name,
    level: e.lv, levelName: LEVELS[e.lv], minutes: e.t || null, grading: e.in ? {} : null };
}
// 頻出度（難易度とは別。高・中・低）。問題データに付いているときだけ出す
function valueChip(e) {
  if (!e || !e.v) return '';
  return `<span class="chip val v${e.v}" title="頻出度：入試で繰り返し問われる内容ほど高い">頻出度 ${VALUES[e.v] || ''}</span>`;
}
// この問題を出した理由（今日の演習・類題）。理由のない演習（単元・検索など）では出さない
function whyHTML() {
  const st = session && session.steps[session.i];
  return st && st.why ? `<p class="why">${esc(st.why)}</p>` : '';
}
function metaChips(e) {
  const p = problemView(e);
  return `<span class="chip">${esc(p.subject)}・${esc(p.fieldName)}</span>
        <span class="chip lv${p.level}">${p.levelName}</span>
        ${valueChip(e)}
        ${p.minutes ? `<span class="chip plain">目安 ${p.minutes}分</span>` : ''}
        ${p.kind === 'written' ? `<span class="chip plain">${gradingLabel(e, p)}</span>` : ''}`;
}
// ---------- ブックマーク：{ 問題ID: { on, t } }（外したときも on:false で残し、同期でよみがえらないようにする）
function isMarked(id) { const b = (state.bookmarks || {})[id]; return !!(b && b.on); }
function markedIds() {
  return Object.entries(state.bookmarks || {}).filter(([id, b]) => b && b.on && entryById.has(id)).sort((a, b) => b[1].t - a[1].t).map(([id]) => id);
}
function bmButton(id) {
  const on = isMarked(id);
  return `<button class="chip bm ${on ? 'on' : ''}" data-bm="${esc(id)}" aria-pressed="${on}" title="ブックマーク">${on ? '★' : '☆'}<span> ブックマーク</span></button>`;
}
document.addEventListener('click', (ev) => {
  const b = ev.target.closest && ev.target.closest('[data-bm]');
  if (!b) return;
  const id = b.dataset.bm, on = !isMarked(id);
  state.bookmarks = { ...(state.bookmarks || {}), [id]: { on, t: Date.now() } };
  save(); play('tap');
  b.classList.toggle('on', on); b.setAttribute('aria-pressed', on); b.firstChild.textContent = on ? '★' : '☆';
  toast(on ? 'ブックマークしました' : 'ブックマークを外しました');
});

// ---------- 選択肢つきの確認（confirm の代わり。ボタンの文言で何が起きるかを書けるように）
// choose({ title, body, options: [{ key, label, cls }] }) → 押したボタンの key（閉じたら null）
function choose({ title, body, options }) {
  return new Promise((resolve) => {
    const wrap = document.createElement('div');
    wrap.className = 'chooser';
    wrap.innerHTML = `<div class="chooser-box" role="dialog" aria-modal="true" aria-labelledby="ch-title">
        <h3 id="ch-title">${esc(title)}</h3>
        ${body ? `<p>${body}</p>` : ''}
        <div class="chooser-btns">${options.map(o => `<button class="btn block ${o.cls || ''}" data-k="${esc(o.key)}">${esc(o.label)}</button>`).join('')}
          <button class="btn ghost block" data-k="">キャンセル</button></div>
      </div>`;
    const done = (k) => { wrap.remove(); document.removeEventListener('keydown', onKey); resolve(k || null); };
    const onKey = (e) => { if (e.key === 'Escape') done(null); };
    wrap.addEventListener('click', (e) => {
      if (e.target === wrap) return done(null);
      const btn = e.target.closest('button[data-k]');
      if (btn) done(btn.dataset.k);
    });
    document.addEventListener('keydown', onKey);
    document.body.appendChild(wrap);
    wrap.querySelector('button')?.focus();
  });
}
const hasRecords = (st) => !!(st && ((st.log || []).length || Object.keys(st.pstat || {}).length || Object.keys(st.cards || {}).length));

// 解答の形式（設定）：input＝記述式（答えを入力して自動採点）・choice＝選択式・self＝自己採点（紙に書いて答え合わせ）
const ANSWER_MODES = [['input', '記述式'], ['choice', '選択式'], ['self', '自己採点']];
function answerMode() { return ['input', 'choice', 'self'].includes(state.settings.answer) ? state.settings.answer : 'input'; }
function gradingLabel(e, p) {
  const m = answerMode();
  if (m === 'self' || !p.grading) return '自己採点';
  if (m === 'choice' && e.ch) return '選択式';
  return '入力して採点';
}

function recordResult(pid, grade, extra) {
  // 初めて解いたとき：ブラウザに「この記録を勝手に消さないで」と頼む（対応していない環境では何もしない）
  if (!persistAsked) { persistAsked = true; try { navigator.storage?.persist?.(); } catch { /* なし */ } }
  const now = Date.now();
  state.pstat[pid] = review(state.pstat[pid], grade, now);
  state.problems[pid] = { last: now, n: ((state.problems[pid] || {}).n || 0) + 1 };
  state.log.push({ p: pid, t: now, g: grade, ...extra });
  // 5000件での切り捨てはしない（古い記録が消えないように）。保存できなくなったら save() が知らせる
  const today = dayKey(now);
  if (state.streak.last !== today) {
    const y = dayKey(now - 86400000);
    state.streak = { last: today, count: state.streak.last === y ? state.streak.count + 1 : 1 };
  }
  save();
  session.results.push({ pid, grade });
  saveSession();
  return state.log[state.log.length - 1];
}

function renderProblem(step) {
  const entry = entryById.get(step.pid);
  if (!entry) { nextStep(); return; }
  return renderWritten(step, entry);
}

const SELF = [
  { id: 'full', label: '完答', icon: 'full' },
  { id: 'part', label: '一部正解', icon: 'part' },
  { id: 'none', label: '不正解', icon: 'none' },
];
const CONF_ICON = { low: 'low', mid: 'mid', high: 'high' };
const KEYS = ['^', '_', '(', ')', '/', '√', 'π', 'sin', 'cos', 'tan', 'log', ',', '<', '>', '≦', '≧', '±'];
// 化学式・数値を入れる欄のキー
const CHEM_KEYS = ['(', ')', '[', ']', '+', '-', '^', '·', '=', '#', '×10^', '→'];

// 選択肢・並べかえの札：数式（$...$）はそのまま、化学式だけの札は \ce で
function optHTML(o) {
  if (/\$/.test(o)) return o;
  if (/^[A-Z\[(][A-Za-z0-9()\[\]+\-^·.=#]*$/.test(o)) return `$\\ce{${ceFix(o)}}$`;
  return rich(o);
}
function labelHTML(label) {
  if (!label) return '';
  // 物理量の名前（t_1、v_{max}、θ など）は数式として表示する。$...$ や文章はそのまま
  if (/^[A-Za-zα-ωΑ-Ω](_(\{[A-Za-z0-9]+\}|[A-Za-z0-9]+))?$/.test(label)) {
    const t = exprTeX(label);
    if (t) return `$${t}$`;
  }
  return esc(label);
}

// 解説の下：使う法則（法則のページへ）・確かめる概念・よくある誤解
function conceptBlock(e) {
  const K = knowledgeCache;
  if (!K || !(e.c || e.m)) return '';
  const cs = (e.c || []).map(id => K.concept[id]).filter(Boolean);
  const ms = (e.m || []).map(id => K.misconception[id]).filter(Boolean);
  return `<div class="concept-box">
    ${(getBody(e.id) || {}).laws ? `<div class="methods"><b class="cb-label">使う法則（導出つき）</b>${getBody(e.id).laws.map(id => K.law[id]).filter(Boolean).map(l => `<a class="chip plain" href="#law=${esc(l.id)}">${esc(l.name)}</a>`).join('')}</div>` : ''}
    ${cs.length ? `<div class="methods"><b class="cb-label">確かめる概念</b>${cs.map(c => `<span class="chip plain">${esc(c.name)}</span>`).join('')}</div>` : ''}
    ${ms.length ? `<div class="misc-list"><b class="cb-label">よくある誤解</b>${ms.map(m => `<div class="misc"><span class="wrong">✗ ${esc(m.wrong)}</span><span class="right">→ ${esc(m.right)}</span></div>`).join('')}</div>` : ''}
  </div>`;
}

function renderWritten(step, e) {
  const b = getBody(e.id);
  if (!b) { showBanner('<b>問題本文がありません</b>' + esc(e.id)); nextStep(); return; }
  let spec = null;
  try { spec = b.in ? parseSpec(b.in) : null; } catch (err) { spec = null; }
  const mode = answerMode();
  if (mode === 'self') spec = null;                       // 自己採点：入力用の答えがあっても使わない
  const useChoice = mode === 'choice' && !!spec && Array.isArray(b.ch) && b.ch.some(Boolean);
  const t0 = Date.now();
  const a = { conf: null, hint: false, revealed: false, sg: null, feel: null, log: null, min: null, attempts: 0, results: null, firstAllOk: null, disputes: [] };
  // 入力欄の状態：fields[k] = { f, part, mode: 'plain'|'cases', text, rows }
  const fields = [];
  if (spec) spec.parts.forEach((p, pi) => p.fields.forEach(f => fields.push({ f, part: p.n, pi, mode: 'plain', text: '', rows: [{ cond: '', val: '' }, { cond: '', val: '' }], sel: new Set(), seq: [] })));
  const isText = (f) => !['self', 'pick', 'order'].includes(f.kind);
  // 選択式：選択肢のある欄はボタンで選ぶ（b.ch[k] の先頭が正解。並びは毎回混ぜる）。選ぶ・並べるの欄と、選択肢のない欄はそのまま
  if (useChoice) fields.forEach((x, k) => {
    const opts = b.ch[k];
    if (opts && !['self', 'pick', 'order'].includes(x.f.kind)) { x.mode = 'choice'; x.opts = opts; x.order = shuffle(opts.map((_, i) => i)); x.pick = null; }
  });
  const inputNote = mode === 'choice' && spec && !useChoice && fields.some(x => isText(x.f)) ? '<p class="muted pen">この問題には選択肢がないので、答えを入力します。</p>' : '';
  const inputOf = (x) => x.mode === 'cases' ? { mode: 'cases', rows: x.rows } : { mode: 'plain', text: x.text, sel: x.sel, seq: x.seq };
  // 係数など、答えがただの数の欄（数字キーボードで入れる・場合分けは使わない）
  const isInt = (f) => f.kind === 'expr' && /^-?\d+$/.test(String(f.value.src).trim());
  const autoFields = fields.filter(x => x.f.kind !== 'self');
  const selfParts = spec ? [...new Set(fields.filter(x => x.f.kind === 'self').map(x => x.part))] : [];
  let focused = null;

  view.innerHTML = `
    <div class="fade-in written">
      ${topBar()}
      <div class="meta">
        ${metaChips(e)}
        <span class="chip plain timer" id="clock">0:00</span>${bmButton(e.id)}
        <span class="count">${session.i + 1} / ${session.steps.length}</span>
      </div>
      ${whyHTML()}
      <div class="question qhtml">${bigOps(b.q)}</div>

      ${spec ? inputNote + '<div id="inputs" class="inputs"></div>' : '<p class="muted pen">この問題は紙に書いてから、答え合わせで自己採点します。</p>'}
      <div id="conf"><div class="label">${spec ? '判定の前に' : '答え合わせの前に'}：どのくらい自信がある？</div>
        <div class="segment conf">${CONF.map(c => `<button data-c="${c.id}">${ICON[CONF_ICON[c.id]]}<span>${c.label}</span></button>`).join('')}</div></div>
      <div class="btn-row" id="actions"></div>
      <div class="hint-list" id="hints"></div>
      <div id="feedback"></div>
      <div class="pid muted">${esc(e.id)}</div>
    </div>`;
  math($('.question'));
  mountFigures($('.question'), b.figs, { onDone: (box) => math(box) });
  bindQuit();
  const clock = $('#clock');
  const tick = () => {
    const sec = Math.floor((Date.now() - t0) / 1000);
    clock.textContent = `${Math.floor(sec / 60)}:${String(sec % 60).padStart(2, '0')}`;
    if (e.t && sec > e.t * 60) clock.classList.add('over');
    mascotCheckBreak();
  };
  timer = setInterval(tick, 1000);

  // ---------- 入力欄
  function drawInputs() {
    if (!spec) return;
    const box = $('#inputs');
    const byPart = {};
    fields.forEach((x, k) => { (byPart[x.part] = byPart[x.part] || []).push(k); });
    box.innerHTML = Object.entries(byPart).map(([n, ks]) => `
      <div class="part">
        ${n ? `<div class="pn">(${esc(n)})</div>` : ''}
        <div class="pf">${ks.map(k => fieldHTML(k)).join('')}</div>
      </div>`).join('') + keypadHTML();
    math(box);
    bindInputs();
  }
  function keypadHTML() {
    const txt = fields.filter(x => isText(x.f) && !isInt(x.f) && x.mode !== 'choice');
    const kinds = new Set(txt.map(x => x.f.kind));
    if (!kinds.size) return '';
    const mathy = [...kinds].some(k => !CHEM_KINDS.has(k));
    if (!mathy) return `<div class="keypad sticky" id="kp">${CHEM_KEYS.map(k => `<button type="button" data-k="${esc(k)}">${esc(k)}</button>`).join('')}</div>`;
    return `<div class="keypad sticky" id="kp">${KEYS.map(k => `<button type="button" data-k="${esc(k)}">${esc(k)}</button>`).join('')}${extraKeys()}</div>`;
  }
  function extraKeys() {
    const letters = new Set();
    for (const x of autoFields) for (const ch of (x.f.src || '').match(/[α-ω]/g) || []) letters.add(ch);
    return [...letters].map(k => `<button type="button" data-k="${k}">${k}</button>`).join('') + '<button type="button" data-k="または" class="wide">または</button>';
  }
  function fieldHTML(k) {
    const x = fields[k], f = x.f;
    const res = a.results && a.results[k];
    const mark = res ? (res.ok === true ? `<span class="fmark ok">${ICON.ok}</span>` : res.ok === false ? `<span class="fmark ng">${ICON.ng}</span>` : '') : '';
    const lab = f.label ? `<span class="flabel">${labelHTML(f.label)}</span>` : '';
    if (f.kind === 'self') {
      return `<div class="field self">${lab}<span class="selfnote">${ICON.self} ${esc(f.sub === '証明' ? '証明' : f.sub === '図示' ? '図示' : '記述')}：紙に書いて、答え合わせで自己採点</span></div>`;
    }
    if (x.mode === 'choice') {
      const done = !!res || a.revealed;
      return `<div class="field pick ${res ? (res.ok ? 'is-ok' : 'is-ng') : ''}" data-k="${k}">
        <div class="frow">${lab}<span class="muted small">1つ選ぶ</span>${mark}</div>
        <div class="opts ${x.opts.every(o => o.replace(/\$\\ce\{|\}\$|\$/g, '').length <= 14) ? 'compact' : ''}">${x.order.map(i => {
          const cls = ['opt', 'ch-opt', done ? (i === 0 ? 'right' : i === x.pick ? 'wrong' : '') : (i === x.pick ? 'on' : '')].join(' ');
          return `<button type="button" class="${cls}" data-k="${k}" data-i="${i}" ${done ? 'disabled' : ''}><span class="ob">${i === x.pick ? '●' : '○'}</span><span class="ot qhtml">${optHTML(x.opts[i])}</span></button>`;
        }).join('')}</div>
      </div>`;
    }
    const lockedC = a.revealed || (res && res.ok === true);
    const errC = res && res.error ? `<div class="inline-err">${esc(res.error)}</div>` : '';
    if (f.kind === 'pick') {
      const show = a.revealed;
      return `<div class="field pick ${res ? (res.ok ? 'is-ok' : 'is-ng') : ''}" data-k="${k}">
        <div class="frow">${lab}<span class="muted small">${f.multi ? '当てはまるものをすべて選ぶ' : '1つ選ぶ'}</span>${mark}</div>
        <div class="opts ${f.options.every(o => o.replace(/\$\\ce\{|\}\$|\$/g, '').length <= 14) ? 'compact' : ''}">${f.options.map((o, i) => {
          const on = x.sel.has(i), right = f.answer.includes(i);
          const cls = ['opt', on ? 'on' : '', show && right ? 'right' : '', show && on && !right ? 'wrong' : ''].join(' ');
          return `<button type="button" class="${cls}" data-pick="${k}" data-i="${i}" ${lockedC ? 'disabled' : ''}><span class="ob">${f.multi ? (on ? '☑' : '☐') : (on ? '●' : '○')}</span><span class="ot qhtml">${optHTML(o)}</span></button>`;
        }).join('')}</div>${errC}</div>`;
    }
    if (f.kind === 'order') {
      const rest = f.options.map((o, i) => i).filter(i => !x.seq.includes(i));
      return `<div class="field order ${res ? (res.ok ? 'is-ok' : 'is-ng') : ''}" data-k="${k}">
        <div class="frow">${lab}<span class="muted small">順にタップして並べる（並べた札をタップすると戻す）</span>${mark}</div>
        <div class="seq">${x.seq.length ? x.seq.map((i, n) => `<button type="button" class="chipbtn placed" data-unseq="${k}" data-i="${i}" ${lockedC ? 'disabled' : ''}><b>${n + 1}</b><span class="qhtml">${optHTML(f.options[i])}</span></button>`).join('<span class="gt">＞</span>') : '<span class="muted small">ここに並びます</span>'}</div>
        ${rest.length && !lockedC ? `<div class="opts pool">${rest.map(i => `<button type="button" class="chipbtn" data-seq="${k}" data-i="${i}"><span class="qhtml">${optHTML(f.options[i])}</span></button>`).join('')}</div>` : ''}
        ${x.seq.length && !lockedC ? `<button class="mini" data-act="reset" data-k="${k}">並べ直す</button>` : ''}
        ${errC}</div>`;
    }
    const canCases = false;   // 化学では場合分けの入力を使わない（書式としては読める）
    const unit = fieldUnit(f);
    const locked = a.revealed || (res && res.ok === true);
    if (x.mode === 'cases') {
      return `<div class="field cases ${res ? (res.ok ? 'is-ok' : 'is-ng') : ''}" data-k="${k}">
        <div class="frow">${lab}<span class="muted">場合分け</span>${mark}${canCases && !locked ? `<button class="mini" data-act="plain" data-k="${k}">場合分けをやめる</button>` : ''}</div>
        ${x.rows.map((r, i) => `<div class="caserow">
          <input class="answer-input cond" data-k="${k}" data-row="${i}" data-part="cond" value="${esc(r.cond)}" placeholder="条件　例：a<0" ${locked ? 'readonly' : ''} autocomplete="off" autocorrect="off" autocapitalize="off" spellcheck="false">
          <span class="arrow">→</span>
          <input class="answer-input val" data-k="${k}" data-row="${i}" data-part="val" value="${esc(r.val)}" placeholder="値" ${locked ? 'readonly' : ''} autocomplete="off" autocorrect="off" autocapitalize="off" spellcheck="false">
          ${!locked && x.rows.length > 1 ? `<button class="mini del" data-act="delrow" data-k="${k}" data-row="${i}" aria-label="この場合を消す">×</button>` : ''}
        </div>`).join('')}
        ${!locked ? `<button class="mini add" data-act="addrow" data-k="${k}">＋ 場合を追加</button>` : ''}
        ${res && res.error ? `<div class="inline-err">${esc(res.error)}</div>` : ''}
      </div>`;
    }
    const short = isInt(f), ph = short ? '' : f.kind === 'word' ? '語句・化学式' : f.kind === 'approx' ? '数値（例 2.4、1.2×10^-3）' : (fieldHint(f) || (f.kind === 'expr' ? '答えを入力' : ''));
    return `<div class="field ${short ? 'short' : ''} ${res ? (res.ok ? 'is-ok' : 'is-ng') : ''}" data-k="${k}">
      <div class="answer-row">
        ${lab}
        <input class="answer-input" data-k="${k}" value="${esc(x.text)}" placeholder="${esc(ph)}" ${locked ? 'readonly' : ''} autocomplete="off" autocorrect="off" autocapitalize="off" spellcheck="false" inputmode="${short ? 'numeric' : 'text'}">
        ${unit ? `<span class="unit">${esc(unit)}</span>` : ''}
        ${mark}
      </div>
      <div class="preview" data-pv="${k}"></div>
      ${canCases && !locked ? `<button class="mini" data-act="cases" data-k="${k}">${ICON.cases} 場合分けする</button>` : ''}
      ${res && res.error ? `<div class="inline-err">${esc(res.error)}</div>` : ''}
    </div>`;
  }
  function preview(k) {
    const x = fields[k], pv = $(`[data-pv="${k}"]`);
    if (!pv) return;
    const v = x.text.trim();
    if (!v || /^-?\d+(\.\d+)?$/.test(v)) { pv.textContent = ''; return; }
    if (x.f.kind === 'word') {
      // 化学式らしいものだけ下付き・上付きで見せる（日本語の語句はそのまま）
      const w = v.normalize('NFKC').replace(/\s+/g, '');
      if (/^[A-Z\[(][A-Za-z0-9()\[\]+\-^·.=#{}]*$/.test(w)) { tex(pv, `\\ce{${ceFix(w)}}`); pv.className = 'preview'; } else pv.textContent = '';
      return;
    }
    if (x.f.kind === 'approx') { pv.textContent = ''; return; }
    try {
      const t = x.f.kind === 'range' ? condTeX(v) : x.f.kind === 'set' ? v.split(/[,，、]/).map(s => exprTeX(s.trim()) ?? s).join(',\\ ') : (exprTeX(v.replace(/^[^=<>]*=(?!=)/, '')) ?? null);
      if (t == null) throw new Error();
      tex(pv, t); pv.className = 'preview';
    } catch { pv.textContent = '…'; pv.className = 'preview err'; }
  }
  function bindInputs() {
    $$('#inputs .ch-opt').forEach(btn => btn.addEventListener('click', () => {
      const x = fields[+btn.dataset.k];
      if (a.results || a.revealed) return;
      x.pick = +btn.dataset.i;
      play('tap'); drawInputs(); actions();
    }));
    $$('#inputs input').forEach(inp => {
      const k = +inp.dataset.k, x = fields[k];
      inp.addEventListener('focus', () => { focused = inp; });
      inp.addEventListener('input', () => {
        if (inp.dataset.row != null) x.rows[+inp.dataset.row][inp.dataset.part] = inp.value;
        else { x.text = inp.value; preview(k); }
        if (a.results && a.results[k]) { a.results[k] = null; }
        actions();
      });
      inp.addEventListener('keydown', (ev) => { if (ev.key === 'Enter') { ev.preventDefault(); $('#judge')?.click(); } });
      if (inp.dataset.row == null) preview(k);
    });
    $$('#inputs [data-pick]').forEach(btn => btn.addEventListener('click', () => {
      const k = +btn.dataset.pick, x = fields[k], i = +btn.dataset.i;
      if (x.f.multi) { x.sel.has(i) ? x.sel.delete(i) : x.sel.add(i); } else x.sel = new Set(x.sel.has(i) ? [] : [i]);
      if (a.results) a.results[k] = null;
      play('tap'); drawInputs(); actions();
    }));
    $$('#inputs [data-seq]').forEach(btn => btn.addEventListener('click', () => {
      const k = +btn.dataset.seq, x = fields[k];
      x.seq.push(+btn.dataset.i); if (a.results) a.results[k] = null;
      play('tap'); drawInputs(); actions();
    }));
    $$('#inputs [data-unseq]').forEach(btn => btn.addEventListener('click', () => {
      const k = +btn.dataset.unseq, x = fields[k];
      x.seq = x.seq.filter(i => i !== +btn.dataset.i); if (a.results) a.results[k] = null;
      drawInputs(); actions();
    }));
    $$('#inputs [data-act]').forEach(btn => btn.addEventListener('click', () => {
      const k = +btn.dataset.k, x = fields[k];
      if (btn.dataset.act === 'reset') { x.seq = []; if (a.results) a.results[k] = null; drawInputs(); actions(); return; }
      if (btn.dataset.act === 'cases') { x.mode = 'cases'; if (x.text) x.rows[0].val = x.text; }
      if (btn.dataset.act === 'plain') x.mode = 'plain';
      if (btn.dataset.act === 'addrow') x.rows.push({ cond: '', val: '' });
      if (btn.dataset.act === 'delrow') x.rows.splice(+btn.dataset.row, 1);
      if (a.results) a.results[k] = null;
      drawInputs(); actions();
      const target = $(`#inputs [data-k="${k}"] input:not([readonly])`); target && target.focus();
    }));
    $$('#kp button').forEach(btn => btn.addEventListener('mousedown', (ev) => ev.preventDefault()));
    $$('#kp button').forEach(btn => btn.addEventListener('click', () => {
      const inp = focused && document.body.contains(focused) ? focused : $('#inputs input:not([readonly])');
      if (!inp || inp.readOnly) return;
      const key = btn.dataset.k;
      const ins = { 'sin': 'sin(', 'cos': 'cos(', 'tan': 'tan(', 'log': 'log(', '√': 'sqrt(', '≦': '<=', '≧': '>=', 'または': ' または ', '×10^': '×10^' }[key] ?? key;
      const st = inp.selectionStart ?? inp.value.length, en = inp.selectionEnd ?? st;
      inp.value = inp.value.slice(0, st) + ins + inp.value.slice(en);
      inp.focus(); inp.setSelectionRange(st + ins.length, st + ins.length);
      inp.dispatchEvent(new Event('input'));
    }));
  }

  $$('#conf button').forEach(btn => btn.addEventListener('click', () => {
    if (a.revealed || a.attempts) return;
    a.conf = btn.dataset.c;
    $$('#conf button').forEach(x => x.classList.toggle('on', x === btn));
    play('tap');
    actions();
  }));

  function filled() {
    return autoFields.every(x => x.mode === 'choice' ? x.pick !== null : x.f.kind === 'pick' ? x.sel.size > 0 : x.f.kind === 'order' ? x.seq.length === x.f.options.length
      : x.mode === 'cases' ? x.rows.some(r => r.cond.trim() && r.val.trim()) : x.text.trim());
  }
  function actions() {
    const box = $('#actions');
    if (a.revealed) { box.innerHTML = ''; return; }
    const hintBtn = b.k && !a.hint ? `<button class="btn ghost" id="hint">${ICON.hint} ヒント：考え方</button>` : '';
    if (!spec) {
      box.innerHTML = `${hintBtn}<button class="btn primary" id="reveal" ${a.conf ? '' : 'disabled'}>答え合わせ</button>`;
      $('#reveal').addEventListener('click', () => reveal());
    } else {
      const ready = a.conf && filled();
      box.innerHTML = `${hintBtn}<button class="btn primary" id="judge" ${ready ? '' : 'disabled'}>${a.attempts ? 'もう一度判定' : '判定'}</button>
        ${a.attempts ? '<button class="btn" id="giveup">答えを見る</button>' : ''}`;
      $('#judge').addEventListener('click', judge);
      $('#giveup')?.addEventListener('click', () => reveal());
    }
    $('#hint')?.addEventListener('click', () => {
      a.hint = true;
      const h = document.createElement('div');
      h.className = 'hint fade-in';
      h.innerHTML = `<b>考え方</b>${b.k}`;
      $('#hints').appendChild(h);
      math(h);
      actions();
    });
  }

  function checkOne(x) {
    if (x.f.kind === 'self') return { ok: null };
    if (x.mode === 'choice') return { ok: x.pick === 0, choice: true };
    return checkField(x.f, inputOf(x));
  }
  function judge() {
    const res = fields.map(checkOne);
    // 書き方の誤りだけなら、試行に数えずに知らせる
    if (res.some(r => r.syntax)) {
      a.results = res.map(r => r.syntax ? r : null);
      drawInputs(); actions(); play('wrong');
      return;
    }
    a.attempts++;
    a.results = res;
    const auto = res.filter(r => r.ok !== null);
    const allOk = auto.every(r => r.ok);
    if (a.firstAllOk === null) a.firstAllOk = allOk;
    drawInputs();
    if (allOk) { play('correct'); reveal(); return; }
    play('wrong');
    if (useChoice) { reveal(); return; }                    // 選択式はやり直しなし（当てずっぽうで正解にならないように）
    const fb = $('#feedback');
    const nOk = auto.filter(r => r.ok).length;
    fb.innerHTML = `<div class="result ng fade-in"><span class="mark">${ICON.ng}</span><div>${nOk ? `${auto.length}つのうち${nOk}つ正解` : '違うようです'}<small>×の欄を見直して、もう一度判定できます。</small></div></div>`;
    actions();
  }

  function reveal() {
    a.revealed = true;
    a.min = Math.max(1, Math.round((Date.now() - t0) / 60000));
    stopTimer();
    $$('#conf button').forEach(x => x.disabled = true);
    if (spec) {
      // 判定していない欄は、いま判定する（答えを見る前の最終状態）
      a.results = fields.map((x, k) => (a.results && a.results[k]) || checkOne(x));
      drawInputs();
    }
    actions();
    const fb = $('#feedback');
    const autoRes = spec ? a.results.filter(r => r.ok !== null) : [];
    const autoAll = autoRes.length && autoRes.every(r => r.ok);
    const answerKey = spec ? specAnswerHTML() : `<div class="ans qhtml">${b.a}</div>`;
    const needSelf = !spec || selfParts.length > 0;   // 自己採点がある＝先に解答・解説を見せてから聞く
    fb.innerHTML = `
      ${spec && autoRes.length ? `<div id="autores">${autoBanner()}</div>` : ''}
      <div class="solution fade-in">
        <div class="title">答</div>
        ${answerKey}
        ${spec ? `<details class="orig"><summary>答えの表記（問題データ）</summary><div class="qhtml">${b.a}</div></details>` : ''}
      </div>
      ${spec && autoRes.some(r => !r.ok && !r.choice) ? '<p class="muted dispute-note">自分の答えも正しいはずなのに×のときは、欄の横の「判定がおかしい？」を押してください。記録に残り、テスト記録のCSVにも出ます。</p>' : ''}
      ${needSelf ? solutionHTML() : ''}
      ${needSelf ? `<div class="label">${spec ? `解答と見比べて、自己採点の小問（${selfParts.map(n => n ? `(${n})` : '').join('')}）はどうだった？` : '解答と見比べて、自分の答案はどうだった？'}</div>
        <div class="segment self" id="self">${SELF.map(x => `<button data-g="${x.id}">${ICON[x.icon]}<span>${x.label}</span></button>`).join('')}</div>` : ''}
      <div id="after"></div>`;
    math(fb);
    mountFigures(fb, b.figs, { allowAnswer: true, onDone: (box) => math(box) });
    if ($('#solbox')) zukanButton($('#solbox'), b.zk);
    addDisputeButtons();
    $$('#self button').forEach(btn => btn.addEventListener('click', () => grade(btn.dataset.g)));
    if (spec && !selfParts.length) grade(null);
  }

  function solutionHTML() {
    return `<div class="solution fade-in" id="solbox">
        <div class="title">解答・解説</div>
        <div class="qhtml sol">${formatSolution(b.s) || '（解説なし）'}</div>
        ${b.k ? `<div class="keypoint"><b>核となる考え方</b>${b.k}</div>` : ''}
        ${conceptBlock(e)}
        ${e.th || e.tg ? `<div class="methods">${e.th ? `<span class="chip plain">${e.th}</span>` : ''}${(e.tg || []).map(t => `<span class="chip plain">${esc(t)}</span>`).join('')}</div>` : ''}

      </div>`;
  }

  function specAnswerHTML() {
    return `<div class="anskey">${spec.parts.map(p => `<div class="akrow">${p.n ? `<span class="pn">(${esc(p.n)})</span>` : ''}<span>${p.fields.map(f => {
      if (f.kind === 'self') return `<span class="muted">${esc(f.sub)}${f.text ? '：' + esc(f.text) : '（解説を参照）'}</span>`;
      const lab = f.label ? `${labelHTML(f.label)} ＝ ` : '';
      const h = fieldAnswerHTML(f);
      if (h != null) return `<span class="akf">${f.label ? `${labelHTML(f.label)}：` : ''}${h}</span>`;
      return `<span class="akf">${lab}$\\displaystyle ${esc(fieldTeX(f))}$</span>`;
    }).join('')}</span></div>`).join('')}</div>`;
  }

  function addDisputeButtons() {
    if (!spec || !a.results) return;
    fields.forEach((x, k) => {
      const r = a.results[k];
      if (!r || r.ok !== false || x.mode === 'choice') return;
      const el = $(`#inputs .field[data-k="${k}"]`);
      if (!el || el.querySelector('.dispute')) return;
      const btn = document.createElement('button');
      btn.className = 'mini dispute';
      btn.textContent = '判定がおかしい？';
      btn.addEventListener('click', () => {
        a.results[k] = { ok: true, disputed: true };
        const said = x.f.kind === 'pick' ? [...x.sel].map(i => x.f.options[i]).join(',') : x.f.kind === 'order' ? x.seq.map(i => x.f.options[i]).join('>')
          : x.mode === 'cases' ? x.rows.map(r2 => `${r2.cond}→${r2.val}`).join(' | ') : x.text;
        a.disputes.push(`(${x.part || '-'})${x.f.label ? ' ' + x.f.label : ''}：${said}`);
        btn.replaceWith(Object.assign(document.createElement('span'), { className: 'muted small', textContent: '正解として記録しました' }));
        if (a.log) { a.log.disp = a.disputes.slice(); save(); }
        el.classList.remove('is-ng'); el.classList.add('is-ok');
        const m = el.querySelector('.fmark'); if (m) { m.className = 'fmark ok'; m.innerHTML = ICON.ok; }
        if (a.sg) regrade();
      });
      el.appendChild(btn);
    });
  }

  function autoScore() {
    const auto = (a.results || []).filter(r => r && r.ok !== null);
    return { n: auto.length, ok: auto.filter(r => r.ok).length };
  }
  function combine(selfSg) {
    const { n, ok } = autoScore();
    if (!spec) return selfSg;
    const autoSg = !n ? null : ok === n ? 'full' : ok === 0 ? 'none' : 'part';
    if (!selfSg) return autoSg;
    if (!autoSg) return selfSg;
    if (autoSg === selfSg) return autoSg;
    return (autoSg === 'none' && selfSg === 'none') ? 'none' : 'part';
  }
  let selfChoice = null;
  function autoBanner() {
    const autoRes = (a.results || []).filter(r => r && r.ok !== null);
    const autoAll = autoRes.every(r => r.ok);
    return `<div class="result ${autoAll ? 'ok' : 'ng'} fade-in"><span class="mark">${autoAll ? ICON.ok : ICON.ng}</span><div>${autoAll ? (a.attempts > 1 ? '正解（2回目以降）' : '正解') : `入力した答え：${autoRes.filter(r => r.ok).length} / ${autoRes.length} 正解`}${a.hint ? '<small>ヒントあり</small>' : ''}</div></div>`;
  }
  // 「判定がおかしい？」で採点が変わったとき：記録と表示を直す（復習の予定は、この回の前の状態から計算し直す）
  function regrade() {
    const sg = combine(selfChoice);
    const ar = $('#autores'); if (ar) ar.innerHTML = autoBanner();
    if (sg === a.sg) return;
    a.sg = sg;
    const v = $('.verdict');
    if (v) { v.className = `verdict v-${sg}`; v.innerHTML = `${ICON[sg]}<b>${SELF.find(x => x.id === sg).label}</b>`; }
    if (a.log) {
      a.log.sg = sg; a.log.ok = sg === 'full';
      a.log.g = gradeOf(sg);
      state.pstat[e.id] = review(a.prevStat, a.log.g, a.log.t);
      const r = [...session.results].reverse().find(x => x.pid === e.id);
      if (r) r.grade = a.log.g;
      save();
    }
  }
  function gradeOf(sg) {
    // 選択式で正解したときは、思い出して書けたときより定着が弱いとみて「確実」でも 3 まで。
    // 自信なしで正解したときは「ヒントあり」と同じ 2（早めに復習）
    return sg === 'full' ? (a.hint || a.attempts > 1 || a.conf === 'low' ? 2 : a.conf === 'high' && !useChoice ? 4 : 3) : sg === 'part' ? 2 : 1;
  }

  function grade(selfSg) {
    if (a.sg) return;
    selfChoice = selfSg;
    const sg = combine(selfSg);
    a.sg = sg;
    a.prevStat = state.pstat[e.id];
    $$('#self button').forEach(x => { x.disabled = true; x.classList.toggle('on', x.dataset.g === selfSg); });
    play(sg === 'full' ? 'full' : sg === 'part' ? 'part' : 'wrong');
    a.log = recordResult(e.id, gradeOf(sg), { h: a.hint ? 1 : 0, c: a.conf, ok: sg === 'full', sg, min: a.min, in: spec ? 1 : 0, md: useChoice ? 'choice' : spec ? 'input' : 'self', disp: a.disputes.length ? a.disputes.slice() : undefined });
    const after = $('#after');
    after.innerHTML = `
      <div class="verdict v-${sg} fade-in">${ICON[sg]}<b>${SELF.find(x => x.id === sg).label}</b></div>
      ${confNoteHTML(sg, a.conf, $('#solbox') ? '' : b.k)}
      ${$('#solbox') ? '' : solutionHTML()}
      <div class="label">手応え（任意）</div>
      <div class="feel" id="feel">${FEEL.map(f => `<button data-k="${f.k}"><b>${f.k}</b><span>${f.label}</span></button>`).join('')}</div>
      <div class="label">${sg === 'full' ? 'この問題で最初に気づくべきことを1行で（任意）' : 'どこで詰まった？（任意）'}</div>
      <textarea class="explain" id="explain" placeholder="${sg === 'full' ? '例：Fe³⁺ は KSCN で血赤色、だから Fe³⁺ と決まる' : '例：酸化数の変化を数え間違えた'}"></textarea>
      <div class="btn-row">
        ${similarTo(state, e.id, session.steps.map(st => st.pid)) ? '<button class="btn" id="similar">類題を1問追加</button>' : ''}
        <button class="btn primary" id="next">${session.i + 1 < session.steps.length ? '次へ' : '結果を見る'}</button>
      </div>`;
    math(after);
    $('#similar')?.addEventListener('click', (ev) => addSimilar(e.id, ev.currentTarget));
    mountFigures(after, b.figs, { allowAnswer: true, onDone: (box) => math(box) });
    if ($('#solbox', after)) zukanButton($('#solbox', after), b.zk);
    $$('#feel button').forEach(btn => btn.addEventListener('click', () => {
      a.feel = +btn.dataset.k;
      $$('#feel button').forEach(x => x.classList.toggle('on', x === btn));
    }));
    $('#next').addEventListener('click', () => {
      const t = $('#explain').value.trim();
      if (a.log) { if (t) a.log.note = t; if (a.feel) a.log.feel = a.feel; if (a.disputes.length) a.log.disp = a.disputes.slice(); }
      if (t) state.explains[e.id] = [...(state.explains[e.id] || []), { t: Date.now(), text: t }].slice(-5);
      save();
      nextStep();
    });
  }

  drawInputs();
  actions();
}

// 自信度と結果のずれへの一言：正解だが自信なし → 早めに復習／確実なのに解けず → 思い込みの確認
function confNoteHTML(sg, conf, key) {
  if (sg === 'full' && conf === 'low') return `<p class="conf-note fade-in">自信がないまま正解しました。定着しきっていないので、早めに復習に出します。</p>`;
  if (sg !== 'full' && conf === 'high') return `<div class="conf-note warn fade-in"><b>「確実」だったのに解けませんでした。</b>覚え違いや見落としがなかったか、解説の「核となる考え方」と照らし合わせてください。復習でも優先して出します。${key ? `<div class="keypoint">${key}</div>` : ''}</div>`;
  return '';
}

// 類題を今の問題の次に差し込む
async function addSimilar(pid, btn) {
  const sim = similarTo(state, pid, session.steps.map(st => st.pid));
  if (!sim) { btn.disabled = true; btn.textContent = '類題が見つかりません'; return; }
  btn.disabled = true; btn.textContent = '読み込んでいます…';
  try { await ensureBodies([sim.id]); } catch (err) { btn.disabled = false; btn.textContent = '類題を1問追加'; toast('類題を読み込めませんでした'); return; }
  session.steps.splice(session.i + 1, 0, { kind: 'problem', pid: sim.id, why: '類題（同じ単元・共通のタグ）' });
  saveSession();
  btn.textContent = `類題を追加しました（${LEVELS[sim.lv]}）`;
  const nx = $('#next'); if (nx) nx.textContent = '次へ（類題）';
}

function renderSummary() {
  const r = session.results;
  const solo = r.filter(x => x.grade >= 3).length;
  const withHint = r.filter(x => x.grade === 2).length;
  const missed = r.filter(x => x.grade === 1).length;
  const next = nextReviewTime(state);
  const today = dayKey();
  const doneN = state.log.filter(l => dayKey(l.t) === today).length, goal = state.practice.size;
  const first = state.log.length === r.length;
  view.innerHTML = `
    <div class="fade-in">
      <div class="kicker">おつかれさま</div>
      <h1>${r.length}問 解きました</h1>
      <div class="stat-row">
        <div class="stat"><b>${solo}</b><span>自力で完答</span></div>
        <div class="stat"><b>${withHint}</b><span>一部・ヒントあり</span></div>
        <div class="stat"><b>${missed}</b><span>解けず</span></div>
      </div>
      <div class="panel">
        <p style="margin:0">解けなかった問題・自信なく正解した問題は近いうちに、解けた問題は間をあけて復習に出します。${next ? `次の復習は <strong>${relTime(next)}</strong>。` : ''}</p>
        <p class="muted" style="margin:8px 0 0">難しく感じた日の方が、あとで解ける力は伸びています。</p>
      </div>
      ${first ? `<div class="panel"><p style="margin:0"><b>記録の使われ方</b></p><p class="small" style="margin:6px 0 0;line-height:1.8">いま解いた結果から、問題ごとに次に解き直す時期を決めました。明日からの「今日の演習」は、解き直しの時期が来た問題・苦手な単元・新しい問題を組み合わせて出します。様子は<a href="#map">定着度</a>で見られます。</p></div>` : ''}
      <p class="sub" style="margin:14px 2px 0">今日 ${doneN} / ${goal}問${doneN >= goal ? '　今日の目標を達成しました。' : `　あと ${goal - doneN}問で今日の目標です。`}</p>
      <div class="btn-row">
        <button class="btn${doneN < goal ? ' primary' : ''}" id="more">${session.cfg && session.cfg.today ? (doneN < goal ? '続けて今日の演習' : 'もう一度今日の演習') : '同じ条件でもう一度'}</button>
        <button class="btn${doneN >= goal ? ' primary' : ''}" id="home">ホームへ</button>
      </div>
    </div>`;
  math(view);
  play('finish');
  setTimeout(() => mascotSay(mascotMessage('summary', { solo, missed, n: r.length }), { sound: false }), 600);
  $('#home').addEventListener('click', () => { session = null; route(); });
  $('#more').addEventListener('click', () => {
    const cfg = session.cfg || { size: 8 };
    if (cfg.collection) { openCollection(cfg.collection); return; }
    if (cfg.ids) { openIds(cfg.ids); return; }
    const steps = cfg.today ? todayPlan({ levels: cfg.levels, size: cfg.size }, weakestField()).steps : cfg.review ? pickReview(state, cfg) : pickPractice(state, cfg);
    if (!steps.length) { toast('今は出題できる問題がありません'); return; }
    startSession(steps, cfg);
  });
}

function shuffle(a) {
  for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(Math.random() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; }
  return a;
}

// ---------- 図鑑（js/zukan.js）----------
// 解いたあとに「図鑑で確かめる」：問題が指す図鑑のページ（イオンの反応表・有機化合物の地図など）へ
const ZUKAN_NAME = { ions: '沈殿・錯イオン・呈色', reactions: '製法・反応式', organic: '有機化合物の反応経路', tests: '検出反応', polymers: '高分子', elements: '周期表・電子配置', laws: '定義・法則' };
const ICON_ZUKAN = '<svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M9 3h6M10 3v6l-5 9a2 2 0 0 0 1.8 3h10.4a2 2 0 0 0 1.8-3l-5-9V3"/><path d="M7.5 14h9"/></svg>';
function zukanButton(host, zk) {
  if (!host || !zk || host.querySelector('.zk-open')) return;
  const a = document.createElement('a');
  a.className = 'btn ghost block zk-open';
  a.href = zk.page === 'laws' ? (zk.key ? `#law=${encodeURIComponent(zk.key)}` : '#concepts') : `#zukan=${zk.page}${zk.key ? '/' + encodeURIComponent(zk.key) : ''}`;
  a.innerHTML = `${ICON_ZUKAN} 図鑑で確かめる：${esc(zk.label || ZUKAN_NAME[zk.page] || '')}`;
  a.addEventListener('click', () => { session = null; });
  host.appendChild(a);
}
async function renderZukanTab(arg) {
  view.innerHTML = '<p class="muted" style="padding:40px 4px">図鑑を読み込んでいます…</p>';
  try {
    const [M, Z] = await Promise.all([import('./zukan.js'), loadZukan()]);
    if (M.FILE_VERSION !== VERSION) showBanner(`<b>古いファイルが混ざっています（v${VERSION}）</b>js/zukan.js が更新されていません。`);
    dataWarn();
    const g = M.renderZukan(view, Z, arg, {
      esc, math, toast, prefs: state.zukanPrefs || (state.zukanPrefs = {}), save,
      problemsFor: (page, key) => ENTRIES.filter(e => e.zk && (key ? e.zk === `${page}/${key}` : e.zk.split('/')[0] === page)),
      openIds, entryById, LEVELS,
    });
    if (g) activeGraphs.push(g);
  } catch (e) { view.innerHTML = `<p class="muted">図鑑を読めませんでした（${esc(e.message)}）。</p>`; }
}

// ---------- 暗記カード（js/cards.js）----------
let cardsData = null;
let cardIdSet = null;   // 記録の点検用（暗記カードのID）
async function loadCards() {
  if (!cardsData) cardsData = fetch('data/cards.json').then(r => { if (!r.ok) throw new Error(`data/cards.json を読めません（${r.status}）`); return r.json(); })
    .then(d => { cardIdSet = new Set(d.decks.flatMap(dk => dk.cards.map(c => c.id))); return d; }).catch(err => { cardsData = null; throw err; });
  return cardsData;
}
async function renderCardsTab(tab, arg) {
  view.innerHTML = '<p class="muted" style="padding:40px 4px">カードを読み込んでいます…</p>';
  try {
    const [M, D] = await Promise.all([import('./cards.js'), loadCards()]);
    if (M.FILE_VERSION !== VERSION) showBanner(`<b>古いファイルが混ざっています（v${VERSION}）</b>js/cards.js が更新されていません。`);
    if (!state.cards) state.cards = {};
    M.renderCards(view, D, tab, arg, { esc, math, state, save, play, prefs: state.cardPrefs || (state.cardPrefs = {}) });
  } catch (e) { view.innerHTML = `<p class="muted">暗記カードを読めませんでした（${esc(e.message)}）。</p>`; }
}
// 今日の画面に、暗記カードの今日の復習の枚数を出す
async function cardsTeaser(host) {
  if (!host) return;
  try {
    const [M, D] = await Promise.all([import('./cards.js'), loadCards()]);
    const n = M.dueCount(D, state), total = D.decks.reduce((a, d) => a + d.cards.length, 0);
    const learned = Object.values(state.cards || {}).filter(r => r.box >= 3).length;
    host.innerHTML = `<a class="cards-teaser" href="${n ? '#study=all' : '#cards'}"><b>暗記カード</b><span>${n ? `今日の復習 ${n}枚` : `覚えた ${learned} / ${total}枚`}</span><i>›</i></a>`;
  } catch { host.remove(); }
}

// ---------- 法則の一覧（暗記カードの「定義・法則」から開く）：分野ごとに法則だけを並べる ----------
async function renderConcepts() {
  view.innerHTML = '<p class="muted" style="padding:40px 4px">読み込んでいます…</p>';
  let K;
  try { K = await loadKnowledge(); knowledgeCache = K; } catch (e) { view.innerHTML = `<p class="muted">${esc(e.message)}</p>`; return; }
  // 法則の分野は、その法則が結びつく概念の分野から決める
  const subjOf = (l) => { const c = (l.concepts || []).map(id => K.concept[id]).find(Boolean); const f = c && FIELDS.find(x => x.fieldKey === c.field); return f ? f.subject : 'その他'; };
  const groups = new Map();
  for (const l of K.laws) { const s = subjOf(l); if (!groups.has(s)) groups.set(s, []); groups.get(s).push(l); }
  const order = [...groups.keys()].sort((x, y) => (SUBJECT_ORDER.indexOf(x) + 1 || 99) - (SUBJECT_ORDER.indexOf(y) + 1 || 99));
  view.innerHTML = `
    <div class="fade-in">
      <div class="kicker">法則の一覧</div>
      <h1>使う法則と導出</h1>
      <p class="sub">法則をタップすると、成り立つ条件・導出（または根拠）・使う問題が出ます。</p>
      ${order.map(sj => `<h2>${esc(sj)}</h2><div class="laws">${groups.get(sj).map(lawRow).join('')}</div>`).join('')}
    </div>`;
  math(view);
}

// 法則の1行（概念ページ・一覧）。タップで法則のページへ
const ORIGIN = { fundamental: '出発点となる法則', empirical: '実験から得た法則', derived: '導出できる', definition: '定義' };
function lawRow(l) {
  return `<a class="law" href="#law=${esc(l.id)}"><b>${esc(l.name)}</b><span class="qhtml">${l.html}</span>
    <small class="muted"><span class="origin o-${esc(l.origin)}">${ORIGIN[l.origin] || ''}</span>${l.condHtml || esc(l.conditions)}</small></a>`;
}

// ---------- 法則のページ（式・適用条件・導出または根拠・関連する概念・問題） ----------
async function renderLaw(id) {
  view.innerHTML = '<p class="muted" style="padding:40px 4px">読み込んでいます…</p>';
  let K;
  try { K = await loadKnowledge(); knowledgeCache = K; } catch (e) { view.innerHTML = `<p class="muted">${esc(e.message)}</p>`; return; }
  const l = K.law[id];
  if (!l) { location.hash = '#concepts'; return; }
  const ps = ENTRIES.filter(e => (e.l || []).includes(l.id));
  const uses = K.laws.filter(x => x.id !== l.id && (x.steps || []).join(' ').includes(l.name));
  view.innerHTML = `
    <div class="fade-in law-page">
      <a class="back" href="#concepts">‹ 法則の一覧</a>
      <div class="kicker">法則 · <span class="origin o-${esc(l.origin)}">${ORIGIN[l.origin] || ''}</span></div>
      <h1>${esc(l.name)}</h1>
      <div class="law-formula qhtml">\\[${l.tex}\\]</div>
      <div class="panel law-cond"><b>成り立つ条件</b><div class="qhtml">${l.condHtml || esc(l.conditions)}</div></div>
      ${l.basisHtml ? `<h2>${l.origin === 'derived' ? 'もとになる事実' : 'この法則の根拠'}</h2><div class="panel qhtml">${l.basisHtml}</div>` : ''}
      ${l.steps ? `<h2>導出</h2>
        <p class="sub">1段ずつ開いて、次の段を自分で考えてから確かめよう。</p>
        <ol class="derivation">${l.steps.map((st, i) => `<li class="${i === 0 ? 'open' : ''}"><div class="qhtml">${st}</div></li>`).join('')}</ol>
        <div class="btn-row"><button class="btn" id="dnext">次の段を見る</button><button class="btn ghost" id="dall">すべて表示</button></div>` : ''}
      ${uses.length ? `<h2>この法則から導かれるもの</h2><div class="laws">${uses.map(lawRow).join('')}</div>` : ''}
      ${ps.length ? `<button class="btn primary block" id="lprob" style="margin-top:18px">この法則を使う問題を解く（${ps.length}問）</button>` : ''}
    </div>`;
  math(view);
  const items = $$('.derivation li');
  const reveal = () => { const n = items.find(li => !li.classList.contains('open')); if (n) n.classList.add('open'); if (!items.some(li => !li.classList.contains('open'))) { $('#dnext')?.remove(); $('#dall')?.remove(); } };
  $('#dnext')?.addEventListener('click', reveal);
  $('#dall')?.addEventListener('click', () => { items.forEach(li => li.classList.add('open')); $('#dnext')?.remove(); $('#dall')?.remove(); });
  if (items.length === 1) { $('#dnext')?.remove(); $('#dall')?.remove(); }
  $('#lprob')?.addEventListener('click', () => openIds(ps.map(e => e.id)));
}

// ---------- 定着度 ----------
function fieldStats(keys) {
  const set = new Set(keys);
  const list = ENTRIES.filter(e => set.has(e.f));
  const logs = state.log.filter(l => { const e = entryById.get(l.p); return e && set.has(e.f); });
  // 各問題の最新の結果
  const last = new Map();
  for (const l of logs) last.set(l.p, l);
  const latest = [...last.values()];
  return {
    total: list.length,
    solved: list.filter(e => state.problems[e.id]).length,
    full: latest.filter(l => l.g >= 3).length,
    tried: latest.length,
    due: dueEntries(state, { fields: keys }).length,
  };
}

// 得意度：各問題の最新の結果（完答 1・一部やヒント 0.5・解けず 0）を難易度で重みづけした平均。
// 少ない回数で極端にならないよう 0.5 に2問分寄せる。3問未満は「データ不足」。
function strength(keys) {
  const set = new Set(keys);
  const last = new Map();
  for (const l of state.log) { const e = entryById.get(l.p); if (e && set.has(e.f)) last.set(l.p, { l, e }); }
  let sw = 0, ss = 0;
  for (const { l, e } of last.values()) {
    const w = 1 + 0.25 * e.lv;
    const sc = l.g >= 3 ? 1 : l.g === 2 ? 0.5 : 0;
    sw += w; ss += w * sc;
  }
  const n = last.size;
  return { n, score: n >= 3 ? (ss + 1) / (sw + 2) : null };
}
const SUBJ_SHORT = { '分野横断': '融合' };
// 得意度 → タイルの色（苦手＝赤系、得意＝緑系、データ不足＝無地）
// 苦手な単元：得意度が 0.5 未満の中で最も低いもの（3問以上の結果がある単元だけ）
function weakestField() {
  let w = null;
  for (const f of FIELDS) {
    const st = strength([f.key]);
    if (st.score != null && st.score < 0.5 && (!w || st.score < w.score)) w = { f, ...st };
  }
  return w;
}
function strengthStyle(sc) {
  if (sc == null) return '';
  const t = Math.round(Math.min(1, Math.abs(sc - 0.5) * 2) * 55 + 8);
  return `background:color-mix(in srgb, var(${sc >= 0.5 ? '--good' : '--bad'}) ${t}%, var(--card));`;
}
// 科目のレーダー図（SVG）：外側ほど得意。データ不足の科目は中心に点線
function radarSVG(items) {
  const N = items.length, R = 104, cx = 150, cy = 136;
  const pt = (i, r) => { const a = -Math.PI / 2 + i * 2 * Math.PI / N; return [cx + r * Math.cos(a), cy + r * Math.sin(a)]; };
  const ring = (f) => items.map((_, i) => pt(i, R * f).map(v => v.toFixed(1)).join(',')).join(' ');
  const val = items.map((it, i) => pt(i, R * (it.score ?? 0)).map(v => v.toFixed(1)).join(',')).join(' ');
  const has = items.some(it => it.score != null);
  return `<svg class="radar" viewBox="0 0 300 278" role="img" aria-label="科目ごとの得意度">
    ${[0.25, 0.5, 0.75, 1].map(f => `<polygon points="${ring(f)}" class="rg${f === 0.5 ? ' mid' : ''}"/>`).join('')}
    ${items.map((_, i) => { const [x, y] = pt(i, R); return `<line x1="${cx}" y1="${cy}" x2="${x.toFixed(1)}" y2="${y.toFixed(1)}" class="rg"/>`; }).join('')}
    ${has ? `<polygon points="${val}" class="rv"/>` : ''}
    ${items.map((it, i) => it.score == null ? '' : (() => { const [x, y] = pt(i, R * it.score); return `<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="3.5" class="rp"/>`; })()).join('')}
    ${items.map((it, i) => { const [x, y] = pt(i, R + 20); return `<text x="${x.toFixed(1)}" y="${(y + 4).toFixed(1)}" text-anchor="middle" class="rl${it.score == null ? ' none' : ''}">${esc(it.label)}<tspan class="rs" x="${x.toFixed(1)}" dy="13">${it.score == null ? '—' : Math.round(it.score * 100)}</tspan></text>`; }).join('')}
    <text x="${cx + 4}" y="${cy - R * 0.5 + 3}" class="rt">50</text>
  </svg>`;
}
function strengthHTML() {
  const subs = subjectsWithFields();
  const fs = FIELDS.map(f => ({ f, ...strength([f.key]) }));
  const judged = fs.filter(x => x.score != null).sort((a, b) => b.score - a.score);
  const good = judged.filter(x => x.score >= 0.6).slice(0, 3);
  const weak = judged.filter(x => x.score < 0.5).slice(-3).reverse();
  const radar = radarSVG(subs.map(s => ({ label: s.name, ...strength(s.fields.map(f => f.key)) })));
  const name = (x) => `${esc(SUBJ_SHORT[x.f.subject] || x.f.subject)}・${esc(x.f.name)}`;
  return `
    <h2>得意・苦手</h2>
    <div class="panel strength">
      ${judged.length ? '' : '<p class="muted" style="margin:0 0 8px">1つの分野で3問以上解くと、その分野の得意・苦手が出ます。</p>'}
      ${radar}
      <div class="sw-cols">
        <div><b class="sw-good">得意</b>${good.length ? good.map(x => `<span>${name(x)} <em>${Math.round(x.score * 100)}</em></span>`).join('') : '<span class="muted">まだありません</span>'}</div>
        <div><b class="sw-weak">苦手</b>${weak.length ? weak.map(x => `<span>${name(x)} <em>${Math.round(x.score * 100)}</em></span>`).join('') : '<span class="muted">まだありません</span>'}</div>
      </div>
      <div class="tiles">${subs.map(s => `<div class="tile-row"><span class="tile-sub">${esc(SUBJ_SHORT[s.name] || s.name)}</span><div class="tile-grid">${s.fields.map(f => {
        const x = fs.find(y => y.f.key === f.key);
        return `<button class="tile" data-f="${esc(f.key)}" style="${strengthStyle(x.score)}"><span>${esc(f.name)}</span><b>${x.score == null ? (x.n ? `${x.n}問` : '—') : Math.round(x.score * 100)}</b></button>`;
      }).join('')}</div></div>`).join('')}</div>
      <p class="muted legend"><i class="lg-weak"></i>苦手<i class="lg-mid"></i>ふつう<i class="lg-good"></i>得意　数字は得意度（0〜100）。最後に解いた結果を難しい問題ほど重く数えます。タップでその分野を演習。</p>
    </div>`;
}

function renderMap() {
  const withConf = state.log.filter(l => l.c);
  const high = withConf.filter(l => l.c === 'high');
  const low = withConf.filter(l => l.c === 'low');
  const over = high.filter(l => !l.ok).length;
  const under = low.filter(l => l.ok).length;
  const total = state.log.length;
  const solo = state.log.filter(l => l.g >= 3).length;
  const pct = (a, b) => b ? Math.round(a / b * 100) + '%' : '—';
  view.innerHTML = `
    <div class="fade-in">
      <div class="kicker">${FIELDS.length}分野・${ENTRIES.length}問</div>
      <h1>定着度</h1>
      <p class="sub">分野ごとの進み具合と、最後に解いたときに完答できた割合です。</p>
      ${strengthHTML()}
      <h2>自信と結果のズレ</h2>
      <div class="calib">
        <div class="stat"><b>${pct(over, high.length)}</b><span>「確実」で完答できず（過信）</span></div>
        <div class="stat"><b>${pct(under, low.length)}</b><span>「自信なし」で完答</span></div>
      </div>
      <p class="muted" style="margin:8px 2px 0">${total ? `これまで ${total}回解答、うち自力で完答 ${solo}回。` : 'まだ記録がありません。'}過信が多いときは、解き終えたら検算する習慣を。</p>
      ${subjectsWithFields().map(sub => `
        <h2>${esc(sub.name)}</h2>
        ${sub.fields.map(f => {
          const st = fieldStats([f.key]);
          const rate = st.tried ? st.full / st.tried : 0;
          return `<button class="card-row" data-f="${esc(f.key)}">
            <span><span class="nm">${esc(f.name)}</span><br>
              <span class="muted">解いた ${st.solved} / ${st.total}問・完答率 ${pct(st.full, st.tried)}${st.due ? `・復習待ち ${st.due}問` : ''}</span></span>
            <span class="chip ${st.due ? 'st-due' : st.tried ? 'st-learning' : 'st-new'}">${st.due ? '復習' : st.tried ? '学習中' : '未着手'}</span>
            <span class="bar"><div style="width:${Math.round(st.solved / Math.max(1, st.total) * 100)}%;background:${rate >= 0.7 ? 'var(--good)' : 'var(--accent)'}"></div></span>
          </button>`;
        }).join('')}`).join('')}
    </div>`;
  $$('.card-row, .tile').forEach(b => b.addEventListener('click', () => {
    const f = b.dataset.f;
    const rev = pickReview(state, { fields: [f], size: state.practice.size });
    if (rev.length) { startSession(rev, { review: true, fields: [f], size: state.practice.size }); return; }
    const cfg = { fields: [f], levels: state.practice.levels.slice(), size: state.practice.size };
    const st = pickPractice(state, cfg);
    if (st.length) startSession(st, cfg); else toast('この条件の問題はありません');
  }));
}

// ---------- 分野を選ぶ ----------
function renderPractice() {
  const pr = state.practice;
  // 保存された選択のうち、今はない単元は外す
  const saved = (pr.fields || []).filter(k => fieldByKey[k]);
  let fields = new Set(saved.length ? saved : FIELDS.map(f => f.key));
  const draw = () => {
    const n = pool({ fields: [...fields], levels: pr.levels }).length;
    const bar = levelToolbar(pr, draw);
    view.innerHTML = `
      <div class="fade-in">
        <div class="kicker">演習モード</div>
        <h1>分野を選ぶ</h1>
        <p class="sub">選んだ範囲から、まだ解いていない問題・久しぶりの問題を優先して、分野を混ぜて出します。</p>
        ${bar.html}
        <div class="btn-row" style="margin-top:4px"><button class="btn ghost" id="all">すべて選択</button><button class="btn ghost" id="none">すべて解除</button></div>
        ${subjectsWithFields().map(sub => {
          const allOn = sub.fields.every(f => fields.has(f.key));
          return `<div class="subject">
            <button class="subj ${allOn ? 'on' : ''}" data-s="${esc(sub.name)}">${esc(sub.name)}</button>
            <div class="toggles">${sub.fields.map(f => `<button class="tog ${fields.has(f.key) ? 'on' : ''}" data-f="${esc(f.key)}">${esc(f.name)}<small>${pool({ fields: [f.key], levels: pr.levels }).length}</small></button>`).join('')}</div>
          </div>`;
        }).join('')}
        <button class="btn primary block" id="go" style="margin-top:20px" ${n ? '' : 'disabled'}>${n ? `始める（該当 ${n}問から ${Math.min(n, pr.size)}問）` : '該当する問題がありません'}</button>
      </div>`;
    bar.bind();
    const persist = () => { pr.fields = [...fields]; save(); draw(); };
    $$('.tog[data-f]').forEach(b => b.addEventListener('click', () => { const k = b.dataset.f; fields.has(k) ? fields.delete(k) : fields.add(k); persist(); }));
    $$('.subj').forEach(b => b.addEventListener('click', () => {
      const ks = FIELDS.filter(f => f.subject === b.dataset.s).map(f => f.key);
      const allOn = ks.every(k => fields.has(k));
      ks.forEach(k => allOn ? fields.delete(k) : fields.add(k));
      persist();
    }));
    $('#all').addEventListener('click', () => { fields = new Set(FIELDS.map(f => f.key)); persist(); });
    $('#none').addEventListener('click', () => { fields = new Set(); persist(); });
    $('#go').addEventListener('click', () => {
      const cfg = { fields: [...fields], levels: pr.levels.slice(), size: pr.size };
      const st = pickPractice(state, cfg);
      if (st.length) startSession(st, cfg);
    });
  };
  draw();
}

// ---------- 設定 ----------
// 問題管理表の「弟テスト記録」シートと同じ列で書き出す
function testLogCSV(since = 0) {
  const head = ['日付', '問題ID', '科目', '分野', '難易度', '正誤', '所要時間(分)', '手応え(1-5)', 'どこで詰まったか', '感想・改善点'];
  const q = (v) => { const t = v == null ? '' : String(v); return /[",\n]/.test(t) ? `"${t.replace(/"/g, '""')}"` : t; };
  const mark = (l) => l.sg ? { full: '○', part: '△', none: '×' }[l.sg] : (l.g >= 3 ? '○' : l.g === 2 ? '△' : '×');
  const rows = state.log.filter(l => l.t > since).map(l => {
    const e = entryById.get(l.p), f = e && fieldByKey[e.f];
    const d = new Date(l.t);
    const date = `${d.getFullYear()}/${d.getMonth() + 1}/${d.getDate()}`;
    const disp = l.disp && l.disp.length ? `判定に異議：${l.disp.join(' / ')}` : '';
    const comment = [l.sg === 'full' ? (l.note || '') : '', disp].filter(Boolean).join('　');
    return [date, l.p, f ? f.subject : '', f ? f.name : '', e ? e.lv : '', mark(l), l.min || '', l.feel || '', l.sg === 'full' ? '' : (l.note || ''), comment].map(q).join(',');
  });
  return '﻿' + [head.join(','), ...rows].join('\r\n');
}
function download(name, text, type) {
  const blob = new Blob([text], { type });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

// 記録の点検に使う ID の一覧（今の問題データ・暗記カード。カードを読み込む前はカードの ID は調べない）
function idContext() {
  return { problemIds: new Set(ENTRIES.map(e => e.id)), ...(cardIdSet ? { cardIds: cardIdSet } : {}) };
}
function backupFile() {
  download(`alchemia-${dayKey()}.json`, JSON.stringify(makeBackup(state, VERSION), null, 1), 'application/json');
}
// 読み込み・控えへの復帰の共通部分：今の記録の控えを残す → 新しい記録を保存 → 点検。失敗したら元のまま
async function applyRestore(incoming, mode, reason) {
  const before = counts(state);
  const next = restoreState(state, incoming, mode);
  const rep = audit(next, idContext());
  if (rep.errors.length) { toast('読み込んだ結果に正しくない値があるため、取りやめました（今の記録はそのままです）'); return false; }
  const r = store.replace(state, next, reason);
  if (!r.ok) { toast(`${SAVE_MSG[r.reason] || SAVE_MSG.write}（読み込みは取りやめました。今の記録はそのままです）`); return false; }
  if (!r.snapshotted) backupFile();   // 端末に控えを残せなかったときは、置き換える前の記録をファイルで残す（この時点では保存前の state）
  state = r.state; lwwSig = null; session = null;
  const after = counts(state);
  const warn = rep.warnings.reduce((a, w) => a + w.count, 0);
  await choose({ title: '読み込みました',
    body: `解答の記録 ${before.log} → ${after.log}件、復習状態 ${before.pstat} → ${after.pstat}問、暗記カード ${before.cards} → ${after.cards}枚、ブックマーク ${before.bookmarks} → ${after.bookmarks}件。${warn ? `注意 ${warn}件（今のデータに無い問題IDの記録など。記録は残してあります）。` : '点検で問題はありませんでした。'}${r.snapshotted ? '読み込む前の記録は、設定の「自動の控え」から戻せます。' : '読み込む前の記録はファイルに書き出しました。'}`,
    options: [{ key: 'ok', label: 'OK', cls: 'primary' }] });
  route();
  if (sync && sync.loggedIn()) sync.syncNow();
  return true;
}

// 設定の「学習記録のバックアップ」：点検の結果・自動の控え・読めなかった記録
function integrityHTML() {
  const rep = audit(state, idContext());
  const nw = rep.warnings.reduce((a, w) => a + w.count, 0), ne = rep.errors.reduce((a, w) => a + w.count, 0);
  const items = [...rep.errors.map(x => ['bad', x]), ...rep.warnings.map(x => ['', x])];
  const snaps = store.snapshots(), broken = store.brokenKeys(), use = store.usage();
  const kb = (n) => `${Math.round(n / 1024).toLocaleString('ja-JP')}KB`;
  return `<details class="integrity" style="margin-top:10px"><summary class="muted small">記録の点検：${ne ? `<b style="color:var(--bad)">正しくない値 ${ne}件</b>` : '正しくない値はありません'}${nw ? `・注意 ${nw}件` : ''}（解答の記録 ${rep.stats.log}件）</summary>
      ${items.length ? `<ul class="small" style="margin:6px 0;padding-left:1.1em">${items.map(([c, x]) => `<li${c ? ' style="color:var(--bad)"' : ''}>${esc(x.msg)}：${x.count}件（${esc(x.ids.slice(0, 5).join(', '))}${x.count > 5 ? ' ほか' : ''}）</li>`).join('')}</ul>
        <p class="muted small" style="margin:0">記録は自動では直しません（消える記録が出ないように）。今のデータに無い問題の記録は、そのまま残してあります。</p>` : ''}
      <p class="muted small" style="margin:4px 0 0">この端末に保存している量：記録 ${kb(use.main)}・控え ${kb(use.snaps)}${use.broken ? `・退避 ${kb(use.broken)}` : ''}（多くのブラウザの上限は約5MB）</p>
    </details>
    ${snaps.length || broken.length ? `<div class="small" style="margin-top:10px"><b>自動の控え</b>（この端末の中だけ。移行・読み込み・消去の前に残します）
      <ul style="margin:6px 0;padding-left:1.1em">${snaps.map((sn, i) => `<li>${new Date(sn.at).toLocaleString('ja-JP')}：${esc(sn.reason || '')}（解答 ${sn.n}件） <button class="btn ghost snap-back" data-i="${i}" style="padding:2px 10px">戻す</button></li>`).join('')}
      ${broken.map(k => `<li>読めなかった記録（${new Date(+k.slice(BROKEN_PREFIX.length)).toLocaleString('ja-JP')}） <button class="btn ghost broken-dl" data-k="${esc(k)}" style="padding:2px 10px">書き出す</button></li>`).join('')}</ul>
      <button class="btn ghost" id="snap-clear" style="padding:2px 10px">控えを消す</button></div>` : ''}`;
}

function renderSettings() {
  const since = state.csvExportedAt || 0;
  const newCount = state.log.filter(l => l.t > since).length;
  view.innerHTML = `
    <div class="fade-in">
      <div class="kicker">v${VERSION}</div>
      <h1>設定</h1>
      ${sync ? '<h2>アカウントと同期</h2><div class="panel" id="acct"></div>' : ''}
      <h2>解答の形式</h2>
      <div class="panel">
        <div class="setting"><span>答え方</span>
          <div class="segment amode" id="amode">${ANSWER_MODES.map(([k, v]) => `<button data-m="${k}" class="${answerMode() === k ? 'on' : ''}">${v}</button>`).join('')}</div></div>
        <p class="muted" style="margin:8px 0 0">${{
          input: '<b>記述式</b>：語句・化学式・数値を入力して自動で採点します。「選ぶ」「並べる」の問題はそのまま選びます。説明・図の小問と、入力用の答えがない問題は自己採点です。',
          choice: '<b>選択式</b>：語句・化学式・数値の欄を4つの選択肢から選びます（やり直しはなし）。選択肢のない欄は入力、説明の小問は自己採点です。',
          self: '<b>自己採点</b>：紙に書いてから、答えと解説を見て自分で採点します。',
        }[answerMode()]}</p>
      </div>
      <h2>表示と音</h2>
      <div class="panel">
        <div class="setting"><span>画面のテーマ</span>
          <div class="segment theme" id="theme">${[['auto', '自動'], ['light', 'ライト'], ['dark', 'ダーク']].map(([k, v]) => `<button data-t="${k}" class="${(state.settings.theme || 'auto') === k ? 'on' : ''}">${v}</button>`).join('')}</div></div>
        <div class="setting"><span>効果音</span><button class="switch ${state.settings.sound !== false ? 'on' : ''}" id="snd" role="switch" aria-checked="${state.settings.sound !== false}"><i></i></button></div>
        <div class="setting"><span>${MASCOT_NAME}（サラマンダー）の声かけ</span><button class="switch ${state.settings.mascot !== false ? 'on' : ''}" id="brd" role="switch" aria-checked="${state.settings.mascot !== false}"><i></i></button></div>
      </div>
      <h2>テスト記録の書き出し</h2>
      <div class="panel">
        <div class="btn-row" style="margin-top:0">
          <button class="btn primary" id="csvnew" ${newCount ? '' : 'disabled'}>新しい${newCount}件</button>
          <button class="btn" id="csvall" ${state.log.length ? '' : 'disabled'}>すべて（${state.log.length}件）</button>
        </div>
        ${since ? `<p class="muted" style="margin:8px 0 0">前回の書き出し：${new Date(since).toLocaleString('ja-JP')}</p>` : ''}
      </div>
      <h2>学習記録のバックアップ</h2>
      <div class="panel">
        <p class="muted" style="margin:0 0 10px">${sync && sync.loggedIn() ? '学習記録はこの端末とアカウント（サーバー）に保存されています。手元にも残したいときは書き出してください。' : '学習記録はこの端末のブラウザ内だけに保存されます。機種変更の前に書き出してください。Safari で使っていて7日以上開かないと、記録が消えることがあります（ホーム画面に追加すると消えにくくなります）。'}</p>
        <div class="btn-row" style="margin-top:0">
          <button class="btn" id="exp">書き出す</button>
          <label class="btn" for="imp">読み込む</label>
          <input type="file" id="imp" accept="application/json" hidden>
        </div>
        ${integrityHTML()}
        <button class="btn ghost block" id="reset" style="margin-top:10px;color:var(--bad)">記録をすべて消す</button>
      </div>
      <h2>設計の根拠</h2>
      <div class="panel">
        <ul class="refs" style="margin:0;padding-left:1.1em">
          <li><strong>想起練習・分散学習</strong>：学習法10種の比較で最も評価が高い2つ。<a href="https://www.psychologicalscience.org/news/releases/which-study-strategies-make-the-grade.html" target="_blank" rel="noopener">Dunlosky ほか 2013</a></li>
          <li><strong>交互練習</strong>：中学1年787人の無作為化試験で、1か月後のテストが 61% 対 38%。<a href="https://www.structural-learning.com/research-briefings/mix-it-up" target="_blank" rel="noopener">Rohrer ほか 2020</a></li>
          <li><strong>答えを出さないヒント</strong>：制限なしのAIで練習した生徒は、AIなしの試験で成績が下がった。<a href="https://papers.ssrn.com/abstract=4895486" target="_blank" rel="noopener">Bastani ほか 2025</a></li>
        </ul>
        <p class="muted" style="margin:10px 0 0">復習間隔は FSRS の考え方を簡略化したモデルで推定しています。難関・最難関レベルでの効果は未検証です。</p>
        <p class="muted small" style="margin:8px 0 0"><a href="privacy.html">プライバシーポリシー</a>・<a href="terms.html">利用規約</a></p>
      </div>
    </div>`;
  drawAccount();
  $$('#amode button').forEach(b => b.addEventListener('click', () => {
    state.settings.answer = b.dataset.m; save(); play('tap'); renderSettings();
  }));
  $$('#theme button').forEach(b => b.addEventListener('click', () => {
    state.settings.theme = b.dataset.t; save(); applyTheme(b.dataset.t); renderSettings();
  }));
  $('#snd').addEventListener('click', () => {
    state.settings.sound = state.settings.sound === false; save(); setSound(state.settings.sound); if (state.settings.sound) play('correct'); renderSettings();
  });
  $('#brd').addEventListener('click', () => {
    state.settings.mascot = state.settings.mascot === false; save();
    if (state.settings.mascot) { mountMascot(); mascotSay('よろしくね。'); } else $('#mascot')?.remove();
    renderSettings();
  });
  const csv = (from) => {
    download(`alchemia-記録-${dayKey()}.csv`, testLogCSV(from), 'text/csv');
    state.csvExportedAt = Date.now(); save(); renderSettings();
  };
  $('#csvnew').addEventListener('click', () => csv(since));
  $('#csvall').addEventListener('click', () => csv(0));
  $('#exp').addEventListener('click', backupFile);
  $('#imp').addEventListener('change', async (e) => {
    const f = e.target.files[0]; e.target.value = ''; if (!f) return;
    let data;
    try {
      if (f.size > 50e6) throw new Error('大きすぎるファイルです');
      let obj;
      try { obj = JSON.parse(await f.text()); } catch { throw new Error('ファイルが壊れていて読めません（JSON ではありません）'); }
      data = readBackup(obj, idContext());
    } catch (err) { toast(err && err.message ? err.message : '読み込めないファイルです'); return; }
    const k = hasRecords(state) ? await choose({ title: 'バックアップを読み込む',
      body: `ファイルには解答の記録が ${data.state.log.length}件あります${data.meta.createdAt ? `（${new Date(data.meta.createdAt).toLocaleString('ja-JP')} に書き出し）` : ''}。今の記録（${state.log.length}件）をどうしますか？`,
      options: [{ key: 'merge', label: '今の記録と合わせる', cls: 'primary' }, { key: 'replace', label: '今の記録を、ファイルの記録に置き換える' }] }) : 'replace';
    if (!k) return;
    await applyRestore(data.state, k, k === 'merge' ? 'バックアップと合わせる前' : 'バックアップで置き換える前');
  });
  $$('.snap-back').forEach(b => b.addEventListener('click', async () => {
    const sn = store.snapshots()[+b.dataset.i];
    if (!sn) return;
    let data;
    try { data = readBackup(JSON.parse(sn.raw), idContext()); } catch (err) { toast(`この控えは読めません：${err.message}`); return; }
    const k = await choose({ title: '控えの時点に戻す', body: `${new Date(sn.at).toLocaleString('ja-JP')}（${esc(sn.reason || '')}）の記録 ${data.state.log.length}件に戻します。今の記録（${state.log.length}件）は新しい控えとして残します。`,
      options: [{ key: 'replace', label: 'この時点に戻す', cls: 'primary' }, { key: 'merge', label: '今の記録と合わせる' }] });
    if (!k) return;
    await applyRestore(data.state, k, '控えに戻す前');
  }));
  $('#snap-clear')?.addEventListener('click', async () => {
    const k = await choose({ title: '自動の控えを消す', body: 'この端末に残している控え（と、読めなかった記録の退避）を消します。学習記録の本体は消えません。',
      options: [{ key: 'y', label: '消す', cls: 'danger' }] });
    if (!k) return;
    const r = store.clearSnapshots();
    if (!r.ok) toast(`${r.failed.length}件の控えを消せませんでした（学習記録の本体は変わっていません）`);
    renderSettings();
  });
  $$('.broken-dl').forEach(b => b.addEventListener('click', () => {
    const v = store.get(b.dataset.k);
    if (v) download(`alchemia-読めなかった記録-${b.dataset.k.slice(BROKEN_PREFIX.length)}.json`, v, 'application/json');
  }));
  $('#reset').addEventListener('click', async () => {
    const logged = sync && sync.loggedIn();
    const k = await choose({ title: '記録をすべて消す',
      body: (logged ? '解答の記録・復習の予定・暗記カードの記録・ブックマーク・メモを消します（設定は残ります）。' : '学習記録をすべて消します（設定は残ります）。')
        + '消す前の記録は、この端末に控えとして残ります（設定の「自動の控え」から戻せます。控えも消すには「控えを消す」）。',
      options: logged
        ? [{ key: 'all', label: 'この端末とアカウントの記録を消す', cls: 'danger' }, { key: 'device', label: 'ログアウトして、この端末の記録だけ消す' }]
        : [{ key: 'all', label: '消す', cls: 'danger' }] });
    if (!k) return;
    const keep = keepPrefs();
    if (k === 'device') {
      // この端末の分だけ消す：先にサーバーへ送れたことを確かめる（送れていない記録を消さない）
      if (!(await sync.syncNow())) { toast('サーバーと同期できなかったため、消していません'); return; }
      await sync.logout();
      const r = store.replace(state, { ...freshState(), ...keep }, null, { fence: true });
      if (!r.ok) { toast(`${SAVE_MSG[r.reason] || SAVE_MSG.write}（消していません）`); return; }
      state = r.state;
    } else {
      const next = { ...freshState(), ...keep, resetAt: Date.now() };   // 消した印：同期で、ほかの端末・サーバーの古い記録を持ち込まない
      const r = store.replace(state, next, '消去の前');
      if (!r.ok) { toast(`${SAVE_MSG[r.reason] || SAVE_MSG.write}（消していません）`); return; }
      if (!r.snapshotted) backupFile();
      state = r.state;
    }
    lwwSig = null; session = null; clearSavedSession(); toast('消去しました'); route();
    if (k === 'all' && logged) sync.syncNow();
  });
}

// テスト・確認用：指定した問題IDでセットを始める（例：__alchemia.open(['MEC-INC-001'])）
window.__alchemia = { open: (ids) => startSession(ids.map(pid => ({ kind: 'problem', pid })), { size: ids.length }) };

// ---------- 起動 ----------
applyTheme(state.settings.theme || 'auto');
setSound(state.settings.sound !== false);
window.matchMedia('(prefers-color-scheme: dark)').addEventListener?.('change', () => applyTheme(state.settings.theme || 'auto'));
mountMascot();
route();
bootWarn();
loadCatalog().then((meta) => {
  ready = true;
  loadKnowledge().then(k => { knowledgeCache = k; }).catch(() => {});
  route();
  initSync();
  if (meta.ok) dataWarn();
  if (!meta.ok) showBanner(`<b>問題データを読み込めませんでした</b>${esc(meta.error || '')}。リポジトリに data フォルダがアップロードされているか確認してください。`);
});
// アップロード漏れの検出：各ファイルの版が app.js と一致するか
function staleFiles() {
  const cssV = getComputedStyle(document.documentElement).getPropertyValue('--alchemia-version').replace(/['"\s]/g, '');
  const list = { 'js/answer.js': M_ANS.FILE_VERSION, 'js/ui.js': M_UI.FILE_VERSION, 'js/store.js': M_STORE.FILE_VERSION, 'js/catalog.js': M_CAT.FILE_VERSION, 'js/scheduler.js': M_SCHED.FILE_VERSION, 'js/parser.js': M_PARSER.FILE_VERSION, 'js/problem/figures.js': M_FIG.FILE_VERSION, 'css/style.css': cssV };
  return Object.entries(list).filter(([, v]) => v !== VERSION).map(([k]) => k);
}
// 起動時の記録の読み込みで起きたこと（読めない・一部壊れていた・新しい版の記録）
function bootWarn() {
  const where = '設定の「学習記録のバックアップ」から書き出せます。';
  const msg = {
    broken: boot.quarantined ? `<b>保存されていた記録を読めませんでした</b>元の記録は消さずに退避してあります（${where}）。` : '<b>保存されていた記録を読めませんでした</b>元の記録を退避できなかったため、この端末への保存を止めています。新しく解いた記録は、この画面を閉じると残りません。',
    salvaged: `<b>記録の一部が壊れていました</b>壊れていた部分（${esc(boot.issues.map(i => i.key || i.code).join('・'))}）だけを初期値に戻し、ほかの記録はそのまま使います。元の記録は退避してあります（${where}）。`,
    newer: '<b>新しい版の ALCHEMIA で保存された記録です</b>この版のアプリは古いため、記録を守るよう<b>保存も同期もしません</b>（解いた結果は残りません）。ページを再読み込みしてアプリを更新してください。',
    migrate_failed: `<b>記録を新しい形に移せませんでした</b>元の記録を守るため、この端末への保存を止めています（解いた結果は残りません）。アプリを更新してから再読み込みすると、もう一度移します。元の記録は退避もしてあります（${where}）。`,
    unavailable: '<b>この端末に保存した記録を読めませんでした</b>記録を上書きしないよう、保存を止めています。しばらくしてから再読み込みしてください（保存できるようになれば、保存されている記録と合わせてから保存します）。',
  }[boot.status] || (noStorage ? '<b>この環境では記録を保存できません</b>プライベートモードなどで、ブラウザの保存領域が使えません。解いた結果は画面を閉じると残りません。' : '');
  if (msg) showBanner(msg);
}
// 同じ端末の別のタブが記録を保存した：この画面の記録と合わせる（解いている途中なら、画面は描き直さない）
window.addEventListener('storage', (e) => {
  if (e.key !== M_STORE.KEY || e.storageArea !== LS) return;
  save();
  if (!session && !document.activeElement?.matches('input, textarea')) route();
});
function showBanner(html) {
  let b = $('.banner');
  if (!b) { b = document.createElement('div'); b.className = 'banner'; document.body.prepend(b); }
  b.innerHTML = html;
}
// 問題データの版の食い違い（一部の data/ ファイルだけ古い、形式が新しすぎる）
function dataWarn() {
  const files = M_CAT.staleDataFiles ? M_CAT.staleDataFiles() : [];
  const msgs = [...(META.warnings || [])];
  if (files.length) msgs.push(`${files.map(esc).join('、')} が data/index.json と違う版です。data フォルダをすべて上書きアップロードしてください。`);
  if (msgs.length) showBanner(`<b>問題データの版が揃っていません</b>${msgs.join('<br>')}`);
}
setTimeout(() => {
  const stale = staleFiles();
  if (stale.length) showBanner(`<b>古いファイルが混ざっています（v${VERSION}）</b>${stale.map(esc).join('、')} が更新されていません。リポジトリにzipの中身をすべて上書きアップロードしてください。`);
}, 400);
// 予期しないエラーを画面に出す（報告用）
window.addEventListener('error', (e) => showBanner(`<b>エラー</b>${esc(e.message || '')} <small>${esc((e.filename || '').split('/').pop())}:${e.lineno || ''}</small>`));
window.addEventListener('unhandledrejection', (e) => showBanner(`<b>エラー</b>${esc(String(e.reason && e.reason.message || e.reason))}`));

if ('serviceWorker' in navigator && (location.protocol === 'https:' || location.hostname === 'localhost') && !/claude|claudeusercontent/.test(location.hostname)) {
  navigator.serviceWorker.register('./sw.js', { updateViaCache: 'none' }).then(r => r.update()).catch(() => {});
}
