// 暗記カード：図鑑・法則・誤概念から作ったカード（data/cards.json）を、間隔をあけてくり返す
//   #cards            デッキの一覧と、今日の暗記
//   #deck=ID          デッキのカード一覧（覚えた・まだ）
//   #study=ID|all     学習（表 → 答えを見る → もう一度・あいまい・覚えた）
// 記録は state.cards[カードID] = { box, due, last, n }（box 0〜6。覚えるたびに次に出るまでの間隔がのびる）
export const FILE_VERSION = '0.3.0';

const FIELD_NAME = { THE: '理論化学', INO: '無機化学', ORG: '有機化学', POL: '高分子化合物' };
const DAY = 86400000;
const INTERVAL = [0, 1, 3, 7, 14, 30, 60];   // 日（box ごと）
const NEW_PER_SESSION = 20;
const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];

let D = null, C = null;
const e = (s) => C.esc(s ?? '');
const rec = (id) => C.state.cards[id];
const learned = (id) => { const r = rec(id); return r && r.box >= 3; };
const isDue = (id, now) => { const r = rec(id); return r && r.due <= now; };

export function renderCards(view, data, tab, arg, ctx) {
  D = data; C = ctx;
  if (!C.state.cards) C.state.cards = {};
  if (tab === 'deck') return deckPage(view, arg);
  if (tab === 'study') return study(view, arg);
  return home(view);
}

export function dueCount(data, state, now = Date.now()) {
  return data.decks.reduce((s, d) => s + d.cards.filter(c => state.cards && state.cards[c.id] && state.cards[c.id].due <= now).length, 0);
}

function stats(deck, now) {
  let nw = 0, due = 0, ok = 0;
  for (const c of deck.cards) {
    const r = rec(c.id);
    if (!r) nw++; else { if (r.due <= now) due++; if (r.box >= 3) ok++; }
  }
  return { nw, due, ok, n: deck.cards.length };
}

function home(view) {
  const now = Date.now();
  const all = D.decks.map(d => ({ d, s: stats(d, now) }));
  const due = all.reduce((a, x) => a + x.s.due, 0);
  const ok = all.reduce((a, x) => a + x.s.ok, 0), n = all.reduce((a, x) => a + x.s.n, 0);
  const groups = [...new Set(D.decks.map(d => d.field))];
  view.innerHTML = `<div class="fade-in cards-home">
    <div class="kicker">ALCHEMIA · 暗記</div>
    <h1>暗記カード</h1>
    <p class="sub">表を見て思い出してから、答えを見て「もう一度・あいまい・覚えた」で答えます。覚えたカードほど、間をあけて出ます（間隔反復）。</p>
    <div class="card-hero">
      <div><b>${due}</b><span>今日の復習</span></div>
      <div><b>${ok}</b><span>覚えた / ${n}</span></div>
      <a class="btn primary" href="#study=all">${due ? `今日の復習を始める（${due}枚）` : '新しいカードを始める'}</a>
    </div>
    <label class="sw-toggle"><input type="checkbox" id="crev" ${C.prefs.reverse ? 'checked' : ''}> 裏から出す（答え → 名前）</label>
    ${groups.map(g => `<h2>${e(FIELD_NAME[g] || g)}</h2><div class="deck-list">${all.filter(x => x.d.field === g).map(({ d, s }) => `
      <div class="deck-card">
        <a class="deck-main" href="#deck=${e(d.id)}"><b>${e(d.name)}</b><span>${e(d.desc)}</span>
          <div class="deck-bar"><i style="width:${Math.round(100 * s.ok / s.n)}%"></i></div>
          <small>${s.ok} / ${s.n} 覚えた${s.due ? `・<em>復習 ${s.due}</em>` : ''}${s.nw ? `・未学習 ${s.nw}` : ''}</small></a>
        <a class="btn ghost deck-go" href="#study=${e(d.id)}">学習</a>
      </div>`).join('')}</div>`).join('')}
  </div>`;
  $('#crev').addEventListener('change', (ev) => { C.prefs.reverse = ev.target.checked; C.save(); });
}

function deckPage(view, id) {
  const d = D.decks.find(x => x.id === id);
  if (!d) { location.hash = '#cards'; return; }
  const now = Date.now(), s = stats(d, now);
  const tag = (c) => { const r = rec(c.id); return !r ? '<span class="ctag new">未学習</span>' : r.box >= 3 ? '<span class="ctag ok">覚えた</span>' : r.due <= now ? '<span class="ctag due">復習</span>' : '<span class="ctag mid">学習中</span>'; };
  view.innerHTML = `<div class="fade-in">
    <a class="back" href="#cards">‹ 暗記カード</a>
    <div class="kicker">${e(FIELD_NAME[d.field] || '')}</div>
    <h1>${e(d.name)}</h1>
    <p class="sub">${e(d.desc)}　${s.ok} / ${s.n} 覚えた</p>
    <div class="row-btns"><a class="btn primary" href="#study=${e(d.id)}">学習する</a><button class="btn ghost" id="dreset">このデッキの記録を消す</button></div>
    <input type="search" id="dq" class="deck-search" placeholder="カードを探す" autocomplete="off">
    <div class="card-rows">${d.cards.map(c => `<details class="card-row" data-t="${e((c.f + ' ' + c.b).replace(/<[^>]+>|\$|\\ce|[{}]/g, ''))}">
      <summary><span class="qhtml">${c.f}</span>${tag(c)}</summary>
      <div class="qhtml cb">${c.b}${c.n ? `<p class="muted small">${e(c.n)}</p>` : ''}${c.link ? `<a class="linkbtn" href="${e(c.link)}">${c.link.startsWith('#law') ? '導出・根拠を見る' : '図鑑で見る'} ›</a>` : ''}</div></details>`).join('')}</div>
  </div>`;
  C.math(view);
  $('#dq').addEventListener('input', (ev) => { const q = ev.target.value.trim(); $$('.card-row', view).forEach(r => { r.hidden = q && !r.dataset.t.includes(q); }); });
  $('#dreset').addEventListener('click', () => {
    if (!confirm(`「${d.name}」の暗記の記録を消しますか？`)) return;
    for (const c of d.cards) delete C.state.cards[c.id];
    C.save(); deckPage(view, id);
  });
}

// 出すカード：期限の来たもの（古い順）＋ 新しいカード（最大 NEW_PER_SESSION 枚）
function queueFor(arg) {
  const now = Date.now();
  const decks = arg === 'all' ? D.decks : D.decks.filter(d => d.id === arg);
  const cards = decks.flatMap(d => d.cards.map(c => ({ ...c, deck: d })));
  const due = cards.filter(c => isDue(c.id, now)).sort((a, b) => rec(a.id).due - rec(b.id).due);
  const nw = cards.filter(c => !rec(c.id)).slice(0, NEW_PER_SESSION);
  return [...due, ...nw];
}

function study(view, arg) {
  const q = queueFor(arg || 'all');
  const deck = D.decks.find(d => d.id === arg);
  const back = deck ? `#deck=${deck.id}` : '#cards';
  if (!q.length) {
    view.innerHTML = `<div class="fade-in"><a class="back" href="${back}">‹ もどる</a><h1>今は出すカードがありません</h1>
      <p class="sub">期限の来たカードも、新しいカードもありません。時間をあけると、覚えたカードが復習として出てきます。</p>
      <a class="btn primary" href="#cards">暗記カードの一覧へ</a></div>`;
    return;
  }
  const total = q.length;
  const res = { again: 0, hard: 0, good: 0 };
  let cur = null, flipped = false;
  const rev = !!C.prefs.reverse;
  const show = () => {
    if (!q.length) return finish();
    cur = q[0]; flipped = false;
    const done = total - q.length;
    const r = rec(cur.id);
    view.innerHTML = `<div class="fade-in study">
      <div class="study-top"><a class="close-x" href="${back}" aria-label="やめる">×</a><div class="progress"><i style="width:${Math.round(100 * done / total)}%"></i></div><span class="muted small">${done} / ${total}</span></div>
      <div class="chip plain">${e(cur.deck.name)}${r ? '' : '・新しいカード'}</div>
      <div class="flash" id="flash" role="button" tabindex="0" aria-label="答えを見る">
        <div class="fface qhtml">${rev ? cur.b : cur.f}</div>
        <div class="fback" hidden>
          ${cur.c ? `<span class="sw big" style="--c:${e(cur.c)}"></span>` : ''}
          <div class="qhtml">${rev ? cur.f : cur.b}</div>
          ${cur.n ? `<p class="muted small">${e(cur.n)}</p>` : ''}
          ${cur.link ? `<a class="linkbtn" href="${e(cur.link)}">${cur.link.startsWith('#law') ? '導出・根拠を見る' : '図鑑で見る'} ›</a>` : ''}
        </div>
        <div class="fhint muted small">タップして答えを見る</div>
      </div>
      <div class="grade" id="grade" hidden>
        <button class="g-again" data-g="again">もう一度<small>すぐにまた出す</small></button>
        <button class="g-hard" data-g="hard">あいまい<small>${days(nextBox(r, 'hard'))}</small></button>
        <button class="g-good" data-g="good">覚えた<small>${days(nextBox(r, 'good'))}</small></button>
      </div>
    </div>`;
    C.math(view);
    const flip = () => {
      if (flipped) return;
      flipped = true;
      $('.fback', view).hidden = false; $('.fhint', view).hidden = true; $('#grade').hidden = false;
      $('#flash').classList.add('open');
      C.play && C.play('tap');
    };
    $('#flash').addEventListener('click', (ev) => { if (!ev.target.closest('a')) flip(); });
    $$('#grade button').forEach(b => b.addEventListener('click', () => grade(b.dataset.g)));
    keyHandler = (ev) => {
      if (ev.target.closest && ev.target.closest('input,textarea')) return;
      if (!flipped && (ev.key === ' ' || ev.key === 'Enter')) { ev.preventDefault(); flip(); }
      else if (flipped && ['1', '2', '3'].includes(ev.key)) grade(['again', 'hard', 'good'][+ev.key - 1]);
    };
  };
  const grade = (g) => {
    const now = Date.now();
    const r = rec(cur.id);
    const box = nextBox(r, g);
    C.state.cards[cur.id] = { box, due: now + INTERVAL[box] * DAY, last: now, n: (r ? r.n : 0) + 1 };
    res[g]++;
    q.shift();
    if (g === 'again') q.splice(Math.min(q.length, 3), 0, cur);   // 少しあとでもう一度
    C.save();
    C.play && C.play(g === 'good' ? 'ok' : 'tap');
    show();
  };
  const finish = () => {
    keyHandler = null;
    view.innerHTML = `<div class="fade-in"><h1>おつかれさま</h1>
      <p class="sub">${total}枚を学習しました（覚えた ${res.good}・あいまい ${res.hard}・もう一度 ${res.again}）。</p>
      <div class="row-btns"><a class="btn primary" href="#cards">暗記カードの一覧へ</a>${deck ? `<a class="btn ghost" href="#deck=${e(deck.id)}">このデッキを見る</a>` : ''}</div></div>`;
  };
  show();
}
let keyHandler = null;
document.addEventListener('keydown', (ev) => { if (keyHandler && location.hash.startsWith('#study')) keyHandler(ev); });

function nextBox(r, g) {
  const b = r ? r.box : 0;
  if (g === 'again') return 0;
  if (g === 'hard') return Math.max(1, b);
  return Math.min(INTERVAL.length - 1, b + 1);
}
function days(box) {
  const d = INTERVAL[box];
  return d === 0 ? 'すぐ' : d < 30 ? `${d}日後` : `${Math.round(d / 30)}か月後`;
}
