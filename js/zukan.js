// 図鑑：CELESTIA のシミュレーションの代わりに、化学の「事実の地図」を操作して確かめるページ
//   #zukan                 トップ（7つのページ）
//   #zukan=ions[/イオン]    沈殿・錯イオン・呈色（反応表・試験管・系統分析・炎色反応）
//   #zukan=reactions[/ID]  製法・反応式（気体の製法・工業的製法・元素ごとの反応）
//   #zukan=organic[/ID]    有機化合物の反応経路（地図・経路さがし）
//   #zukan=tests[/ID]      有機化合物の検出反応（反応ごと・化合物×反応の表）
//   #zukan=polymers[/ID]   高分子の構造・製法・性質
//   #zukan=elements[/記号]  周期表・電子配置・周期性のグラフ
// データは data/zukan.json（tools/build_data.py が content/chem/ から作る）。構造式は data/struct/*.svg
import { sanitizeSVG } from './problem/figures.js';

export const FILE_VERSION = '0.3.0';

let Z = null, C = null;   // 図鑑のデータ、アプリから渡される道具（esc・math・toast・prefs・save・problemsFor・openIds）
const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const e = (s) => C.esc(s ?? '');
const ce = (f) => `$\\ce{${f}}$`;

export const PAGES = [
  { id: 'ions', name: '沈殿・錯イオン・呈色反応', desc: '金属イオンに試薬を加えたときの沈殿の色、過剰で溶けるか、錯イオン、系統分析、炎色反応', icon: 'tube' },
  { id: 'reactions', name: '無機化学の製法・反応式', desc: '気体の製法（捕集・乾燥）、工業的製法、元素ごとの反応式', icon: 'flask' },
  { id: 'organic', name: '有機化学の反応経路', desc: '脂肪族・芳香族の反応の地図。2つの化合物を選ぶと経路をさがす', icon: 'map' },
  { id: 'tests', name: '有機化合物の検出反応', desc: 'ヨードホルム・銀鏡・塩化鉄(III)…。化合物ごとに陽性・陰性を一覧', icon: 'drop' },
  { id: 'polymers', name: '高分子の構造・製法・性質', desc: '単量体と繰り返し単位、重合の種類、結合、性質と用途', icon: 'chain' },
  { id: 'laws', name: '理論化学の定義・法則', desc: '定義・実験から得た法則・導出できる法則（導出を1段ずつ）', icon: 'law', href: '#deck=laws' },
  { id: 'elements', name: '元素の周期性・電子配置', desc: '周期表を電気陰性度・イオン化エネルギーなどで塗り分け、電子配置を表示', icon: 'table' },
];
const ICONS = {
  tube: '<path d="M9 3h6M10 3v13a2 2 0 0 0 4 0V3"/><path d="M10 11h4"/>',
  flask: '<path d="M9 3h6M10 3v6l-5 9a2 2 0 0 0 1.8 3h10.4a2 2 0 0 0 1.8-3l-5-9V3"/><path d="M7.5 14h9"/>',
  map: '<circle cx="6" cy="6" r="2.5"/><circle cx="18" cy="8" r="2.5"/><circle cx="9" cy="18" r="2.5"/><path d="M8.3 7l7.3 1M17 10.3l-6.5 6M6.6 8.4l1.8 7.2"/>',
  drop: '<path d="M12 3c3 4 6 7.5 6 11a6 6 0 0 1-12 0c0-3.5 3-7 6-11z"/><path d="M9.5 14.5a2.5 2.5 0 0 0 2.5 2.5"/>',
  chain: '<rect x="2.5" y="9" width="7" height="6" rx="3"/><rect x="8.5" y="9" width="7" height="6" rx="3"/><rect x="14.5" y="9" width="7" height="6" rx="3"/>',
  law: '<path d="M5 4h14v16H5z"/><path d="M8 8h8M8 12h8M8 16h5"/>',
  table: '<path d="M3 5h4v4H3zM17 5h4v4h-4zM3 10h4v4H3zM8 10h4v4H8zM12 10h4v4h-4zM17 10h4v4h-4zM3 15h18v4H3z"/>',
};
const icon = (k) => `<svg viewBox="0 0 24 24" width="22" height="22" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${ICONS[k]}</svg>`;

export function renderZukan(view, data, arg, ctx) {
  Z = data; C = ctx;
  const [page, ...rest] = String(arg || '').split('/');
  const key = rest.join('/') || null;
  const fn = { ions: ionsPage, reactions: reactionsPage, organic: organicPage, tests: testsPage, polymers: polymersPage, elements: elementsPage }[page];
  if (!fn) return topPage(view);
  C.prefs.page = page; C.save();
  return fn(view, key);
}

// ---------------------------------------------------------------- 共通
function head(title, kicker, sub) {
  return `<a class="back" href="#zukan">‹ 図鑑</a><div class="kicker">図鑑 · ${e(kicker)}</div><h1>${e(title)}</h1>${sub ? `<p class="sub">${sub}</p>` : ''}`;
}
function tabs(id, items, cur) {
  return `<div class="segment zk-tabs" id="${id}">${items.map(([k, label]) => `<button data-t="${k}" class="${k === cur ? 'on' : ''}">${e(label)}</button>`).join('')}</div>`;
}
function bindTabs(id, fn) {
  $$(`#${id} button`).forEach(b => b.addEventListener('click', () => { $$(`#${id} button`).forEach(x => x.classList.toggle('on', x === b)); fn(b.dataset.t); }));
}
function swatch(color, big = false) {
  const hex = Z.ions.colors[color];
  if (color === '無色' || hex === null) return `<span class="sw clear ${big ? 'big' : ''}" title="無色"></span>`;
  return `<span class="sw ${big ? 'big' : ''}" style="--c:${hex || '#999'}" title="${e(color)}"></span>`;
}
function probButton(page, key, label) {
  const ps = C.problemsFor(page, key);
  if (!ps.length) return '';
  return `<button class="btn primary block zk-prob" data-ids="${ps.map(p => p.id).join(',')}">${e(label || 'この内容の問題を解く')}（${ps.length}問）</button>`;
}
function bindProb(root) {
  $$('.zk-prob', root).forEach(b => b.addEventListener('click', () => C.openIds(b.dataset.ids.split(','))));
}
const structCache = new Map();
async function structSVG(name) {
  const src = Z.struct[name];
  if (!src) return null;
  if (!structCache.has(src)) structCache.set(src, fetch(src).then(r => (r.ok ? r.text() : null)).catch(() => null));
  return structCache.get(src);
}
async function mountStructs(root) {
  for (const box of $$('[data-struct]', root)) {
    if (box.dataset.done) continue;
    box.dataset.done = '1';
    const t = await structSVG(box.dataset.struct);
    if (!t) { box.remove(); continue; }
    try { const svg = sanitizeSVG(t); svg.setAttribute('role', 'img'); svg.setAttribute('aria-label', box.dataset.alt || '構造式'); box.appendChild(document.importNode(svg, true)); }
    catch { box.remove(); }
  }
}
function done(root) { C.math(root); mountStructs(root); bindProb(root); }

// ---------------------------------------------------------------- トップ
function topPage(view) {
  view.innerHTML = `<div class="fade-in">
    <div class="kicker">ALCHEMIA · 図鑑</div>
    <h1>図鑑</h1>
    <p class="sub">反応表・地図・周期表を動かして、覚えたことを確かめる場所です。解いた問題の「図鑑で確かめる」からもここに来られます。</p>
    <div class="zk-grid">${PAGES.map(p => `<a class="zk-card" href="${p.href || '#zukan=' + p.id}">${icon(p.icon)}<b>${e(p.name)}</b><span>${e(p.desc)}</span></a>`).join('')}</div>
  </div>`;
}

// ---------------------------------------------------------------- 沈殿・錯イオン・呈色
const SHORT = { Cl: 'HCl', SO4: 'H₂SO₄', CrO4: 'K₂CrO₄', CO3: '(NH₄)₂CO₃', 'H2S-a': 'H₂S 酸性', 'H2S-b': 'H₂S 塩基性', NaOH: 'NaOH 少量', 'NaOH-ex': 'NaOH 過剰',
  NH3: 'NH₃ 少量', 'NH3-ex': 'NH₃ 過剰', SCN: 'KSCN', 'Fe(CN)6-4': 'K₄[Fe(CN)₆]', 'Fe(CN)6-3': 'K₃[Fe(CN)₆]' };
const KIND = { ppt: '沈殿', dissolve: '溶ける', remain: '溶けない', none: '変化なし', solution: '溶液の色が変わる', na: '—' };

function ionsPage(view, key) {
  const I = Z.ions;
  const mode = (C.prefs.ionsMode) || 'table';
  view.innerHTML = `<div class="fade-in zk-ions">${head('沈殿・錯イオン・呈色反応', '無機化学', '金属イオンの水溶液に試薬を加えたときに起こることを、表・試験管・系統分析の3通りで確かめられます。')}
    ${tabs('imode', [['table', '反応表'], ['tube', '試験管'], ['sys', '系統分析'], ['flame', '炎色反応']], mode)}
    <div id="ibody"></div>
    <div id="iprob">${probButton('ions', null, '沈殿・錯イオンの問題を解く')}</div>
  </div>`;
  const draw = (m) => { C.prefs.ionsMode = m; C.save(); ({ table: ionTable, tube: ionTube, sys: ionSystematic, flame: ionFlame })[m]($('#ibody'), key); };
  bindTabs('imode', draw);
  draw(mode);
  done(view);
}

function ionTable(box, key) {
  const I = Z.ions;
  const hide = !!C.prefs.ionsHide;
  box.innerHTML = `
    <div class="zk-toolbar"><label class="sw-toggle"><input type="checkbox" id="ihide" ${hide ? 'checked' : ''}> かくして覚える（タップで表示）</label></div>
    <div class="ion-table-wrap"><table class="ion-table ${hide ? 'hide' : ''}">
      <thead><tr><th class="corner">イオン ＼ 試薬</th>${I.reagents.map(r => `<th title="${e(r.name)}">${e(SHORT[r.id] || r.id)}</th>`).join('')}</tr></thead>
      <tbody>${I.cations.map(c => `<tr class="${key === c.id ? 'sel' : ''}"><th>${ce(c.id)}<small>${e(c.name)}</small></th>${I.reagents.map(r => cell(c, r)).join('')}</tr>`).join('')}</tbody>
    </table></div>
    <div id="idetail" class="panel zk-detail"><span class="muted">マスをタップすると、沈殿の名前・色・反応式が出ます。</span></div>`;
  C.math(box);
  $('#ihide').addEventListener('change', (ev) => { C.prefs.ionsHide = ev.target.checked; C.save(); ionTable(box, key); });
  $$('.ion-table td[data-c]', box).forEach(td => td.addEventListener('click', () => {
    td.classList.add('shown');
    $$('.ion-table td.cur', box).forEach(x => x.classList.remove('cur')); td.classList.add('cur');
    showIon(td.dataset.c, td.dataset.r);
  }));
  if (key && I.results[key]) setTimeout(() => $(`.ion-table tr.sel`, box)?.scrollIntoView({ block: 'center', behavior: 'smooth' }), 50);
}
function cell(c, r) {
  const v = Z.ions.results[c.id][r.id];
  if (v.kind === 'na') return `<td class="na">—</td>`;
  const inner = v.kind === 'ppt' ? `${swatch(v.color)}<span class="f">${v.f ? ce(v.f) : e(v.color)}</span>`
    : v.kind === 'dissolve' ? `<span class="tag t-dis">溶</span><span class="f">${ce(v.f)}</span>`
      : v.kind === 'remain' ? `<span class="tag t-rem">不溶</span>`
        : v.kind === 'solution' ? `${swatch(v.color)}<span class="f">${e(v.color)}色</span>` : `<span class="muted">×</span>`;
  return `<td data-c="${e(c.id)}" data-r="${e(r.id)}" class="k-${v.kind}"><span class="cv">${inner}</span></td>`;
}
function ionResultHTML(cid, rid) {
  const I = Z.ions, c = I.cations.find(x => x.id === cid), r = I.reagents.find(x => x.id === rid), v = I.results[cid][rid];
  let what;
  if (v.kind === 'ppt') what = `${swatch(v.color, true)} <b>${e(v.color)}色の沈殿</b>${v.f ? `：${ce(v.f)}` : ''}${v.name ? `（${e(v.name)}）` : ''}`;
  else if (v.kind === 'dissolve') what = `<b>沈殿が溶ける</b> → ${ce(v.f)}${v.name ? `（${e(v.name)}）` : ''} ${swatch(v.color)} ${e(v.color)}`;
  else if (v.kind === 'remain') what = '<b>沈殿は溶けない</b>';
  else if (v.kind === 'solution') what = `${swatch(v.color, true)} <b>${e(v.color)}色の溶液</b>${v.f ? `：${ce(v.f)}` : ''}（沈殿しない）`;
  else if (v.kind === 'none') what = '<b>沈殿しない</b>（変化なし）';
  else what = '高校では扱わない組み合わせ';
  return `<div class="zk-dhead">${ce(c.id)}（${e(c.name)}）＋ ${e(r.name)}</div>
    <div class="zk-what">${what}</div>
    ${v.eq ? `<div class="zk-eq">${ce(v.eq)}</div>` : ''}
    ${v.note ? `<p class="muted small">${e(v.note)}</p>` : ''}`;
}
function showIon(cid, rid) {
  const d = $('#idetail');
  d.innerHTML = ionResultHTML(cid, rid);
  C.math(d);
}

// 試験管：イオンを選んで試薬を加える。少量で沈殿 → 「さらに過剰に」
function ionTube(box, key) {
  const I = Z.ions;
  let ion = key && I.results[key] ? key : (C.prefs.tubeIon || 'Cu^2+');
  let st = null;   // { reagent, v, ex }
  const draw = () => {
    const c = I.cations.find(x => x.id === ion);
    const base = I.reagents.filter(r => !r.id.endsWith('-ex'));
    const liquid = st ? (st.ex && st.ex.kind === 'dissolve' ? st.ex.color : st.v.kind === 'solution' ? st.v.color : c.aq) : c.aq;
    const ppt = st ? (st.ex && st.ex.kind === 'dissolve' ? null : st.v.kind === 'ppt' ? st.v.color : null) : null;
    const exId = st && (st.reagent === 'NaOH' || st.reagent === 'NH3') && st.v.kind === 'ppt' ? st.reagent + '-ex' : null;
    box.innerHTML = `<div class="tube-lab">
      <div class="tube-left">${tubeSVG(liquid, ppt)}<div class="muted small tube-cap">${ce(c.id)} の水溶液（${e(c.aq)}）</div></div>
      <div class="tube-right">
        <label class="zk-label">イオン</label>
        <select id="tion">${I.cations.map(x => `<option value="${e(x.id)}" ${x.id === ion ? 'selected' : ''}>${e(x.name)}</option>`).join('')}</select>
        <label class="zk-label">加える試薬</label>
        <div class="reag">${base.map(r => `<button class="chipbtn ${st && st.reagent === r.id ? 'on' : ''}" data-r="${e(r.id)}">${e(SHORT[r.id] || r.name)}</button>`).join('')}</div>
        ${exId ? `<button class="btn block" id="tex">さらに ${e(SHORT[exId])}</button>` : ''}
        <button class="linkbtn" id="treset">試験管を新しくする</button>
        <div id="tres" class="zk-detail">${st ? ionResultHTML(ion, st.reagent) + (st.ex ? `<hr>${ionResultHTML(ion, st.reagent + '-ex')}` : '') : '<span class="muted">試薬を選ぶと、試験管の中のようすが変わります。</span>'}</div>
      </div></div>`;
    C.math(box);
    $('#tion').addEventListener('change', (ev) => { ion = ev.target.value; C.prefs.tubeIon = ion; C.save(); st = null; draw(); });
    $$('.reag [data-r]', box).forEach(b => b.addEventListener('click', () => { st = { reagent: b.dataset.r, v: I.results[ion][b.dataset.r] }; draw(); }));
    $('#tex')?.addEventListener('click', () => { st.ex = I.results[ion][exId]; draw(); });
    $('#treset').addEventListener('click', () => { st = null; draw(); });
  };
  draw();
}
function tubeSVG(liquid, ppt) {
  const hex = (c) => (c && Z.ions.colors[c]) || null;
  const lq = hex(liquid);
  const pp = ppt ? hex(ppt) || '#ccc' : null;
  return `<svg class="tube" viewBox="0 0 80 200" aria-hidden="true">
    <defs><clipPath id="tclip"><path d="M18 20 V170 a22 22 0 0 0 44 0 V20 Z"/></clipPath></defs>
    <g clip-path="url(#tclip)">
      <rect x="0" y="80" width="80" height="140" fill="${lq || 'var(--tube-water)'}" opacity="${lq ? 0.55 : 1}"/>
      ${pp ? `<ellipse cx="40" cy="186" rx="26" ry="16" fill="${pp}" stroke="var(--ink-3)" stroke-width="0.6"/>
        ${[...Array(14)].map((_, i) => `<circle cx="${22 + (i * 37) % 36}" cy="${150 + (i * 23) % 30}" r="1.6" fill="${pp}" stroke="var(--ink-3)" stroke-width="0.3"/>`).join('')}` : ''}
    </g>
    <path d="M18 20 V170 a22 22 0 0 0 44 0 V20" fill="none" stroke="currentColor" stroke-width="2"/>
    <path d="M12 20 h56" stroke="currentColor" stroke-width="2" stroke-linecap="round"/>
  </svg>`;
}

// 系統分析：混ぜるイオンを選んで、1段ずつ
const SYS_REAGENT = { 1: 'Cl', 2: 'H2S-a', 3: 'NH3', 4: 'H2S-b', 5: 'CO3' };
function ionSystematic(box) {
  const I = Z.ions;
  let pick = new Set(C.prefs.sysIons || ['Ag+', 'Cu^2+', 'Fe^3+', 'Al^3+', 'Zn^2+', 'Ba^2+', 'Na+']);
  let step = 0;
  const draw = () => {
    const ions = I.cations.filter(c => pick.has(c.id));
    const shown = I.systematic.slice(0, step);
    let remain = ions.map(c => c.id);
    const rows = shown.map(st => {
      const here = remain.filter(i => st.ions.includes(i));
      remain = remain.filter(i => !st.ions.includes(i));
      const rid = SYS_REAGENT[st.step];
      const ppts = here.map(i => {
        const src = i === 'Fe^2+' ? 'Fe^3+' : i;    // Fe²⁺ は希硝酸で Fe³⁺ にしてから
        const v = rid ? I.results[src][rid] : null;
        return v && v.kind === 'ppt' ? `<span class="pchip">${swatch(v.color)}${ce(v.f)}<small>${e(v.color)}</small></span>` : `<span class="pchip">${ce(i)}<small>${st.step === 6 ? '炎色反応で確認' : ''}</small></span>`;
      });
      return `<li><div class="sys-head"><b>${e(st.group)}</b><span>${e(st.reagent)}</span></div>
        <div class="sys-ppt">${ppts.length ? (st.step === 6 ? '<span class="muted small">ろ液に残るイオン：</span>' : '<span class="muted small">沈殿：</span>') + ppts.join('') : '<span class="muted small">この段では何も沈殿しない</span>'}</div>
        <p class="muted small">${e(st.note)}</p></li>`;
    });
    box.innerHTML = `
      <p class="muted small" style="margin:4px 0 8px">混ぜるイオンを選び、「次の操作」で1段ずつ進めます。沈殿はろ過して取り除き、ろ液を次の操作に使います。</p>
      <div class="toggles sys-pick">${I.cations.map(c => `<button class="tog ${pick.has(c.id) ? 'on' : ''}" data-i="${e(c.id)}">${ce(c.id)}</button>`).join('')}</div>
      <ol class="sys-steps">${rows.join('')}</ol>
      ${step < I.systematic.length ? `<button class="btn primary block" id="snext">${step ? '次の操作' : '分析を始める'}</button>` : `<button class="btn block" id="sagain">はじめから</button>`}
      ${step ? `<p class="muted small">ろ液に残っているイオン：${remain.length ? remain.map(ce).join('，') : 'なし'}</p>` : ''}`;
    C.math(box);
    $$('.sys-pick [data-i]', box).forEach(b => b.addEventListener('click', () => {
      const i = b.dataset.i; pick.has(i) ? pick.delete(i) : pick.add(i); C.prefs.sysIons = [...pick]; C.save(); step = 0; draw();
    }));
    $('#snext')?.addEventListener('click', () => { step++; draw(); });
    $('#sagain')?.addEventListener('click', () => { step = 0; draw(); });
  };
  draw();
}
function ionFlame(box) {
  const I = Z.ions;
  const list = I.cations.filter(c => c.flame);
  box.innerHTML = `<div class="flame-grid">${list.map(c => `<div class="flame"><svg viewBox="0 0 40 60" aria-hidden="true"><path d="M20 4c8 12 14 20 14 32a14 14 0 0 1-28 0c0-12 6-20 14-32z" fill="${Z.ions.flameColors[c.flame]}"/><path d="M20 26c4 6 6 10 6 15a6 6 0 0 1-12 0c0-5 2-9 6-15z" fill="#fff" opacity=".45"/></svg>
    <b>${ce(c.id.replace(/\^?\d*[+-]$/, ''))}</b><span>${e(c.flame)}</span></div>`).join('')}</div>
    <p class="muted small">覚え方の例：「リアカー（Li 赤）なき（Na 黄）K 村（K 赤紫）、動力（Cu 青緑）に借りる（Ca 橙赤）と するも（Sr 紅）くれない、馬力（Ba 黄緑）で」</p>`;
  C.math(box);
}

// ---------------------------------------------------------------- 製法・反応式
function reactionsPage(view, key) {
  const R = Z.reactions;
  const all = R.reactions;
  const cur = key && all.find(r => r.id === key);
  const mode = cur ? (cur.cat === 'gas' ? 'gas' : cur.cat === 'industry' ? 'industry' : 'reaction') : (C.prefs.rxMode || 'gas');
  view.innerHTML = `<div class="fade-in zk-rx">${head('無機化学の製法・反応式', '無機化学', '反応式の右辺をかくして、思い出してからタップで確かめられます。')}
    ${tabs('rmode', [['gas', '気体の製法'], ['industry', '工業的製法'], ['reaction', '元素の反応'], ['gtable', '気体の性質']], mode)}
    <div class="zk-toolbar"><label class="sw-toggle"><input type="checkbox" id="rhide" ${C.prefs.rxHide ? 'checked' : ''}> 右辺をかくす</label></div>
    <div id="rbody"></div>
    ${probButton('reactions', null, '製法・反応式の問題を解く')}
  </div>`;
  const draw = (m) => {
    C.prefs.rxMode = m; C.save();
    const body = $('#rbody');
    const hide = !!C.prefs.rxHide;
    if (m === 'gtable') { body.innerHTML = gasTable(); C.math(body); return; }
    const list = all.filter(r => (m === 'gas' ? r.cat === 'gas' : m === 'industry' ? r.cat === 'industry' : r.cat === 'reaction'));
    const units = [...new Set(list.map(r => r.unit))];
    body.innerHTML = units.map(u => `${m === 'reaction' ? `<h2>${e(unitName(u))}</h2>` : ''}<div class="rx-list">${list.filter(r => r.unit === u).map(r => rxCard(r, hide, key)).join('')}</div>`).join('');
    C.math(body);
    $$('.rhs', body).forEach(x => x.addEventListener('click', () => x.classList.add('shown')));
    if (cur) setTimeout(() => $(`[data-rx="${CSS.escape(cur.id)}"]`, body)?.scrollIntoView({ block: 'center', behavior: 'smooth' }), 60);
  };
  bindTabs('rmode', draw);
  $('#rhide').addEventListener('change', (ev) => { C.prefs.rxHide = ev.target.checked; C.save(); draw(C.prefs.rxMode); });
  draw(mode);
  done(view);
}
function unitName(u) {
  return { 'INO-01': 'ハロゲン', 'INO-02': '酸素・硫黄', 'INO-03': '窒素・リン', 'INO-04': '炭素・ケイ素', 'INO-05': 'アルカリ金属・アルカリ土類金属',
    'INO-06': 'アルミニウム・亜鉛・スズ', 'INO-07': '遷移元素', 'INO-09': '気体の製法', 'ORG-03': '有機化合物' }[u] || u;
}
function splitEq(eq) {
  const m = eq.match(/\s*(<=>|->|<-)\s*/);
  if (!m) return [eq, '', ''];
  return [eq.slice(0, m.index), m[1], eq.slice(m.index + m[0].length)];
}
function eqHTML(eq, hide) {
  const [l, arrow, r] = splitEq(eq);
  if (!arrow) return ce(eq);
  return `<span class="eqline">${ce(l + ' ' + arrow)}<span class="rhs ${hide ? 'hid' : ''}">${ce(r)}</span></span>`;
}
function rxCard(r, hide, key) {
  const T = Z.reactions.types;
  const g = r.gas && Z.reactions.gases[r.gas];
  const steps = r.steps ? `<ol class="rx-steps">${r.steps.map(s => `<li>${eqHTML(s.eq, hide)}${s.note ? `<small class="muted">${e(s.note)}</small>` : ''}</li>`).join('')}</ol>${r.total ? `<div class="rx-total"><span class="muted small">全体</span>${eqHTML(r.total, hide)}</div>` : ''}` : `<div class="zk-eq">${eqHTML(r.eq, hide)}</div>`;
  return `<div class="rx-card ${key === r.id ? 'sel' : ''}" data-rx="${e(r.id)}">
    <div class="rx-head"><b>${e(r.title)}</b>${r.type ? `<span class="chip plain">${e(T[r.type] || r.type)}</span>` : ''}</div>
    ${r.cond ? `<div class="muted small">条件：${e(r.cond)}</div>` : ''}
    ${steps}
    ${g ? `<div class="gas-line">${g.collect ? `<span class="chip plain">捕集：${e(g.collect)}</span>` : ''}${g.dry.length ? `<span class="chip plain">乾燥：${e(g.dry.join('・'))}</span>` : ''}<span class="chip plain">${e(g.color)}・${e(g.smell)}・${e(g.nature)}</span></div>` : ''}
    ${r.note ? `<p class="muted small">${e(r.note)}</p>` : ''}
  </div>`;
}
function gasTable() {
  const G = Z.reactions.gases;
  return `<div class="ion-table-wrap"><table class="gas-table">
    <thead><tr><th>気体</th><th>色</th><th>におい</th><th>水溶液</th><th>水への溶け方</th><th>捕集</th><th>乾燥剤</th></tr></thead>
    <tbody>${Object.entries(G).map(([f, g]) => `<tr><th>${ce(f)}<small>${e(g.name)}</small></th><td>${e(g.color)}</td><td>${e(g.smell)}</td><td>${e(g.nature)}</td><td>${e(g.sol)}</td><td>${e(g.collect || '—')}</td><td>${e(g.dry.join('・') || '—')}</td></tr>`).join('')}</tbody>
  </table></div><p class="muted small">乾燥剤：酸性の気体に塩基性のソーダ石灰は使えない。NH₃ は濃硫酸（酸性）・塩化カルシウム（反応する）が使えない。H₂S は還元性があるので濃硫酸が使えない。</p>`;
}

// ---------------------------------------------------------------- 有機化学の反応経路
function organicPage(view, key) {
  const O = Z.organic;
  const comp = Object.fromEntries(O.compounds.map(c => [c.id, c]));
  let map = C.prefs.orgMap || 'aliphatic';
  if (key && comp[key] && !O.maps[map].pos[key]) map = O.maps.aromatic.pos[key] ? 'aromatic' : 'aliphatic';
  let sel = key && comp[key] ? key : null, from = null, to = null, path = [];
  view.innerHTML = `<div class="fade-in zk-org">${head('有機化学の反応経路', '有機化学', '化合物をタップすると、構造式とその化合物に出入りする反応が出ます。「経路さがし」で2つを選ぶと、つなぐ反応を順にたどります。')}
    ${tabs('omap', [['aliphatic', '脂肪族'], ['aromatic', '芳香族']], map)}
    <div class="zk-toolbar"><label class="sw-toggle"><input type="checkbox" id="opath"> 経路さがし</label><span class="muted small" id="ohint"></span></div>
    <div class="org-map" id="omapbox"></div>
    <div id="odetail" class="panel zk-detail"><span class="muted">化合物をタップしてください。</span></div>
    ${probButton('organic', null, '反応経路の問題を解く')}
  </div>`;
  const edgesOf = (m) => O.reactions.filter(r => (r.map === m || r.map === 'both') && O.maps[m].pos[r.from] && O.maps[m].pos[r.to]);
  const drawMap = () => {
    const P = O.maps[map].pos, E = edgesOf(map);
    const xs = Object.values(P).map(p => p[0]), ys = Object.values(P).map(p => p[1]);
    const W = 120, H = 74, pad = 60;
    const X = (x) => pad + x * W / 2 * 1.25, Y = (y) => pad + y * H;
    const w = X(Math.max(...xs)) + pad, h = Y(Math.max(...ys)) + pad;
    const onPath = new Set(path.map(r => r.id));
    const pathNodes = new Set(path.flatMap(r => [r.from, r.to]));
    const edgeSVG = E.map(r => {
      const [x1, y1] = [X(P[r.from][0]), Y(P[r.from][1])], [x2, y2] = [X(P[r.to][0]), Y(P[r.to][1])];
      const dx = x2 - x1, dy = y2 - y1, L = Math.hypot(dx, dy) || 1;
      const rev = E.some(q => q.from === r.to && q.to === r.from);
      const off = rev ? 6 : 0, nx = -dy / L * off, ny = dx / L * off;
      // 箱の端で止める（箱 104×36）
      const t = (bw, bh) => Math.min(Math.abs((bw / 2) / (dx || 1e-9)), Math.abs((bh / 2) / (dy || 1e-9)));
      const s1 = t(108, 40), s2 = t(112, 44);
      const ax = x1 + dx * s1 + nx, ay = y1 + dy * s1 + ny, bx = x2 - dx * s2 + nx, by = y2 - dy * s2 + ny;
      const on = onPath.has(r.id), dim = (sel && r.from !== sel && r.to !== sel && !on) || (path.length && !on);
      return `<g class="oedge ${on ? 'on' : ''} ${dim ? 'dim' : ''}" data-r="${r.id}"><line x1="${ax}" y1="${ay}" x2="${bx}" y2="${by}" marker-end="url(#oarr${on ? 'on' : ''})"/>
        <text x="${(ax + bx) / 2 + nx * 1.5}" y="${(ay + by) / 2 + ny * 1.5 - 3}">${e(shortReagent(r.reagent))}</text></g>`;
    }).join('');
    const nodeSVG = Object.entries(P).map(([id, [x, y]]) => {
      const c = comp[id];
      const on = id === sel || id === from || id === to || pathNodes.has(id);
      const dim = (sel && id !== sel && !E.some(r => (r.from === sel && r.to === id) || (r.to === sel && r.from === id))) || (path.length && !pathNodes.has(id));
      return `<g class="onode ${on ? 'on' : ''} ${dim ? 'dim' : ''} cls-${clsKey(c.cls)}" data-c="${id}" tabindex="0" role="button" aria-label="${e(c.name)}">
        <rect x="${X(x) - 54}" y="${Y(y) - 19}" width="108" height="38" rx="9"/>
        <text x="${X(x)}" y="${Y(y) - 2}" class="nm">${e(c.name.length > 9 ? c.name.slice(0, 9) + '…' : c.name)}</text>
        <text x="${X(x)}" y="${Y(y) + 12}" class="cf">${e(plainFormula(c.cf))}</text></g>`;
    }).join('');
    $('#omapbox').innerHTML = `<svg viewBox="0 0 ${w} ${h}" style="min-width:${Math.round(w * 0.82)}px">
      <defs>
        <marker id="oarr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" class="ah"/></marker>
        <marker id="oarron" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" class="ah on"/></marker>
      </defs>${edgeSVG}${nodeSVG}</svg>`;
    $$('.onode', $('#omapbox')).forEach(g => {
      const act = () => clickNode(g.dataset.c);
      g.addEventListener('click', act);
      g.addEventListener('keydown', (ev) => { if (ev.key === 'Enter') act(); });
    });
    $$('.oedge', $('#omapbox')).forEach(g => g.addEventListener('click', () => showEdge(O.reactions.find(r => r.id === g.dataset.r))));
  };
  const pathMode = () => $('#opath').checked;
  const hint = () => { $('#ohint').textContent = pathMode() ? (from ? (to ? '' : '目的の化合物を選ぶ') : '出発の化合物を選ぶ') : ''; };
  function clickNode(id) {
    if (pathMode()) {
      if (!from || (from && to)) { from = id; to = null; path = []; sel = null; }
      else { to = id; path = findPath(edgesOf(map), from, to); }
      hint(); drawMap(); showPath(); return;
    }
    sel = sel === id ? null : id; drawMap(); if (sel) showCompound(sel); else $('#odetail').innerHTML = '<span class="muted">化合物をタップしてください。</span>';
  }
  function showCompound(id) {
    const c = comp[id];
    const out = O.reactions.filter(r => r.from === id), inn = O.reactions.filter(r => r.to === id);
    const d = $('#odetail');
    d.innerHTML = `<div class="org-card">
      <div class="struct" data-struct="${e(id)}" data-alt="${e(c.name)}の構造式"></div>
      <div><div class="zk-dhead">${e(c.name)}</div>${c.cf ? `<div>${ce(c.cf)}</div>` : ''}<span class="chip plain">${e(c.cls)}</span>${c.note ? `<p class="muted small">${e(c.note)}</p>` : ''}</div></div>
      ${inn.length ? `<h3>この化合物ができる反応</h3>${inn.map(r => edgeRow(r, comp)).join('')}` : ''}
      ${out.length ? `<h3>この化合物から進む反応</h3>${out.map(r => edgeRow(r, comp)).join('')}` : ''}
      ${probButton('organic', id, `${c.name}の問題`)}`;
    done(d);
    $$('[data-goto]', d).forEach(a => a.addEventListener('click', () => { const t = a.dataset.goto; if (O.maps[map].pos[t]) { sel = t; drawMap(); showCompound(t); } }));
  }
  function showEdge(r) {
    const d = $('#odetail');
    d.innerHTML = edgeRow(r, comp, true);
    done(d);
  }
  function showPath() {
    const d = $('#odetail');
    if (!from) { d.innerHTML = '<span class="muted">出発の化合物をタップしてください。</span>'; return; }
    if (!to) { d.innerHTML = `出発：<b>${e(comp[from].name)}</b>　→　目的の化合物をタップ`; return; }
    if (!path.length) { d.innerHTML = `<b>${e(comp[from].name)}</b> から <b>${e(comp[to].name)}</b> へは、この地図の反応だけではたどれません。`; return; }
    d.innerHTML = `<div class="zk-dhead">${e(comp[from].name)} → ${e(comp[to].name)}（${path.length}段階）</div><ol class="path-steps">${path.map(r => `<li>${edgeRow(r, comp)}</li>`).join('')}</ol>`;
    done(d);
  }
  bindTabs('omap', (m) => { map = m; C.prefs.orgMap = m; C.save(); sel = null; from = to = null; path = []; drawMap(); $('#odetail').innerHTML = '<span class="muted">化合物をタップしてください。</span>'; });
  $('#opath').addEventListener('change', () => { from = to = null; path = []; sel = null; hint(); drawMap(); showPath(); if (!pathMode()) $('#odetail').innerHTML = '<span class="muted">化合物をタップしてください。</span>'; });
  drawMap();
  if (sel) { showCompound(sel); setTimeout(() => $(`.onode[data-c="${sel}"]`)?.scrollIntoView({ block: 'nearest', inline: 'center' }), 50); }
  done(view);
}
function edgeRow(r, comp, full = false) {
  return `<div class="edge-row"><div><a class="linkbtn" data-goto="${e(r.from)}">${e(comp[r.from].name)}</a> <span class="arrow">→</span> <a class="linkbtn" data-goto="${e(r.to)}">${e(comp[r.to].name)}</a></div>
    <div class="muted small"><span class="chip plain">${e(r.type)}</span> ${e(r.reagent)}</div>
    ${r.eq ? `<div class="zk-eq">${ce(r.eq)}</div>` : ''}${r.note ? `<p class="muted small">${e(r.note)}</p>` : ''}</div>`;
}
function findPath(E, a, b) {
  const prev = new Map([[a, null]]), q = [a];
  while (q.length) {
    const u = q.shift();
    if (u === b) break;
    for (const r of E.filter(x => x.from === u)) if (!prev.has(r.to)) { prev.set(r.to, r); q.push(r.to); }
  }
  if (!prev.has(b) || a === b) return [];
  const out = [];
  for (let v = b; prev.get(v); v = prev.get(v).from) out.unshift(prev.get(v));
  return out;
}
const shortReagent = (s) => { const t = String(s).replace(/（.*?）/g, '').replace(/、.*$/, ''); return t.length > 12 ? t.slice(0, 11) + '…' : t; };
const plainFormula = (f) => {
  if (!f) return '';
  const sub = '₀₁₂₃₄₅₆₇₈₉';
  return f.replace(/_\{n\}/g, 'ₙ').replace(/#/g, '≡').replace(/(?<=[A-Za-z)\]])(\d+)/g, (m) => [...m].map(d => sub[d]).join(''));
};
const clsKey = (c) => ({ 'アルカン': 'hc', 'アルケン': 'hc', 'アルキン': 'hc', '芳香族炭化水素': 'hc', 'シクロアルカン': 'hc', 'アルコール': 'oh', 'フェノール類': 'oh', 'アルデヒド': 'co', 'ケトン': 'co',
  'カルボン酸': 'cooh', 'エステル': 'est', '高分子': 'poly' }[c] || 'etc');

// ---------------------------------------------------------------- 検出反応
function testsPage(view, key) {
  const T = Z.tests.tests;
  const comp = Object.fromEntries(Z.organic.compounds.map(c => [c.id, c]));
  const mode = key ? 'list' : (C.prefs.testMode || 'list');
  view.innerHTML = `<div class="fade-in zk-tests">${head('有機化合物の検出反応', '有機化学', '反応ごとに、陽性になる化合物と、その理由（どの構造に反応するか）を確かめられます。')}
    ${tabs('tmode', [['list', '反応ごと'], ['matrix', '化合物 × 反応']], mode)}
    <div id="tbody"></div>
    ${probButton('tests', null, '検出反応の問題を解く')}
  </div>`;
  const draw = (m) => {
    C.prefs.testMode = m; C.save();
    const body = $('#tbody');
    if (m === 'list') {
      body.innerHTML = T.map(t => `<div class="test-card ${key === t.id ? 'sel' : ''}" data-t="${e(t.id)}">
        <div class="rx-head"><b>${e(t.name)}</b></div>
        <div class="muted small">操作：${e(t.reagent)}</div>
        <div><span class="muted small">陽性：</span><b>${e(t.obs)}</b></div>
        <div class="muted small">反応する構造：${e(t.target)}</div>
        <div class="tchips">${Object.entries(t.pos).map(([k, n]) => `<button class="chipbtn pos" data-c="${e(k)}" data-n="${e(n)}" data-t="${e(t.id)}">＋ ${e(comp[k].name)}</button>`).join('')}
          ${Object.entries(t.tri).map(([k, n]) => `<button class="chipbtn tri" data-c="${e(k)}" data-n="${e(n)}" data-t="${e(t.id)}">△ ${e(comp[k].name)}</button>`).join('')}
          ${Object.entries(t.neg_note).map(([k, n]) => `<button class="chipbtn neg" data-c="${e(k)}" data-n="${e(n)}" data-t="${e(t.id)}">− ${e(comp[k].name)}</button>`).join('')}</div>
        <div class="tdet"></div></div>`).join('');
      $$('.tchips .chipbtn', body).forEach(b => b.addEventListener('click', () => {
        const box = b.closest('.test-card').querySelector('.tdet'), c = comp[b.dataset.c];
        box.innerHTML = `<div class="org-card"><div class="struct" data-struct="${e(c.id)}"></div><div><b>${e(c.name)}</b>${c.cf ? `<div>${ce(c.cf)}</div>` : ''}${b.dataset.n ? `<p class="small">${e(b.dataset.n)}</p>` : ''}</div></div>`;
        done(box);
      }));
      if (key) setTimeout(() => $(`.test-card[data-t="${CSS.escape(key)}"]`, body)?.scrollIntoView({ block: 'center', behavior: 'smooth' }), 60);
    } else {
      const cs = [...new Set(T.flatMap(t => t.scope))];
      body.innerHTML = `<div class="zk-toolbar"><label class="sw-toggle"><input type="checkbox" id="thide" ${C.prefs.testHide ? 'checked' : ''}> かくして覚える</label></div>
        <div class="ion-table-wrap"><table class="test-table ${C.prefs.testHide ? 'hide' : ''}"><thead><tr><th class="corner">化合物 ＼ 反応</th>${T.map(t => `<th>${e(t.name.replace(/水溶液|による呈色|との反応/g, ''))}</th>`).join('')}</tr></thead>
        <tbody>${cs.map(k => `<tr><th>${e(comp[k].name)}</th>${T.map(t => {
          if (!t.scope.includes(k)) return '<td class="na">·</td>';
          const v = t.pos[k] !== undefined ? ['pos', '＋', t.pos[k]] : t.tri[k] !== undefined ? ['tri', '△', t.tri[k]] : ['neg', '−', t.neg_note[k] || ''];
          return `<td class="tv ${v[0]}" title="${e(v[2])}" data-n="${e(v[2])}" data-c="${e(k)}" data-t="${e(t.name)}"><span class="cv">${v[1]}</span></td>`;
        }).join('')}</tr>`).join('')}</tbody></table></div>
        <div id="tdetail" class="panel zk-detail"><span class="muted">＋ 陽性　− 陰性　△ 条件による・注意（タップで理由）</span></div>`;
      $('#thide').addEventListener('change', (ev) => { C.prefs.testHide = ev.target.checked; C.save(); draw('matrix'); });
      $$('.test-table td.tv', body).forEach(td => td.addEventListener('click', () => {
        td.classList.add('shown');
        const d = $('#tdetail');
        d.innerHTML = `<b>${e(comp[td.dataset.c].name)}</b> × ${e(td.dataset.t)}：<b>${e(td.querySelector('.cv').textContent)}</b>${td.dataset.n ? `<p class="small">${e(td.dataset.n)}</p>` : ''}`;
      }));
    }
    C.math(body);
  };
  bindTabs('tmode', draw);
  draw(mode);
  done(view);
}

// ---------------------------------------------------------------- 高分子
function polymersPage(view, key) {
  const P = Z.polymers.polymers;
  const masses = Object.fromEntries(Z.elements.elements.filter(x => x.exam).map(x => [x.sym, x.exam]));
  const groups = [...new Set(P.map(p => p.cls))];
  view.innerHTML = `<div class="fade-in zk-poly">${head('高分子の構造・製法・性質', '高分子化合物', '単量体から繰り返し単位ができる様子と、結合・性質・用途をまとめています。繰り返し単位の式量（高校の原子量で計算）も出します。')}
    <div class="toggles" id="pgrp">${groups.map(g => `<button class="tog" data-g="${e(g)}">${e(g)}</button>`).join('')}</div>
    <div id="pbody"></div>
    ${probButton('polymers', null, '高分子の問題を解く')}
  </div>`;
  const draw = (g) => {
    const list = P.filter(p => !g || p.cls === g);
    $('#pbody').innerHTML = list.map(p => {
      const fm = p.rf ? formulaMass(p.rf, masses) : null;
      return `<div class="poly-card ${key === p.id ? 'sel' : ''}" data-p="${e(p.id)}">
        <div class="rx-head"><b>${e(p.name)}</b>${p.abbr ? `<span class="chip plain">${e(p.abbr)}</span>` : ''}<span class="chip">${e(p.poly)}</span>${p.heat ? `<span class="chip plain">${e(p.heat)}</span>` : ''}</div>
        <div class="poly-flow">
          <div class="mono">${p.monomers.map((m, i) => `<div class="mono-item">${m.smiles ? `<div class="struct sm" data-struct="mono-${e(p.id)}-${i + 1}" data-alt="${e(m.name)}"></div>` : ''}<small>${e(m.name)}</small></div>`).join('<span class="plus">＋</span>')}</div>
          <span class="arrow">→</span>
          <div class="rep">${p.repeat ? `<div class="struct" data-struct="poly-${e(p.id)}" data-alt="${e(p.name)}の繰り返し単位"></div>` : '<span class="muted small">（網目状・多様な構造）</span>'}
            ${p.rf ? `<small>繰り返し単位 ${ce(p.rf)}${fm ? `　式量 ${fm}` : ''}</small>` : ''}</div>
        </div>
        ${p.bond ? `<div class="small">結合：<b>${e(p.bond)}</b></div>` : ''}
        ${p.props ? `<p class="small">${e(p.props)}</p>` : ''}
        ${p.uses ? `<div class="muted small">用途：${e(p.uses)}</div>` : ''}
      </div>`;
    }).join('');
    done($('#pbody'));
    if (key) setTimeout(() => $(`.poly-card[data-p="${CSS.escape(key)}"]`)?.scrollIntoView({ block: 'center', behavior: 'smooth' }), 80);
  };
  let curG = key ? (P.find(p => p.id === key) || {}).cls : null;
  $$('#pgrp .tog').forEach(b => {
    b.classList.toggle('on', b.dataset.g === curG);
    b.addEventListener('click', () => { curG = curG === b.dataset.g ? null : b.dataset.g; $$('#pgrp .tog').forEach(x => x.classList.toggle('on', x.dataset.g === curG)); draw(curG); });
  });
  draw(curG);
  done(view);
}
function formulaMass(f, masses) {
  // 簡単な化学式（かっこなし・水和なし）の式量
  let sum = 0;
  for (const [, el, n] of f.matchAll(/([A-Z][a-z]?)(\d*)/g)) {
    if (!(el in masses)) return null;
    sum += masses[el] * (n ? +n : 1);
  }
  return Math.round(sum * 10) / 10;
}

// ---------------------------------------------------------------- 周期表・電子配置
const PROPS = {
  cat: { name: '分類' },
  en: { name: '電気陰性度', unit: '', fmt: (v) => v.toFixed(2) },
  ie1: { name: '第1イオン化エネルギー', unit: 'kJ/mol', fmt: (v) => String(v) },
  radius: { name: '原子半径', unit: 'pm', fmt: (v) => String(Math.round(v)) },
  state: { name: '常温の状態' },
  block: { name: 'ブロック' },
};
const CAT_COLOR = { 'アルカリ金属': '#e8836b', 'アルカリ土類金属': '#ecb067', '遷移元素': '#d8c36a', '典型金属': '#9cc2a6', '半金属': '#88c0c4', '非金属': '#8fb2e2',
  'ハロゲン': '#b7a1e0', '貴ガス': '#d59ad0', 'ランタノイド': '#c9b79c', 'アクチノイド': '#b9a28d' };
function elementsPage(view, key) {
  const E = Z.elements.elements;
  const by = Object.fromEntries(E.map(x => [x.sym, x]));
  let prop = C.prefs.elProp || 'cat';
  let sel = key && by[key] ? key : (C.prefs.elSel || 'Fe');
  view.innerHTML = `<div class="fade-in zk-el">${head('元素の周期性・電子配置', '理論化学', '元素をタップすると電子配置が出ます。性質を選ぶと周期表が塗り分けられ、周期性（右上ほど電気陰性度が大きい、など）が見えます。')}
    <div class="zk-toolbar"><label class="zk-label" for="eprop">塗り分け</label><select id="eprop">${Object.entries(PROPS).map(([k, p]) => `<option value="${k}" ${k === prop ? 'selected' : ''}>${e(p.name)}</option>`).join('')}</select></div>
    <div class="ptable-wrap"><div class="ptable" id="ptable"></div></div>
    <div id="elegend" class="muted small"></div>
    <div id="edetail" class="panel zk-detail"></div>
    <h2>原子番号と第1イオン化エネルギー</h2>
    <div class="ie-chart" id="iechart"></div>
    ${probButton('elements', null, '周期性・電子配置の問題を解く')}
  </div>`;
  const vals = (k) => E.map(x => x[k]).filter(v => typeof v === 'number');
  const colorOf = (x) => {
    if (prop === 'cat') return CAT_COLOR[x.cat] || '#ccc';
    if (prop === 'state') return { '気体': '#8fb2e2', '液体': '#88c0c4', '固体': '#d9d4c7' }[x.state] || '#eee';
    if (prop === 'block') return { s: '#e8836b', p: '#8fb2e2', d: '#d8c36a', f: '#c9b79c' }[x.block];
    const v = x[prop];
    if (typeof v !== 'number') return 'var(--line)';
    const vs = vals(prop), lo = Math.min(...vs), hi = Math.max(...vs);
    const t = (v - lo) / (hi - lo || 1);
    return `hsl(${220 - 200 * t} 70% ${78 - 30 * t}%)`;
  };
  const drawTable = () => {
    const cells = [];
    for (const x of E) {
      let col = x.group, row = x.period;
      if (x.cat === 'ランタノイド' || x.cat === 'アクチノイド') { row = x.cat === 'ランタノイド' ? 9 : 10; col = 3 + (x.z - (x.cat === 'ランタノイド' ? 57 : 89)); }
      if (!col) continue;
      const v = PROPS[prop].fmt && typeof x[prop] === 'number' ? PROPS[prop].fmt(x[prop]) : '';
      cells.push(`<button class="el ${x.sym === sel ? 'on' : ''}" data-s="${x.sym}" style="grid-column:${col};grid-row:${row};--bg:${colorOf(x)}" title="${e(x.ja)}">
        <small>${x.z}</small><b>${x.sym}</b><i>${v || e(x.ja)}</i></button>`);
    }
    cells.push('<div class="el ph" style="grid-column:3;grid-row:6">57–71</div><div class="el ph" style="grid-column:3;grid-row:7">89–103</div>');
    for (let g = 1; g <= 18; g++) cells.push(`<div class="gnum" style="grid-column:${g};grid-row:1 / span 1;align-self:start" aria-hidden="true"></div>`);
    $('#ptable').innerHTML = cells.join('');
    $$('#ptable .el[data-s]').forEach(b => b.addEventListener('click', () => { sel = b.dataset.s; C.prefs.elSel = sel; C.save(); drawTable(); showEl(); drawChart(); }));
    const lg = $('#elegend');
    if (prop === 'cat') lg.innerHTML = Object.entries(CAT_COLOR).map(([k, c]) => `<span class="lg"><i style="background:${c}"></i>${e(k)}</span>`).join('');
    else if (prop === 'state' || prop === 'block') lg.innerHTML = '';
    else { const vs = vals(prop); lg.innerHTML = `<span class="lg grad"><i></i>${Math.min(...vs)} 〜 ${Math.max(...vs)} ${e(PROPS[prop].unit || '')}（色が赤いほど大きい）</span>`; }
  };
  const showEl = () => {
    const x = by[sel];
    const valence = x.group && x.group <= 2 ? x.group : x.group >= 13 ? (x.group === 18 ? 0 : x.group - 10) : null;
    const d = $('#edetail');
    d.innerHTML = `<div class="el-head"><div class="el-big" style="--bg:${CAT_COLOR[x.cat] || '#ccc'}"><small>${x.z}</small><b>${x.sym}</b></div>
      <div><div class="zk-dhead">${e(x.ja)}</div><div class="muted small">${e(x.cat)}・${x.period}周期${x.group ? `・${x.group}族` : ''}・${x.block}ブロック${x.state ? `・常温で${x.state}` : ''}</div>
      <div class="small">原子量 ${x.mass}${x.exam ? `<span class="muted">（計算では ${x.exam} を使う）</span>` : ''}</div></div></div>
      <div class="econf"><span class="muted small">電子配置</span> <b>${confHTML(x.conf)}</b></div>
      <div class="shells">${x.shells.map((n, i) => `<span class="shell"><b>${'KLMNOPQ'[i]}</b>${n}</span>`).join('')}${valence !== null ? `<span class="muted small">価電子 ${valence}個</span>` : ''}</div>
      ${orbitalBoxes(x)}
      <div class="el-vals">
        <span>電気陰性度 <b>${x.en ?? '—'}</b></span><span>第1イオン化エネルギー <b>${x.ie1 ?? '—'}</b> kJ/mol</span><span>原子半径 <b>${x.radius ?? '—'}</b> pm</span>
      </div>
      ${['Cr', 'Cu'].includes(x.sym) ? '<p class="muted small">Cr・Cu は 4s 軌道の電子が1個になる（3d 軌道が半分・全部うまると安定）。</p>' : ''}`;
  };
  const drawChart = () => {
    const pts = E.filter(x => x.z <= 54 && x.ie1);
    const W = 640, H = 220, L = 40, B = 24, maxv = 2500;
    const X = (z) => L + (z - 1) / 53 * (W - L - 10), Y = (v) => H - B - v / maxv * (H - B - 10);
    const s = by[sel];
    $('#iechart').innerHTML = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="原子番号 1〜54 の第1イオン化エネルギー">
      ${[0, 500, 1000, 1500, 2000, 2500].map(v => `<line x1="${L}" x2="${W - 10}" y1="${Y(v)}" y2="${Y(v)}" class="grid"/><text x="${L - 4}" y="${Y(v) + 4}" class="ax" text-anchor="end">${v}</text>`).join('')}
      ${[2, 10, 18, 36, 54].map(z => `<line x1="${X(z)}" x2="${X(z)}" y1="10" y2="${H - B}" class="grid noble"/><text x="${X(z)}" y="${H - 6}" class="ax" text-anchor="middle">${by[E[z - 1].sym].sym}</text>`).join('')}
      <polyline points="${pts.map(x => `${X(x.z)},${Y(x.ie1)}`).join(' ')}" class="line"/>
      ${pts.map(x => `<circle cx="${X(x.z)}" cy="${Y(x.ie1)}" r="${x.sym === sel ? 5 : 2.6}" class="pt ${x.sym === sel ? 'on' : ''} ${x.group === 1 ? 'alk' : ''}" data-s="${x.sym}"><title>${x.ja} ${x.ie1} kJ/mol</title></circle>`).join('')}
      ${s && s.z <= 54 && s.ie1 ? `<text x="${X(s.z) + 7}" y="${Y(s.ie1) - 7}" class="lab">${s.sym} ${s.ie1}</text>` : ''}
    </svg><p class="muted small">貴ガス（縦線）で極大、アルカリ金属で極小。同じ周期でも 2族→13族、15族→16族で小さくなる。</p>`;
    $$('#iechart .pt').forEach(c => c.addEventListener('click', () => { sel = c.dataset.s; C.prefs.elSel = sel; C.save(); drawTable(); showEl(); drawChart(); }));
  };
  $('#eprop').addEventListener('change', (ev) => { prop = ev.target.value; C.prefs.elProp = prop; C.save(); drawTable(); });
  drawTable(); showEl(); drawChart();
  done(view);
}
function confHTML(conf) {
  return e(conf).replace(/(\d)([spdf])(\d+)/g, '$1$2<sup>$3</sup>');
}
function orbitalBoxes(x) {
  // 最外殻と、その内側の d 軌道（遷移元素）の箱。フントの規則で上向きから入れる
  const orb = x.orbitals;
  const last = Math.max(...orb.map(([o]) => +o[0]));
  const show = orb.filter(([o]) => +o[0] === last || (o[1] === 'd' && +o[0] === last - 1 && x.block === 'd') || (o[1] === 'f' && x.block === 'f'));
  if (!show.length) return '';
  const boxes = (type, n) => {
    const k = { s: 1, p: 3, d: 5, f: 7 }[type];
    const cnt = Array(k).fill(0);
    for (let i = 0; i < n; i++) cnt[i % k]++;
    return cnt.map(c => `<span class="ob">${c >= 1 ? '↑' : ''}${c >= 2 ? '↓' : ''}</span>`).join('');
  };
  return `<div class="orbitals">${show.map(([o, n]) => `<div class="orb"><div class="boxes">${boxes(o[1], n)}</div><small>${o}</small></div>`).join('')}</div>`;
}
