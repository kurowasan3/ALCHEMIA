// ブラウザでの通し操作（Playwright・Chromium、iPhone 幅 390×844）
//   node tests/e2e.mjs
// 5つの出題形式の採点、図鑑の7ページ、一部のデータだけ古いときの警告、オフライン、同期、
// Lemma v0.27 に合わせた学習の流れ（今日の学習・前回の続き・ブックマーク・解答の形式・類題・バックアップ・消去）を確かめる。
import fs from 'fs';
import http from 'http';
import path from 'path';
import { createRequire } from 'module';
import { fileURLToPath } from 'url';
import { createMockSupabase } from './mock-supabase.mjs';

const require = createRequire(import.meta.url);
const { chromium } = require('playwright');
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.json': 'application/json',
  '.png': 'image/png', '.webmanifest': 'application/manifest+json', '.woff2': 'font/woff2', '.svg': 'image/svg+xml' };
let fakeDataVersion = null;
const mock = createMockSupabase({ key: 'sb_publishable_e2e' });   // 同期のテスト用（このサーバーの /mock/ で応える）   // Service Worker のテスト用：index.json の dataVersion を差し替える
const server = http.createServer((req, res) => {
  let p = decodeURIComponent(new URL(req.url, 'http://x').pathname);
  if (p.endsWith('/')) p += 'index.html';
  if (p.startsWith('/mock/')) {
    let body = ''; req.on('data', c => { body += c; });
    req.on('end', async () => {
      const r = await mock.handle({ method: req.method, url: req.url, headers: req.headers, body: body || null });
      res.writeHead(r.status, { 'Content-Type': 'application/json' }); res.end(r.body);
    });
    return;
  }
  if (p === '/js/config.js') {
    res.writeHead(200, { 'Content-Type': 'text/javascript' });
    const v = fs.readFileSync(path.join(root, 'js/config.js'), 'utf8').match(/FILE_VERSION = '([^']+)'/)[1];
    return res.end(`export const FILE_VERSION = '${v}';
export const SUPABASE_URL = 'http://localhost:${server.address().port}/mock';
export const SUPABASE_KEY = 'sb_publishable_e2e';`);
  }
  if (fakeDataVersion && p === '/data/index.json') {
    const d = JSON.parse(fs.readFileSync(path.join(root, p), 'utf8')); d.dataVersion = fakeDataVersion;
    res.writeHead(200, { 'Content-Type': 'application/json' }); return res.end(JSON.stringify(d));
  }
  const f = path.join(root, p);
  if (!f.startsWith(root) || !fs.existsSync(f)) { res.writeHead(404); return res.end(); }
  res.writeHead(200, { 'Content-Type': TYPES[path.extname(f)] || 'application/octet-stream' });
  fs.createReadStream(f).pipe(res);
});
await new Promise(r => server.listen(0, r));
const U = `http://localhost:${server.address().port}/`;

let fail = 0, pass = 0;
const ok = (cond, msg) => { if (cond) pass++; else { fail++; console.log('FAIL', msg); } };
const browser = await chromium.launch();

async function newPage({ state = null, route = null, sw = false } = {}) {
  const ctx = await browser.newContext({ viewport: { width: 390, height: 844 }, acceptDownloads: true, serviceWorkers: sw ? 'allow' : 'block' });
  if (state) await ctx.addInitScript((s) => { if (!sessionStorage.getItem('seeded')) { sessionStorage.setItem('seeded', '1'); localStorage.setItem('alchemia.v1', s); } }, JSON.stringify(state));
  const page = await ctx.newPage();
  page.errors = [];
  page.on('pageerror', e => page.errors.push(String(e)));
  page.on('console', m => { if (m.type() === 'error') page.errors.push(m.text()); });
  page.on('dialog', d => d.accept());
  page.requested = [];
  page.on('request', r => page.requested.push(new URL(r.url()).pathname));
  if (route) await route(page);
  return page;
}
const next = (page) => page.locator('button:has-text("次")').first().click();

try {
  const fixture = JSON.parse(fs.readFileSync(path.join(root, 'tests/fixtures/state-v1.json'), 'utf8'));
  const LS = (page) => page.evaluate(() => JSON.parse(localStorage.getItem('alchemia.v1')));
  const conf = (page) => page.locator('#conf button').nth(1).click();
  // 1. 起動・タブ・5つの出題形式を入力して採点
  {
    const page = await newPage();
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    const tabs = await page.locator('.tabbar a span').allInnerTexts();
    ok(tabs.join() === '今日,図鑑,暗記,定着度,設定', `タブが化学版の構成 ${tabs}`);
    ok(await page.locator('#lv button').count() === 6, '難易度ボタンが6つ');
    ok(await page.locator('.fmt-card').count() === 6 && await page.locator('.fmt-card:disabled').count() === 0, '「形式ごとに演習」に6つの形式があり、すべて問題がある');
    ok(await page.locator('.banner').count() === 0, '起動時に警告が出ない');
    await page.evaluate("__alchemia.open(['INO-PPT-001', 'INO-PPT-007', 'THE-PER-001', 'INO-GPR-015', 'THE-MOL-001', 'ORG-ARO-001'])"); await page.waitForTimeout(1200);
    // 一問一答（語句・化学式）
    await conf(page);
    const inp = page.locator('#inputs input');
    ok(await inp.count() === 3, '一問一答：小問の欄が3つ');
    await inp.nth(0).fill('Cu(OH)₂'); await inp.nth(1).fill('[Cu(NH3)4]2+'); await inp.nth(2).fill('てとらあんみん銅（Ⅱ）イオン');
    ok(await page.locator('.preview .katex').count() >= 2, '化学式の入力は下付き・上付きで見える');
    await page.click('#judge'); await page.waitForTimeout(800);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '一問一答：表記ゆれ（全角・下付き文字・ひらがな）を吸収して完答');
    ok(await page.locator('.zk-open').count() === 1, '解いたあとに「図鑑で確かめる」');
    await next(page); await page.waitForTimeout(900);
    // 条件判断（すべて選ぶ）：選びすぎると不正解、正しく選ぶと正解
    await conf(page);
    ok(await page.locator('.field.pick .opt').count() === 6, '条件判断：選択肢が6つ');
    for (const i of [0, 1, 2]) await page.locator('.field.pick .opt').nth(i).click();
    await page.click('#judge'); await page.waitForTimeout(500);
    ok((await page.locator('#feedback').innerText()).includes('違う'), '条件判断：選びすぎは不正解');
    await page.locator('.field.pick .opt').nth(2).click(); await page.locator('.field.pick .opt').nth(4).click();
    await page.click('#judge'); await page.waitForTimeout(800);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '条件判断：選び直して正解（2回目）');
    await next(page); await page.waitForTimeout(900);
    // 比較（並べかえ）
    await conf(page);
    for (const t of ['Si', 'Mg', 'Al', 'Na']) await page.locator(`.field.order .pool button:has-text("${t}")`).click();
    ok(await page.locator('.field.order .placed').count() === 4, '比較：タップした順に並ぶ');
    await page.click('#judge'); await page.waitForTimeout(800);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '比較：Si > Mg > Al > Na で完答');
    await next(page); await page.waitForTimeout(900);
    // 反応式補完（係数）
    await conf(page);
    ok(await page.locator('.field.short').count() === 5 && await page.locator('#kp').count() === 0, '反応式補完：係数の欄は小さく、数式キーボードは出ない');
    for (const [i, v] of ['1', '4', '1', '2', '1'].entries()) await page.locator('#inputs input').nth(i).fill(v);
    await page.click('#judge'); await page.waitForTimeout(800);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '反応式補完：係数を答えて完答');
    await next(page); await page.waitForTimeout(900);
    // 計算（有効数字つきの数値）
    await conf(page);
    for (const [i, v] of ['0.10', '6.0×10^22', '2.24'].entries()) await page.locator('#inputs input').nth(i).fill(v);
    await page.click('#judge'); await page.waitForTimeout(800);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '計算：6.0×10^22 のような書き方も採点できる');
    await next(page); await page.waitForTimeout(900);
    // 反応経路
    await conf(page);
    ok(await page.locator('.question .katex').count() >= 1, '反応経路：経路の図（数式）が出る');
    await page.locator('#inputs input').nth(0).fill('ニトロベンゼン'); await page.locator('#inputs input').nth(1).fill('C6H5NH2');
    await page.locator('.field.pick .opt').first().click();
    await page.click('#judge'); await page.waitForTimeout(800);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '反応経路：名称でも示性式でも正解');
    await page.locator('.zk-open').click(); await page.waitForTimeout(1500);
    ok(page.url().includes('#zukan=organic') && await page.locator('.org-map svg .onode').count() > 20, '「図鑑で確かめる」で反応経路の地図が開く');
    const st = await LS(page);
    ok(['INO-PPT-001', 'INO-PPT-007', 'THE-PER-001', 'INO-GPR-015', 'THE-MOL-001', 'ORG-ARO-001'].every(id => st.pstat[id]), '結果が問題ID で記録される');
    ok(st.log.find(l => l.p === 'INO-PPT-007').g === 2, '2回目で正解した問題は「一部・ヒントあり」と同じ扱い');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 2. 図鑑の7つのページ
  {
    const page = await newPage();
    await page.goto(U + '#zukan'); await page.waitForTimeout(1200);
    ok(await page.locator('.zk-card').count() === 7, '図鑑のトップに7つのページ');
    await page.goto(U + '#zukan=ions'); await page.waitForTimeout(1200);
    ok(await page.locator('.ion-table tbody tr').count() === 17 && await page.locator('.ion-table thead th').count() === 14, '反応表：17イオン × 13試薬');
    await page.locator('.ion-table td[data-c="Cu^2+"][data-r="NH3-ex"]').click(); await page.waitForTimeout(300);
    ok((await page.locator('#idetail').innerText()).includes('テトラアンミン銅') && await page.locator('#idetail .katex').count() >= 2, '反応表：マスをタップすると錯イオンと反応式');
    await page.check('#ihide'); await page.waitForTimeout(300);
    ok(await page.locator('.ion-table.hide').count() === 1, '反応表：かくして覚える');
    await page.locator('#imode button[data-t="tube"]').click(); await page.waitForTimeout(300);
    await page.selectOption('#tion', 'Al^3+');
    await page.locator('.reag button[data-r="NaOH"]').click(); await page.waitForTimeout(200);
    ok((await page.locator('#tres').innerText()).includes('白色の沈殿') && await page.locator('#tex').count() === 1, '試験管：Al³⁺ に NaOH 少量で白色沈殿、「さらに過剰」ボタン');
    await page.click('#tex'); await page.waitForTimeout(200);
    ok((await page.locator('#tres').innerText()).includes('沈殿が溶ける'), '試験管：過剰で溶ける');
    await page.locator('#imode button[data-t="sys"]').click(); await page.waitForTimeout(300);
    for (let i = 0; i < 6; i++) { await page.click('#snext'); await page.waitForTimeout(80); }
    ok(await page.locator('.sys-steps li').count() === 6 && (await page.locator('.sys-steps').innerText()).includes('AgCl'), '系統分析：6段階をたどる');
    await page.locator('#imode button[data-t="flame"]').click(); await page.waitForTimeout(300);
    ok(await page.locator('.flame').count() === 7, '炎色反応：7元素');
    await page.goto(U + '#zukan=reactions'); await page.waitForTimeout(1000);
    ok(await page.locator('.rx-card').count() >= 18, '製法：気体の製法の一覧');
    await page.check('#rhide'); await page.waitForTimeout(300);
    const hid = page.locator('.rhs.hid').first();
    await hid.click();
    ok(await page.locator('.rhs.hid.shown').count() === 1, '製法：右辺をかくして、タップで表示');
    await page.locator('#rmode button[data-t="industry"]').click(); await page.waitForTimeout(300);
    ok((await page.locator('#rbody').innerText()).includes('オストワルト法'), '製法：工業的製法');
    await page.locator('#rmode button[data-t="gtable"]').click(); await page.waitForTimeout(300);
    ok(await page.locator('.gas-table tbody tr').count() >= 15, '製法：気体の性質の表');
    await page.goto(U + '#zukan=organic'); await page.waitForTimeout(1000);
    await page.locator('#omap button[data-t="aromatic"]').click(); await page.waitForTimeout(300);
    await page.locator('.onode[data-c="aniline"]').click(); await page.waitForTimeout(800);
    ok((await page.locator('#odetail').innerText()).includes('アニリン') && await page.locator('#odetail .struct svg').count() === 1, '反応経路：化合物をタップすると構造式と反応');
    await page.check('#opath'); await page.waitForTimeout(200);
    await page.locator('.onode[data-c="benzene"]').click(); await page.waitForTimeout(200);
    await page.locator('.onode[data-c="aspirin"]').click(); await page.waitForTimeout(500);
    ok(await page.locator('.path-steps > li').count() >= 4 && await page.locator('.oedge.on').count() >= 4, '反応経路：経路さがし（ベンゼン → アセチルサリチル酸）');
    await page.goto(U + '#zukan=tests'); await page.waitForTimeout(1000);
    ok(await page.locator('.test-card').count() >= 13, '検出反応：反応ごとの一覧');
    await page.locator('#tmode button[data-t="matrix"]').click(); await page.waitForTimeout(300);
    ok(await page.locator('.test-table td.tv.pos').count() >= 30, '検出反応：化合物 × 反応の表');
    await page.goto(U + '#zukan=polymers/nylon66'); await page.waitForTimeout(1500);
    ok(await page.locator('.poly-card[data-p="nylon66"] .struct svg').count() >= 2 && (await page.locator('.poly-card[data-p="nylon66"]').innerText()).includes('226'), '高分子：単量体と繰り返し単位の構造式、式量');
    await page.goto(U + '#zukan=elements/Cu'); await page.waitForTimeout(1200);
    ok(await page.locator('#ptable .el[data-s]').count() === 118, '周期表：118元素');
    ok((await page.locator('#edetail').innerText()).includes('3d10 4s1'), '電子配置：Cu は 3d¹⁰ 4s¹');
    await page.selectOption('#eprop', 'en'); await page.waitForTimeout(200);
    ok((await page.locator('#elegend').innerText()).includes('〜'), '塗り分け：電気陰性度');
    await page.locator('#iechart .pt[data-s="Na"]').click(); await page.waitForTimeout(200);
    ok((await page.locator('#edetail').innerText()).includes('ナトリウム'), 'イオン化エネルギーのグラフの点をタップすると元素が選ばれる');
    await page.goto(U + '#zukan=ions'); await page.waitForTimeout(800);
    await page.locator('.zk-prob').click(); await page.waitForTimeout(1200);
    ok(await page.locator('#conf button').count() === 3, '図鑑から関連する問題を解ける');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 3. 自己採点を含む問題・法則・設定・定着度
  {
    const page = await newPage();
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    await page.evaluate("__alchemia.open(['INO-PPT-016'])"); await page.waitForTimeout(1000);
    await conf(page);
    ok(await page.locator('.field.self').count() === 1, '記述の小問は自己採点');
    const inp = page.locator('#inputs input');
    for (const [i, v] of ['AgCl', 'CuS', 'Fe(OH)3', '赤褐色', 'ZnS', 'BaCO3', '炎色反応'].entries()) await inp.nth(i).fill(v);
    await page.click('#judge'); await page.waitForTimeout(600);
    ok(await page.locator('#self button').count() === 3 && await page.locator('#solbox').count() === 1, '記述の小問があると、解説を見てから自己採点');
    await page.locator('#self button').first().click(); await page.waitForTimeout(400);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '自動採点と自己採点を合わせて完答');
    await page.goto(U + '#concepts'); await page.waitForTimeout(800);
    ok(await page.locator('a.law').count() === 41, '「法則」タブに理論化学の定義・法則（41件）');
    await page.goto(U + '#law=ideal-gas-law'); await page.waitForTimeout(1000);
    ok(await page.locator('ol.derivation li.open').count() === 1, '法則のページ：導出は1段目だけ開いている');
    await page.click('#dall'); await page.waitForTimeout(200);
    ok(await page.locator('ol.derivation li:not(.open)').count() === 0 && await page.locator('.law-page .katex').count() > 5, 'すべて表示で導出が全部出る');
    await page.goto(U + '#law=boyle'); await page.waitForTimeout(800);
    ok((await page.locator('.law-page').innerText()).includes('この法則の根拠'), '実験から得た法則は根拠を出す');
    await page.goto(U + '#settings'); await page.waitForTimeout(600);
    await page.locator('.segment.theme button:has-text("ダーク")').click();
    ok(await page.evaluate(() => document.documentElement.dataset.theme) === 'dark', 'ダークモード');
    await page.goto(U + '#zukan=elements'); await page.waitForTimeout(1000);
    ok(await page.locator('#ptable .el').count() > 100, 'ダークモードでも図鑑が出る');
    await page.goto(U + '#map'); await page.waitForTimeout(600);
    ok((await page.locator('#view').innerText()).includes('定着度'), '定着度画面');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 3'. 暗記カード
  {
    const page = await newPage();
    await page.goto(U + '#cards'); await page.waitForTimeout(1200);
    ok(await page.locator('.deck-card').count() >= 15, `暗記：デッキが並ぶ（${await page.locator('.deck-card').count()}）`);
    await page.goto(U + '#study=flame'); await page.waitForTimeout(1000);
    ok(await page.locator('#grade').isHidden(), '暗記：答えを見るまで評価ボタンは出ない');
    await page.click('#flash'); await page.waitForTimeout(200);
    ok(await page.locator('.fback').isVisible() && await page.locator('#grade button').count() === 3, '暗記：タップで裏が出て、3つの評価ボタン');
    await page.click('#grade .g-again'); await page.waitForTimeout(200);
    for (let i = 0; i < 8; i++) { if (!(await page.locator('#flash').count())) break; await page.click('#flash'); await page.click('#grade .g-good'); await page.waitForTimeout(120); }
    ok((await page.locator('#view').innerText()).includes('おつかれさま'), '暗記：「もう一度」のカードももう一度出て、最後まで学習できる');
    const st = await LS(page);
    const ks = Object.keys(st.cards || {}).filter(k => k.startsWith('flame/'));
    ok(ks.length === 7 && ks.every(k => st.cards[k].box === 1 && st.cards[k].due > Date.now()), '暗記：記録がカードごとに残り、次は1日後');
    await page.goto(U + '#deck=flame'); await page.waitForTimeout(800);
    ok(await page.locator('.ctag.mid').count() === 7, '暗記：デッキの一覧で「学習中」になる');
    await page.goto(U + '#deck=laws'); await page.waitForTimeout(800);
    await page.locator('.card-row summary').first().click(); await page.locator('.card-row .linkbtn').first().click(); await page.waitForTimeout(1000);
    ok(page.url().includes('#law=') && await page.locator('.law-page').count() === 1, '暗記：「定義・法則」のカードから導出・根拠のページへ');
    await page.goto(U + '#cards'); await page.waitForTimeout(600);
    await page.check('#crev');
    await page.goto(U + '#study=constants'); await page.waitForTimeout(800);
    ok((await page.locator('.fface').innerText()).includes('10'), '暗記：裏から出す（値 → 名前）');
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    ok(await page.locator('.cards-teaser').count() === 1, '今日の画面に暗記カードの入口');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 4. 学習記録のある端末：CSV・続きに記録
  {
    const page = await newPage({ state: fixture });
    await page.goto(U + '#settings'); await page.waitForTimeout(1200);
    const dl = page.waitForEvent('download');
    await page.click('#csvall');
    const csv = fs.readFileSync(await (await dl).path(), 'utf8');
    ok(csv.split('\n').filter(l => l.trim()).length === fixture.log.length + 1, `CSV に記録 ${fixture.log.length}件がすべて出る`);
    await page.evaluate("__alchemia.open(['THE-MOL-001'])"); await page.waitForTimeout(1000);
    await conf(page);
    for (const [i, v] of ['0.1', '6e22', '2.24'].entries()) await page.locator('#inputs input').nth(i).fill(v);
    await page.click('#judge'); await page.waitForTimeout(600);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '記録のある端末で解ける');
    const st = await LS(page);
    ok(Object.keys(fixture.pstat).every(k => st.pstat[k]), '前の復習状態がすべて残る');
    ok(st.log.length === fixture.log.length + 1, '学習記録が続きに足される');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 5. 一部の data/ だけ古い
  {
    const page = await newPage({ route: (p) => p.route('**/data/f-INO-08.json', async (r) => {
      const res = await r.fetch(); const d = await res.json(); d._meta = { dataVersion: 'old000000000' };
      await r.fulfill({ response: res, body: JSON.stringify(d) });
    }) });
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    await page.evaluate("__alchemia.open(['INO-PPT-001'])"); await page.waitForTimeout(800);
    const b = page.locator('.banner');
    ok(await b.count() === 1 && (await b.innerText()).includes('data/f-INO-08.json'), '古い本文ファイルを読むと警告が出る');
    ok(await page.locator('#inputs input').count() > 0, '警告が出ても問題は解ける');
  }
  // 6. 検索・問題セット・共有リンク
  {
    const page = await newPage({ state: fixture });
    await page.goto(U + '#today'); await page.waitForTimeout(1500);
    ok(await page.locator('#sets .set-card').count() === 3, '今日の画面に問題セット3つ');
    await page.fill('#sq', '錯イオン'); await page.press('#sq', 'Enter'); await page.waitForTimeout(1200);
    ok(await page.locator('.res-row').count() >= 2, '検索できる');
    await page.goto(U + '#set=col-ions-start'); await page.waitForTimeout(1500);
    ok(page.url().endsWith('#today') && await page.locator('#conf button').count() === 3, '#set= のリンクで問題セットが始まる');
    await page.goto(U + '#ids=ORG-ARO-001'); await page.waitForTimeout(1500);
    ok((await page.locator('.chip').allInnerTexts()).some(t => t.includes('有機化学・芳香族化合物')), '#ids= のリンクで指定の問題が出る');
    await page.goto(U + '#today'); await page.waitForTimeout(800);
    ok(await page.locator('#conf button').count() === 3 || await page.locator('.hero').count() === 1, '今日の画面へ戻れる');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 7'. Service Worker：自分のキャッシュだけを管理し、同じドメインの他のアプリのキャッシュは消さない。オフラインでも問題・図鑑が使える
  {
    const page = await newPage({ sw: true });
    await page.goto(U + '#today');
    await page.evaluate(async () => { await (await caches.open('lemma-v0.9.0')).put('/x', new Response('lemma')); await (await caches.open('celestia-v0.5.0')).put('/y', new Response('c')); });
    await page.waitForFunction(() => navigator.serviceWorker.controller, null, { timeout: 15000 }).catch(() => {});
    await page.reload(); await page.waitForTimeout(1200);
    await page.evaluate("__alchemia.open(['INO-PPT-001'])"); await page.waitForTimeout(1500);
    await page.goto(U + '#zukan=organic'); await page.waitForTimeout(1500);
    await page.locator('.onode[data-c="ethanol"]').click(); await page.waitForTimeout(800);
    const v = JSON.parse(fs.readFileSync(path.join(root, 'data/index.json'), 'utf8')).dataVersion;
    const keys = await page.evaluate(async () => caches.keys());
    ok(keys.includes('alchemia-data-' + v) && keys.includes('alchemia-v0.3.1'), `コード用とデータ用のキャッシュ ${keys}`);
    ok(keys.includes('lemma-v0.9.0') && keys.includes('celestia-v0.5.0'), '同じドメインの他のアプリのキャッシュを消さない');
    await page.context().setOffline(true);
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    await page.evaluate("__alchemia.open(['INO-PPT-001'])"); await page.waitForTimeout(1500);
    ok(await page.locator('#inputs input').count() === 3, 'オフラインでも開いた問題が解ける');
    await page.goto(U + '#zukan=organic'); await page.waitForTimeout(1500);
    await page.locator('.onode[data-c="ethanol"]').click(); await page.waitForTimeout(800);
    ok(await page.locator('#odetail .struct svg').count() === 1, 'オフラインでも図鑑（構造式）が使える');
    await page.context().setOffline(false);
    fakeDataVersion = 'e2enewversion';
    await page.goto(U + '#today'); await page.waitForTimeout(1500);
    await page.evaluate(async () => fetch('data/index.json', { cache: 'no-cache' }));
    await page.waitForTimeout(500);
    const keys2 = await page.evaluate(async () => caches.keys());
    ok(keys2.includes('alchemia-data-e2enewversion') && !keys2.includes('alchemia-data-' + v) && keys2.includes('lemma-v0.9.0'), `データの版が変わると自分の古いデータ用キャッシュだけを消す ${keys2}`);
    fakeDataVersion = null;
  }
  // 7. アカウントと同期：2台の端末でログインすると記録が合わさる（にせ Supabase）
  {
    const login = async (page) => {
      await page.goto(U + '#settings'); await page.waitForTimeout(1500);
      await page.fill('#acct-email', 'me@example.com'); await page.click('#acct-send'); await page.waitForTimeout(500);
      await page.fill('#acct-code', '123456'); await page.click('#acct-verify'); await page.waitForTimeout(1500);
      // 端末とアカウントの両方に記録があるときは、合わせるかを聞かれる（Lemma v0.27 と同じ）
      page.asked = await page.locator('.chooser').count() > 0;
      if (page.asked) { await page.click('.chooser [data-k="merge"]'); await page.waitForTimeout(1500); }
    };
    const A = await newPage({ state: { ...fixture, log: fixture.log.slice(0, 25) } });
    await login(A);
    ok((await A.locator('#acct').innerText()).includes('Lemma と共通'), 'アカウントが Lemma と共通であることを表示');
    const B = await newPage({ state: { ...fixture, log: fixture.log.slice(20) } });
    await login(B);
    ok(!A.asked && B.asked, 'ログイン時：アカウントにも端末にも記録があるときだけ、合わせるかを聞く');
    const nB = (await LS(B)).log.length;
    ok(nB === fixture.log.length, `端末B：両方の端末の記録が合わさる（${nB}件）`);
    ok(await B.evaluate(() => !!localStorage.getItem('alchemia.auth') && !localStorage.getItem('lemma.auth')), 'ログイン情報は alchemia.auth に保存（Lemma と別）');
    ok([...mock.rows.values()][0] && mock.rows.size === 1, 'サーバーの表 alchemia_state に1行');
    const real = (es) => es.filter(e => !/status of 403/.test(e));
    ok(real(A.errors).length === 0 && real(B.errors).length === 0, `ページのエラーなし ${A.errors} ${B.errors}`);
  }
  // 8. Lemma v0.27 に合わせた学習の流れ：今日の学習・前回の続き・ブックマーク・解答の形式・類題・バックアップ・消去
  {
    mock.users.clear(); mock.rows.clear();
    const P = await newPage();
    await P.goto(U + '#today'); await P.waitForTimeout(1500);
    ok(await P.locator('.first-run #firstgo').count() === 1 && (await P.locator('h1').innerText()) === '今日の学習', '初めて：「今日の学習」と「まず8問」');
    await P.click('#firstgo'); await P.waitForTimeout(1500);
    ok((await P.locator('.why').innerText()).includes('はじめての演習'), '出題の理由を表示する');
    // 1問解いて中断 → 前回の続き
    await P.evaluate("__alchemia.open(['INO-GPR-007', 'INO-HAL-001', 'THE-MOL-001'])"); await P.waitForTimeout(1200);
    await P.click('[data-bm="INO-GPR-007"]'); await P.waitForTimeout(200);
    ok(await P.locator('[data-bm="INO-GPR-007"].on').count() === 1, 'ブックマーク：問題画面で付けられる');
    await P.locator('#conf button').first().click();
    await P.fill('#inputs input', 'NH3'); await P.click('#judge'); await P.waitForTimeout(800);
    ok((await P.locator('.conf-note').innerText()).includes('自信がないまま正解'), '自信なしで正解：早めに復習に出すと知らせる');
    const g = (await LS(P)).log.at(-1);
    ok(g.g === 2 && g.md === 'input', `自信なしの正解は「ヒントあり」と同じ評価（g=${g.g}）・答え方を記録`);
    ok(await P.locator('#similar').count() === 1, '解いたあとに「類題を1問追加」');
    await P.click('#similar'); await P.waitForTimeout(800);
    ok((await P.locator('#similar').innerText()).includes('追加しました'), '類題を次に差し込む');
    await P.click('#next'); await P.waitForTimeout(800);
    ok((await P.locator('.why').innerText()).includes('類題'), '類題には理由「類題」');
    await P.click('#quit'); await P.waitForTimeout(800);
    ok((await P.locator('#resume').innerText()).includes('4問の 2問目'), `中断すると「前回の続き」（${await P.locator('#resume').innerText()}）`);
    ok((await P.locator('#bmset').innerText()).includes('1問') && await P.locator('.goal').count() === 1, 'ホームに今日の目標とブックマーク');
    await P.click('#resume'); await P.waitForTimeout(1200);
    ok((await P.locator('.count').innerText()).startsWith('2 /'), '前回の続きから再開できる');
    await P.click('#quit'); await P.waitForTimeout(500);
    await P.click('#todaygo'); await P.waitForTimeout(1500);
    ok(await P.locator('.why').count() === 1, '今日の演習：出題の理由つき');
    await P.click('#quit'); await P.waitForTimeout(500);
    // 解答の形式：選択式（語句の欄が4択に）・自己採点
    await P.goto(U + '#settings'); await P.waitForTimeout(800);
    await P.click('#amode [data-m="choice"]'); await P.waitForTimeout(300);
    await P.evaluate("__alchemia.open(['INO-GPR-007'])"); await P.waitForTimeout(1200);
    ok(await P.locator('#inputs .ch-opt').count() === 4 && await P.locator('#inputs input').count() === 0, '選択式：語句・化学式の欄が4つの選択肢になる');
    await P.locator('#conf button').nth(2).click();
    await P.locator('#inputs .ch-opt').filter({ hasNotText: 'NH' }).first().click(); await P.click('#judge'); await P.waitForTimeout(800);
    ok(await P.locator('#inputs .ch-opt.right').count() === 1 && await P.locator('.verdict').count() === 1, '選択式：まちがえるとやり直しなしで答えを見せる');
    ok((await P.locator('.conf-note.warn').innerText()).includes('確実'), '「確実」で不正解：思い込みの確認を促す');
    ok((await LS(P)).log.at(-1).md === 'choice', '答え方（選択式）を記録');
    await P.goto(U + '#settings'); await P.waitForTimeout(800);
    await P.click('#amode [data-m="self"]'); await P.waitForTimeout(300);
    await P.evaluate("__alchemia.open(['INO-GPR-007'])"); await P.waitForTimeout(1200);
    ok(await P.locator('#inputs').count() === 0 && await P.locator('#reveal').count() === 1, '自己採点：入力欄なしで答え合わせ');
    await P.goto(U + '#settings'); await P.waitForTimeout(800);
    await P.click('#amode [data-m="input"]'); await P.waitForTimeout(300);
    // バックアップ：形式つきで書き出し → 別のファイルは断る → 置き換え → 控えから戻す
    const [exp] = await Promise.all([P.waitForEvent('download'), P.click('#exp')]);
    const bk = JSON.parse(fs.readFileSync(await exp.path(), 'utf8'));
    ok(bk.format === 'alchemia-backup' && bk.version === 1 && bk.schema === 2 && /^fnv1a32:/.test(bk.hash) && bk.state.log.length === 2, '書き出し：形式・版・確認用の値つきのバックアップ');
    await P.setInputFiles('#imp', { name: 'x.json', mimeType: 'application/json', buffer: Buffer.from('{"hello":1}') });
    await P.waitForTimeout(400);
    ok((await P.locator('.toast').innerText()).includes('ALCHEMIA の学習記録のファイルではありません'), '読み込み：別のファイルは断る');
    await P.setInputFiles('#imp', { name: 'b.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(fixture)) });
    await P.waitForTimeout(400);
    await P.click('.chooser [data-k="replace"]'); await P.waitForTimeout(800);
    ok((await P.locator('.chooser').innerText()).includes('読み込みました'), '読み込み：結果（件数・点検）を表示する');
    await P.click('.chooser [data-k="ok"]'); await P.waitForTimeout(400);
    const st = await LS(P);
    ok(st.log.length === fixture.log.length && st.resetAt > 0 && st.schema === 2, `読み込み：v0.2 の記録を移行して置き換え、消した印が付く（${st.log.length}件）`);
    await P.goto(U + '#settings'); await P.waitForTimeout(800);
    ok(await P.locator('.snap-back').count() >= 1 && (await P.locator('.integrity summary').innerText()).includes('記録の点検'), '設定：記録の点検と自動の控えが見える');
    await P.click('.snap-back >> nth=0'); await P.waitForTimeout(300);
    await P.click('.chooser [data-k="replace"]'); await P.waitForTimeout(800);
    await P.click('.chooser [data-k="ok"]'); await P.waitForTimeout(400);
    ok((await LS(P)).log.length === bk.state.log.length, '控えの時点に戻せる');
    // 消去：控えを残してから消す
    await P.goto(U + '#settings'); await P.waitForTimeout(800);
    await P.click('#reset'); await P.waitForTimeout(300);
    await P.click('.chooser [data-k="all"]'); await P.waitForTimeout(800);
    const st2 = await LS(P);
    ok(st2.log.length === 0 && st2.resetAt > 0 && (await P.evaluate(() => JSON.parse(localStorage.getItem('alchemia.snap.v1')).length)) >= 1, '記録をすべて消す：消す前の記録は控えに残る');
    ok(P.errors.length === 0, `ページのエラーなし ${P.errors}`);
    // 壊れた記録は上書きせずに退避する
    const Q = await browser.newContext({ viewport: { width: 390, height: 844 }, serviceWorkers: 'block' });
    await Q.addInitScript(() => { if (!sessionStorage.getItem('seeded')) { sessionStorage.setItem('seeded', '1'); localStorage.setItem('alchemia.v1', '{壊れた'); } });
    const q = await Q.newPage();
    await q.goto(U + '#today'); await q.waitForTimeout(1500);
    ok((await q.locator('.banner').innerText()).includes('読めませんでした') && await q.evaluate(() => Object.keys(localStorage).some(k => k.startsWith('alchemia.v1.broken.'))), '読めない記録は退避して知らせる');
  }
} finally {
  await browser.close();
  server.close();
}
console.log(`e2e: ${pass} 件 OK${fail ? `、${fail} 件失敗` : ''}`);
process.exitCode = fail ? 1 : 0;
