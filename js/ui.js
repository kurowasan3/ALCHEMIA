// アイコン・効果音・テーマ・マスコットのサラマンダー「ラヴォ」
export const FILE_VERSION = '0.3.0';

// ---------------------------------------------------------------- アイコン（線画、currentColor）
const S = (body, extra = '') => `<svg class="ico" viewBox="0 0 24 24" aria-hidden="true" ${extra}>${body}</svg>`;
export const ICON = {
  // 自信度：自信なし（?）・たぶん（半分の円）・確実（チェック付きの円）
  low: S('<circle cx="12" cy="12" r="9" fill="none" stroke="currentColor" stroke-width="1.8" stroke-dasharray="3 3"/><path d="M9.6 9.6a2.5 2.5 0 1 1 3.4 2.3c-.7.3-1 .8-1 1.5v.4" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"/><circle cx="12" cy="16.9" r="1.1" fill="currentColor"/>'),
  mid: S('<circle cx="12" cy="12" r="9" fill="none" stroke="currentColor" stroke-width="1.8"/><path d="M12 3a9 9 0 0 1 0 18z" fill="currentColor" opacity=".55"/>'),
  high: S('<circle cx="12" cy="12" r="9" fill="currentColor" opacity=".18"/><circle cx="12" cy="12" r="9" fill="none" stroke="currentColor" stroke-width="1.8"/><path d="M7.8 12.4l2.8 2.8 5.6-6" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>'),
  // 採点：完答（二重丸）・一部正解（三角）・不正解（×）
  full: S('<circle cx="12" cy="12" r="9" fill="none" stroke="currentColor" stroke-width="2"/><circle cx="12" cy="12" r="4.6" fill="none" stroke="currentColor" stroke-width="2"/>'),
  part: S('<path d="M12 3.8l8.6 15.4H3.4z" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>'),
  none: S('<path d="M6 6l12 12M18 6L6 18" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round"/>'),
  ok: S('<circle cx="12" cy="12" r="8.5" fill="none" stroke="currentColor" stroke-width="2.2"/>'),
  ng: S('<path d="M6.5 6.5l11 11M17.5 6.5l-11 11" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round"/>'),
  self: S('<path d="M4 20l4-1 10.5-10.5a2.1 2.1 0 0 0-3-3L5 16z" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linejoin="round"/><path d="M13.5 7.5l3 3" stroke="currentColor" stroke-width="1.8"/>'),
  hint: S('<path d="M9 18h6M10 21h4M12 3a6 6 0 0 0-3.6 10.8c.6.5 1 1.2 1 2V16h5.2v-.2c0-.8.4-1.5 1-2A6 6 0 0 0 12 3z" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linejoin="round"/>'),
  cases: S('<path d="M5 4v16M5 7h4M5 17h4" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"/><path d="M12 7h7M12 17h7" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"/>'),
  sound: S('<path d="M4 9h4l5-4v14l-5-4H4z" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linejoin="round"/><path d="M16 9a4 4 0 0 1 0 6M18.5 6.5a7.5 7.5 0 0 1 0 11" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"/>'),
};

// ---------------------------------------------------------------- 効果音（音声ファイルを使わず合成）
let ctx = null;
let soundOn = true;
export function setSound(on) { soundOn = on; }
function ac() {
  if (!ctx) {
    const C = window.AudioContext || window.webkitAudioContext;
    if (!C) return null;
    ctx = new C();
  }
  if (ctx.state === 'suspended') ctx.resume();
  return ctx;
}
function tone(freq, start, dur, { type = 'sine', gain = 0.16, glide = null } = {}) {
  const a = ac(); if (!a) return;
  const t = a.currentTime + start;
  const o = a.createOscillator(), g = a.createGain();
  o.type = type;
  o.frequency.setValueAtTime(freq, t);
  if (glide) o.frequency.exponentialRampToValueAtTime(glide, t + dur);
  g.gain.setValueAtTime(0.0001, t);
  g.gain.exponentialRampToValueAtTime(gain, t + 0.012);
  g.gain.exponentialRampToValueAtTime(0.0001, t + dur);
  o.connect(g).connect(a.destination);
  o.start(t); o.stop(t + dur + 0.02);
}
export function play(name) {
  if (!soundOn) return;
  try {
    switch (name) {
      case 'correct': tone(1046.5, 0, 0.16); tone(1318.5, 0.09, 0.32); break;            // ド→ミ
      case 'full': tone(784, 0, 0.14); tone(1046.5, 0.08, 0.16); tone(1318.5, 0.16, 0.42); break;
      case 'part': tone(659.3, 0, 0.18, { gain: 0.12 }); tone(784, 0.1, 0.3, { gain: 0.12 }); break;
      case 'wrong': tone(233, 0, 0.26, { type: 'triangle', gain: 0.2, glide: 175 }); break;
      case 'tap': tone(1600, 0, 0.04, { gain: 0.05 }); break;
      case 'finish': [523.3, 659.3, 784, 1046.5].forEach((f, i) => tone(f, i * 0.09, 0.4, { gain: 0.12 })); break;
      case 'mascot': tone(1568, 0, 0.06, { gain: 0.05, glide: 2093 }); tone(2093, 0.08, 0.08, { gain: 0.04, glide: 2637 }); break;   // ぴょこっ
    }
  } catch { /* 音が出せない環境 */ }
}

// ---------------------------------------------------------------- テーマ（自動・ライト・ダーク）
export function applyTheme(mode) {
  const root = document.documentElement;
  if (mode === 'light' || mode === 'dark') root.dataset.theme = mode; else delete root.dataset.theme;
  const dark = mode === 'dark' || (mode !== 'light' && window.matchMedia('(prefers-color-scheme: dark)').matches);
  document.querySelectorAll('meta[name="theme-color"]').forEach(m => m.setAttribute('content', dark ? '#16171b' : '#f6f4ee'));
}

// ---------------------------------------------------------------- サラマンダー「ラヴォ」
// 錬金術で火の精とされたサラマンダー（オリジナル）。名前は質量保存の法則のラボアジエから。尾の先に小さな炎。色はテーマに合わせて CSS 変数から
export const MASCOT_SVG = `<svg viewBox="0 0 64 64" class="mascot-svg" aria-hidden="true">
  <path d="M45 41c9-1 14-6 12-12-1.4-4.2-6.6-4.6-7.6-1.4-.7 2.2 1.4 3.6 3 2.8" fill="none" stroke="var(--mascot-line)" stroke-width="7.6" stroke-linecap="round"/>
  <path d="M45 41c9-1 14-6 12-12-1.4-4.2-6.6-4.6-7.6-1.4-.7 2.2 1.4 3.6 3 2.8" fill="none" stroke="var(--mascot-body)" stroke-width="5.4" stroke-linecap="round"/>
  <path d="M52.6 27.6c-1.8-2.4-1.2-4.8.6-6.6.2 1.6 1 2.4 2 3 .2-1.6 1-2.8 2.2-3.6-.4 2.4.6 3.6.6 5.2 0 1.8-1.6 3.2-3.4 3.2-.8 0-1.6-.4-2-1.2z" fill="var(--mascot-star)"/>
  <ellipse cx="25" cy="54" rx="4.6" ry="2.8" fill="var(--mascot-body)" stroke="var(--mascot-line)" stroke-width="1.1"/>
  <ellipse cx="42" cy="54" rx="4.6" ry="2.8" fill="var(--mascot-body)" stroke="var(--mascot-line)" stroke-width="1.1"/>
  <ellipse cx="34" cy="45" rx="16" ry="10" fill="var(--mascot-body)" stroke="var(--mascot-line)" stroke-width="1.1"/>
  <ellipse cx="33" cy="49" rx="10" ry="4.6" fill="var(--mascot-belly)"/>
  <circle cx="38" cy="39.5" r="1.7" fill="var(--mascot-star)"/><circle cx="44" cy="42" r="1.3" fill="var(--mascot-star)"/><circle cx="31.5" cy="38.5" r="1.2" fill="var(--mascot-star)"/>
  <ellipse cx="19.5" cy="34" rx="13" ry="10.5" fill="var(--mascot-body)" stroke="var(--mascot-line)" stroke-width="1.1"/>
  <circle cx="14" cy="31" r="2.6" fill="var(--mascot-eye)"/><circle cx="14.8" cy="30.1" r="0.95" fill="#fff"/>
  <circle cx="23.5" cy="30.5" r="2.6" fill="var(--mascot-eye)"/><circle cx="24.3" cy="29.6" r="0.95" fill="#fff"/>
  <path d="M15.5 37.5c1.6 1.6 4.6 1.6 6.2 0" fill="none" stroke="var(--mascot-line)" stroke-width="1" stroke-linecap="round"/>
  <ellipse cx="10.5" cy="35.6" rx="2.3" ry="1.4" fill="var(--mascot-cheek)" opacity=".7"/><ellipse cx="27.5" cy="35" rx="2.3" ry="1.4" fill="var(--mascot-cheek)" opacity=".7"/>
</svg>`;

export const MASCOT_NAME = 'ラヴォ';

const HEALTH = [
  'お身体は大丈夫ですか？ 肩に力が入っていたら、一度深呼吸しましょう。',
  '水分はとっていますか？ のどが渇く前に一口どうぞ。',
  '20分ごとに、20秒ほど遠くを見ると目が楽になります。',
  '背中が丸まっていませんか？ 画面は目の高さに近づけて。',
  '眠気が強いときは、5〜15分の短い仮眠も手です。',
];
const ADVICE = [
  '反応式は「何が何に変わるか」を言葉にしてから書こう。酸化剤・還元剤なら、まず半反応式から。',
  '沈殿の色は、表を眺めるより「試薬を加えたら何が起こる？」と自分に問う方が残ります（想起練習）。',
  '有機の反応経路は、官能基が何に変わったかだけを追うと見通しがよくなります。',
  '計算の前に、単位（mol・g・L）を式の中で一緒に動かしてみよう。答えの単位が合えば式の形もたいてい正しい。',
  '解説を読み返すより、閉じてもう一度解いてみる方が記憶に残ります。',
  '同じ単元を一気にやるより、日をあけて解き直す方が長持ちします（分散学習）。',
  '無機・有機・理論を混ぜて解くと、「どの知識を使うか」を選ぶ力がつきます。',
  '間違えた問題は宝です。どこで考えがずれたか、1行で書いておこう。',
  '例外（Cr・Cu の電子配置、Al(OH)₃ は NH₃ 水に溶けない…）は、理由と一緒に覚えると忘れにくい。',
  '眠っている間に記憶は整理されます。夜ふかしより、早めに寝る方が得です。',
];
const QUOTES = [
  ['自然においては、何ものも失われず、何ものも創られない。すべては変化するだけである。', 'ラボアジエ（趣意）'],
  ['人生に恐れるべきものは何もない。ただ理解すべきものがあるだけだ。', 'マリー・キュリー'],
  ['偶然は、準備のできた心にだけ味方する。', 'ルイ・パスツール'],
  ['自分で作れないものは、理解していない。', 'リチャード・ファインマン'],
  ['よいアイデアを得る最良の方法は、たくさんのアイデアを得ることだ。', 'ライナス・ポーリング'],
  ['千里の行も足下に始まる。', '老子'],
  ['学びて思わざれば則ち罔し、思いて学ばざれば則ち殆し。', '孔子『論語』'],
];

export function mascotMessage(kind, ctxInfo = {}) {
  const pick = (a) => a[Math.floor(Math.random() * a.length)];
  const h = new Date().getHours();
  if (kind === 'late') return 'もう遅い時間です。今日はここまでにして、ぐっすり眠ろう。睡眠中に記憶が整理されます。';
  if (kind === 'break') return ctxInfo.min >= 50
    ? `${ctxInfo.min}分続けて勉強しています。お身体は大丈夫ですか？ 5〜10分、立って体を伸ばしましょう。`
    : `${ctxInfo.min}分たちました。少し目を閉じて休憩しませんか？`;
  if (kind === 'summary') {
    const { solo = 0, missed = 0, n = 0 } = ctxInfo;
    if (n && solo === n) return '全部自力で解けました！ 次は難易度を1つ上げてみては？';
    if (n && missed >= n / 2) return '難しかったですね。でも、難しく感じた日の方があとで伸びます。間違えた問題は近いうちにまた出します。';
    return 'おつかれさま。間違えたところは、復習で必ず取り返せます。';
  }
  if (kind === 'greet') {
    const g = h < 5 ? 'こんばんは。夜遅くまでおつかれさま。' : h < 11 ? 'おはようございます。' : h < 18 ? 'こんにちは。' : 'こんばんは。';
    const r = Math.random();
    if (r < 0.3) return g + pick(HEALTH);
    if (r < 0.7) return g + pick(ADVICE);
    const q = pick(QUOTES); return `${g}今日のことば：「${q[0]}」（${q[1]}）`;
  }
  const r = Math.random();
  if (r < 0.3) return pick(HEALTH);
  if (r < 0.65) return pick(ADVICE);
  const q = pick(QUOTES); return `「${q[0]}」\n― ${q[1]}`;
}
