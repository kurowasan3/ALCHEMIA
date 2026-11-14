# 検証スクリプト THE-ION-087  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/THE-ION-087/item.json'
import sys, json, re, subprocess, os
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
FIELDS = ['原子の構造と電子配置','元素の周期表と周期性','化学結合と結晶','物質量と化学反応式','酸と塩基','酸化還元反応','電池と電気分解','気体の性質','溶液の性質','化学反応とエネルギー','反応速度と化学平衡','電離平衡と溶解平衡','ハロゲン・貴ガス・水素','酸素・硫黄','窒素・リン','炭素・ケイ素','アルカリ金属・アルカリ土類金属','アルミニウム・亜鉛・スズ・鉛','遷移元素（鉄・銅・銀・クロム・マンガン）','金属イオンの沈殿・錯イオン・系統分析','気体の製法と性質','有機化合物の構造と元素分析','脂肪族炭化水素','アルコール・アルデヒド・ケトン','カルボン酸・エステル・油脂','芳香族化合物','有機化合物の検出反応と分離','高分子の基礎と合成樹脂','合成繊維','ゴム・機能性高分子','糖類','アミノ酸・タンパク質・核酸']
FORMATS = {'qa','equation','compare','condition','pathway','calc','structure'}

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0:
            return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    arrows = re.findall(r'<=>|->|<-', eq)
    return len(arrows) == 1 and (' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None)

def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if should_check(b):
                r = check_equation(clean_ce(b))
                assert r is None, (b, r)
                n += 1
    print('反応式', n, '本 つり合いOK')

def check_structure():
    assert it['field'] in FIELDS, it['field']
    assert set(it['format'].split(',')) <= FORMATS and 1 <= len(it['format'].split(',')) <= 2
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items:
        assert [int(x) for x in n_ans] == list(range(1, n_items + 1)), (n_items, n_ans)
    else:
        assert not n_ans
    for k in ('problem_fig', 'solution_fig'):
        if it[k]:
            assert os.path.exists(os.path.join(os.path.dirname(ITEM), it[k])), it[k]
    print('小問数', n_items, '答え', len(n_ans), '単元・形式OK')

KJS = r'''
const path='/home/claude/kurowasan3/alchemia/vendor/';
const Module=require('module');const r0=Module._resolveFilename;
Module._resolveFilename=function(r,...a){return r==='katex'?path+'katex.min.js':r0.call(this,r,...a)};
const katex=require(path+'katex.min.js');require(path+'mhchem.min.js');
const it=JSON.parse(require('fs').readFileSync(process.argv[1],'utf8'));let bad=0,n=0;
for(const k of ['problem','answer','solution']){let s=it[k];const segs=[];
s=s.replace(/\\\[([\s\S]*?)\\\]/g,(_,m)=>{segs.push([m,true]);return ' '});
let i=0;while(i<s.length){if(s[i]==='$'){let j=i+1,d=0;while(j<s.length){if(s[j]==='\\'){j+=2;continue}if(s[j]==='{')d++;else if(s[j]==='}')d--;else if(s[j]==='$'&&d===0)break;j++}
if(j>=s.length){console.log('閉じていない $',k,s.slice(i,i+50));bad++;break}segs.push([s.slice(i+1,j),false]);i=j+1}else i++}
for(const [m,dsp] of segs){n++;try{katex.renderToString(m,{throwOnError:true,displayMode:dsp,strict:false})}catch(e){bad++;console.log(k,e.message.slice(0,150),'|',m.slice(0,80))}}}
console.log('KaTeX',n,'式',bad?'NG '+bad:'OK');process.exit(bad?1:0);
'''

def check_katex():
    r = subprocess.run(['node', '-e', KJS, ITEM], capture_output=True, text=True)
    print(r.stdout.strip())
    assert r.returncode == 0, r.stderr

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

check_structure()
check_all_equations()
check_katex()
# ---- 以下 問題ごとの検証 ----
import math
pk = 7.21; K1 = 7.1e-3; K2 = 10**-pk
X0 = pk; X1 = pk+math.log10(0.048/0.052)
Y0 = pk+math.log10(0.01/0.09); Y1 = pk+math.log10(0.008/0.092)
assert (round(X0, 2), round(X1, 2), round(Y0, 2), round(Y1, 2)) == (7.21, 7.18, 6.26, 6.15)
assert 2.8 < (Y0-Y1)/(X0-X1) < 3.3
# Z を電荷・物質収支で厳密に解く：P 全量 1.0e-3，Na+ 1.5e-3（=5e-4+2*5e-4），Cl- 2.0e-3
def charge(h):
    D = h**3 + K1*h*h + K1*K2*h + K1*K2*4.5e-13
    a1 = K1*h*h/D; a2 = K1*K2*h/D
    return h + 1.5e-3 - (2.0e-3 + 1.0e-3*(a1 + 2*a2) + 1e-14/h)
lo, hi = 1e-6, 1e-1
for _ in range(200):
    m = (lo*hi)**0.5
    if charge(lo)*charge(m) <= 0: hi = m
    else: lo = m
pHZ = -math.log10(m)
b = 9.6e-3; x = (b-math.sqrt(b*b-6e-6))/2; h = 1.5e-3-x
assert round(pHZ, 2) == round(-math.log10(h), 2) == 2.87, (pHZ, h)
r = 10**-0.05; C = 2*2.0e-3*(1+r)/(1-r)
assert sig(C, 2) == 0.070, C
assert max(range(1, 100), key=lambda i: (i/100)*(1-i/100)) == 50
print(X1, Y0, Y1, pHZ, C)
print("ALL OK")
