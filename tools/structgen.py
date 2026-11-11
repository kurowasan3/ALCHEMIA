#!/usr/bin/env python3
"""構造式の SVG を SMILES から作る（RDKit が必要：pip3 install rdkit）

    python3 tools/structgen.py           図鑑（content/chem/struct/）と、問題の図（"smiles" を持つ図）の SVG を作る
    python3 tools/structgen.py --check   作り直したものと今の SVG が同じか（作り忘れの検査。RDKit が無ければ飛ばす）

色は currentColor（ライト・ダークに合わせて変わる）。原子記号は線（path）で描くのでフォントに依存しない。
高分子の繰り返し単位は、SMILES の * を「鎖が続く」印（波線）として描く。
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from content_lib import CONTENT, load_chem, load_problems  # noqa: E402

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    Chem = None


def svg_of(smiles, w=260, h=180, polymer=False):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f'SMILES を読めません：{smiles}')
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    o = d.drawOptions()
    o.useBWAtomPalette()
    o.clearBackground = False
    o.bondLineWidth = 1.6
    o.padding = 0.08
    o.minFontSize = 11
    o.fixedBondLength = 28
    if polymer:
        o.dummiesAreAttachments = True
    d.DrawMolecule(m)
    d.FinishDrawing()
    s = d.GetDrawingText()
    s = re.sub(r'<\?xml[^>]*\?>\s*', '', s)
    s = re.sub(r"<!-- END OF HEADER -->\s*", '', s)
    s = s.replace('#000000', 'currentColor').replace('#000', 'currentColor')
    s = re.sub(r"xmlns:rdkit='[^']*'\s*", '', s)
    s = re.sub(r"xmlns:xlink='[^']*'\s*", '', s)
    s = re.sub(r" class='[^']*'", '', s)
    s = re.sub(r"<rect[^>]*style='opacity:1.0;fill:#FFFFFF[^>]*>\s*</rect>", '', s)
    s = re.sub(r"\s+width='\d+px'\s+height='\d+px'", '', s, count=1)
    s = s.replace('<svg ', f'<svg viewBox="0 0 {w} {h}" ', 1) if 'viewBox' not in s else s
    s = re.sub(r'\n\s*\n', '\n', s)
    return s.strip() + '\n'


def targets(content=CONTENT):
    """[(出力先, SMILES, 高分子か)]"""
    C = load_chem(content)
    st = Path(content) / 'chem' / 'struct'
    out = []
    for c in C['organic']['compounds']:
        if c.get('smiles'):
            out.append((st / f'{c["id"]}.svg', c['smiles'], False))
    for p in C['polymers']['polymers']:
        if p.get('repeat'):
            out.append((st / f'poly-{p["id"]}.svg', p['repeat'], True))
        for i, mo in enumerate(p['monomers']):
            if mo.get('smiles'):
                out.append((st / f'mono-{p["id"]}-{i + 1}.svg', mo['smiles'], False))
    for d, pr in load_problems(content):
        for fg in pr['figures']:
            if fg.get('smiles'):
                out.append((d / fg['src'], fg['smiles'], '*' in fg['smiles']))
    return out


def main():
    check = '--check' in sys.argv
    if Chem is None:
        if check:
            print('RDKit が無いので構造式の検査を飛ばしました（pip3 install rdkit）')
            return
        sys.exit('RDKit が必要です：pip3 install rdkit')
    bad, n = [], 0
    for path, smi, poly in targets():
        s = svg_of(smi, polymer=poly)
        n += 1
        if check:
            if not path.exists() or path.read_text(encoding='utf-8') != s:
                bad.append(str(path.relative_to(CONTENT)))
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists() or path.read_text(encoding='utf-8') != s:
                path.write_text(s, encoding='utf-8')
    if check:
        print(f'構造式 {n}枚を検査、作り直しが必要なもの {len(bad)}枚')
        for b in bad[:20]:
            print('  ', b)
        sys.exit(1 if bad else 0)
    print(f'構造式 {n}枚を書き出しました')


if __name__ == '__main__':
    main()
