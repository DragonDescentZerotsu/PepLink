#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Peptide-like SMILES → amino-acid sequence 解析器（区分 L / D）。

更新要点
--------
1. N-端模板兼容中性 NH₂ 和质子化 NH₃⁺：
      [*:1]N  →  [N,N+;H1,H2]
2. 仍保留哑原子 * 参与匹配；GetMolFrags(sanitizeFrags=False)。
3. 按模板原子数降序匹配，避免 Asp → Ala 错配。

适用范围：线性、无环、无复杂修饰的 20 标准氨基酸肽链；其余情况需扩展模板。
"""

from rdkit import Chem
from typing import Dict, List, Tuple, Optional
import re
try:
    import selfies as sf
except ImportError:  # pragma: no cover - optional at runtime
    sf = None

# ------------------------------------------------------------
# 1. 常用工具
# ------------------------------------------------------------

def flip_stereo(smarts: str) -> str:
    """把 SMARTS 中 '@@' ↔ '@'，生成 D 对映体模板。"""
    if '@' not in smarts:
        return smarts
    tmp = smarts.replace('@@', '§§')
    tmp = tmp.replace('@', '@@')
    return tmp.replace('§§', '@')

# ------------------------------------------------------------
# 2. 肽键模式：C(=O)–N–Cα
# ------------------------------------------------------------

PATT_PEPTIDE = Chem.MolFromSmarts(
    '[C;X3,X4:1](=O)[N;!a;X2,X3,X4:2][C;X4;H1,H2:3]'
)

# def peptide_bond_indices(mol: Chem.Mol) -> List[int]:
#     """返回需断开的所有 C–N 单键编号（肽键）。"""
#     idxs = []
#     for c_idx, _, n_idx, _ in mol.GetSubstructMatches(PATT_PEPTIDE, useChirality=False):
#         bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
#         if bond and bond.GetBondType() == Chem.BondType.SINGLE:
#             idxs.append(bond.GetIdx())
#     return idxs

def peptide_bonds(mol: Chem.Mol):
    """
    返回 [(bondIdx, c_idx, n_idx), ...]，
    其中 c_idx 是羰基碳，n_idx 是脊椎 N。
    """
    out = []
    for c_idx, _, n_idx, _ in mol.GetSubstructMatches(PATT_PEPTIDE, useChirality=False):
        bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
        if bond and bond.GetBondType() == Chem.BondType.SINGLE:
            out.append((bond.GetIdx(), c_idx, n_idx))
    return out

# ------------------------------------------------------------
# 2b. 识别 Cα 并判定 L / D  ------------------  # --- MODIFIED ---
# ------------------------------------------------------------

# PATT_CA = Chem.MolFromSmarts('[N;!H0][C;X4:1][C](=O)')   # 捕获残基自身的 Cα
PATT_CA = Chem.MolFromSmarts('[N][C;X4:1][C](=O)')

def _find_ca_idx(frag: Chem.Mol) -> Optional[int]:
    """返回片段中 Cα 的原子索引（若未找到则 None）。"""
    matches = frag.GetSubstructMatches(PATT_CA, useChirality=False)
    return matches[0][1] if matches else None

def _chirality_of_ca(frag: Chem.Mol, ca_idx: int, aa: str) -> Optional[str]:
    """
    根据 CIP ('R'/'S') → 'L'/'D'。
    对绝大多数氨基酸:  L = S, D = R
    对 Cys (含硫):     L = R, D = S
    """
    Chem.AssignStereochemistry(frag, cleanIt=True, force=True)
    atom = frag.GetAtomWithIdx(ca_idx)
    if not atom.HasProp('_CIPCode'):
        return None                # Gly 或无手性信息
    # cip = atom.GetProp('_CIPCode') # 'R' or 'S'
    cip = atom.GetProp('_CIPCode').upper()
    if aa in {'C', 'U'}:                  # Cys 例外
        cip = 'S' if cip == 'R' else 'R'
    return 'L' if cip == 'S' else 'D'

# ------------------------------------------------------------
# 3. 20 种氨基酸内残基 (L)
# ------------------------------------------------------------

AA_INNER_L = {
    'A': '[*:1]N[C@@H](C)C(=O)[*:2]',
    'R': '[*:1]N[C@@H](CCCNC(N)=N)C(=O)[*:2]',
    'B': '[*:1]N[C@@H](CCCN=C(N)N)C(=O)[*:2]',
    'N': '[*:1]N[C@@H](CC(=O)N)C(=O)[*:2]',
    'D': '[*:1]N[C@@H](CC(=O)O)C(=O)[*:2]',
    'C': '[*:1]N[C@@H](CS)C(=O)[*:2]',
    'E': '[*:1]N[C@@H](CCC(=O)O)C(=O)[*:2]',
    'Q': '[*:1]N[C@@H](CCC(=O)N)C(=O)[*:2]',
    'G': '[*:1]NC(C(=O)[*:2])',
    'H': '[*:1]N[C@@H](Cc1c[nH]cn1)C(=O)[*:2]',
    'L': '[*:1]N[C@@H](CC(C)C)C(=O)[*:2]',
    'I': '[*:1]N[C@@H](C(C)CC)C(=O)[*:2]',
    'K': '[*:1]N[C@@H](CCCCN)C(=O)[*:2]',
    'M': '[*:1]N[C@@H](CCSC)C(=O)[*:2]',
    'F': '[*:1]N[C@@H](Cc1ccccc1)C(=O)[*:2]',
    'P': '[*:1]N1CCC[C@@H]1C(=O)[*:2]',
    'S': '[*:1]N[C@@H](CO)C(=O)[*:2]',
    'T': '[*:1]N[C@@H](C(O)C)C(=O)[*:2]',
    'W': '[*:1]N[C@@H](Cc1c[nH]c2ccccc12)C(=O)[*:2]',
    'Y': '[*:1]N[C@@H](Cc1ccc(O)cc1)C(=O)[*:2]',
    'V': '[*:1]N[C@@H](C(C)C)C(=O)[*:2]',
}

# ------------------------------------------------------------
# 4. N-端 / C-端 模板（L）
# ------------------------------------------------------------

def make_nterm(smarts: str) -> str:
    """内残基 SMARTS → N 端：允许中性 N 或 N⁺，带 1–2 H。"""
    return re.sub(r'\[\*\:1\]N', '[N,N+;H1,H2]', smarts, count=1)

def make_cterm(smarts: str) -> str:
    """内残基 SMARTS → C 端：羧基自由。"""
    return smarts.replace('[*:2]', '[O;H1,H0-]')

AA_NTERM_L = {aa: make_nterm(s) for aa, s in AA_INNER_L.items()}
AA_NTERM_L['G'] = '[N,N+;H1,H2]C(C(=O)[*:2])'      # Gly N 端特例

AA_CTERM_L = {aa: make_cterm(s) for aa, s in AA_INNER_L.items()}
AA_CTERM_L['G'] = '[*:1]NC(C(=O)[O;H1,H0-])'            # Gly C 端特例

# ------------------------------------------------------------
# 5. D-系列模板（小写字母）
# ------------------------------------------------------------

def make_D_dict(src: Dict[str, str]) -> Dict[str, str]:
    return {aa.lower(): flip_stereo(s) for aa, s in src.items()}

# ------------------------------------------------------------
# 6. 编译模板并按原子数降序存放
# ------------------------------------------------------------

def compile_patts(base: Dict[str, str]):
    return [
        (aa, patt := Chem.MolFromSmarts(smarts), patt.GetNumAtoms())
        for aa, smarts in base.items()
    ]

# PATT_ORDER: List[Tuple[List[Tuple[str, Chem.Mol, int]], str]] = [
#     (compile_patts(AA_INNER_L),               'L'),
#     (compile_patts(make_D_dict(AA_INNER_L)),  'D'),
#     (compile_patts(AA_NTERM_L),               'L'),
#     (compile_patts(make_D_dict(AA_NTERM_L)),  'D'),
#     (compile_patts(AA_CTERM_L),               'L'),
#     (compile_patts(make_D_dict(AA_CTERM_L)),  'D'),
# ]

PATT_ORDER: List[List[Tuple[str, Chem.Mol, int]]] = [
    compile_patts(AA_INNER_L),
    compile_patts(make_D_dict(AA_INNER_L)),
    compile_patts(AA_NTERM_L),
    compile_patts(make_D_dict(AA_NTERM_L)),
    compile_patts(AA_CTERM_L),
    compile_patts(make_D_dict(AA_CTERM_L)),
]

# ------------------------------------------------------------
# 7. 片段 → 氨基酸字母
# ------------------------------------------------------------

# def frag_to_aa(frag: Chem.Mol) -> str:
#     for patt_list, chir in PATT_ORDER:
#         for aa, patt, natoms in sorted(patt_list, key=lambda x: -x[2]):
#             if frag.HasSubstructMatch(patt, useChirality=True) and frag.GetNumAtoms() == natoms:
#                 return aa if chir == 'L' else aa.lower()
#     return 'X'


def frag_to_aa(frag: Chem.Mol) -> str:                      # --- CHANGED ---
    """
    逻辑：
    1. 先忽略手性，按模板找出是哪一种氨基酸（侧链指纹）
    2. 再用 CIP (R/S) 判定 L / D → 大写 / 小写
    """
    aa_base: Optional[str] = None

    # 1. 侧链识别（忽略手性）
    for patt_list in PATT_ORDER:
        for aa, patt, natoms in sorted(patt_list, key=lambda x: -x[2]):
            if frag.GetNumAtoms() == natoms and frag.HasSubstructMatch(patt, useChirality=False):
                aa_base = aa.upper()      # 先统一成大写
                break
        if aa_base:
            break

    if aa_base is None:
        return 'X'                       # 未识别

    # 2. 判定 L / D
    ca_idx = _find_ca_idx(frag)
    if ca_idx is None:                   # Gly 或非标准结构
        return aa_base

    chir = _chirality_of_ca(frag, ca_idx, aa_base)
    if chir is None:
        return aa_base                   # 找不到 CIP，默认大写

    return aa_base if chir == 'L' else aa_base.lower()

# ------------------------------------------------------------
# 8. SMILES → 序列
# ------------------------------------------------------------

def smiles_to_pepseq(smiles: str) -> Optional[str]:
    if '[C]' in smiles:
        if sf is None:
            return smiles, None
        try:
            smiles = sf.decoder(smiles)
        except:
            return smiles, None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, None
    # bonds = peptide_bond_indices(mol)
    # if not bonds:
    #     return smiles, None
    # frag_mol = Chem.FragmentOnBonds(mol, set(bonds), addDummies=True)
    # frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False)
    # # print([Chem.MolToSmiles(f) for f in frags])
    # smiles = Chem.MolToSmiles(mol)
    # return smiles, ''.join(frag_to_aa(f) for f in frags)

    pbonds = list(set(peptide_bonds(mol)))
    if not pbonds:                        # 单残基
        return smiles, frag_to_aa(mol)

    bond_ids = {b[0] for b in pbonds}
    frag_mol = Chem.FragmentOnBonds(mol, bond_ids, addDummies=True)

    # 1) 片段列表 & 原子→片段映射
    # frags, atom_groups = Chem.GetMolFrags(
    #     frag_mol, asMols=True, sanitizeFrags=False, returnMapping=True
    # )
  # 先拿片段（asMols=True）
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False)
  # 再拿“原子索引列表”——同顺序
    atom_groups = Chem.GetMolFrags(frag_mol, asMols=False, sanitizeFrags=False)
    atom2frag = {}
    for fid, atoms in enumerate(atom_groups):
        for a in atoms:
            atom2frag[a] = fid

    # 2) 构建有向图 frag_c → frag_n
    succ, indeg = {}, {}
    for _, c_idx, n_idx in pbonds:
        f_c = atom2frag[c_idx]
        f_n = atom2frag[n_idx]
        succ[f_c] = f_n
        indeg[f_n] = indeg.get(f_n, 0) + 1

    # # 3) 找 N-端（入度 0）
    # start = next(f for f in range(len(frags)) if indeg.get(f, 0) == 0)
    #
    # # 4) 按链遍历
    # ordered = []
    # cur = start
    # while True:
    #     ordered.append(cur)
    #     cur = succ.get(cur)
    #     if cur is None:
    #         break
    # if len(ordered) != len(frags):        # 检查是否有环或断链
    #     raise ValueError('Peptide not linear or fragmented.')
    #
    # seq = ''.join(frag_to_aa(frags[i]) for i in ordered)
    # return Chem.MolToSmiles(mol), seq

    # 3) 判断是否线性还是头尾环
    starts = [f for f in range(len(frags)) if indeg.get(f, 0) == 0]

    if starts:                           # ----- 线性肽 -----
        start = starts[0]
        cyclic = False
    else:  # 环肽
        cyclic = True
        # 检查是否所有片段都有且只有一个入度
        if len(indeg) != len(frags) or any(count != 1 for count in indeg.values()):
            # raise ValueError('Peptide graph contains branches / multiple cycles.')
            return smiles, None
        # 找到起始片段（原子序号最小）
        start = min(range(len(frags)), key=lambda f: min(atom_groups[f]))

        # 4) 按链遍历
    ordered, visited = [], set()
    cur = start

    # 修复2：处理环肽的连接关系
    if cyclic:
        # 找到连接回起始片段的片段
        last_to_first = None
        for _, c_idx, n_idx in pbonds:
            f_c = atom2frag[c_idx]
            f_n = atom2frag[n_idx]
            if f_n == start:  # 找到连接到起始片段的肽键
                last_to_first = f_c
                break

        if last_to_first is not None:
            succ[last_to_first] = start  # 设置闭环连接

    while cur not in visited:
        ordered.append(cur)
        visited.add(cur)
        next_cur = succ.get(cur)
        if next_cur is None:
            break
        cur = next_cur

    # 修复3：正确处理环肽的完整遍历
    if cyclic and len(ordered) == len(frags) - 1 and last_to_first is not None:
        # 确保最后一个片段被包含
        if last_to_first not in visited:
            ordered.append(last_to_first)

    if len(ordered) != len(frags):
        # raise ValueError('Peptide graph fragmented or branched.')
        return smiles, None

    seq = ''.join(frag_to_aa(frags[i]) for i in ordered)
    if cyclic:
        seq = 'cyclo-' + seq
    return Chem.MolToSmiles(mol), seq

# ------------------------------------------------------------
# 9. 简单 CLI
# ------------------------------------------------------------

if __name__ == '__main__':
    tests = [
        'CC[C@H1](C)[C@H1](NC(=O)[C@H1](CC(C)C)NC(=O)[C@H1](CCCN=C(N)N)NC(=O)[C@H1](C)NC(=O)[C@H1](CC(C)C)NC(=O)[C@H1](CCCCN)NC(=O)[C@@H1]1CCCN1C(=O)[C@H1](CCCN=C(N)N)NC(=O)[C@H1](C)NC(=O)[C@H1](CCCN=C(N)N)NC(=O)[C@@H1](N)CCSC)C(=O)N[C@@H1](CCCN=C(N)N)C(=O)N[C@H1](C(=O)N[C@@H1](CO)C(=O)NCC(=O)N[C@@H1](CC(C)C)C(=O)N[C@H1](C(=O)N[C@@H1](CC(C)C)C(=O)N[C@@H1](CO)C(=O)N[C@@H1](CCCN=C(N)N)C(=O)N[C@@H1](CC(C)C)C(=O)N2CCC[C@H1]2C(=O)N[C@@H1](CCSC)C(=O)N[C@@H1](CCCCN)C(=O)NCC(=O)N3CCC[C@H1]3C(=O)N[C@@H1](C)C(=O)N[C@@H1](C)C(=O)[OH1])C(C)C)C(C)C',          # Gly-Ala-Ala
        'N[C@H](CC9=CNC=N9)C(=O)N5[C@@H](CCC5)C(=O)N[C@@H](CC4=CC=CC=C4)C(=O)O',                # dAla-lAsp
        '[N][C@@H1][Branch1][=Branch2][C][C][=C][N][C][=N][Ring1][Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][Branch2][C][C][C][=Branch1][C][=O][O][C][=Branch1][C][=O][N][C@@H1][Branch1][#Branch1][C][C][=Branch1][C][=O][N][C][=Branch1][C][=O][N][C@@H1][Branch1][Ring1][C][O][C][=Branch1][C][=O][N][C@@H1][Branch1][P][C][C][=Branch1][Ring1][=C][N][C][=C][Ring1][Ring1][C][=C][C][=C][Ring1][=Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][C][C][C][N][C][=Branch1][C][=O][N][C@@H1][Branch1][=N][C][C][=C][C][=C][Branch1][C][O][C][=C][Ring1][#Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][#Branch2][C][C][=C][C][=C][C][=C][Ring1][=Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][Branch2][C][C][C][=Branch1][C][=O][N][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][C][C][Ring1][Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][Branch1][C][C][C][C][=Branch1][C][=O][N][C@@H1][Branch1][#Branch1][C][C][Branch1][C][C][C][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][C][C][C][N][C][=Branch1][C][=O][N][C@@H1][Branch1][Branch1][C][C][S][C][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][C][C][C][N][C][=Branch1][C][=O][N][C@@H1][Branch1][#Branch1][C@H1][Branch1][Ring1][C][C][C][C][=Branch1][C][=O][N][C@@H1][Branch1][#Branch2][C][C][C][N][C][=Branch1][C][=N][N][C][=Branch1][C][=O][O]',               # lVal-dVal
        '[N][C@H1][Branch1][=Branch1][C][C][C][Ring1][Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][C][C][C][N][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][C][C][Ring1][Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][C][C][Ring1][Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][P][C][C][=Branch1][Ring1][=C][N][C][=C][Ring1][Ring1][C][=C][C][=C][Ring1][=Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][Branch2][C][C][C][=Branch1][C][=O][N][C][=Branch1][C][=O][N][C@@H1][Branch1][P][C][C][=Branch1][Ring1][=C][N][C][=C][Ring1][Ring1][C][=C][C][=C][Ring1][=Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][#Branch2][C][C][C][N][C][=Branch1][C][=N][N][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][Branch1][C][C][C][C][=Branch1][C][=O][N][C@@H1][Branch1][#Branch1][C][Branch1][C][C][C][C][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch2][C][C][=C][N][C][=N][Ring1][Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch2][C][C][=C][N][C][=N][Ring1][Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch1][C][C][C][C][N][C][=Branch1][C][=O][N][C][C][=Branch1][C][=O][N][C@@H1][Branch1][Branch2][C][C][C][=Branch1][C][=O][N][C][=Branch1][C][=O][N][C@@H1][Branch1][#C][C][C][=C][NH1][C][=C][Ring1][Branch1][C][=C][C][=C][Ring1][=Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][=Branch2][C][C][=C][N][C][=N][Ring1][Branch1][C][=Branch1][C][=O][N][C@@H1][Branch1][C][C][C][=Branch1][C][=O][N][C@@H1][Branch1][#Branch2][C][C][C][N][C][=Branch1][C][=N][N][C][=Branch1][C][=O][N][C@@H1][Branch1][#Branch2][C][C][C][N][C][=Branch1][C][=N][N][C][Ring2][=N][=Branch2][=O]'                                                   # non-peptide
    ]
    for s in tests:
        smiles, aa_seq = smiles_to_pepseq(s)
        print(f'{smiles:55} → {aa_seq.replace("B", "R`")}')
