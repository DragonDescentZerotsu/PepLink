![PepLink 展示图](assets/Peplink.png)

<p align="center">
  <a href="https://huggingface.co/spaces/Kiria-Nozan/PepLink">
    <img alt="Hugging Face logo" src="https://huggingface.co/datasets/huggingface/brand-assets/resolve/main/hf-logo.svg" height="44">
  </a>
</p>

<p align="center">
  <a href="https://huggingface.co/spaces/Kiria-Nozan/PepLink"><strong>在线体验 PepLink Hugging Face Space</strong></a><br>
  如果你想直接在线使用，可以一键打开部署在 Hugging Face Spaces 上的 Gradio 应用。
</p>

# PepLink

PepLink 是一个面向 peptide 与分子字符串表示互转的 Python 包 (`SMILES`/`SELFIES`)。

## 一眼看懂支持范围

在进入 API 细节之前，可以先用下面这组数字快速理解 PepLink v1 的能力边界：

- 420 个内置 unusual amino acid，也就是 non-canonical residue 映射
- 296 种内置 terminal modification
- 241 种 N 端修饰
- 55 种 C 端修饰
- 3 类环肽连接拓扑
- 11 种已实现的 intrachain bond 化学类型

这 3 类环肽连接拓扑分别是：

- `SSB`：sidechain-sidechain cyclization
- `SMB`：sidechain-mainchain cyclization
- `MMB`：mainchain-mainchain cyclization，其中包括 head-to-tail macrocyclization

当前 v1 聚焦在一个边界清晰、可靠的范围内，也就是 peptide 与 `SMILES`/`SELFIES` 的互转：

- `aa_seqs_to_smiles(...)`：单链 peptide 定义 -> `SMILES` 或 `SELFIES`
- `smiles_to_aa_seqs(...)`：标准氨基酸 peptide 的 `SMILES` 或 `SELFIES` -> peptide sequence
- `list_supported_noncanonical_aas(...)`：查看内置和用户注册的 non-canonical amino acid 映射
- `register_noncanonical_aa(...)` / `register_noncanonical_aas(...)`：为当前 Python 进程注册自定义 non-canonical amino acid
- `load_noncanonical_aas_from_csv(...)` / `register_noncanonical_aas_from_csv(...)`：从用户 CSV 文件读取自定义 non-canonical amino acid

## 安装

```bash
pip install PepLink
```

运行时依赖：

- `rdkit`
- `selfies`

## 快速开始

如果你想交互式运行本 README 里的示例，可以直接打开 [`examples/quick_start.ipynb`](examples/quick_start.ipynb)。

### `aa_seqs_to_smiles(...)`

```python
from PepLink import aa_seqs_to_smiles

smiles = aa_seqs_to_smiles(
    "RRXXRF",
    unusual_amino_acids=[
        {"position": 3, "name": "1-NAL"},
        {"position": 4, "name": "1-NAL"},
    ],
    n_terminal="ACT",
    c_terminal="AMD",
)

print(smiles)
```

### `smiles_to_aa_seqs(...)`

```python
from PepLink import smiles_to_aa_seqs

result = smiles_to_aa_seqs("C[C@H](N)C(=O)N[C@@H](CS)C(=O)O")

print(result.sequence)            # AC
print(result.is_cyclic)           # False
print(result.cyclization)         # linear
print(result.unsupported_reason)  # None
```

### non-canonical amino acid 注册表

```python
from PepLink import (
    aa_seqs_to_smiles,
    list_supported_noncanonical_aas,
    register_noncanonical_aa,
)

supported = list_supported_noncanonical_aas()
print(supported["1-NAL"])

register_noncanonical_aa("MyAA", "N[C@@H](CC)C(=O)O")

smiles = aa_seqs_to_smiles(
    "AXA",
    unusual_amino_acids=[{"position": 2, "name": "MyAA"}],
)
print(smiles)
```

```python
from PepLink import register_noncanonical_aas_from_csv

register_noncanonical_aas_from_csv("examples/example_custom_noncanonical_aas.csv")
```

## 支持范围

### `aa_seqs_to_smiles(...)`

PepLink v1 的 `aa_seqs_to_smiles(...)` 支持单链 `monomer` peptide，且可包含：

- 20 个标准氨基酸，以及用小写单字母表示的 D-构型标准氨基酸
- 420 个内置 non-canonical amino acid 映射
- `all_peptides_data.json` 中出现的全部 241 种 N 端修饰
- `all_peptides_data.json` 中出现的全部 55 种 C 端修饰
- 3 类环肽连接拓扑：`SSB`、`SMB` 和 `MMB`
- 11 种已实现的链内成环或交联类型

支持的 intrachain bond 类型：

- `DSB`
- `AMD`
- `TIE`
- `DCB`
- `EST`
- `AMN`
- `p-XylB`
- `TRZB`
- `(E)-but-2-enyl-B`
- `BisMeBn-B`
- `but-2-ynyl-B`

这些支持的 intrachain bond 缩写含义如下：

| Bond | 全称 | 含义 |
| --- | --- | --- |
| `DSB` | Disulfide Bond | 两个半胱氨酸硫原子之间形成的共价 `S-S` 键。 |
| `AMD` | Amide Bond | 羧基与氮形成的酰胺键；在 peptide 里这类键带有部分双键性质，所以这里的 `C-N` 键不能自由旋转。 |
| `TIE` | Thioether Bond | 一般结构形式为 `R-S-R'` 的硫醚键。 |
| `DCB` | Dicarbon Bond (C=C) | 一个碳碳双键 `C=C` 形式的交联。 |
| `EST` | Ester Bond | 由羧基和羟基形成的酯键。 |
| `AMN` | Amine Bond | 泛指含有氨基或胺基，比如 `-NH2`、`-NH-`、`-N-` 的化学键。 |
| `p-XylB` | para-Xylene thioether bridge | 以对二甲苯为桥、并通过硫原子连接两侧残基的硫醚桥。 |
| `TRZB` | Triazole bridge | 通过三唑环桥实现的侧链-侧链连接。 |
| `(E)-but-2-enyl-B` | (E)-but-2-enyl bridge | 由 `(E)-but-2-enyl` 基团桥接的侧链-侧链交联，包含一个 `C=C` 单元。 |
| `BisMeBn-B` | Bismethylenebenzene bridge | 通过带两个亚甲基连接位点的苯环实现的侧链-侧链交联。 |
| `but-2-ynyl-B` | but-2-ynyl bridge | 由 `but-2-ynyl` 基团桥接的侧链-侧链交联，包含一个碳碳三键。 |

示例里常见的 `chain_participating` 缩写：

- `SSB`：Sidechain-Sidechain Bond
- `MMB`：Mainchain-Mainchain Bond
- `SMB`：Sidechain-Mainchain Bond

### `smiles_to_aa_seqs(...)`

PepLink v1 的逆向解析是刻意保守的。

正式支持的范围：

- 仅标准氨基酸
- 支持 L/D 构型
- 线性肽
- 头尾成环肽
- `SMILES` 输入
- `SELFIES` 输入

以下情况不承诺支持逆向解析：

- non-canonical amino acids
- 侧链交联型环肽
- 带 terminal modification 的 peptide
- 配位复合物

当输入超出上述可靠范围时，`smiles_to_aa_seqs(...)` 会返回带有 `unsupported_reason` 的 `PeptideParseResult`。

## 公开 API

### `aa_seqs_to_smiles(...)`

```python
aa_seqs_to_smiles(
    sequence,
    *,
    unusual_amino_acids=None,
    intrachain_bonds=None,
    n_terminal=None,
    c_terminal=None,
    output_format="smiles",
    aa_overrides=None,
    n_terminal_overrides=None,
    c_terminal_overrides=None,
) -> str
```

关键约定：

- `sequence` 使用单字母氨基酸编码
- non-canonical residue 在序列中用 `X` 或 `x` 占位
- `unusual_amino_acids` 的位置必须与 `X/x` 占位位置完全一致
- `intrachain_bonds` 可以使用简化 dict，也可以直接传 DBAASP 风格嵌套 dict
- `output_format` 只能是 `"smiles"` 或 `"selfies"`

最小直接调用示例：

### 线性 peptide

数据集示例：`id=11`

```python
from PepLink import aa_seqs_to_smiles

smiles = aa_seqs_to_smiles("RVKRVWPLVIRTVIAGYNLYRAIKKK")
```

### 单个 non-canonical residue

数据集示例：`id=151`

```python
smiles = aa_seqs_to_smiles(
    "GIKEXKRIVQRIKDFLRNLV",
    unusual_amino_acids=[
        {"position": 5, "name": "Phg"},
    ],
)
```

### 多个 non-canonical residues

数据集示例：`id=157`

```python
smiles = aa_seqs_to_smiles(
    "GRFKRXRKKXKKLFKKIS",
    unusual_amino_acids=[
        {"position": 6, "name": "Phg"},
        {"position": 10, "name": "Phg"},
    ],
)
```

### 端基修饰

数据集示例：`id=10360`

```python
smiles = aa_seqs_to_smiles(
    "K",
    n_terminal="C16",
    c_terminal="AMD",
)
```

另一个包含 D-氨基酸的小真实例子是 `id=8`：

```python
smiles = aa_seqs_to_smiles(
    "KVvvKWVvKvVK",
    n_terminal="C16",
    c_terminal="AMD",
)
```

### Intrachain bond 示例

下面每一种 bond 类型都对应 `all_peptides_data.json` 里的真实记录。

#### `DSB`

数据集示例：`id=57`

```python
smiles = aa_seqs_to_smiles(
    "VTCDILSVEAKGVKLNDAACAAHCLFRGRSGGYCNGKRVCVCR",
    intrachain_bonds=[
        {"position1": 3, "position2": 34, "type": "DSB", "chain_participating": "SSB"},
        {"position1": 20, "position2": 40, "type": "DSB", "chain_participating": "SSB"},
        {"position1": 24, "position2": 42, "type": "DSB", "chain_participating": "SSB"},
    ],
)
```

#### `AMD` 头尾成环

数据集示例：`id=105`

```python
smiles = aa_seqs_to_smiles(
    "SwFkTkSk",
    intrachain_bonds=[
        {"position1": 1, "position2": 8, "type": "AMD", "chain_participating": "MMB"},
    ],
)
```

#### `TIE`

数据集示例：`id=1079`

```python
smiles = aa_seqs_to_smiles(
    "IXSIXLCTPGCKTGALMGCNMKTATCHCSIHVXK",
    unusual_amino_acids=[
        {"position": 2, "name": "DHB"},
        {"position": 5, "name": "DHA"},
        {"position": 33, "name": "DHA"},
    ],
    intrachain_bonds=[
        {"position1": 3, "position2": 7, "type": "TIE", "chain_participating": "SSB"},
        {"position1": 8, "position2": 11, "type": "TIE", "chain_participating": "SSB"},
        {"position1": 13, "position2": 19, "type": "TIE", "chain_participating": "SSB"},
        {"position1": 23, "position2": 26, "type": "TIE", "chain_participating": "SSB"},
        {"position1": 25, "position2": 28, "type": "TIE", "chain_participating": "SSB"},
    ],
)
```

#### `DCB`

数据集示例：`id=4419`

```python
smiles = aa_seqs_to_smiles(
    "FLPILASLAAKFGPKLFXLVTKKX",
    unusual_amino_acids=[
        {"position": 18, "name": "AGL"},
        {"position": 24, "name": "AGL"},
    ],
    intrachain_bonds=[
        {"position1": 18, "position2": 24, "type": "DCB", "chain_participating": "SSB"},
    ],
)
```

#### `EST`

数据集示例：`id=6917`

```python
smiles = aa_seqs_to_smiles(
    "SadAssX",
    unusual_amino_acids=[
        {"position": 7, "name": "D-Allo-Thr"},
    ],
    n_terminal="3,4-OH-4-Me-C16",
    intrachain_bonds=[
        {"position1": 0, "position2": 7, "type": "EST", "chain_participating": "MMB"},
    ],
)
```

#### `AMN`

数据集示例：`id=19104`

```python
smiles = aa_seqs_to_smiles(
    "CANSCXYGPLTWSCXGNTK",
    unusual_amino_acids=[
        {"position": 6, "name": "DHA"},
        {"position": 15, "name": "3-OH-Asp"},
    ],
    intrachain_bonds=[
        {"position1": 1, "position2": 18, "type": "TIE", "chain_participating": "SSB"},
        {"position1": 5, "position2": 11, "type": "TIE", "chain_participating": "SSB"},
        {"position1": 4, "position2": 14, "type": "TIE", "chain_participating": "SSB"},
        {"position1": 6, "position2": 19, "type": "AMN", "chain_participating": "SSB"},
    ],
)
```

#### `p-XylB`

数据集示例：`id=11913`

```python
smiles = aa_seqs_to_smiles(
    "cWkKkC",
    c_terminal="AMD",
    intrachain_bonds=[
        {"position1": 1, "position2": 6, "type": "p-XylB", "chain_participating": "SSB"},
    ],
)
```

#### `TRZB`

数据集示例：`id=14660`

```python
smiles = aa_seqs_to_smiles(
    "FKXRRWQWRMKKLGAPSITXVRRAF",
    unusual_amino_acids=[
        {"position": 3, "name": "BisHomo-Pra"},
        {"position": 20, "name": "Lys(N3)"},
    ],
    intrachain_bonds=[
        {"position1": 3, "position2": 20, "type": "TRZB", "chain_participating": "SSB"},
    ],
)
```

#### `(E)-but-2-enyl-B`

数据集示例：`id=17263`

```python
smiles = aa_seqs_to_smiles(
    "KFFKKLKKAVKKGFKKFAKV",
    intrachain_bonds=[
        {"position1": 4, "position2": 8, "type": "(E)-but-2-enyl-B", "chain_participating": "SSB"},
    ],
)
```

#### `BisMeBn-B`

数据集示例：`id=17273`

```python
smiles = aa_seqs_to_smiles(
    "KFFKKLKKAVKKGFKKFAKV",
    intrachain_bonds=[
        {"position1": 12, "position2": 16, "type": "BisMeBn-B", "chain_participating": "SSB"},
    ],
)
```

#### `but-2-ynyl-B`

数据集示例：`id=19191`

```python
smiles = aa_seqs_to_smiles(
    "VKRFKKFFRKFKKFV",
    c_terminal="AMD",
    intrachain_bonds=[
        {"position1": 6, "position2": 10, "type": "but-2-ynyl-B", "chain_participating": "SSB"},
    ],
)
```

### `smiles_to_aa_seqs(...)`

```python
smiles_to_aa_seqs(text, *, input_format="auto") -> PeptideParseResult
```

返回字段：

- `sequence`
- `is_cyclic`
- `cyclization`
- `normalized_smiles`
- `input_format`
- `unsupported_reason`

示例：

```python
from PepLink import aa_seqs_to_smiles, smiles_to_aa_seqs

linear_smiles = aa_seqs_to_smiles("AC")
print(smiles_to_aa_seqs(linear_smiles))
```

```python
head_to_tail_smiles = aa_seqs_to_smiles(
    "SwFkTkSk",
    intrachain_bonds=[
        {"position1": 1, "position2": 8, "type": "AMD", "chain_participating": "MMB"},
    ],
)
print(smiles_to_aa_seqs(head_to_tail_smiles))
```

对于头尾环肽，返回的序列会被规范化为一个固定 rotation，因为闭环本身没有唯一的起始残基。

### non-canonical amino acid 注册表

```python
list_supported_noncanonical_aas(*, include_custom=True) -> dict[str, str]
load_noncanonical_aas_from_csv(csv_path) -> dict[str, str]
register_noncanonical_aa(name, smiles) -> str
register_noncanonical_aas(mapping) -> dict[str, str]
register_noncanonical_aas_from_csv(csv_path) -> dict[str, str]
clear_registered_noncanonical_aas() -> None
```

关键约定：

- `list_supported_noncanonical_aas(...)` 只返回 non-canonical residue 的 `name -> SMILES` 映射
- 默认内置映射里包含 420 个 non-canonical residue 名称
- CSV 辅助函数要求列名包含 `name` 或 `aa`，以及 `SMILES`
- `register_noncanonical_aa(...)` 会校验并规范化输入的 `SMILES`
- 注册后的映射只在当前 Python 进程内生效，`aa_seqs_to_smiles(...)` 会自动读取
- 如果你只想做单次调用覆盖，而不是修改进程级注册表，依然可以使用 `aa_overrides`

## DBAASP 辅助函数

如果你的输入数据已经是 `all_peptides_data.json` 这种 DBAASP 风格结构，可以直接使用 `from_dbaasp_record(...)`。

```python
import json
from pathlib import Path

from PepLink import aa_seqs_to_smiles, from_dbaasp_record

records = json.loads(Path("all_peptides_data.json").read_text())
record = next(item for item in records if item["id"] == 57)

inputs = from_dbaasp_record(record)
smiles = aa_seqs_to_smiles(**inputs.to_api_kwargs())
```

## 数据集兼容性

`all_peptides_data.json` 是这个仓库当前实现能力的参考数据集。

当前覆盖情况：

- 数据集中的 N 端修饰：`241 / 241` 已内置
- 数据集中的 C 端修饰：`55 / 55` 已内置
- 数据集中的 unusual amino acid 名称：`420 / 545` 已内置
- 尚未收录的 unusual amino acid 名称：`125`

## 不支持的情况

PepLink v1 会明确拒绝以下类别，而不是静默给出不可靠结果。

- 多链 peptide 与 interchain bond
- coordination bond
- 带 non-canonical residue、terminal modification、侧链交联的 peptide 的逆向解析
- 尚未实现的 intrachain bond 类型：`ETH`、`CAR`、`IMN`

真实数据示例：

- multimer / interchain bond: `id=1`
- coordination bond: `id=15`
- unsupported bond types appear in records such as `id=17389` and `id=21130`
- 一个目前在 v1 中仍会失败的已知正向边界案例：`id=5779`

## 扩展映射

你可以在不修改 PepLink 源码的情况下补充内置映射。

### 为当前进程注册自定义 unusual amino acids

```python
from PepLink import register_noncanonical_aas

register_noncanonical_aas(
    {
        "MyAA": "N[C@@H](CC)C(=O)O",
        "MyAA2": "N[C@@H](CO)C(=O)O",
    }
)
```

### 从 CSV 文件为当前进程注册自定义 unusual amino acids

示例文件：[`examples/example_custom_noncanonical_aas.csv`](examples/example_custom_noncanonical_aas.csv)

```python
from PepLink import register_noncanonical_aas_from_csv

register_noncanonical_aas_from_csv("examples/example_custom_noncanonical_aas.csv")
```

### 为单次调用补充缺失的 unusual amino acids

```python
smiles = aa_seqs_to_smiles(
    "AXA",
    unusual_amino_acids=[{"position": 2, "name": "MyAA"}],
    aa_overrides={"MyAA": "N[C@@H](CC)C(=O)O"},
)
```

### 添加 terminal modifications

```python
smiles = aa_seqs_to_smiles(
    "AK",
    n_terminal="MyNCap",
    c_terminal="MyCTail",
    n_terminal_overrides={"MyNCap": "CC(=O)O"},
    c_terminal_overrides={"MyCTail": "N"},
)
```

## 说明

- 正向 `SELFIES` 输出现在已经通过公开 API 支持。
- 逆向解析的能力边界刻意比正向生成功能更窄。
- 现在真正的运行时实现已经全部放进 `PepLink/` 包内部。
- 自定义 non-canonical amino acid 注册是当前进程内的运行时状态。

## 引用

如果你觉得这个项目有用，请引用：

```bibtex
@article{leng2025predicting,
  title={Predicting and generating antibiotics against future pathogens with ApexOracle},
  author={Leng, Tianang and Wan, Fangping and Torres, Marcelo Der Torossian and de la Fuente-Nunez, Cesar},
  journal={arXiv preprint arXiv:2507.07862},
  year={2025}
}
```
