<div align="center">

# 🔄 Flip Vacuum

**一键翻转 POSCAR / CONTCAR / CIF 模型的真空层方向**

把 `c` 轴真空层从一侧「换边」到另一侧 —— 面向第一性原理计算（VASP / DFT）的小工具，带原生图形界面。

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey.svg)](#-安装)
[![Python](https://img.shields.io/badge/python-3.8%2B-3776AB.svg?logo=python&logoColor=white)](#-安装)
[![GUI](https://img.shields.io/badge/GUI-Tkinter-FF9800.svg)](#-图形界面)
[![Release](https://img.shields.io/github/v/release/moyulyy/flip-vacuum?label=latest&color=green)](https://github.com/moyulyy/flip-vacuum/releases)

</div>

---

## 📖 这是什么？

在做**表面 / 二维材料 / 异质结**等第一性原理计算时，我们经常需要调整 `POSCAR` 中真空层的位置。例如把底部的真空层换到顶部，以便搭建新的界面结构。

**Flip Vacuum** 做的事情非常简单：

```
z → 1 - z
```

也就是对分数坐标做一次**镜像翻转**，将真空层从一侧搬到另一侧。支持 `POSCAR`、`CONTCAR`、`*.vasp` 与 `*.cif` 格式，可选同时翻转 `a`、`b` 方向。

> ✨ 无需安装任何依赖即可处理 VASP 文件；处理 CIF 时可选安装 [ASE](https://wiki.fysik.dtu.dk/ase/)。

---

## ✨ 特性

- 🖱️ **原生图形界面** —— 鼠标选文件即可，无需记忆命令
- 🧮 **零依赖处理 VASP 格式** —— POSCAR / CONTCAR 读写不依赖 `numpy`、`ase`
- 🧬 **支持 CIF 格式** —— 借助 ASE 展开对称性后按 `P1` 写出
- 🎛️ **可选翻转 a / b 方向** —— `x → 1-x`、`y → 1-y`
- 📝 **保留原始格式信息** —— 元素名、选择性动力学标记（Selective dynamics）等
- 📦 **可打包为单个 exe** —— 双击即用，见 [Releases](https://github.com/moyulyy/flip-vacuum/releases)
- 🖥️ **同时支持命令行** —— 适合脚本批处理

---

## 🖼️ 界面预览

图形界面采用 iOS 风格的圆角弹窗设计：

```
┌─────────────────────────────────────────────┐
│             翻转真空层方向                    ✕ │
│        POSCAR / CIF · c 轴真空层换边          │
├─────────────────────────────────────────────┤
│  结构文件                              0 个文件 │
│  ┌─────────────────────────────────────┐   │
│  │ 尚未添加文件                          │   │
│  │ 点击「添加文件」选择 POSCAR / CIF      │   │
│  └─────────────────────────────────────┘   │
│  [ ＋ 添加文件 ]  清空列表                   │
├─────────────────────────────────────────────┤
│  选项                                        │
│  翻转 c 方向（z → 1-z）                [ON]  │
│  同时翻转 a、b 方向                    [OFF] │
├─────────────────────────────────────────────┤
│  运行日志                                    │
│  就绪：点击「添加文件」选择 POSCAR / CIF。     │
├─────────────────────────────────────────────┤
│           [        开始处理        ]          │
└─────────────────────────────────────────────┘
```

---

## 📦 安装

### 方式一：直接下载 exe（Windows，推荐）

前往 [Releases](https://github.com/moyulyy/flip-vacuum/releases) 下载 `FlipVacuum.exe`，双击即可运行，**无需安装 Python**。

### 方式二：Python 源码运行

```bash
# 克隆仓库
git clone https://github.com/moyulyy/flip-vacuum.git
cd flip-vacuum

# 仅处理 VASP 文件（POSCAR / CONTCAR）：无需任何第三方依赖
python flip_vacuum.py

# 如需处理 CIF 文件，安装 ASE
pip install ase
```

> 双击 `run.bat`（Windows）也可以直接启动图形界面。

---

## 🚀 使用方法

### 🖱️ 图形界面

1. 运行程序，点击 **「＋ 添加文件」** 选择 `POSCAR` / `CONTCAR` / `*.vasp` / `*.cif`；
2. 按需打开 / 关闭 **「翻转 c 方向」** 与 **「同时翻转 a、b 方向」** 开关；
3. 点击 **「开始处理」**。

输出文件保存在**输入文件同一目录**，文件名自动加后缀 `_flip`：

| 输入文件 | 输出文件 |
| --- | --- |
| `POSCAR` | `POSCAR_flip` |
| `slab.cif` | `slab_flip.cif` |
| `CONTCAR` | `CONTCAR_flip` |

### 🖥️ 命令行

```bash
# 翻转一个或多个文件（默认翻转 c 方向）
python flip_vacuum.py POSCAR
python flip_vacuum.py POSCAR slab.cif CONTCAR

# 同时翻转 a、b、c 三个方向
python flip_vacuum.py --flip-ab POSCAR

# 指定输出文件（仅单个输入时有效）
python flip_vacuum.py POSCAR -o POSCAR_new
```

| 参数 | 说明 |
| --- | --- |
| `files` | 一个或多个输入文件（不带参数时启动 GUI） |
| `--flip-ab` | 同时翻转 a、b 方向（`x → 1-x`、`y → 1-y`） |
| `-o, --output` | 指定输出文件路径（仅单个输入时有效） |

---

## 🔬 工作原理

对晶胞分数坐标做如下变换（默认仅翻转 `c` 方向）：

```
(x, y, z)  →  (x, y, 1 - z)
```

翻转后再将坐标回卷到 `[0, 1)` 区间。由于是**分数坐标**的翻转，晶格矢量保持不变，因此：

- ✅ 原子相对位置、键长、键角**完全不变**
- ✅ 晶胞体积、对称性不变
- ✅ 只有原子沿 `c` 方向「倒过来」，真空层随之换边

若开启 `--flip-ab`：

```
(x, y, z)  →  (1 - x, 1 - y, 1 - z)
```

对于笛卡尔坐标格式的 POSCAR（`Cartesian`），程序会先用晶格矩阵的逆矩阵转换到分数坐标，翻转后再写回。

---

## 🗂️ 支持的文件格式

| 格式 | 说明 | 依赖 |
| --- | --- | --- |
| `POSCAR` / `CONTCAR` | VASP 结构文件（VASP 4 / 5） | ✅ 无 |
| `*.vasp` | 同 POSCAR 格式 | ✅ 无 |
| `*.cif` | 晶体学信息文件（多帧支持） | 需要 [ASE](https://wiki.fysik.dtu.dk/ase/) |

---

## 🛠️ 打包为 exe

```bash
pip install pyinstaller
pyinstaller --onefile --windowed --name FlipVacuum flip_vacuum.py
```

生成的 `dist/FlipVacuum.exe` 即为可独立分发的单文件程序。

---

## 📄 许可证

本项目使用 [MIT License](LICENSE)。© 2026 [moyulyy](https://github.com/moyulyy)
