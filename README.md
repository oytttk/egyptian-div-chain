# Egyptian Divisibility-Chain Egyptian Fractions

本仓库提供一个 Python 脚本，用于在**分母上界**约束下，为任意有理数 `a/b` 构造埃及分数表示（单位分数和）：

\[
\frac{a}{b} = k + \sum_{i=1}^{n}\frac{1}{d_i},\quad 0\le r=\frac{a}{b}-k <1
\]

其中 `d_i` 为互异正整数分母（脚本输出为严格递增序列），并且目标函数偏向让相邻分母满足**整除关系**：

\[
d_{i+1}\bmod d_i = 0.
\]

---

## 1. 核心约束与优化目标

### 1.1 硬约束：分母上界
脚本强制：

- `max(d_i) <= 2 * b_orig`

其中 `b_orig` 指你输入形式为 `a/b` 时的 `b`（原始分母）；若输入不是 `a/b`（例如小数），则使用 Python `Fraction` 规约后的分母作为 `b_orig`。

> 注意：在该硬约束下，**并非所有分数一定有解**。如果在给定搜索预算内找不到可行解，会明确输出 “No exact Egyptian expansion found ...”。

### 1.2 目标：相邻整除链
脚本通过搜索（回溯 DFS）优化“相邻整除链”的结构，主要指标：

- **longest_divisible_run**：最长连续段，使得段内每一对相邻满足 `d_{i+1} % d_i == 0`
- **divisible_edges**：总共有多少对相邻满足整除（最多是 `n-1`）

并提供两种优化模式：

- `--mode shortest`（默认）：**先最少项数**，在同项数下尽量让整除链更强
- `--mode longest`：**先最大化最长整除段**（可能会产生更多项，因为长链可通过拆分不断延长）

以及一个更强的硬约束：

- `--force-chain`：要求整个分母序列都是整除链（即对所有相邻都整除）。

---

## 2. 环境与安装

### 2.1 Python 版本
- 推荐：Python 3.9+（使用 `fractions.Fraction`，无第三方依赖）

### 2.2 获取脚本
将 `egyptian_div_chain.py` 放在任意目录即可运行。

---

## 3. 快速开始

### 3.1 基本用法
```bash
python egyptian_div_chain.py 13/32
