# 进度 / 断点

中断后只看本文「当前断点」和下一小段。目标仍以 [PLAN.md](PLAN.md) 为准。

## 当前断点

| 项 | 值 |
| --- | --- |
| 阶段 | G（可认真调用的 1.x） |
| 小段 | **G2 已完成**；下一刀 **G3** |
| 不要做 | 对 \(x_j\) 退火、搜索宇宙学 \(z\)、`backend=`、用真星系当测试真理、未完成 G 就开 H、再长 `easyppxf` |

## 小段清单

### v0.1 基线

- [x] A 框架、前向模型、固定 \(v,\sigma=0\) 时收回 \(x,A_V\)；LOSVD 在均匀 lnλ
- [x] **B1** 吸收线合成器；告诉真值 \(v,\sigma\) 时仍收回 \(x,A_V\)
- [x] **B2** 运动学网格搜索
- [x] **B3** 真值不落网格点时收回
- [x] **C 跳过**（不对 \(x_j\) 退火；非线性加密改 G5）
- [x] D clip + EX0
- [x] E `.cxt` / gzip / 批量 / 可选 SDSS FITS
- [x] F 冻结 `easyppxf`

### 阶段 G（1.x）

- [x] **G1** `FitResult` + 防错
- [x] G2 已知 \(z\) + 真空/空气
- [ ] G3 LSF + 发射线 mask 表
- [ ] G4 缺误差时 RMS 估计
- [ ] G5 非线性局部加密
- [ ] G6 `light_to_mass` / 光加权产品
- [ ] G7 README 用法清单

### 阶段 H

未开始。完成 G 之前不准做。

## G1 验收（先写再做，2026-09-02）

走 `fit_spectrum`，断言在 `tests/starlightpy/test_starlight_g1.py`。现有合成收回测试必须仍绿（字段只增不删）。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 字段 | 无噪声、少模板、固定 \(v,\sigma=0\) 的现有收回路径 | `FitResult` 含 `good`、`obs_scale`、`config`、`dropped`；`good` 长度等于波长；默认 `dropped` 全 False |
| 坏输入 | flux 含 NaN；或 `norm_window` 与波长无交 | `ValueError`（或同类明确错误），不要给出 χ² |
| 警告 | `search_kinematics=False` 且 `sigma_kms=0` | 发出 `warning`（`UserWarning` 即可），拟合仍可进行 |
| 回归 | `pytest tests/starlightpy` | 原先收回 / clip / IO 测试全绿 |

失败则记本段「未通过」原因，不开始 G2。

**G1 结果（2026-09-02）：通过。** 字段、NaN/空窗口、σ=0 警告均达验收。

## G2 验收（先写再做，2026-09-02）

`preprocess.py`：`to_rest_frame`、`air_to_vacuum` / `vacuum_to_air`。`FitConfig.redshift` 只应用已知 \(z\)，不搜索。`wave_frame` 为 `air` 时把观测波长转到真空。模板视为静止系、真空。

走公开 API。测试用假红移 \(z=0.02\)，无噪声吸收线合成谱。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 静止系 | 把静止系合成谱移到观测系，再 `to_rest_frame` + `resample_to` 后 `fit_spectrum` | \(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12 |
| 对照 | 同一观测系谱**不**移回静止系就拟合 | \(x\) 或 \(A_V\) 超出上表，或 χ² 明显更差 |
| 空气/真空 | 真空 λ 来回转换 | 相对误差 \< \(10^{-6}\)；`wave_frame="air"` 的谱能收回 |

未开始 G2 代码。失败不开始 G3。

**G2 结果（2026-09-02）：通过。** 静止系对照、空气/真空来回、`redshift`/`wave_frame` 均达验收。

## G3 验收（先写再做，2026-09-02）

`fwhm_data` / `fwhm_template` 单位为 **Å**（仪器 FWHM，在波长方向高斯抹平更锐的一方）。仅当两者都给出且 `fwhm_data > fwhm_template` 时把模板展宽到数据；不能反向锐化。光学发射线 mask 表返回 STARLIGHT 风格窗口，给 `apply_mask` 用。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| LSF | 无噪声合成后把观测再展宽 \(\mathrm{FWHM}=4\)Å，模板 2Å；打开对齐 | \(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12 |
| 对照 | 同一展宽观测不对齐 | χ² 大于对齐后 |
| 线表 | `optical_emission_mask_regions` + `apply_mask` | Hα 6563 附近为 False |

未开始 G3 代码。

## B2 验收（历史，2026-09-01）

实现：`search_kinematics=True` 时外层网格扫 \(v,\sigma\)，每个点内层仍 \(A_V\)+NNLS；`False` 时行为与 B1 相同。

合成：`default_absorption_bases`，真值 \(x=(0.50,0.35,0.15)\)，\(A_V=0.30\)，\(v=100\)，\(\sigma=150\)（真值落在测试网格上）。走 `fit_spectrum`，不告诉真值 \(v,\sigma\)。

| 路径 | 网格（测试用，比默认密/小） | 断言 |
| --- | --- | --- |
| 无噪声 | \(v\in\{0,50,100,150,200\}\)，\(\sigma\in\{50,100,150,200\}\)，\(A_V\) 步长 0.1、范围 [0,0.6] | \(|\Delta v|\le 50\)，\(|\Delta\sigma|\le 50\)，\(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12 |
| SNR=30 | 同上 | \(|\Delta v|\le 100\)，\(|\Delta\sigma|\le 100\)，\(|\Delta A_V|\le 0.20\)，\(x\) atol 0.20 |
| 门禁 | 1×1 运动学点，确认不再 `NotImplementedError` | 当时 clip 未实现；D 之后 clip 已实现 |

**B2 结果（2026-09-01）：通过。**

## B3 验收（历史，2026-09-01）

| 路径 | 真值 | 测试网格 | 断言 |
| --- | --- | --- | --- |
| 离网 1 | \(x=(0.50,0.35,0.15)\)，\(A_V=0.30\)，\(v=80\)，\(\sigma=130\) | \(v\in[50,125]\) 步长 25；\(\sigma\in[100,175]\) 步长 25；\(A_V\in[0,0.6]\) 步长 0.1 | \(|\Delta v|\le 25\)，\(|\Delta\sigma|\le 25\)，\(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12 |
| 离网 2 | \(x=(0.20,0.55,0.25)\)，\(A_V=0.35\)，\(v=120\)，\(\sigma=90\) | \(v\in[50,175]\) 步长 25；\(\sigma\in[50,150]\) 步长 25；\(A_V\in[0,0.6]\) 步长 0.1 | 同上 |

**B3 结果（2026-09-01）：通过。** C 跳过；非线性加密改 G5。

## D 验收（历史，2026-09-01）

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| clip | 无噪声合成谱打一个 10σ 尖峰（不在归一窗口） | `1 ≤ n_clipped ≤ 5`；\(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12 |
| EX0 | 真值 \(x=(0.55,0.40,0.05)\)，`x_min_keep=0.10` | 第三分量 \(x=0\)；其余分数接近丢掉后的归一 |
| 门禁 | `metropolis_anneal` 仍 `NotImplementedError` | clip 不再抛未实现 |

**D 结果（2026-09-01）：通过。** v0.1 拟合器门槛达到。

## E 验收（历史，2026-09-01）

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| `.cxt` | 首行 Npix，四列，`flag=2` 一像素 | `load_spectrum` 读出；`good_from_flags` 该像素 False |
| gzip | `.cxt.gz` 与 `.spec.gz` 模板 | 与明文相同 |
| 批量 | 目录两份 `.cxt` | `iter_ascii_spectra` 产出 2 条 |
| FITS | 自造 SDSS 风格 table（需 astropy） | 波长/流量/误差/mask；两包各自 loader，互不 import |

**E 结果（2026-09-01）：通过。** v0.1 基线（拟合器 + 读文件）完成。

## 日志

| 日期 | 小段 | 做了什么 |
| --- | --- | --- |
| 2026-08-31 | A | 框架、合同、LOSVD 修正 |
| 2026-09-01 | B1 | `simulate.py` + 固定真值 \(v,\sigma\) 收回 \(x,A_V\) |
| 2026-09-01 | B2 | 外层 \(v,\sigma\) 网格；无噪声 + SNR=30 测试 |
| 2026-09-01 | B3 | 离网 \(v,\sigma\) 收回；C 跳过 |
| 2026-09-01 | D | NSIGMA clip + EX0；v0.1 拟合器 |
| 2026-09-01 | E | `.cxt` / gzip / 批量 / SDSS FITS；v0.1 基线 |
| 2026-09-02 | 合同 | PLAN 改为全面好的 1.x；断点 G1 |
| 2026-09-02 | G1 | FitResult 字段、NaN/空窗口、σ=0 警告 |
| 2026-09-02 | G2 | 已知 z 静止系、空气/真空 |
