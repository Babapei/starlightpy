# 进度 / 断点

中断后只看本文「当前断点」和下一小段。目标仍以 [PLAN.md](PLAN.md) 为准。

## 当前断点

| 项 | 值 |
| --- | --- |
| 阶段 | H（研究可用） |
| 小段 | **H1 已完成**；下一刀 **H2** |
| 不要做 | 对 \(x_j\) 退火、搜索宇宙学 \(z\)、`backend=`、用真星系当测试真理、再长 `easyppxf`、未写验收就开下一 H 项 |

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
- [x] G3 LSF + 发射线 mask 表
- [x] G4 缺误差时 RMS 估计
- [x] G5 非线性局部加密
- [x] G6 `light_to_mass` / 光加权产品
- [x] G7 README 用法清单

### 阶段 H

- [x] **H1** 多模板非负正则或按年龄箱
- [ ] H2 χ² 切片 / 粗误差
- [ ] H3 LOSVD 垫边与设计矩阵预计算
- [ ] H4 可选 AYV
- [ ] H5 可选 `dust_extinction`
- [ ] H6 npz/json 存盘

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

**G3 结果（2026-09-02）：通过。** 4Å vs 2Å 对齐收回 \(x,A_V\)；不对齐 χ² 更大；Hα 6563 被 mask。

## G4 验收（先写再做，2026-09-02）

`error=None` 且 `estimate_error=True` 时，用观测在 `good` 像素上的 RMS（相对连续谱或相对流量中位数）构造常数误差，并发出警告：这不是真 χ²。`estimate_error=False`（默认）且未给误差则 `ValueError`。有真实误差时不要覆盖。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 估计 | 无噪声合成后丢掉误差数组，`estimate_error=True` | 收回 \(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12；`UserWarning` 含「not a true」或「不是真」χ² |
| 拒绝 | 同一调用 `estimate_error=False` 且 `error=None` | `ValueError` |
| 真误差 | 给出误差数组 | 不因 `estimate_error` 改写该数组 |

未开始 G4 代码。

**G4 结果（2026-09-02）：通过。** 缺误差 + `estimate_error` 收回；默认拒绝；真误差不被覆盖。

## G5 验收（先写再做，2026-09-02）

粗网格之后对 \((A_V,v,\sigma)\) 局部加密（可再连续优化）。每个试探点 \(x_j\) 仍 NNLS。禁止走 `optimize.py` 对 \(x\) 的 Metropolis。`refine_kinematics=True` 才开；默认关，B3 粗网格行为不变。

与 B3 同一类离网真值：\(v=80\)，\(\sigma=130\)，\(A_V=0.30\)，\(x=(0.50,0.35,0.15)\)。测试粗网格 \(v\in[50,125]\) 步长 25，\(\sigma\in[100,175]\) 步长 25，\(A_V\in[0,0.6]\) 步长 0.1。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 加密 | `search_kinematics=True` 且 `refine_kinematics=True` | \(|\Delta v|\) 与 \(|\Delta\sigma|\) 均小于未加密；或加密 χ² ≤ 未加密 χ²。\(x\) atol 0.12，\(|\Delta A_V|\le 0.10\) |
| 对照 | 同一数据 `refine_kinematics=False` | 落在粗网格最近节点附近（与 B3 一致） |
| 门禁 | `metropolis_anneal` | 仍 `NotImplementedError` |

未开始 G5 代码。

**G5 结果（2026-09-02）：通过。** 离网 \(v=80,\sigma=130\) 加密后收到真值；粗网格仍落最近节点；Metropolis 仍关。

## G6 验收（先写再做，2026-09-02）

新文件 `products.py`。`light_to_mass(x, mass_to_light)` 把光分数变成质量权重/质量分数。光加权年龄或 Z 必须自带元数据数组。无 \(M/L\) 或无年龄/Z 则 raise。不要写进 `fit()`，不要无元数据输出 \(\mu_j\)。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 质量 | 假 \(x=(0.5,0.5)\)，\(M/L=(1,3)\) | 质量分数为 \((0.25,0.75)\) |
| 拒绝 M/L | `mass_to_light=None` 或省略 | `ValueError` |
| 光加权 | 假年龄 \((1,3)\) Gyr 与上列 \(x\) | 光加权年龄 \(=2\)；无年龄数组则 raise |

未开始 G6 代码。

**G6 结果（2026-09-02）：通过。** `light_to_mass` 分数正确；无 \(M/L\) / 年龄 / Z 均 raise。

## G7 验收（先写再做，2026-09-02）

README 增加最小正确用法（`fit_spectrum` + 预处理开关）和真谱清单：静止系、同一波长、分辨率、发射线 mask、成分简并。不要写成第三份架构文。勾选 PLAN G7；断点改为阶段 G 完成。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| README | 含用法示例与五条清单 | 仓库根 README 能对照 PLAN §2 核心一直要的最后一条 |
| 勾选 | G1–G7 全勾 | PLAN / PROGRESS 清单与断点一致 |
| 回归 | `pytest` | 全绿 |

**G7 结果（2026-09-02）：通过。** README 含最小正确用法与静止系 / 波长 / 分辨率 / 发射线 / 简并清单；G1–G7 全勾。

## H1 验收（先写再做，2026-09-02）

默认 `regularize_x=None`，行为与 G 相同。打开时 \(x_j\) 仍 NNLS（可加差分行或先按年龄箱塌缩），需要 `template_ages`。同龄差分用 \(\exp(-|\Delta\log_{10}t|/0.5\,\mathrm{dex})\) 加权，避免年轻/年老被强行拉平。

合成：6 个吸收线模板，三年轻（几乎一样）三年老；真值 \(x=(0.40,0,0,0.60,0,0)\)，\(A_V=0.30\)，无运动学。走 `fit_spectrum`。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 年龄箱 | `regularize_x="age_bins"`，`age_bin_edges=(0,1e8,1e11)` | 年轻箱 \(\sum x\) 与年老箱 \(|\Delta|\le 0.12\)；箱内分量极差 \(\le 0.05\)；\(x\ge 0\) |
| 对照 | 同一数据默认关正则 | 两箱光和仍 \(|\Delta|\le 0.12\)，但年轻箱内极差 \(> 0.08\) |
| 平滑 | `smooth_age` 且 `regularize_strength=1` | 年轻箱内极差小于对照；\(x\ge 0\) |
| 拒绝 | 开正则且 `template_ages=None` | `ValueError` |

未开始 H1 代码。

**H1 结果（2026-09-02）：通过。** 年龄箱均分光和；未正则箱内稀疏；smooth_age 缩小同龄极差；缺年龄 raise。

## H2 验收（先写再做，2026-09-02）

默认关。`error_method=None` 时 `FitResult` 不含误差产品。打开后用 χ² 切片（固定其它参数，扫 \(A_V\) 或 \(v,\sigma\)）给出粗 1σ 半宽：χ²(θ) ≤ χ²_min + 1。不声称协方差或贝叶斯区间。可选对合成噪声做少量重复拟合，报告 \(x\) 的标准差。\(x_j\) 仍 NNLS。

合成：无噪声 3 模板收回路径，真值 \(A_V=0.30\)。另：SNR=30 重复拟合。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 切片 | `error_method="chi2_slice"` | `FitResult.errors` 含 `a_v` 半宽；真值落在 \(\hat A_V\pm\) 半宽（无噪声半宽可很小但有限） |
| 对照 | `error_method=None` | `errors is None` |
| 重复 | SNR=30，`error_method="repeat"`，`n_repeat=5`，固定种子 | `errors["x"]` 长度 = n_comp，非负 |

未开始 H2 代码。

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
| 2026-09-02 | G3 | LSF/FWHM 对齐 + 光学发射线 mask |
| 2026-09-02 | G4 | 缺误差 RMS 估计 + 非真 χ² 警告 |
| 2026-09-02 | G5 | (A_V,v,σ) 局部加密，x_j 仍 NNLS |
| 2026-09-02 | G6 | light_to_mass 与光加权年龄/Z |
| 2026-09-02 | G7 | README 用法 + 真谱清单；阶段 G 完成 |
| 2026-09-02 | H1 | 年龄箱 / smooth_age 非负正则 |
