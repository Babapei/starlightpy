# 进度 / 断点

中断后只看本文「当前断点」和下一小段。目标仍以 [PLAN.md](PLAN.md) 为准。

## 当前断点

| 项 | 值 |
| --- | --- |
| 阶段 | O4-PyPI 已完成（1.0.0 已上架） |
| 小段 | README 主安装句为 `pip install starlightpy`。不要自动开工 O1–O3、O5–O10；Zenodo 不做 |
| 不要做 | 对 \(x_j\) 退火、`fit()` 内搜宇宙学 \(z\)、`backend=`、用真星系当测试真理、发行包捆绑 SSP、把 `easyppxf` 做成对等产品、无证据改 `pad_losvd` 默认、把发射线当非线性网格、把 pPXF 多项式搬进 `starlightpy` |

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
- [x] H2 χ² 切片 / 粗误差
- [x] H3 LOSVD 垫边与设计矩阵预计算
- [x] H4 可选 AYV
- [x] H5 可选 `dust_extinction`
- [x] H6 npz/json 存盘

### 阶段 I（1.0 工作流）

- [x] **I1** `FitResult.wavelength` 与拟合误差
- [x] I2 存盘读回新字段
- [x] I3 合成谱端到端
- [x] I4 CI 与 `1.0.0`
- [x] I5 README 可选开关与存盘

### 阶段 J（存盘可画图）

- [x] **J1** `FitResult.flux`
- [x] J2 存盘读回 `flux`
- [x] J3 README 画图示例

### 真实使用补丁 U（不是新阶段）

- [x] **U1** README 从文件到拟合 / `redshift=` 约定
- [x] **U2** `easyppxf` 同波长不再 pPXF 断言崩溃
- [x] **U3** 合成谱：`redshift=` 与 `to_rest_frame` 两条路径收回

### 产品文档 D

- [x] **D1** 根 README 产品页；工作版 → `docs/DEVELOP.md`；真星系示例图

### 合同再审查 R

- [x] **R1** 再审查永不做；以后可做 O1–O10 写入 PLAN（本段无代码）

### O4 发行

- [x] **O4-PyPI** 1.0.0 上 PyPI；README 改为 `pip install starlightpy`

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

**H2 结果（2026-09-02）：通过。** χ² 切片半宽盖住真 \(A_V\)；默认 `errors is None`；重复拟合给出非负 \(x\) 标准差。

## H3 验收（先写再做，2026-09-02）

默认 `pad_losvd=False`、`losvd_oversample=1`，与 G 的 LOSVD 相同。打开垫边：在均匀 lnλ 上向两端延拓再卷积，压低谱端伪结构。`losvd_oversample≥2` 把 lnλ 网格加密。同一 \((v,\sigma)\) 的设计矩阵可预计算，数值与逐列 LOSVD 一致。

合成：吸收线靠近蓝端（~3850Å），真值 \(v=100\)，\(\sigma=150\)，告诉真值运动学。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 垫边 | `pad_losvd=True`，`losvd_oversample=2` | \(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12；χ² ≤ 不垫边 |
| 对照 | 同一谱默认 LOSVD | χ² 不小于垫边 |
| 预计算 | 同一 \((v,\sigma)\) 两次 `losvd_design` | `allclose` |

未开始 H3 代码。

**H3 结果（2026-09-02）：通过。** 蓝端吸收线垫边+加密收回；χ² 不差于默认 LOSVD；`losvd_design` 两次一致。

## H4 验收（先写再做，2026-09-02）

默认 `fit_ayv=False`，全体模板共用 \(A_V\)。打开后年轻旗标为真的模板用 \(A_V+A_{YV}\)，年老只用 \(A_V\)。需要 `young_flags`。两套真值：有 AYV / 无 AYV。

合成：3 模板，第一年轻；无运动学。有 AYV 真值 \(A_V=0.20\)，\(A_{YV}=0.40\)，\(x=(0.50,0.35,0.15)\)。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 有 AYV | `fit_ayv=True` + `young_flags` | \(|\Delta A_V|\le 0.10\)，\(|\Delta A_{YV}|\le 0.15\)，\(x\) atol 0.12 |
| 对照 | 同一有 AYV 的谱但 `fit_ayv=False` | χ² 更大，或 \(A_V\) 偏到年轻尘埃上 |
| 无 AYV | 真值 \(A_{YV}=0\)，`fit_ayv=True` | \(A_{YV}\) 落在 0 附近（≤ 一步网格），\(x,A_V\) 仍收回 |
| 拒绝 | `fit_ayv=True` 无旗标 | `ValueError` |

**H4 结果（2026-09-02）：通过。** 打开 AYV 收回 \(A_V,A_{YV},x\)；共用 \(A_V\) 更差；真值 \(A_{YV}=0\) 时估计落在 0 附近；无旗标 raise。

## H5 验收（先写再做，2026-09-02）

默认仍用自带 CCM / CAL / Gordon。`law="dust:<Model>"` 时才对接 `dust_extinction`（例如 `dust:F99`）。不改默认 `law="CCM"` 的数值。缺包时明确报错，不静默回退到自带曲线。

合成：3 模板，真值 \(x=(0.50,0.35,0.15)\)，\(A_V=0.30\)，无运动学。真消光用 F99。走 `fit_spectrum`。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| dust F99 | `law="dust:F99"` | \(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12 |
| 对照 | 同一 F99 谱但 `law="CCM"` | χ² 更大，或 \(|\Delta A_V|>0.05\) |
| 默认 | 自带 CCM 合成 + 默认 `law` | 仍收回；此路径不 `import dust_extinction` |
| 缺依赖 | 请求 `dust:F99` 但包不可用 | `ImportError` 或 `ValueError`，消息含 `dust_extinction` |
| 未知模型 | `law="dust:NOTALAW"` | `ValueError` |

**H5 结果（2026-09-02）：通过。** `dust:F99` 收回；错用自带 CCM 更差；默认 CCM 不变；缺包与未知模型明确报错。

## H6 验收（先写再做，2026-09-02）

自己的 npz / json 存盘，**不是** Fortran `.out`。默认 `fit_spectrum` 不写盘。保存并读回 `FitResult` 的关键字段：\(x\)、\(x\) 分数、\(A_V\)、\(A_{YV}\)、\(v,\sigma\)、χ²、model、good。config 有则一并保存。走公开 API。

合成：3 模板无噪声收回路径，真值 \(x=(0.50,0.35,0.15)\)，\(A_V=0.30\)。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| npz | `save_fit_result(path.npz)` 再 `load_fit_result` | \(x,A_V,v,\sigma,\chi^2\) 与内存一致；model/good 形状一致 |
| json | 同上，后缀 `.json` | 同上 |
| 默认 | 只 `fit_spectrum` | 不产生 `.out` 文件 |
| 拒绝 | 后缀既不是 npz 也不是 json | `ValueError` |

**H6 结果（2026-09-02）：通过。** npz 与 json 读回 \(x,A_V,v,\sigma,\chi^2\) 与 model/good；`fit_spectrum` 不写 `.out`；其它后缀 raise。

## I1 验收（先写再做，2026-09-03）

`FitResult` 增加 `wavelength`（拟合用的网格）和 `error`（进入 χ² 的误差）。`error` 与返回的 `model` 同长度、同流量单位，不是除过 `obs_scale` 的内部数组。字段只增不删。走 `fit_spectrum`。不新增 `FitConfig` 字段。

合成：3 模板，真值 \(x=(0.50,0.35,0.15)\)，\(A_V=0.30\)，无运动学。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 字段 | 无噪声收回 | `wavelength` 与输入 `wave` allclose；`error` 与 `model`、`good` 同长 |
| 单位 | 给出误差数组 | `error` 与传入误差 allclose（不是 `err/obs_scale`） |
| 估计 | `estimate_error=True` 且 `error=None` | `error` 有限、为正、同长；仍收回 \(A_V,x\) |
| 回归 | G1 字段 | `good` / `obs_scale` / `config` / `dropped` 仍在 |

**I1 结果（2026-09-03）：通过。** 波长与输入网格一致；误差与传入流量单位一致；估误差路径同长为正；G1 字段仍在。

## I2 验收（先写再做，2026-09-03）

存盘写入 I1 新字段。payload `version` 为 2。仍能读 H6 的 v1（缺字段则为 `None`）。不是 Fortran `.out`。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| npz / json | I1 路径的结果存了再读 | `wavelength`、`error`、`model` 形状一致且数值 allclose |
| 旧文件 | 无 `wavelength` 的 v1 payload | `load_fit_result` 不崩；`wavelength` 与 `error` 为 `None` |
| 拒绝 | 后缀 `.out` | `ValueError` |

**I2 结果（2026-09-03）：通过。** npz/json 读回波长与误差；v1 缺字段加载为 `None`；`.out` raise。

## I3 验收（先写再做，2026-09-03）

合成谱走公开工作流，不用真巡天谱。可把合成写成临时 `.cxt` 再 `load_spectrum`。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 端到端 | 写 `.cxt` → `load_spectrum` → `apply_mask` + 线表 → `fit_spectrum` → `light_to_mass` / `light_weighted_age`（自带假 \(M/L\)、年龄）→ `save_fit_result` → `load_fit_result` | \(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12；读回 `wavelength` 与 `model` 同长；无 `.out` |
| 拒绝 | `light_to_mass` 无 \(M/L\) | `ValueError` |

**I3 结果（2026-09-03）：通过。** `.cxt` 读入 + 线 mask + fit 收回；假 \(M/L\) 产品有限；存盘读回波长与 model 同长；无 `.out`。

## I4 验收（先写再做，2026-09-03）

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 版本 | `pyproject.toml` 与 `__version__` | 均为 `1.0.0` |
| CI | `.github/workflows/tests.yml` | `pip install -e ".[dev]"` 后 `pytest` |

**I4 结果（2026-09-03）：通过。** 版本 `1.0.0`；workflow 安装 `.[dev]` 并跑 `pytest`。

## I5 验收（先写再做，2026-09-03）

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 可选 | README 一小段「可选（默认关）」 | 出现 `regularize_x`、`error_method`、`pad_losvd`、`fit_ayv`、`dust:F99`、`save_fit_result`；提到 `.[dust]` / `.[fits]` |
| 当前做到 | 一句 | 阶段 I / 1.0.0 |
| 禁止 | — | 不新增长文、不把真星系当用法真理 |

**I5 结果（2026-09-03）：通过。** README 有可选开关与存盘；提到 `.[dust]` / `.[fits]`；当前做到阶段 I / 1.0.0。

## J1 验收（先写再做，2026-09-03）

`FitResult` 增加 `flux`：拟合用的观测流量，与 `model`、`wavelength` 同长，单位与 `model` 相同（未除 `obs_scale`）。字段只增不删。走 `fit_spectrum`。

合成：3 模板，真值 \(x=(0.50,0.35,0.15)\)，\(A_V=0.30\)，无运动学。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 字段 | 无噪声收回 | `flux` 与 `model`、`wavelength` 同长 |
| 单位 | 流量乘 10 再拟合 | `flux` 与传入观测 allclose；不是 `obs/obs_scale` |
| 回归 | I1 字段 | `wavelength`、`error` 仍在 |

**J1 结果（2026-09-03）：通过。** `flux` 与观测同单位、与 `model` 同长；`wavelength`/`error` 仍在。

## J2 验收（先写再做，2026-09-03）

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| npz / json | J1 结果存了再读 | `flux`、`model`、`wavelength` 同形状且 allclose |
| 旧文件 | v2 payload 无 `flux` | 加载不崩；`flux is None` |
| 拒绝 | `.out` | `ValueError` |

**J2 结果（2026-09-03）：通过。** npz/json 读回 `flux`；v2 缺字段为 `None`；`.out` raise。

## J3 验收（先写再做，2026-09-03）

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| README | 存盘示例 | 出现 `flux` 与 `model`，表明可画观测 vs 模型 |

**J3 结果（2026-09-03）：通过。** README 存盘示例写出 `loaded.wavelength` / `loaded.flux` / `loaded.model`。

## U1 验收（先写再做，2026-09-03）

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 文件入口 | README | 出现 `load_spectrum` 或 `load_sdss_fits`，以及 `to_rest_frame` |
| 约定 | README | 写明先静止系再 `redshift=0`，或 `redshift=` 时 wave/bases 已是静止系网格 |
| 禁止 | — | 不把真星系当用法真理、不新增长文 |

未开始 U1 正文。

**U1 结果（2026-09-03）：通过。** README 含 `load_sdss_fits` / `load_spectrum`、`to_rest_frame`，以及先静止系再 `redshift=0` 与 `redshift=` 网格约定。

## U2 验收（先写再做，2026-09-03）

同一套线性波长调用 `easyppxf.fit_spectrum`（pPXF 默认需要模板比星系更宽约 ±2900 km/s）。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 同波长 | 观测与模板同一 `wave`，`mask_emission=False` | 不抛 `AssertionError`；返回有限的 `velocity`/`sigma`/`bestfit` |
| 过短 | 模板只覆盖星系中段一小段 | `ValueError`，消息含 wavelength 或 2900 / km/s |
| 回归 | 现有「模板更宽」测试 | 仍绿 |

未开始 U2 代码。

**U2 结果（2026-09-03）：通过。** 同波长调用返回有限 \(v,\sigma\)；过短模板 `ValueError`；原「模板更宽」测试仍绿。

## U3 验收（先写再做，2026-09-03）

合成吸收线，真值 \(x=(0.50,0.35,0.15)\)，\(A_V=0.30\)，\(v=80\)，\(\sigma=130\)。测试网格与 G5/B3 同类缩小网格。走 `fit_spectrum`。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| G2 `redshift=` | 观测系 \(F_\lambda\) 采样到静止系网格，`redshift=z` | \(|\Delta v|\le 50\)，\(|\Delta\sigma|\le 50\)，\(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12 |
| `to_rest_frame` | 观测系波长/流量先 `to_rest_frame` + `resample_to`，`redshift=0` | 同上 |
| 禁止 | — | 不用真巡天谱 |

未开始 U3 代码。

**U3 结果（2026-09-03）：通过。** `redshift=` 与 `to_rest_frame` 两条路径均收回 \(x,A_V,v,\sigma\)。

## D1 验收（先写再做，2026-09-03）

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 工作版 | 原根 README 迁走 | `docs/DEVELOP.md` 含「按 PLAN 开发」与 `阶段 I` |
| 产品页 | 根 `README.md` | 含 `git+https` 安装、`default_absorption_bases`、`to_rest_frame`、`loaded.flux`；一段合成谱可直接执行 |
| 示例图 | `docs/figures/` | 有 NGC 3522 图；README 图注含「不是」标准解或「不是」真理 |
| 回归 | G7/I5/J3/U1 用户可见字符串 | 产品 README 仍含静止系清单与可选开关；`阶段 I` 改到 DEVELOP |

未开始 D1 正文。

**D1 结果（2026-09-03）：通过。** 工作版在 `docs/DEVELOP.md`；根 README 可 `git+https` 安装、合成谱片段可执行、含 NGC 3522 示例图（图注写明不是标准解）。

## R1 验收（合同再审查，2026-09-06）

只改文档。不实现 O1–O10，不加 `FitConfig`，不改 `pad_losvd` 默认。

| 路径 | 做法 | 断言 |
| --- | --- | --- |
| 永不做 | 再审查旧八条 + 分析里的不要做 | PLAN §2 表含维持/改口径/新写入；退火、`fit()` 搜 \(z\)、`backend=`、真星系 pytest 真理仍在永不做 |
| 改口径 | 发射线 / IFU / easyppxf | 旧「同一 fit 拟合发射线」不再当禁令原文；O7 写明线性气体模板；F 改为可修现有入口 |
| 以后可做 | 优化写入合同 | 有 O1–O10；写明不是自动开工的阶段 K；多项式/退火不进该表 |
| 代码 | 对比本段 diff | `src/` 无改动 |

**R1 结果（2026-09-06）：通过。** 合同已更新；拟合代码未改。

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
| 2026-09-02 | H2 | χ² 切片与重复拟合粗误差 |
| 2026-09-02 | H3 | LOSVD 垫边、密采样、设计矩阵预计算 |
| 2026-09-03 | U1–U3 | 真实使用：README 文件入口；easyppxf 同波长边界；两条静止系路径合成收回 |
| 2026-09-03 | D1 | 根 README 改为产品页；工作版迁到 docs/DEVELOP.md；NGC 3522 示例图 |
| 2026-09-06 | R1 | 再审查永不做（维持身份禁令；发射线/IFU/easyppxf 改口径）；优化项写入 PLAN §2 以后可做；不实现代码 |
| 2026-09-06 | O4-PyPI | starlightpy 1.0.0 上架 PyPI；README 主安装句改为 pip install |
