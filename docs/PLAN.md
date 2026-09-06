# 开发计划（按这个做，不要改目标）

本文是本仓库的开发合同。**以本文为准，不以聊天记录为准。** 以后加功能前先对这里；和本文冲突的想法默认不做。聊天里说的若要生效，必须改成本文的一节。

最后更新：2026-09-06（再审查「永不做」；1.0 之后的优化写入以后可做。未自动开工那些项）

---

## 0. 一句话

主项目 `starlightpy` 是这条定位上**尽可能好的 Python 库**：用恒星模板的非负线性组合 + 尘埃 + 视向速度/弥散解释星系连续谱。种群和尘埃是第一公民，运动学要正确且够用。它不是 pPXF，不是 Bagpipes/Prospector，不是 `starlight.exe` 的 Python 发行版。

附带 `easyppxf`（线性波长交给 pPXF）。两个包并排，禁止合成一个 `backend=`。

这是**库开发**。没有投稿任务。「引用」只是用了别人算法时文档里该写的出处。

完善 = 工作流、防错、推断品质和科学产品，不是 Fortran 开关数量。v0.1（拟合器 + 读文件）是**基线，不是上限**。

---

## 1. 两个包，永远分开

| 包 | 干什么 | 谁做拟合 | 文档里写清出处 |
| --- | --- | --- | --- |
| `starlightpy` | STARLIGHT 风格连续谱拟合（主业） | 我们 | 算法来自 Cid Fernandes et al. 2005，不是官方 Fortran |
| `easyppxf` | 线性波长 → 调 pPXF | pPXF | Cappellari；不要把本包装写成一种新拟合方法 |

禁止：

- `fit_spectrum(..., backend="starlight"|"ppxf")`
- 在 `starlightpy` 里 `import easyppxf` 或反过来当拟合后端
- 给 Fortran 生成 grid / config / `.cxt` 当主路径（那是包装器，不是复现）
- 再引入 FIREFLY、Bagpipes、PST 等到同一个 `fit()`

目录（不要改这层语义）：

```text
docs/PLAN.md                 ← 本文
src/starlightpy/             ← 主业
src/easyppxf/                ← 附带
tests/starlightpy/
tests/easyppxf/
```

---

## 2. 三档：一直要 / 以后可做 / 永不做

聊天或产品讨论里的「不要做 / 可以优化」**不自动生效**。本节是 2026-09-06 对旧 §2 永不做、以及当时分析里那批「不要做」的再审查结果。审查过仍在「永不做」里的，不能用「优化」「以后再说」绕过去。

### 核心一直要（库的底线；G–J / D1 已完成，下列不降级）

- Python 里 `from starlightpy import fit_spectrum`
- 模型（风格对齐 Cid Fernandes et al. 2005，不是逐字节克隆）：

\[
M_\lambda = \left(\sum_j x_j\, b_{j,\lambda}\, r_\lambda(A_V)\right) \otimes G(v_\star,\sigma_\star)
\]

- \(x_j \ge 0\)，窗口光度分数；验收用**合成谱**，不对 Fortran 二进制
- v0.1 已完成：数组 API、消光、NNLS、运动学网格、clip/EX0、`.cxt` / gzip / FITS 读入
- 拟合**前**的显式管道：已知 \(z\) 静止系、真空/空气、LSF/FWHM 对齐、光学发射线 mask 表、缺误差时的 RMS 估计（并警告不是真 χ²）
- 完整 `FitResult`：波长、流量、模型、最终 good、归一尺度、config 副本、EX0 丢掉哪些；未搜运动学且 \(\sigma=0\) 时警告
- 粗网格后对 \((A_V,v,\sigma)\) 局部加密或连续优化（旧阶段 C 的正确版；**不动** \(x_j\)）
- 有 \(M/L\) 才做 `light_to_mass`；有年龄/Z 元数据才做光加权年龄/Z；没有就拒绝
- README：最小正确用法 + 真谱清单（静止系、同一波长、分辨率、mask、简并）

### 研究可用（阶段 H；已实现；默认关）

- 多模板时非负正则或按年龄箱；χ² 切片或重复拟合的粗误差
- LOSVD 垫边 / 更密速度采样；同一 \((v,\sigma)\) 对设计矩阵预计算
- 可选 AYV（模板要有年轻旗标 + 有/无 AYV 两套真值）
- 可选对接 `dust_extinction`；自己的 npz/json 存盘（**不是**官方 `.out`）

G（1.x）已完成。阶段 H（H1–H6）已完成；各项默认关。阶段 I/J 是工作流与存盘，不是新物理开关。H 已有的开关不要再发明一套并行 API。

### 永不做（审查后仍不是「以后可选」）

| 永不做 | 审查 | 原因 |
| --- | --- | --- |
| 对 \(x_j\) 做 Metropolis / 退火 / MCMC | **维持** | 固定 \((A_V,v,\sigma)\) 后是线性非负最小二乘；NNLS 才是正解。Fortran 退火是搜索器，不是物理。 |
| `fit_spectrum` **搜索**宇宙学红移（`search_redshift`、把 \(z\) 放进非线性网格） | **维持** | 与 \(v_\star\) 共线。已知 \(z\) 只允许预处理。独立粗估见「以后可做」，不进 `fit()`。 |
| 把发射线 EW / 高斯线心线宽与 \((A_V,v,\sigma)\) 一起网格；用 clip 当发射线处理；光电离星云模型 | **维持（改口径）** | 那是另一套非线性物理。默认仍 mask。可选气体模板不是这条禁令，见改口径。 |
| 与 `starlight.exe` 输出 1:1、官方 `.out` / grid 生成器当产品 | **维持** | 暗示克隆；周边已有 `starlight_toolkit`。 |
| `backend=`、把 pPXF/FIREFLY/Bagpipes 接进同一个 `fit()` | **维持** | 别的算法；历史上合成一个入口已经证伪。 |
| 发行包 / wheel 里捆绑 BC03、E-MILES 或其他 SSP 库文件 | **维持** | 版权与体积。运行时下载的示例脚本可以（D1 已如此）。 |
| 本包做成贝叶斯 SFH、测光+光谱联合 SED、IFU 三维建模产品 | **维持（拆开 IFU）** | 那是 Bagpipes / Prospector / Pipe3D 一类软件。对每个 spaxel 调现有 `fit()` 的薄示例见以后可做。 |
| 用真星系谱当单元测试真理 | **维持** | 无法分辨是求解器错还是数据没准备好。示例图可以，assert 不行。 |
| 把 pPXF 的加性/乘性多项式、对 \(x\) 的 pPXF 正则、LOSVD 高阶矩 \(h_3,h_4\)、每个年龄不同的 \(\sigma\) 搬进 `starlightpy` | **新写入** | 会变成另一种拟合器。多项式禁令原在 §4.8，上升到与退火同级。 |
| 把 `easyppxf` 做成第二套恒星种群产品，或加一套自己的拟合器 | **新写入** | 它只是 pPXF 薄封装。现有入口的样板错误可修（U2）；不是「一行代码都不能动」。 |
| 未改本文就加 `FitConfig` 开关、写第三份架构愿景长文、在根 README 发明需求 | **新写入（过程）** | 聊天不生效；产品页不是合同。 |
| 无合成谱端证据就改 `pad_losvd` 默认 | **新写入（默认政策）** | H3 开关已有。默认仍 `False`。以后若改默认，必须先改本节并写验收。 |

旧句「同一 `fit()` 里拟合发射线」「IFU 工作流」「不要再为 easyppxf 加功能」**不再当禁令原文**。见下一小节。审查记录见本节末。

### 改口径（旧「永不」过宽或过死）

1. **发射线。** 旧表写成「同一 `fit()` 里拟合发射线 → 永不」。过宽：会禁止以后把若干**已知静止波长的气体模板**当作额外非负线性列（仍 NNLS），也会禁止残差上的第二次拟合。真正有害的是：把线当成与尘埃、运动学绑在一起的非线性未知数，或靠 clip 当发射线。**默认行为不变：光学 mask，不拟合气体。** 气体模板若做，必须先写以后可做对应小段，默认关，合成验收。
2. **宇宙学 \(z\)。** `fit_spectrum` 永不搜 \(z\)。仓库可以以后有一个**独立**的交叉相关粗估函数，返回建议 \(z\)，由调用者交给 `to_rest_frame`。禁止自动灌进拟合网格。
3. **IFU。** 本包不做立方体建模、空间正则、二维运动学场。以后可做：文档或 `examples/` 里对每个 spaxel 调现有 `fit()`。
4. **`easyppxf`。** 阶段 F「不要再加功能」过死（U2 已经是正确的例外）。改为：不加新拟合方法、不做对等产品；可修现有入口的错误与文档。

### 1.0 之后可做（优化改进；不是自动开工的阶段 K）

一次只点名一项：先在 [docs/PROGRESS.md](PROGRESS.md) 写验收，再改代码。不把本表一次性做成一堆开关。真星系仍不当 pytest 真理。未点名的项不要「顺便做」。

**拟合品质与速度（优先；不改前向模型身份）**

| 编号 | 项 | 怎么做才算还是这个库 | 不要做成 |
| --- | --- | --- | --- |
| O1 | 波长相关 LSF | 预处理：每个像素的 FWHM（Å）或巡天 `wdisp` 一类数组，只把模板抹到数据，不能反向锐化。常数 `fwhm_data`/`fwhm_template` 仍可用。合成：λ 相关展宽对齐后能收回少模板 \(x,A_V\)；不对齐 χ² 更大。 | 把逐像素 FWHM 放进 \((A_V,v,\sigma)\) 网格；反向锐化 |
| O2 | 网格 / LOSVD 性能 | H3 的 `losvd_design` 已在垫边/加密路径。默认路径也可缓存同一 \((v,\sigma)\) 的设计矩阵、避免重复卷积。现有合成收回测试数值不能无故变差。 | 为了快改回线性 Å 上的固定 σ 平滑；悄悄改 χ² 定义 |
| O3 | 光学 mask 可扩展 | `optical_emission_mask_regions` 现为写死表；年老星族 Balmer 吸收会被挡住。允许调用者传入额外窗口或选用「少挡吸收」预设。默认表保持光学星云线。 | 默默删掉 Hβ/Hα 却不改文档；用 clip 代替 mask |

**使用、误差、发行**

| 编号 | 项 | 怎么做才算还是这个库 | 不要做成 |
| --- | --- | --- | --- |
| O4 | PyPI / DOI / 版本节奏 | 包名 `starlightpy`；安装说明与仓库一致。不改拟合物理。 | 借发行改模型或捆绑 SSP |
| O5 | 更多真实使用示例 | 与 `examples/fit_ngc3522_emiles.py` 同类：运行时下载、图注写明不是标准种群解。 | 把真巡天谱写进 pytest 真理 |
| O6 | 误差表述 | H2 的 χ² 切片 / 重复拟合已有。文档写清：这是粗误差，不是 \(x_j\) 的协方差，也不是官方 Fortran 误差公式。 | 对简并的 \(x\) 输出假协方差矩阵当科学误差 |

**可选扩展（必须默认关；先改 §7 白名单）**

| 编号 | 项 | 怎么做才算还是这个库 | 不要做成 |
| --- | --- | --- | --- |
| O7 | 可选气体模板 | 额外非负线性列（已知静止系波长的线或窄模板），与恒星一起 NNLS。第一版与恒星共用 \((v,\sigma)\)。默认仍 mask、此项关。合成：有线时打开能收回恒星 \(x,A_V\)；关掉则 χ² 更差或恒星 \(x\) 偏。 | `fit_emission` 把 EW/线宽与 \(A_V\) 一起网格；CLOUDY / 光电离 |
| O8 | 独立红移粗估 | 新函数（名称不要叫 `fit_spectrum`）。交叉相关或同类；返回建议 \(z\)。合成：已知移动能收回到够 `to_rest_frame` 的精度。 | `FitConfig.search_redshift`；在 `fit()` 里扫 \(z\) |
| O9 | loader 可选读出 SPECOBJ \(z\) | 只返回，不自动写进 `FitConfig.redshift`。文档写明仍要自己做静止系。 | 读 FITS 就偷偷搜 \(z\) 或改流量 |
| O10 | 评估 `pad_losvd` 默认 | 仅当谱端合成证据表明默认垫边明显更好、且旧测试仍绿，才改本文把默认改为 `True`。 | 无验收改默认 |

审查时**明确不当优化、也不进上表**的：对 \(x_j\) 退火或 MCMC；加性/乘性多项式；逐成分不同 \(\sigma\)；测光联合；贝叶斯 SFH；捆绑 SSP；Fortran `.out`；`backend=`。

### 审查记录（旧永不做 + 当时分析里的「不要做」）

| 来源 | 原表述 | 结论 |
| --- | --- | --- |
| 旧 §2 | 对 \(x_j\) 退火 | 维持永不做。阶段 C 跳过是对的。 |
| 旧 §2 | 库内搜索宇宙学红移 | `fit()` 内维持永不做。独立粗估降为 O8，不自动开工。 |
| 旧 §2 | 同一 `fit()` 里拟合发射线 | **过宽。** 禁令收窄为非线性气体 + clip 当发射线 + 星云模型。线性气体模板 → O7。当前库仍不拟合发射线。 |
| 旧 §2 | 官方 `.out` / 1:1 exe | 维持永不做。 |
| 旧 §2 | `backend=` | 维持永不做。 |
| 旧 §2 | 捆绑 BC03 / E-MILES | 维持：不进发行包。运行时下载示例已允许。 |
| 旧 §2 | 贝叶斯 SFH、测光联合、IFU 工作流 | 前两项维持永不做，且**不要当下一步**。IFU 产品维持永不做；逐 spaxel 调 `fit()` → 以后可做示例，未编号开工。 |
| 旧 §2 | 真星系当单元测试真理 | 维持永不做。 |
| 分析 | 不要把 pPXF 多项式搬进 `starlightpy` | 写入永不做（原只在 §4.8）。 |
| 分析 | 不要把 `easyppxf` 做成对等产品 | 写入永不做。F 的「不再加功能」改为可修现有入口。 |
| 分析 | 不要未写合同加开关；不要第三份架构文 | 写入永不做（过程）。 |
| 分析 | 不要改 `pad_losvd` 默认 | 写入默认政策；评估放 O10，不自动改。 |
| 分析 | 网格性能、波长相关 LSF、PyPI/DOI、更多真实使用、误差表述 | 写入以后可做 O1–O6。**本审查不实现它们。** |
| 分析 | 不要赌贝叶斯 SFH / 搜 \(z\) / 同 fit 发射线 / 捆绑 SSP 当创新 | 同意：那不是这个库的下一步。发射线仅允许 O7 那种线性模板，且须另开小段。 |

### 别人写过的：该用就用

- **有、直接当依赖：** NumPy / SciPy。消光阶段 H 可改用 `dust_extinction`。
- **有、但是另一件事：** pPXF → 只放在 `easyppxf`。
- **有、但是不拟合：** `starlight_toolkit` 只读 Fortran 输出。
- **没有：** 可维护的 STARLIGHT **拟合器** Python 上游。主业只能自己写求解器。

---

## 3. `starlightpy` 模块与职责

只在这些文件里长代码，不要再复制一份 `fit_model_v2`。标「阶段 G 才建」的文件，未到 G 对应小段不要创建。

| 文件 | 职责 | 现在 | 以后 |
| --- | --- | --- | --- |
| `config.py` | `FitConfig` 白名单 | 有 v0.1 字段 | G 预告字段见 §7；不要加退火温度表 |
| `io.py` | `.cxt` / mask / base / gzip / FITS / `iter_ascii_spectra` / `resample_to`；npz/json 存盘（含波长与误差） | 有 | 不要作业调度；不要 Fortran `.out` |
| `extinction.py` | \(q_\lambda=A_\lambda/A_V\)：CCM、CAL、Gordon；H5 可选 `dust:<Model>` | 有 | — |
| `model.py` | 红化后的线性组合；H4 可选按年轻旗标加 \(A_{YV}\) | 有 | — |
| `kinematics.py` | 均匀 lnλ 上的高斯 LOSVD | 有 | H3：垫边/密采样/设计矩阵预计算 |
| `fit.py` | NNLS；\(A_V\) 与 \(v,\sigma\) 网格；clip/EX0；`FitResult`（含波长与拟合误差） | 有 | — |
| `clip.py` | NSIGMA clip + EX0 | 有 | — |
| `optimize.py` | 占位 | 阶段 C 已跳过 | **禁止对 \(x_j\) 做 Metropolis**；G5 不要走这条对 \(x\) 的退火 |
| `refine.py` | \((A_V,v,\sigma)\) 局部加密 + Nelder-Mead；每点 NNLS | G5 有 | 禁止对 \(x_j\) 退火 |
| `regularize.py` | 多模板非负正则 / 年龄箱（仍 NNLS） | H1 | 禁止改成对 \(x_j\) 退火 |
| `errors.py` | χ² 切片 / 重复拟合粗误差 | H2 | 不是协方差、不是官方误差公式 |
| `preprocess.py` | 静止系、空气/真空、LSF、线 mask、估误差 | 有 | 不要把 \(z\)/FWHM 当拟合参数 |
| `products.py` | \(M/L\) 后处理；有元数据时的光加权量 | 有 | 无 \(M/L\) 禁止输出 \(\mu_j\) |
| `simulate.py` | 合成谱 | 有 | G 测试继续用 |

对外保证（preprocess / products 已导出）：

```python
from starlightpy import fit_spectrum, FitResult, FitConfig, build_model
```

`easyppxf.fit_spectrum` 保持另一个函数，不要改名去「统一」。

---

## 4. 前向模型约定（写代码必须遵守）

1. **静止系。** `fit_spectrum` 不搜索宇宙学红移。\(v_\star\) 只是剩余视向速度。已知 \(z\) 在拟合前应用（阶段 G2），不是拟合参数。独立粗估若做，见 §2 O8，不得进本函数网格。
2. **波长单位**：埃（Å）。观测与模板必须先到**同一套波长**再拟合；用 `resample_to`，不要在 `build_model` 里插值。真空/空气转换是预处理（G2），不是拟合参数。
3. **归一**：`norm_window` 必须落在数据覆盖范围内（默认 `(4010, 4060)`）。每个模板一列单独除该窗口 median；观测同样。禁止用整块矩阵一个 median。
4. **红化**：\(r_\lambda = 10^{-0.4(q_\lambda-q_{\lambda0})A_V}\)，\(q_{\lambda0}\) 为归一窗口内 \(q\) 的 median。先红化各列再混合。
5. **运动学**：在**均匀 lnλ（速度）网格**上平移并高斯展宽，再插回原波长。禁止在线性 Å 像素上用固定 σ 做 `gaussian_filter1d`。\(v_\star>0\) 表示退行，\(\lambda \to \lambda e^{v/c}\)。第一版所有成分共用同一 \((v,\sigma)\)。同一 \((v,\sigma)\) 下每列先 LOSVD 再 NNLS。
6. **\(\chi^2\)**：\(\sum_\lambda w_\lambda^2 (O_\lambda-M_\lambda)^2\)，\(w=1/e_\lambda\)，mask 处 \(w=0\)。发射线**默认**用 mask，不要靠 clip 当发射线处理。禁止把 EW/高斯线宽与 \((A_V,v,\sigma)\) 一起网格。可选气体模板（额外非负列）见 §2 O7，未写验收前不要实现。
7. **\(x_j\)**：非负，**始终用 NNLS**。报告保留原始权重和归一化分数。
8. **不要加 pPXF 那种加性或乘性多项式，也不要 LOSVD 高阶矩或每个年龄不同的 \(\sigma\)。** 这是永不做，不是以后可做。
9. **合成运动学测试必须用带吸收线的模板**，不要用纯 Planck/黑体连续谱。
10. **仪器分辨率与 \(z\)、发射线、真空/空气是拟合前管道**，不要把 FWHM 或 \(z\) 塞进非线性网格当未知数。合成测试里模板与假观测一致即可；真数据先做 G3。
11. **卷积边缘**（谱两端）：默认 `pad_losvd=False` 时两端仍不可信。H3 已提供垫边与 lnλ 加密，须显式打开。无 §2 O10 验收前**不要改默认值**。

---

## 5. 阶段（必须按顺序，未完成不要跳）

每阶段结束条件：对应测试绿、本文勾选、没有新的「顺便做的」大功能。小段在 [docs/PROGRESS.md](PROGRESS.md) 先写验收再编码。

### v0.1 基线（已完成）

拟合器 + 读文件。不是「库的终点」。

#### 阶段 A — 框架

- [x] 双包目录 + `FitConfig`
- [x] 消光、模型、LOSVD、NNLS+\(A_V\) 网格原型
- [x] 合成谱：无运动学时收回 \(x, A_V\)
- [x] `easyppxf` 现有测试仍通过
- [x] 文档只指向本文，不写第三产品

#### 阶段 B — 运动学进入拟合

- [x] **B1** 吸收线合成器；告诉真值 \(v,\sigma\) 时仍收回 \(x,A_V\)
- [x] **B2** `search_kinematics=True`：外层粗搜 \(v,\sigma\)，内层 \(A_V\)+NNLS
- [x] **B3** 真值不落网格节点时仍收回

#### 阶段 C — 对 \(x_j\) 退火（跳过，且永不做）

- [x] **跳过（2026-09-01）。** 对 \(x_j\) 退火不是物理要求。非线性参数的**局部加密**改放到 **G5**，仍禁止 Metropolis 扫 \(x_j\)。`optimize.py` 保持占位。

#### 阶段 D — clip 与稀疏 \(x\)

- [x] NSIGMA clip 后冻结 \((v,\sigma)\) 再拟合
- [x] EX0（`x_min_keep`）
- [x] 合成谱验收

#### 阶段 E — 输入

- [x] `.cxt` / mask / base / gzip
- [x] 两包各自 SDSS FITS loader（不硬耦合）
- [x] `iter_ascii_spectra`

#### 阶段 F — `easyppxf`（薄封装，不是对等产品）

- 不要捆绑 E-MILES
- [x] SDSS FITS 已在 E 各写一份；拟合引 Cappellari
- **不做第二套恒星种群产品、不加自己的拟合器。** 现有入口的样板错误可修（U2）。不要借机加 pPXF 多项式以外的「新方法」。

---

### 阶段 G — 可认真调用的 1.x（已完成）

把「进拟合器之前会算错」和「结果看不清」补上。G 已完成。H 默认关。一次一个小段。

- [x] **G1** `FitResult` + 防错：返回 `good`、`obs_scale`、`config`、`dropped`；NaN / 空 norm 窗口硬报错；`search_kinematics=False` 且 `sigma_kms=0` 时警告
- [x] **G2** 已知 \(z\) 静止系 + 真空/空气（`redshift=` 只应用，不搜索）
- [x] **G3** LSF/FWHM 对齐 + 光学发射线 mask 表
- [x] **G4** 缺误差谱时 RMS 估计 + 警告（这不是真 χ²）
- [x] **G5** 网格最佳附近对 \((A_V,v,\sigma)\) 局部加密或连续优化；\(x_j\) 仍 NNLS。验收：离网真值比纯粗网格更近，或 χ² 不差
- [x] **G6** `light_to_mass`；可选光加权年龄/Z（无元数据则报错）
- [x] **G7** README：最小正确用法 + 真谱清单（仍不准写成第三份架构文）

**1.x 完成线 = G1–G7。** H 不是 1.x 门槛。**G 已完成。**

---

### 阶段 H — 研究可用（G 全部完成之后）

默认关；每项先改本文对应小节和 PROGRESS 验收再写代码。一次一个小段。

- [x] **H1** 多模板非负正则或按年龄箱（仍 NNLS；需 `template_ages`）
- [x] **H2** χ² 切片或重复拟合的粗误差
- [x] **H3** LOSVD 垫边 / 更密速度采样；同一 \((v,\sigma)\) 预计算设计矩阵
- [x] **H4** 可选 AYV（年轻旗标 + 有/无 AYV 两套真值）
- [x] **H5** 可选 `dust_extinction`
- [x] **H6** 自己的 npz/json 存盘（不是 Fortran `.out`）

---

### 阶段 I — 可交付的 1.0（工作流，不是新物理）

H 全部完成之后。不新增 `FitConfig` 开关，不改 `pad_losvd` 默认，不把真星系当单元测试真理，不再长 `easyppxf`。一次一个小段；先改本文和 PROGRESS 验收再编码。

- [x] **I1** `FitResult` 带上拟合用的波长，以及进入 χ² 的误差（观测流量单位，与 `model` 一致）
- [x] **I2** `save_fit_result` / `load_fit_result` 写入并读回上述字段；旧 payload 仍能加载
- [x] **I3** 合成谱端到端：读入 → mask → `fit_spectrum` → 光加权产品 → 存盘 → 读回
- [x] **I4** pytest CI；版本号 `1.0.0`
- [x] **I5** README 一小段可选（默认关）开关与存盘；仍禁止第三份架构文

---

### 阶段 J — 存盘可画图（不是新物理）

I 完成之后只补这一处：读回的 `FitResult` 要能单独对齐画「观测 vs 模型」。不新增 `FitConfig` 开关，不改 `pad_losvd` 默认，不用真星系当测试真理，不再长 `easyppxf`。

- [x] **J1** `FitResult.flux`：拟合网格上的观测流量，与 `model` 同长度、同单位（已知 \(z\) / 空气转换之后，不是归一后的内部数组）
- [x] **J2** 存盘 payload 写入 `flux`（version 3）；缺该字段的旧文件加载为 `None`
- [x] **J3** README 存盘示例写出可用 `wavelength` / `flux` / `model` 画图

### 真实使用补丁（不是新物理阶段）

J 完成之后按 README 做真实调用，并用公开 SDSS FITS 测读入（**不断言 SFH；禁止把真巡天谱写进 pytest 真理**）。只修实际卡点，不新增长 `FitConfig`，不改 `pad_losvd` 默认，不再把 `easyppxf` 做成对等产品。

- [x] **U1** README 补「从文件到拟合」一小段：`.cxt` / `load_sdss_fits`；loader 不返回宇宙学 \(z\)；更不易错的路径是 `to_rest_frame` + `resample_to` 之后 `redshift=0`；`FitConfig.redshift=` 的前提是 `wave`/`bases` 已是静止系真空网格、`flux` 仍是观测系 \(F_\lambda\) 采样到该网格。FITS `good` 与发射线 mask 要 AND。`clip_nsigma` 不是发射线处理。仍禁止第三份架构文。
- [x] **U2** `easyppxf.fit_spectrum`：观测与模板同一波长时，不要把 pPXF 因默认速度边界（约 ±2900 km/s 需要模板更宽）抛出的 `AssertionError` 交给用户。裁剪 `goodpixels` 到模板能覆盖的范围；重叠太短则 `ValueError` 说明原因。这是现有入口的样板错误，不是新拟合功能。
- [x] **U3** 合成谱测试：G2 的 `redshift=` 路径与显式 `to_rest_frame` 路径都能收回少模板 \(x,A_V,v,\sigma\)。禁止用真巡天谱当 assert 真理。

### 产品文档（不是新物理）

根目录 README 改为给使用者看的产品页。原先那份开发工作 README 原样迁到 [docs/DEVELOP.md](DEVELOP.md)。不新增 `FitConfig`，不捆绑 SSP，不用真星系当 pytest 真理。

- [x] **D1** 根 README：安装、一段能直接跑的合成谱、真星系示例图（NGC 3522 + 用户自备 E-MILES，图注写明不是标准种群解）、FitResult / 预处理 / 可选开关 / 引用。开发合同仍只在本文与 PROGRESS。

### 合同再审查 R（不是新物理）

D1 之后只改本文与进度，不实现以后可做、不加 `FitConfig`、不改 `pad_losvd` 默认。

- [x] **R1** 再审查旧 §2 永不做与分析里的「不要做」：维持身份禁令；发射线 / IFU / `easyppxf` 冻结改口径；可优化项写入 §2 以后可做 O1–O10。本段不写拟合代码。

---

## 6. 每阶段怎么验收（防止「看起来能跑」）

合成流程（拟合器本身）：

1. 用 `build_model` + `apply_losvd` 造无噪声 \(M\)
2. 加高斯噪声（给定 SNR）
3. `fit_spectrum` 只看见 \(O,e,b\)（以及该阶段允许的预处理）
4. 断言写在 `tests/starlightpy/`（混合物：`test_recovery.py`；运动学：`test_kinematics_search.py`；clip：`test_clip.py`；G 段新文件 `test_starlight_g*.py` 或同等唯一名）

G2–G4 额外：**先做错预处理 → 参数应偏；做对 → 收回。**  
G6：测试必须自带假 \(M/L\)（或假年龄/Z）；无元数据必须 raise。  
G5：与 B3 同一类离网真值，加密后更近或 χ² 不差于粗网格。

H1：简并多模板时，未正则箱内分量不平分；年龄箱或平滑后箱内光和收回，且 \(x_j\ge 0\)。  
H4：年轻模板真值带 \(A_{YV}\) 时，打开 `fit_ayv` 能收回；关掉则 χ² 更大或 \(A_V\) 偏。  
I1：`FitResult.wavelength` 与输入网格一致；`error` 与 `model` 同单位、同长度。  
I3：公开读入 + 预处理 + fit + 产品 + 存盘走一遍，合成谱收回；不用真巡天谱当真理。  
J1：`FitResult.flux` 与传入（或对齐后的）观测流量同单位、与 `model` 同长。

U1：README 能对照「文件 → 静止系 → 同一网格 → mask → fit」，并写明 `redshift=` 与 `to_rest_frame` 的分工。  
U2：同一波长数组调用 `easyppxf.fit_spectrum` 不出现 pPXF `AssertionError`；模板过短则明确 `ValueError`。  
U3：两条静止系路径走 `fit_spectrum`，收回容差与 B2 同级（\(|\Delta v|\le 50\)，\(|\Delta\sigma|\le 50\)，\(|\Delta A_V|\le 0.10\)，\(x\) atol 0.12）。  
D1：根 README 含 `git+https` 安装与一段可执行合成谱；`docs/DEVELOP.md` 仍含开发工作页；示例图在 `docs/figures/`，图注不把真星系当真理。  
R1：§2 含审查后的永不做表、改口径、以后可做 O1–O10 与审查记录；代码与默认 `pad_losvd` 不变。

禁止：用和拟合器同一套近似去「验收自己」却不经过 `fit_spectrum`；禁止只画图不 assert；禁止用真巡天谱当单元测试绿灯。

---

## 7. 配置项白名单

只允许往 `FitConfig` 加本节出现过的字段。想加新开关：先改本文，再改代码。

**v0.1 已有（实现了）：**

- `norm_window: tuple[float, float]`
- `law: str`（CCM / CAL / GD1 / GD2 / GD3）
- `r_v: float`
- `a_v_bounds, a_v_step`
- `v0_kms, sigma_kms`（`search_kinematics=False` 时使用的固定值）
- `search_kinematics: bool`
- `v_bounds, v_step`（仅搜索时；单位 km/s）
- `sigma_bounds, sigma_step`（仅搜索时；σ 网格不宜含过大一段 0，默认下限 40 km/s）
- `clip_nsigma: float | None`
- `x_min_keep: float`
- `redshift: float`（G2；只应用已知 \(z\)）
- `wave_frame: str`（G2；`air` / `vacuum` / `as_is`）
- `fwhm_data`, `fwhm_template`（G3；单位 **Å** 仪器 FWHM；仅当数据更宽时展宽模板）
- `estimate_error: bool`（G4；缺 \(e_\lambda\) 时才用）
- `refine_kinematics: bool`（G5）
- `regularize_x: str | None`（H1；`None` / `smooth_age` / `age_bins`；默认 `None`）
- `regularize_strength: float`（H1；仅 `smooth_age`）
- `age_bin_edges: tuple[float, ...] | None`（H1；仅 `age_bins`，也可按年龄唯一值自动分箱）
- `error_method: str | None`（H2；`None` / `chi2_slice` / `repeat`）
- `n_repeat: int`（H2；仅 `repeat`）
- `repeat_seed: int`（H2）
- `pad_losvd: bool`（H3；默认关）
- `losvd_oversample: int`（H3；均匀 lnλ 加密倍数，默认 1）
- `fit_ayv: bool`（H4；默认关）
- `a_yv_bounds`, `a_yv_step`（H4）
- `law` 可取 `dust:<Model>`（H5；例如 `dust:F99`）。未装 `dust_extinction` 时明确报错。默认仍是自带 CCM / CAL / Gordon。

`fit_spectrum(..., template_ages=)` 不是 `FitConfig` 字段：年龄是模板元数据，长度必须等于成分数。未开正则时忽略。

`fit_spectrum(..., young_flags=)` 也不是 `FitConfig` 字段：年轻旗标是模板元数据，长度必须等于成分数。未开 `fit_ayv` 时忽略。

`save_fit_result` / `load_fit_result`（H6）也不是 `FitConfig` 字段：npz 或 json；默认不自动写盘；不是 Fortran `.out`。

`FitResult.wavelength` / `FitResult.error`（I1）也不是 `FitConfig` 字段：波长是拟合用的网格；误差是进入 χ² 的那份，单位与返回的 `model` 相同。旧存盘缺这两项时加载为 `None`。

`FitResult.flux`（J1）也不是 `FitConfig` 字段：拟合网格上的观测流量，与 `model` 同单位。旧存盘缺此项时加载为 `None`。

**阶段 J 预告：** 无（J1–J3 已实现）。

不要加：N_chains、Fortran 同名配置几十条、学习率、CNN 权重路径、`anneal_x`、`search_redshift`（进 `FitConfig` 或进 `fit()` 网格）、把发射线 EW/高斯参数当非线性网格的 `fit_emission`。§2 以后可做的项要进白名单时，先改本节对应字段再写代码；一次一项。

---

## 8. 走偏检查表（每次 PR 过一遍）

- [ ] 有没有新的第二个拟合入口和 `starlightpy.fit_spectrum` 抢名字？
- [ ] 有没有 `sys.path.append` 本机路径？
- [ ] 有没有 `fit_model_v2` 这种并列文件？
- [ ] 归一是不是按列、按窗口？
- [ ] 新功能有没有合成谱测试？
- [ ] 有没有把 pPXF 接到 `starlightpy` 里当「加速后端」？
- [ ] 有没有对 \(x_j\) 做 Metropolis/退火（应保持 NNLS）？
- [ ] LOSVD 是否在线性 Å 像素上用固定 σ 平滑（错误）？
- [ ] 有没有把 \(z\) 或 FWHM 当成拟合未知数（应为预处理）？
- [ ] 有没有无 \(M/L\) 却输出 \(\mu_j\)？
- [ ] 有没有用真星系当单元测试真理？
- [ ] 有没有把 pPXF 多项式 / \(h_3,h_4\) / 逐成分 \(\sigma\) 搬进 `starlightpy`？
- [ ] 有没有把发射线 EW/线宽与 \(A_V\) 一起网格，或用 clip 当发射线？
- [ ] 有没有无证据改 `pad_losvd` 默认？
- [ ] 有没有把 §2「以后可做」整表在一个 PR 里实现？

任一为是：先停，改回本文。

---

## 9. 写 1.x 时仍适用的提醒

1. **一次只做一件。** 阶段 J / D1 已完成。§2 以后可做须点名一项再写验收，不要开阶段 K 堆开关。不要为「更像 Fortran」加开关、不要无证据改 `pad_losvd` 默认。
2. **两个 `fit_spectrum` 不要混。** 主业是 `from starlightpy import fit_spectrum`；pPXF 请 `import easyppxf` 并起别的名字。
3. **先合成谱，不要一上来拟合真星系。** 真谱有红移、真空/空气、仪器展宽、发射线；没做 G2–G3 时锅会甩给拟合器。不要往仓库塞 BC03。
4. **观测和模板必须已经在同一套波长上**（或先走 G 的预处理再 `resample_to`）。`build_model` 不负责插值。
5. **45 个 SSP 收不回「真实 \(x_j\)」。** 成分高度简并。验收用 2～4 个差得开的模板；G6 的年龄/Z 看光加权量，不要要求向量逐元相等。
6. **网格会爆炸。** \(A_V \times v \times \sigma\) 每个点一次 NNLS。测试用缩小网格；G5 是局部加密，不是全空间细网格。
7. **卷积：** 同一 \((v,\sigma)\) 下每列先 LOSVD 再 NNLS。不要给每个年龄不同的 σ（那是另一产品）。
8. **新开关先改本文再改进 `FitConfig`。** 禁止 `*_v2.py` 和本机 `sys.path`。
9. **每次改拟合：合成 → `fit_spectrum`（或该段公开预处理 + fit）→ `assert`。** 只出图不算过。
10. **`easyppxf` 不要做成对等产品。** 可修现有入口；不加自己的拟合器。
11. **不要对 \(x_j\) 退火。** G5 只动 \((A_V,v,\sigma)\)。

---

## 10. 当前下一步

见 [docs/PROGRESS.md](PROGRESS.md)。产品文档 D1 已完成。2026-09-06 已再审查 §2：永不做仍有效（口径见审查记录）；可优化项写在 §2「以后可做」，**不要自动开工**。下一步若做代码，从 PROGRESS 点名的一项开始（建议先 O1 或 O2），先写验收再编码。

---

## 11. 仓库规范（文档怎么管、代码怎么长）

权威顺序：**PLAN > PROGRESS > docs/DEVELOP.md > 根 README > 代码注释 > 聊天。** 目标只认本文；做到哪只认 [docs/PROGRESS.md](PROGRESS.md)。根 README 是产品页，禁止在那里发明新需求或新开关。

不要再开「设计思想 / 架构愿景」第三份长文。阶段勾选、验收数字、禁止项都写在本文。开发抄写片段放 [docs/DEVELOP.md](DEVELOP.md)。

**完成一个阶段 / 小段：**

1. 对应测试绿（合成谱走公开 API + `assert`）。
2. 本文该清单打勾，改「最后更新」。
3. 产品 README（根目录）保持可安装、可抄一段就跑；G7/I5/J3/U1 的用户可见约定也写在产品页。开发阶段用语留在 DEVELOP / 本文。
4. 不新增 `*_v2.py`、不把 pPXF 接进 `starlightpy`。

**代码：**

- 主业只进 `src/starlightpy/`，附带只进 `src/easyppxf/`。
- 测试文件名全局唯一（`test_starlight_*.py` / `test_easyppxf_*.py`）。
- 观测与模板波长不同时，先预处理 / `resample_to`，不要在 `build_model` 里插值。
- 主包装不强制依赖 pPXF；`pip install -e ".[dev]"` 才跑 `easyppxf` 测试。

**出处（不是投稿）：** 算法来自 Cid Fernandes et al. 2005；`easyppxf` 注明 Cappellari；消光曲线各写自己的出处。LICENSE 为 MIT。

**两条完成线：**

- **v0.1（已完成）：** A + B + D + E。合成谱能收回少模板 \(x\)、\(A_V\)、\(v,\sigma\)，并支持 mask + clip + 读文件。C 已跳过。
- **1.x（G1–G7）：** 可认真调用：管道、防错、完整结果、局部加密、有元数据时的质量/年龄产品、README 清单。H 不是 1.x 门槛。
- **1.0.0（I1–I5）：** 结果带波长与拟合误差并可存盘读回；合成谱端到端；CI；README 可选开关。

---

## 12. 历史上会适得其反、已经改掉的部分

旧阶段 C 要把 \(x_j\) 改成退火，那是在模仿 Fortran 的搜索器，不是重写物理模型。会丢掉已经正确的线性求解。已改为：\(x_j\) 保持 NNLS。非线性加密改到 G5。

旧 LOSVD 在线性波长像素上做固定宽度平滑，σ 的单位不是 km/s。已改成均匀 lnλ 网格。

旧合同在 E 之后写「不要做 μ_j / FWHM / 真星系、库已完整」。那是作业最小集，不是上限。μ_j 与 FWHM 对齐改为 G 的后处理/预处理；真星系仍不当测试真理。

---

## 13. 定位（和现在不必拆掉的部分）

若今天空白开这个库，仍会保留：

- 物理拆成 红化 / 线性混合 / LOSVD（均匀 lnλ）/ \(\chi^2\)
- \(x_j\) 用 NNLS，\((A_V,v,\sigma)\) 用搜索再局部加密
- 合成谱 + `assert` 当验收，不对 Fortran 1:1
- 一份合同文档，禁止聊天当需求
- 拟合前管道与完整 `FitResult`，否则只有原型拟合器

从头**不会**做的（现在仓库里有，是历史，不是最优形状）：

- 把 pPXF 封装做成和主库对等的产品（留下当薄可选模块，不再长）
- 把 `.cxt` 当第一版核心（第一版仍是数组进、`FitResult` 出；文件是便利）
- 设「对 \(x_j\) 退火」当里程碑
- 并列两个 `fit_spectrum` 当统一 API（主入口只有 `starlightpy.fit_spectrum`）

三层：

- **G（1.x）**：能被认真调用。
- **H**：研究可用，默认关。
- **I**：1.0 工作流（结果自洽、存盘、CI、版本）。
- **J**：存盘可画图（观测流量进 `FitResult`）。
- **U（真实使用补丁）**：README 文件入口与 `redshift=` 约定；`easyppxf` 同波长不再断言崩溃。
- **D（产品文档）**：根 README 给使用者；`docs/DEVELOP.md` 给开发。
- **1.0 之后**：§2 以后可做（O1–O10）。须点名一项，不是新的总阶段。
- **永不做**：§2 第三档（2026-09-06 再审查后的表）。

一个仓库。`easyppxf` 太薄，不够单独开库；不统一成 `backend=`。
