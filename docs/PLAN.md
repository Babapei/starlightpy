# 开发计划（按这个做，不要改目标）

本文是本仓库的开发合同。**以本文为准，不以聊天记录为准。** 以后加功能前先对这里；和本文冲突的想法默认不做。聊天里说的若要生效，必须改成本文的一节。

最后更新：2026-09-01（B3 与 PROGRESS 对齐：真值不落网格点）

---

## 0. 一句话

主项目是用 **Python 重写 STARLIGHT 那套算法**（`starlightpy`）：模板线性组合 + 尘埃 + 运动学。附带 **pPXF 封装**（`easyppxf`）。两个包并排，禁止合成一个 `backend=`。

这是**库开发**。没有投稿任务。下面「引用」只是用了别人算法时库文档里该写的出处，不是要你写论文。

---

## 1. 两个包，永远分开

| 包 | 干什么 | 谁做拟合 | 文档里写清出处 |
| --- | --- | --- | --- |
| `starlightpy` | Python 重写 STARLIGHT 算法（主业） | 我们 | 算法来自 Cid Fernandes et al. 2005，不是官方 Fortran |
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

## 2. 锁定：做什么 / 不做什么

### 要做

- Python 里 `from starlightpy import fit_spectrum`
- 观测谱 \(O_\lambda\)、误差 \(e_\lambda\)、可选 mask、SSP 矩阵 \(b_{j,\lambda}\)
- 模型（风格对齐 Cid Fernandes et al. 2005，不是逐字节克隆）：

\[
M_\lambda = \left(\sum_j x_j\, b_{j,\lambda}\, r_\lambda(A_V)\right) \otimes G(v_\star,\sigma_\star)
\]

- \(x_j \ge 0\)；在 \(\lambda_0\) 窗口内各模板与观测分别归一，使 \(x_j\) 表示该窗口光度分数
- 验收用**合成谱**（已知 \(x, A_V, v, \sigma\)），不是对上 Fortran 输出

### 明确不做

| 不做 | 原因 |
| --- | --- |
| 与 `starlight.exe` 输出 1:1 | 验收别设成对上二进制；算法仍按 STARLIGHT 那套写 |
| 把 pPXF/FIREFLY 接进同一个 `fit()` | 那是别的算法，不是「完善 STARLIGHT」 |
| Fortran 输入生成器当产品 | 周边已有 `starlight_toolkit` |

### 别人写过的：该用就用，没有拟合器可接着改

- **有、直接当依赖：** NumPy / SciPy（NNLS、插值、卷积）。消光若哪天换成 `dust_extinction` 也可以，公式是公开的。
- **有、但是另一件事：** pPXF → 只放在 `easyppxf` 里调用，不要改成 `starlightpy` 的求解器（退火 / 非负 \(x_j\) 不是 pPXF 那套）。
- **有、但是不拟合：** `starlight_toolkit`、`starlight_tools` 只读 Fortran 输出。可以参考文件格式，不要 fork 成主项目。
- **没有：** 可维护的、把 STARLIGHT **拟合器**用 Python 写完的库。所以主业（`starlightpy`）没有「拿来继续完善」的上游，只能自己写求解器。

1:1 若将来要对着二进制做对照，单开一章，不挡主线。

---

## 3. `starlightpy` 模块与职责

只在这些文件里长代码，不要再复制一份 `fit_model_v2`。

| 文件 | 职责 | 现在 | 以后才做 |
| --- | --- | --- | --- |
| `config.py` | `FitConfig`：归一窗口、\(A_V\) 网格、消光定律、运动学初值 | 有（含 `clip_nsigma` / `x_min_keep` 字段，阶段 D 才启用） | 退火温度表 |
| `io.py` | ASCII 谱 / mask / base master；`resample_to` | 有 | gzip SSP、SDSS FITS、目录循环 |
| `extinction.py` | \(q_\lambda=A_\lambda/A_V\)：CCM（含 Fa/Fb）、CAL、Gordon | 有 CCM/CAL/GD | HyperZ 表 |
| `model.py` | 红化后的线性组合 | 有 | AYV |
| `kinematics.py` | 均匀 lnλ 上的高斯 LOSVD（\(v_\star,\sigma_\star\)） | 有（已按速度空间修正） | 与仪器分辨率匹配 |
| `fit.py` | 优化：\(x_j\) **永远 NNLS**；\(A_V\) 网格（阶段 A）；\(v,\sigma\) 外层搜索（阶段 B） | 有原型 | 阶段 C 只加密非线性参数 |
| `clip.py` | 残差 clip 再拟合 | 空模块 + 文档 | 阶段 D |
| `optimize.py` | 可选：对 \((A_V,v,\sigma)\) 退火 | 占位 | 阶段 C；**禁止对 \(x_j\) 做 Metropolis** |

对外只保证：

```python
from starlightpy import fit_spectrum, FitResult, FitConfig, build_model
```

`easyppxf.fit_spectrum` 保持另一个函数，不要改名去「统一」。

---

## 4. 前向模型约定（写代码必须遵守）

1. **静止系。** 第一版不拟合宇宙学红移。调用方先把谱移到静止系；\(v_\star\) 只是剩余视向速度。
2. **波长单位**：埃（Å）。观测与模板必须先到**同一套波长**再拟合；用 `resample_to`，不要在 `build_model` 里插值。
3. **归一**：`norm_window` 必须落在数据覆盖范围内（默认 `(4010, 4060)`）。每个模板一列单独除该窗口 median；观测同样。禁止用整块矩阵一个 median。
4. **红化**：\(r_\lambda = 10^{-0.4(q_\lambda-q_{\lambda0})A_V}\)，\(q_{\lambda0}\) 为归一窗口内 \(q\) 的 median。先红化各列再混合。
5. **运动学**：在**均匀 lnλ（速度）网格**上平移并高斯展宽，再插回原波长。禁止在线性 Å 像素上用固定 σ 做 `gaussian_filter1d`（那不是 km/s）。\(v_\star>0\) 表示退行，吸收往红端移，\(\lambda \to \lambda e^{v/c}\)。第一版所有成分共用同一 \((v,\sigma)\)。实现上应对**每列模板**做同一 LOSVD 再 NNLS（与「先混合再卷积」线性等价，且 \(x\) 仍是线性参数）。
6. **\(\chi^2\)**：\(\sum_\lambda w_\lambda^2 (O_\lambda-M_\lambda)^2\)，\(w=1/e_\lambda\)，mask 处 \(w=0\)。发射线用 mask，不要靠 clip 当发射线处理（clip 是阶段 D 的野点）。
7. **\(x_j\)**：非负，**始终用 NNLS（或等价非负最小二乘）求解**，不要改成单纯形上的随机游走。报告保留原始权重和归一化分数。
8. **不要加 pPXF 那种加性 Legendre 多项式。** 那会改掉 STARLIGHT 的连续谱模型（尘埃 + 种群）。
9. **合成运动学测试必须用带吸收线的模板**，不要用纯Planck/黑体连续谱（σ 无法约束）。
10. **仪器分辨率：** 合成测试里模板与假观测一致即可。真数据要对齐 FWHM，那是后续，不要塞进阶段 B。

---

## 5. 阶段（必须按顺序，未完成不要跳）

每阶段结束条件：对应测试绿、README/PLAN 勾选、没有新的「顺便做的」大功能。

### 阶段 A — 框架（本提交）

- [x] 双包目录 + `FitConfig`
- [x] 消光、模型、LOSVD、NNLS+\(A_V\) 网格原型
- [x] 合成谱：无运动学时收回 \(x, A_V\)
- [x] `easyppxf` 现有测试仍通过
- [x] 文档只指向本文，不写第三产品

### 阶段 B — 运动学进入拟合

拆成小段，每段在 [docs/PROGRESS.md](PROGRESS.md) 打勾后才能做下一段。

- [x] **B1** 合成器：吸收线模板 + 与拟合器同一套前向模型（红化 → 混合 → LOSVD）造 \(O_\lambda\)。固定真值 \(v,\sigma\) 时，现有 NNLS+\(A_V\) 仍能收回 \(x,A_V\)。
- [x] **B2** `search_kinematics=True`：外层粗搜 \(v,\sigma\)，内层仍 \(A_V\)+NNLS。字段：`v_bounds` / `v_step` / `sigma_bounds` / `sigma_step`。
- [ ] **B3** 真值**不落在网格节点上**时仍收回（例如 \(v=80\)、\(\sigma=130\)），确认是 χ² 最小而不是踩点。无噪声容差不宽于 B2；再加一组不同 \(x,A_V\)。B2 已满足第一轮 \(|\Delta v|\le 50\)、\(|\Delta\sigma|\le 50\)。

未完成 B 前不要退火，不要拟合真星系。

### 阶段 C — 只改进非线性参数的搜索（可选）

Fortran STARLIGHT 对 \(x_j\) 也走 Metropolis，是 2005 年的实现选择，不是物理要求。\(x_j\) 在固定 \((A_V,v,\sigma)\) 下是线性非负最小二乘。**对 \(x_j\) 做退火会更慢、更不稳，阶段 B 的 χ² 还可能变差。**

阶段 C 若做：只对 \((A_V, v_\star, \sigma_\star)\) 做更细的网格或短退火；\(x_j\) 仍 NNLS。3～5 个随机起点即可。不要把 Fortran 配置项搬进来。

验收：同一合成谱，χ² **不差于** 阶段 B；\(A_V\sim 1\) 不系统性崩。若网格已经够，本阶段可以跳过并在本文注明。

### 阶段 D — clip 与稀疏 \(x\)

- `clip.py`：第一次拟合后 \(|O-M|>n\sigma\) 置权重 0，再拟合一次
- 丢掉过小 \(x_j\) 再拟合（EX0 风格，阈值写进 `FitConfig`）
- 仍用合成谱；不要为了像手册而加没测试的分支

### 阶段 E — 输入体验（次要）

- 读 STARLIGHT 的 `.cxt` / mask / base master（`io.py` 扩）
- 可选：SDSS 一维 FITS 读入（给 `starlightpy` 和 `easyppxf` 各写自己的 loader，可复制不可硬耦合）
- 批量：目录里多条谱的循环即可，不要作业调度系统

### 阶段 F — `easyppxf` 补强（永远低于 A–D）

- 文档示例用用户自己的模板
- 不要捆绑下载 E-MILES 进本仓库
- 需要时再加 SDSS FITS 读入；加完必须在文档写「引 Cappellari」

---

## 6. 每阶段怎么验收（防止「看起来能跑」）

合成流程固定：

1. 用 `build_model` + `apply_losvd` 造无噪声 \(M\)
2. 加高斯噪声（给定 SNR）
3. `fit_spectrum` 只看见 \(O,e,b\)
4. 断言写在 `tests/starlightpy/`（混合物/尘埃：`test_recovery.py`；运动学搜索：`test_kinematics_search.py`；clip：`test_clip.py`）

禁止：用和拟合器同一套近似去「验收自己」却不经过 `fit_spectrum`；禁止只画图不 assert。

---

## 7. 配置项白名单

只允许往 `FitConfig` 加计划里出现过的字段。想加新开关：先改本文，再改代码。

当前字段：

- `norm_window: tuple[float, float]`
- `law: str`（CCM / CAL / GD1 / GD2 / GD3）
- `r_v: float`
- `a_v_bounds, a_v_step`
- `v0_kms, sigma_kms`（`search_kinematics=False` 时使用的固定值）
- `search_kinematics: bool`
- `v_bounds, v_step`（仅 `search_kinematics=True`；单位 km/s）
- `sigma_bounds, sigma_step`（仅搜索时；σ 网格不宜含过大一段 0，默认下限 40 km/s）
- `clip_nsigma: float | None`（阶段 D 才用）
- `x_min_keep: float`（阶段 D）

不要加：N_chains、Fortran 同名配置几十条、学习率、CNN 权重路径。

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

任一为是：先停，改回本文。

---

## 9. 开始开发前（读完再写阶段 B）

1. **一次只做一件。** 现在只做阶段 B（\(v,\sigma\) 搜索）。不要同时改 `easyppxf`、不要写退火、不要加 AYV / 发射线 / IFU。
2. **两个 `fit_spectrum` 不要混。** `from starlightpy import fit_spectrum` 才是主业；pPXF 那条请 `import easyppxf` 并起别的名字。
3. **先合成谱，不要一上来拟合真星系。** 真谱有红移、真空/空气波长、仪器展宽、发射线，阶段 B 会把锅甩给拟合器。SSP 模板版权（如 BC03）也不要往仓库里塞。
4. **观测和模板必须已经在同一套波长上。** 现在 `build_model` / `fit_spectrum` 不负责插值。网格不同会直接报错或默默拟合错。
5. **45 个 SSP 收不回「真实 \(x_j\)」。** 成分高度简并。验收用 2～4 个差得开的模板；\(x\) 看大类（年轻/年老）对不对，不要要求向量逐元相等。
6. **阶段 B 网格会爆炸。** \(A_V \times v \times \sigma\) 每个点一次 NNLS。先粗网格、波长短一点的合成谱，跑通测试再加密。
7. **卷积：** 同一 \((v,\sigma)\) 下，每列先 LOSVD 再 NNLS。不要给每个年龄不同的 σ。
8. **新开关先改 PLAN 再改进 `FitConfig`。** 禁止再出现 `fit_model_v2.py` 和本机 `sys.path`。
9. **每次改拟合：合成 → `fit_spectrum` → `assert`。** 只出图不算过。
10. **`easyppxf` 可以几个月不碰。** 它已经能跑；完善它排在阶段 F。
11. **不要为了「更像 Fortran」而对 \(x_j\) 做退火。**

---

## 10. 当前下一步

见 [docs/PROGRESS.md](PROGRESS.md)。B2 已完成；当前小段是 **B3**。

---

## 11. 仓库规范（文档怎么管、代码怎么长）

权威顺序：**PLAN > PROGRESS > README > 代码注释 > 聊天。** 做到哪只认 [docs/PROGRESS.md](PROGRESS.md)；目标仍只认本文。

不要再开「设计思想 / 架构愿景」第三份长文。阶段勾选、验收数字、禁止项都写在本文。

**完成一个阶段：**

1. 对应测试绿（合成谱走 `fit_spectrum` + `assert`）。
2. 本文该阶段清单打勾，改「最后更新」。
3. README 只改「当前做到哪」一句，不另写故事。
4. 不新增 `*_v2.py`、不把 pPXF 接进 `starlightpy`。

**代码：**

- 主业只进 `src/starlightpy/`，附带只进 `src/easyppxf/`。
- 测试文件名全局唯一（`test_starlight_*.py` / `test_easyppxf_*.py`），避免 pytest 撞模块名。
- 观测与模板波长不同时，先 `resample_to`，不要在 `build_model` 里插值。
- 主包装不强制依赖 pPXF；`pip install -e ".[dev]"` 才跑 `easyppxf` 测试。

**出处（不是投稿）：** 文档注明算法来自 Cid Fernandes et al. 2005；`easyppxf` 注明 Cappellari。LICENSE 为 MIT。

**何时算第一版做完（免得永远加功能）：** 阶段 A + B + D，合成谱能收回 \(x\)（少模板）、\(A_V\)、\(v,\sigma\)，并支持 mask + clip。质量权重 \(μ_j\)、AYV、仪器 FWHM、真巡天 FITS 都不是 v0.1 门槛。阶段 C 可跳过。

---

## 12. 本次通读改了什么（原计划会适得其反的部分）

旧阶段 C 要把 \(x_j\) 改成退火，那是在模仿 Fortran 的搜索器，不是重写物理模型。按那样做，会丢掉已经正确的线性求解，开工后拟合变差。已改为：\(x_j\) 保持 NNLS。

旧 LOSVD 在线性波长像素上做固定宽度平滑，σ 的单位不是 km/s，阶段 B 的速度弥散会假。已改成均匀 lnλ 网格。

## 13. 若从头设计（和现在不必拆掉的部分）

若今天空白开一个「Python 重写 STARLIGHT 算法」的库，仍会保留：

- 物理拆成 红化 / 线性混合 / LOSVD（均匀 lnλ）/ \(\chi^2\)
- \(x_j\) 用 NNLS，\((A_V,v,\sigma)\) 用搜索
- 合成谱 + `assert` 当验收，不对 Fortran 1:1
- 一份合同文档，禁止聊天当需求

从头**不会**做的（现在仓库里有，是历史，不是最优形状）：

- **不会把 pPXF 封装做成和主库对等的产品。** 若要放，只作为本仓库里的小可选模块（太薄，不够单独开 GitHub），不要统一 API。
- **不会把 `.cxt` / grid / base master 当第一版核心。** 第一版 API 就是数组进、`FitResult` 出。STARLIGHT 文件格式是以后的便利，不是算法。
- **不会设「退火阶段」当必做里程碑。** Fortran 用 Metropolis 搜 \(x_j\) 不必复刻。网格不够再加，且只动非线性参数。
- **不会并列两个 `fit_spectrum`。** 主入口只有 `starlightpy.fit_spectrum`。

现在的约定（与最早「两个包并排」一致，不是分两个 GitHub）：

- **一个仓库。** `easyppxf` 太薄，不够单独开库；作为本仓库里的**小可选模块**留下，不统一成 `backend=`，也不再为它加功能。
- **`starlightpy` 要做成比较完整的库：** 数组 API、合成谱验收、运动学、mask/clip、以及后来的文件读入（`.cxt` / mask / base）。完善指这条链路能用，不是 1:1 Fortran，也不是必做退火。
