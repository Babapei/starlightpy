# 开发计划（按这个做，不要改目标）

本文是本仓库的开发合同。以后加功能前先对这里；和本文冲突的想法默认不做。

最后更新：2026-08-30

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
| `config.py` | `FitConfig`：归一窗口、\(A_V\) 网格、消光定律、运动学初值 | 有 | 退火温度表、clip 阈值 |
| `io.py` | ASCII 谱 / mask / base master（STARLIGHT 文件格式可读） | 有骨架 | gzip SSP、FITS |
| `extinction.py` | \(q_\lambda=A_\lambda/A_V\)：CCM（含 Fa/Fb）、CAL、Gordon | 有 CCM/CAL/GD | HyperZ 表 |
| `model.py` | 红化后的线性组合 | 有 | AYV |
| `kinematics.py` | 对数波长上的高斯 LOSVD（\(v_\star,\sigma_\star\)） | 有基础卷积 | 与仪器分辨率匹配 |
| `fit.py` | 优化 \(x_j, A_V\)（先网格+\(A_V\)+NNLS） | 有原型 | 见阶段 C |
| `clip.py` | 残差 clip 再拟合 | 空模块 + 文档 | 阶段 D |
| `optimize.py` | 退火 + Metropolis | 占位，禁止提前堆 | 阶段 C |

对外只保证：

```python
from starlightpy import fit_spectrum, FitResult, FitConfig, build_model
```

`easyppxf.fit_spectrum` 保持另一个函数，不要改名去「统一」。

---

## 4. 前向模型约定（写代码必须遵守）

1. **波长单位**：埃（Å）。观测与模板必须先插到**同一套波长**再拟合；插值放在 `io` 或 `fit` 预处理，不要在 `build_model` 里偷偷插。
2. **归一**：`norm_window` 默认 `(4010, 4060)`。每个模板一列单独除该窗口 median；观测同样。禁止用整块 `base_matrix` 的一个 median 去除。
3. **红化**：\(r_\lambda = 10^{-0.4(q_\lambda-q_{\lambda0})A_V}\)，\(q_{\lambda0}\) 取归一窗口内 \(q\) 的 median。先红化各列再 \(\sum x_j\)。
4. **运动学**：在**对数波长**（速度）上卷积高斯，中心 \(v_\star\) km/s，色散 \(\sigma_\star\) km/s。在合成 \(M_\lambda\)（红化之后）做一次，不要对每个模板用不同 \(\sigma\)（第一版）。
5. **\(\chi^2\)**：\(\sum_\lambda w_\lambda^2 (O_\lambda-M_\lambda)^2\)，\(w=1/e_\lambda\)，mask 处 \(w=0\)。
6. **\(x_j\)**：非负。第一版用 NNLS，不强制 \(\sum x_j=1\)（归一窗口已把尺度吃进 \(x\)）；报告前可再除以 \(\sum x\) 得到分数，原始权重也要保留。

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

- 固定 \(v_\star,\sigma_\star\) 可卷积（已有函数）之后：对 \(v,\sigma\) 粗网格或外层一维/二维搜索，内层仍 NNLS+\(A_V\)
- 验收：合成谱 \(v=100\) km/s、\(\sigma=150\) km/s，收回误差先定宽（例如 \(|\Delta v|<50\), \(|\Delta\sigma|<80\)），再收紧
- 未完成前不要退火

### 阶段 C — 优化器换成风格更近的搜索（仍非 1:1）

- 在 `optimize.py` 写 Metropolis + 降温表
- 参数：\(x\)（可在单纯形或非负空间扰动）、\(A_V\)、\(v\)、\(\sigma\)
- 多起点（「链」）只做 3～5 条短链，不要复制 Fortran 配置文件里几十个开关
- 验收：与阶段 B 同一合成谱，\(\chi^2\) 不差于网格法；高消光 \(A_V\sim1\) 不系统性崩

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
4. 断言写在 `tests/starlightpy/test_recovery.py`

禁止：用和拟合器同一套近似去「验收自己」却不经过 `fit_spectrum`；禁止只画图不 assert。

---

## 7. 配置项白名单

只允许往 `FitConfig` 加计划里出现过的字段。想加新开关：先改本文，再改代码。

当前字段：

- `norm_window: tuple[float, float]`
- `law: str`（CCM / CAL / GD1 / GD2 / GD3）
- `r_v: float`
- `a_v_bounds, a_v_step`
- `v0_kms, sigma_kms`（阶段 A 可当固定值；阶段 B 起可搜索）
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

任一为是：先停，改回本文。

---

## 9. 开始开发前（读完再写阶段 B）

1. **一次只做一件。** 现在只做阶段 B（\(v,\sigma\) 搜索）。不要同时改 `easyppxf`、不要写退火、不要加 AYV / 发射线 / IFU。
2. **两个 `fit_spectrum` 不要混。** `from starlightpy import fit_spectrum` 才是主业；pPXF 那条请 `import easyppxf` 并起别的名字。
3. **先合成谱，不要一上来拟合真星系。** 真谱有红移、真空/空气波长、仪器展宽、发射线，阶段 B 会把锅甩给拟合器。SSP 模板版权（如 BC03）也不要往仓库里塞。
4. **观测和模板必须已经在同一套波长上。** 现在 `build_model` / `fit_spectrum` 不负责插值。网格不同会直接报错或默默拟合错。
5. **45 个 SSP 收不回「真实 \(x_j\)」。** 成分高度简并。验收用 2～4 个差得开的模板；\(x\) 看大类（年轻/年老）对不对，不要要求向量逐元相等。
6. **阶段 B 网格会爆炸。** \(A_V \times v \times \sigma\) 每个点一次 NNLS。先粗网格、波长短一点的合成谱，跑通测试再加密。
7. **卷积顺序。** 计划写的是先混合再 LOSVD；线性卷积下与「每列先卷积再 NNLS」等价。改运动学时别再搞一套对每个模板不同的 \(\sigma\)。
8. **新开关先改 PLAN 再改进 `FitConfig`。** 禁止再出现 `fit_model_v2.py` 和本机 `sys.path`。
9. **每次改拟合：合成 → `fit_spectrum` → `assert`。** 只出图不算过。
10. **`easyppxf` 可以几个月不碰。** 它已经能跑；完善它排在阶段 F。

---

## 9. 当前下一步（阶段 A 收尾 → B）

1. 保持测试绿色。
2. 阶段 B：`fit_spectrum` 增加对 \(v_\star,\sigma_\star\) 的搜索（外层网格即可）。
3. 不要先写退火。
