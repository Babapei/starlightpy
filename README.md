# starlightpy

把 STARLIGHT 的**物理模型**做成可在 Python 里直接调用的库。这是一次审视后的重启：旧代码方向混杂且拟合算法对不上原版，已清空重写。

当前仓库是**原型**，不是 Fortran STARLIGHT 的数值克隆，也还不能替代正式科研流程里的 `starlight.exe`。

## STARLIGHT 是干什么的

[STARLIGHT](https://ui.adsabs.harvard.edu/abs/2011ascl.soft08006C/abstract)（Cid Fernandes et al. 2005, MNRAS, 358, 363）是星系光谱的全谱拟合 / 恒星种群合成程序，传统实现是 **Fortran 77 二进制**，不是开源 Python 库。

它做的事可以写成一句话：把一条观测光谱 \(O_\lambda\) 分解成一组事先准备好的简单恒星种群（SSP / base）的线性组合，同时拟合前景尘埃消光和视线恒星运动：

\[
M_\lambda = \sum_j x_j\, b_{j,\lambda}\, r_\lambda(A_V)\, \otimes G(v_\star,\sigma_\star)
\]

输出大家真正关心的是：光度权重向量 \(x_j\)（以及对应的质量权重）、平均年龄 / 金属丰度、星形成历史、\(A_V\)、速度弥散等。历史上它被大量用在 SDSS、CALIFA 一类巡天里，所以文献里经常要和 STARLIGHT 的结果对齐。

原版优化器**不是** `scipy.optimize.minimize`。论文里写的是 **simulated annealing + Metropolis** 的 \(\chi^2\) 搜索，外加 clip-and-refit、多链、丢掉 \(x_j \approx 0\) 的成分再精修（EX0）等一堆工程细节。源码长期不公开，官方发的是编译好的可执行文件；后来有 [v06 发行包](https://minerva.ufsc.br/~ariel/manual_starlight_v06r1.pdf) 和 `starlight_toolkit`，但 toolkit 只处理输入输出，不负责拟合。

## 现在还有没有必要做成 Python 库

分目标看，结论不一样。

**如果目标是「在 Python 里做全谱拟合，别再碰那个古老界面」——这个需求大体上已经被别人填了。**

| 工具 | 角色 | 和本仓库的关系 |
| --- | --- | --- |
| [pPXF](https://pypi.org/project/ppxf/) | Python 全谱拟合事实标准（运动学 + 种群，正则化线性反演） | 绝大多数「我想用 Python 拟合星系光谱」的人应该直接用它 |
| FIREFLY、pyPipe3D | 同样是 Python 全谱拟合 | 覆盖相近科学问题 |
| [starlight_toolkit](https://github.com/arielwrl/starlight_toolkit) | 官方相关的 Python 工具 | **读 STARLIGHT 输出、画图、合成测光**，不重写拟合器 |
| PST / FSPS / DSPS / Bagpipes | 正向种群合成或 SED 拟合 | 解决的是造谱 / 测光拟合，不是 STARLIGHT 那套像素级 SSP 分解 |

2024–2026 的对比论文里 STARLIGHT 仍在被当作对照码使用（例如 [arXiv:2401.12300](https://arxiv.org/abs/2401.12300)），论文里语言一栏仍是 Fortran。也就是说：**没有人把 STARLIGHT 拟合成一个广泛使用的 Python 库**；大家要么继续跑二进制，要么改用 pPXF。

**如果目标是「STARLIGHT 兼容的 Python 拟合器，结果能和旧文献对齐」——这件事仍然空着，但难、而且用户面窄。**

值得做，当且仅当你接受下面几条：

1. 产品是 *STARLIGHT-compatible fitter*，验收标准是：同一条谱、同一套 base，\(x_j\) / \(A_V\) / \(\chi^2\) 能对上 Fortran 输出（至少在合成谱上）。
2. 你明确**不**和 pPXF 抢「通用全谱拟合」这个市场。
3. 优化器要按原版逻辑做（退火 + Metropolis、clip、EX0、运动学卷积），而不是拿 SLSQP 硬凑一个看起来像的 \(\chi^2\)。

**不值得做**的版本是：再包一层 grid 文件生成器去喂 `starlight.exe`。那条路 `starlight_toolkit` 已经占了，而且你并没有摆脱那个二进制。

## 旧代码错在哪（不是「写得乱」那么简单）

仓库里其实同时走了**两条相反的路**：

1. `starlight.py` / `StarlightInputGenerator.py`：给 Fortran STARLIGHT **生成 .cxt / config / grid.in**。这是包装器，不是库。
2. `core/` + `model/` + `fit/`：试图在 Python 里自己拟合。这才是「重构成库」的方向。

两条路没有收成一个 API，模块互相 `sys.path.append("/Users/xuwangweifan/VSC/starlight")`，`fit_model` / `clip_and_refit` 各有 v1/v2，pipeline 调用的参数（`init_result=`）和实现签名对不上。演示脚本在当时的目录结构下就跑不起来。

物理和算法上的硬伤：

- **优化器选错了。** 原版是退火 + Metropolis，用来对付 \(x_j\) 高度简并、再乘上非线性 \(A_V\) 的 \(\chi^2\)。旧代码用 SLSQP 去优化几十维单纯形上的 \(x_j\) 加 \(A_V\)。这既对不上 STARLIGHT，在真实 45-SSP 问题上也很容易卡在差的局部解。
- **运动学被整段跳过。** 笔记里写了「后面一定要加上」，但没有 \(v_\star,\sigma_\star\) 的卷积，科学上还不能叫 STARLIGHT。
- **归一化不对。** STARLIGHT 在 \(\lambda_0\) 窗口里对**每个** base 和观测谱分别归一。旧 `fit_model` 用整块 `base_matrix` 子矩阵的一个 median 去除，成分之间的相对尺度被搅在一起。
- **CCM 紫外段缺 Fa/Fb 修正**（\(x \ge 5.9\,\mu\mathrm{m}^{-1}\)）。光学波段还凑合，近紫外会偏。
- `extinction.py` 里 `get_extinction_curve` 定义了两次；合成测试谱自己又写了一遍残缺的 CCM。
- clip / burn-in / EX0 只是「再调用一次 SLSQP」，没有温度表，没有从上次最佳模型出发的 Metropolis 精修。名字像 STARLIGHT，计算不是。

所以：当时的直觉（IO + 红化 + 线性组合 + clip 流程）有一部分是对的；**把 SLSQP 当成 STARLIGHT、以及同时写 Fortran 包装器，是错的。**

## 这次重启做了什么

只保留「Python 库」这一条路，物理模型对齐 STARLIGHT 的前半段（线性 mix + 前景尘埃），优化器换成对这个问题更诚实的做法：

- 对 \(A_V\) 做一维网格；
- 每个 \(A_V\) 上用 **非负最小二乘（NNLS）** 解 \(x_j\)；
- 每个模板在归一窗口内单独归一。

这仍然**不是** STARLIGHT 算法。它的用处是：API 和前向模型先立住，测试能在合成谱上收回 \(x_j\) 和 \(A_V\)，后面如果要做兼容克隆，再换成退火 + Metropolis，而不是在 SLSQP 上继续堆 clip。

```python
from starlightpy import fit_spectrum

result = fit_spectrum(wavelength, flux, error, base_matrix, mask=good_pixels)
print(result.a_v, result.x, result.chi2)
```

```text
src/starlightpy/     可安装的库
  extinction.py      CCM / Calzetti / Gordon q(λ)
  model.py           M_λ 前向模型
  io.py              .cxt / mask / base master
  fit.py             AV 网格 + NNLS
tests/               pytest
```

尚未实现（有意留空）：\(v_\star,\sigma_\star\) 卷积、AYV、clip-and-refit、EX0、与 Fortran 输出的回归测试、Rapid-χ²。

## 如果你只听一句

作者是计算机 / DL，不是天文学生。老师说的「方便天文学家」指的是**用户**，不是「再实现一遍 STARLIGHT」。

对这种身份，默认也不该去克隆 Fortran，也不该做各家拟合码的超级库。有用的只有两类：

1. **工程（方便）**：不要自己拟合。把天文学家已经信的东西（pPXF + Astropy）接到他们的数据上——SDSS/LAMOST FITS 进、批量出表、笔记本里三行调用。验收是「一个没装过 Fortran 的人十分钟跑通」，不是 \(\chi^2\) 对上 STARLIGHT。
2. **研究（DL）**：在拟合这件事上做你真正有优势的事——摊销推断 / 光谱 emulator（毫秒级出年龄和 SFH，而不是再跑 100 秒退火）。这有方法论文可写（天文方法栏或 ML workshop 都行）。不要把 CNN 包装成「Python 版 STARLIGHT」去跟 pPXF 抢信任。

克隆 STARLIGHT 给天文学家用，他们不会信一个外行重写的 \(\chi^2\)；整合所有码你也打不过现有引用。老师要的「方便」= 少碰古老软件，这条路 pPXF 已经铺了，你能加的是接口、数据和速度，不是再一个求解器。

## 建议你怎么选

1. **课程作业 / 作品集，做出「方便」**：本仓库改成 pPXF（或其它已有拟合器）的薄封装 + 巡天读入 + 批量，而不是继续自研 NNLS。
2. **想发 DL 方法**：另开方向做光谱的摊销推断 / emulator，用合成 SSP 训练，用 pPXF/STARLIGHT 当教师模型，不要在这个拟合器原型上堆网络。
3. **误以为要当天文学生**：用 pPXF 分析光谱即可，不必维护本仓库。

## 这不是我们的方向：把市面上所有拟合码整合成一个「完整标准库」

pPXF、FIREFLY、pyPipe3D、STARLIGHT、STECKMAP、Bagpipes、FSPS 看起来都在「拟合星系光」，但它们不是同一件事的不同实现，不能借来重写成一个大家都能用的标准库。

- **科学假设不同。** pPXF 是正则化线性反演（运动学可以和种群一起拟合）；STARLIGHT 是退火搜非负 \(x_j\)；FIREFLY 是 \(\chi^2\) 加权的 SFH 集成；Bagpipes / Prospector 是测光/SED 贝叶斯。合成一个 `fit()` 会把不可比的后验假装成同一套数。
- **「借用重写」会两边不讨好。** 重写 pPXF 等于重做 Cappellari 多年的工作，论文也不会引你的克隆；薄封装别人的包则变成永久适配层（API、许可证、引用、版本）。天文学里论文必须引用**方法原文**，一个超级包装器成不了新标准。
- **标准层已经有了，而且不是拟合器。** 共享的该是光谱数据结构（Astropy `Spectrum1D` / `specutils`）、消光曲线、SSP 读入。拟合算法必须保持可替换的后端，而不是揉成一种。
- **已经有人做过「统一界面」。** 例如 SPAN 是 GUI，里面 **调用 pPXF**，并不重写各家算法。那种产品是工作流外壳，不是新的科学库，也很快和各后端版本绑死。

所以本仓库的范围保持狭窄：**STARLIGHT 风格的前向模型 + 将来可选的兼容拟合器**。可以 *依赖* `dust_extinction` 或 `specutils` 这类基础设施，但不会去吞并 pPXF / FIREFLY / Bagpipes。若你真正需要的是「同一条谱用多种码跑一遍做对比」，那是一篇方法论文 + 脚本，不是一个 PyPI 标准库。

## 安装与测试

```bash
pip install -e ".[dev]"
pytest
```
