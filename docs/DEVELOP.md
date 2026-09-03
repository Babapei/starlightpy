# 开发工作页（原根目录 README）

本文是给**仓库开发者**看的工作说明：阶段用语、内部约定、最小抄写片段。

- 目标与禁止项只认 [PLAN.md](PLAN.md)
- 做到哪只认 [PROGRESS.md](PROGRESS.md)
- 对外安装、完整例子、示例图：仓库根目录 [README.md](../README.md)

不要把新需求写在根 README 里却不改 PLAN。

---

# starlightpy + easyppxf

按 [PLAN.md](PLAN.md) 开发。以那份文档为准，不要靠聊天记录记目标。

两个包并排，不要合成一个 `backend=`：

| 包 | 用途 |
| --- | --- |
| `starlightpy` | **主业**：STARLIGHT 风格复现（模板 + 消光 + 运动学）。不是 `starlight.exe` 的 1:1。 |
| `easyppxf` | **附带**：把线性波长交给 pPXF。拟合请引用 Cappellari。 |

## 最小正确用法

观测与模板先准备好（静止系、同一波长、分辨率、发射线 mask）。`fit_spectrum` 不搜索宇宙学红移，也不拟合发射线。\(x_j\) 始终是 NNLS。

```python
from starlightpy import (
    FitConfig,
    apply_mask,
    fit_spectrum,
    light_to_mass,
    light_weighted_age,
    optical_emission_mask_regions,
)

config = FitConfig(
    redshift=0.02,            # 已知 z，只应用
    wave_frame="vacuum",      # 观测若是空气波长则用 "air"
    fwhm_data=2.5,            # 仪器 FWHM，单位 Å
    fwhm_template=1.0,
    search_kinematics=True,
    refine_kinematics=True,   # 粗网格后再加密 (A_V, v, σ)
    clip_nsigma=3.0,
)
good = apply_mask(wave, optical_emission_mask_regions())
result = fit_spectrum(wave, flux, error, bases, mask=good, config=config)
print(result.a_v, result.v0_kms, result.sigma_kms, result.x_fraction)

# 有模板 M/L、年龄才做后处理；没有就不要编
mass, mu = light_to_mass(result.x_fraction, mass_to_light)
age_L = light_weighted_age(result.x_fraction, ages)
```

没有误差谱时才设 `estimate_error=True`（会警告这不是真 χ²）。网格不同先 `resample_to`，不要假设 `build_model` 会插值。

## 从文件到拟合

更不易错：先把观测移到静止系，再与模板采到同一网格，然后 `redshift=0`。`load_sdss_fits` 只返回波长/流量/误差/像素 mask，**不返回宇宙学 \(z\)**（SDSS 的 \(z\) 在 SPECOBJ 等元数据里）。

```python
from starlightpy import (
    FitConfig,
    apply_mask,
    fit_spectrum,
    load_sdss_fits,
    load_spectrum,
    optical_emission_mask_regions,
    resample_to,
    to_rest_frame,
)

# .cxt：wave, flux, error, flags = load_spectrum("galaxy.cxt")
wave_obs, flux_obs, err_obs, good_obs = load_sdss_fits("spec.fits")
wave_rest, flux_rest, err_rest = to_rest_frame(wave_obs, flux_obs, z, err_obs)
flux = resample_to(wave_rest, flux_rest, wave)   # wave 与 bases 同一静止系真空网格
error = resample_to(wave_rest, err_rest, wave)
good = resample_to(wave_rest, good_obs.astype(float), wave) > 0.5
good &= apply_mask(wave, optical_emission_mask_regions())
result = fit_spectrum(wave, flux, error, bases, mask=good, config=FitConfig(redshift=0.0))
```

`FitConfig.redshift=` 只在这种约定下使用：`wave` / `bases` **已经是**静止系真空网格，`flux` 仍是观测系 \(F_\lambda\)、只是已经采样到这些波长数字上。不要把 SDSS 的 \(\lambda_\mathrm{obs}\) 既拿去建模板、又设 `redshift=z`——线心会错一截。

发射线用 mask，不要靠 `clip_nsigma`。clip 是拟合后的离群点；模型很差时会误删大量像素。光学 mask 含 Hβ/Hα 等，年老星族的吸收线也会被挡住，需要可改表。模板请自带（仓库不捆绑 SSP）。

## 可选（默认关）

这些开关默认都不开。打开前请看 [PLAN.md](PLAN.md)。

```python
config.regularize_x = "age_bins"       # 或 "smooth_age"；需要 template_ages=
config.error_method = "chi2_slice"     # 或 "repeat"；粗误差，不是协方差
config.pad_losvd = True                # 谱端 LOSVD 垫边；默认仍关
config.fit_ayv = True                  # 年轻模板额外 A_YV；需要 young_flags=
config.law = "dust:F99"                # 需 pip install "starlightpy[dust]"

from starlightpy import save_fit_result, load_fit_result
save_fit_result("fit.npz", result)     # 或 .json；不是 Fortran .out
loaded = load_fit_result("fit.npz")
# 读回后可画观测 vs 模型：loaded.wavelength, loaded.flux, loaded.model
```

`pad_losvd` 默认仍是 `False`（谱两端默认不可信）。消光可选 `.[dust]`，SDSS FITS 用 `.[fits]`。

## 真谱进拟合器之前

1. **静止系。** 已知 \(z\) 用 `redshift=` 或 `to_rest_frame`。不要让库去搜索宇宙学红移。
2. **同一套波长。** 观测与模板都在静止系真空 Å；不同网格用 `resample_to`。
3. **分辨率。** `fwhm_data` / `fwhm_template` 单位 Å。只把更锐的模板抹到数据，不能反向锐化。
4. **发射线 mask。** `optical_emission_mask_regions` + `apply_mask`。不要靠 clip 当发射线处理。
5. **简并。** 几十个 SSP 收不回「真实 \(x_j\)」。光加权年龄/Z 和 \(\mu_j\) 必须自带年龄、Z、\(M/L\)。

```python
from easyppxf import fit_spectrum as fit_ppxf

pp = fit_ppxf(wave, flux, templates, template_wave, error=err)
print(pp.velocity, pp.sigma)
```

`easyppxf.load_sdss_fits` 只读一维谱；拟合请引用 Cappellari，不要把本包装成一种新方法。SDSS FITS 需要 `pip install astropy`（或 `.[fits]` / `.[dev]`）。模板波长必须比星系谱更宽（pPXF 默认速度边界约 ±2900 km/s）；同一网格时包装器会丢掉两端不够的像素，太短则报错。

当前：**阶段 I 与阶段 J 完成；真实使用补丁见 [PLAN.md](PLAN.md)。版本 1.0.0。** 见 [PLAN.md](PLAN.md) 与 [PROGRESS.md](PROGRESS.md)。对外产品页是仓库根目录 [README.md](../README.md)。

算法出处：Cid Fernandes et al. 2005（STARLIGHT）。`easyppxf` 用 pPXF 时请引用 Cappellari。本库是 MIT 许可的软件，不是那两篇论文的官方实现。

```bash
pip install -e ".[dev]"
pytest
```
