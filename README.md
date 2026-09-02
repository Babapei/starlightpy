# starlightpy + easyppxf

按 [docs/PLAN.md](docs/PLAN.md) 开发。以那份文档为准，不要靠聊天记录记目标。

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

`easyppxf.load_sdss_fits` 只读一维谱；拟合请引用 Cappellari，不要把本包装成一种新方法。SDSS FITS 需要 `pip install astropy`（或 `.[fits]` / `.[dev]`）。

当前：**阶段 H，H4 完成；断点 H5。** 见 [docs/PLAN.md](docs/PLAN.md) 与 [docs/PROGRESS.md](docs/PROGRESS.md)。

算法出处：Cid Fernandes et al. 2005（STARLIGHT）。`easyppxf` 用 pPXF 时请引用 Cappellari。本库是 MIT 许可的软件，不是那两篇论文的官方实现。

```bash
pip install -e ".[dev]"
pytest
```
