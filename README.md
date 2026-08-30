# easyppxf

把一条**线性波长**的星系光谱交给 [pPXF](https://pypi.org/project/ppxf/)，省掉官方示例里那段 log 重采样、速度尺度、模板对齐。

这不是 STARLIGHT，也不自己拟合。拟合是 pPXF 做的，论文请引用 Cappellari，不要引这个包装。

```python
from easyppxf import fit_spectrum, load_spectrum

wave, flux, err = load_spectrum("galaxy.txt")
result = fit_spectrum(wave, flux, templates, template_wave, error=err)
print(result.velocity, result.sigma)  # km/s
```

模板要你自己提供（和官方示例一样）。本包装不附带 E-MILES 库。

```bash
pip install -e ".[dev]"
pytest
```

本仓库以前试过重写 STARLIGHT，那条路停了。现在只做这一件：老师说的「方便天文学家在 Python 里调用」，用已经存在的拟合器，而不是再造一个。
