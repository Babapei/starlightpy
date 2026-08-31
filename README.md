# starlightpy + easyppxf

按 [docs/PLAN.md](docs/PLAN.md) 开发。以那份文档为准，不要靠聊天记录记目标。

两个包并排，不要合成一个 `backend=`：

| 包 | 用途 |
| --- | --- |
| `starlightpy` | **主业**：STARLIGHT 风格复现（模板 + 消光 + 运动学）。不是 `starlight.exe` 的 1:1。 |
| `easyppxf` | **附带**：把线性波长交给 pPXF。拟合请引用 Cappellari。 |

```python
from starlightpy import FitConfig, fit_spectrum

result = fit_spectrum(wave, flux, error, bases, config=FitConfig())
print(result.a_v, result.x_fraction)
```

```python
from easyppxf import fit_spectrum as fit_ppxf

pp = fit_ppxf(wave, flux, templates, template_wave, error=err)
print(pp.velocity, pp.sigma)
```

当前优化器是 \(A_V\) 网格 + NNLS（计划阶段 A）。运动学搜索、退火、clip 见 PLAN 阶段 B–D，对应函数会 `NotImplementedError`。

```bash
pip install -e ".[dev]"
pytest
```
