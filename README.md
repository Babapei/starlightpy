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

当前做到 **阶段 B1**（合成器）。断点见 [docs/PROGRESS.md](docs/PROGRESS.md)，合同见 [docs/PLAN.md](docs/PLAN.md)。

算法出处：Cid Fernandes et al. 2005（STARLIGHT）。`easyppxf` 用 pPXF 时请引用 Cappellari。本库是 MIT 许可的软件，不是那两篇论文的官方实现。

```bash
pip install -e ".[dev]"
pytest
```
