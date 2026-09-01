# 进度 / 断点

中断后只看本文「当前断点」和下一小段。目标仍以 [PLAN.md](PLAN.md) 为准。

## 当前断点

| 项 | 值 |
| --- | --- |
| 阶段 | B（运动学） |
| 小段 | **B1 已完成**（2026-09-01） |
| 下一小段 | B2：`search_kinematics` 外层搜 \(v,\sigma\) |
| 不要做 | 退火、真星系、改 `easyppxf`、对 \(x_j\) 做 Metropolis |

## 小段清单

### 阶段 A — 框架

- [x] 双包、前向模型、固定 \(v,\sigma=0\) 时收回 \(x,A_V\)
- [x] LOSVD 在均匀 lnλ；\(v>0\) 吸收往红端

### 阶段 B

- [x] **B1** `starlightpy.simulate`：吸收线模板 + `mock_observation`（与 `fit_spectrum` 同一套红化/混合/LOSVD）。无噪声、**告诉拟合器真实 \(v,\sigma\)** 时仍收回 \(x,A_V\)。
- [ ] **B2** 打开运动学网格搜索（先改 PLAN 白名单：`v_bounds` 等）
- [ ] **B3** 不告诉真值 \(v,\sigma\)，从合成谱收回（容差见 PLAN B3）

### 阶段 C–F

未开始。C 可选且不对 \(x_j\) 退火。

## 日志

| 日期 | 小段 | 做了什么 |
| --- | --- | --- |
| 2026-08-31 | A | 框架、合同、LOSVD 修正 |
| 2026-09-01 | B1 | `simulate.py` + `tests/starlightpy/test_simulate.py` |
