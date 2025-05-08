# io.py（继续添加）

import numpy as np


def load_mask(filename: str):
    """
    读取 STARLIGHT 的 mask 文件，返回一个遮蔽区域列表。

    文件格式：
    第1行：一个整数 Nmasks（遮蔽区域数量）
    第2行起：每行 3 个数字：λ_ini, λ_fin, weight，可附注释（忽略）

    返回：
        mask_regions: List[Tuple[λ_ini, λ_fin, weight]]
    """
    mask_regions = []

    with open(filename, 'r') as f:
        lines = f.readlines()

        if not lines:
            raise ValueError("Mask file is empty.")

        try:
            n_masks = int(lines[0].strip())
        except ValueError:
            raise ValueError("First line of mask file must be an integer Nmasks.")

        if n_masks == 0:
            return []  # 空 mask 文件

        for line in lines[1:n_masks+1]:
            if line.strip() == '':
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue  # 忽略格式不完整行

            λ_ini = float(parts[0])
            λ_fin = float(parts[1])
            weight = float(parts[2])
            mask_regions.append((λ_ini, λ_fin, weight))

    return mask_regions


def apply_mask(wavelengths: np.ndarray, mask_regions: list):
    """
    根据 mask 区域列表，在给定波长数组中标记出遮蔽位置。

    返回：
        mask_flags: 一个与 wavelengths 等长的数组（0: 正常，≥2: 遮蔽）
    """
    mask_flags = np.zeros_like(wavelengths, dtype=int)

    for λ_ini, λ_fin, _ in mask_regions:
        mask_zone = (wavelengths >= λ_ini) & (wavelengths <= λ_fin)
        mask_flags[mask_zone] = 99  # 与 STARLIGHT 中的 flagλ ≥ 2 一致

    return mask_flags


"""
用于测试
from io import load_mask, apply_mask
import numpy as np

# 假设 wl 是你的波长数组
wl = np.arange(3500, 7500, 1)

mask_regions = load_mask("Masks.EmLines.SDSS.gm")
mask_flags = apply_mask(wl, mask_regions)

print(mask_regions[:3])
print(mask_flags[:20])  # 打印前20个点的 mask 状态
"""


"""
可扩展想法（后续可以加）：
权重处理               未来可加入"不同mask区域设置不同权值"的机制
可视化                标记mask区域在图中,辅助调试或者交互界面
自动合并flag和mask     融合flag和mask,生成最终遮蔽位图数组
"""