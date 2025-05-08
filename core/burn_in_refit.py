"""
	在 clip_and_refit() 之后，模型已经收敛得差不多了
    此函数
    •用 T=1（真实误差）进行一次收敛迭代
	•但不再重新初始化（不像 First Fits 那样到处乱跳）
	•尝试再微调一次 x、A_V、A_YV，拿到最终结果

    只保留“重启一次拟合，以 clip 后的光谱 + T=1 设置”的核心动作
"""

import numpy as np
from core import fit_model  #注意引用的是哪一个版本


def burn_in_refit(wavelengths, flux_obs, error_obs, base_matrix, YAV_flags,
                  clipped_mask=None, config=None):
    """
    Burn-In：对 clip 后的谱再做一次精细拟合
    """
    if clipped_mask is None:
        clipped_mask = np.ones_like(wavelengths, dtype=bool)

    print("[Burn-In] Performing final refinement with Temperature = 1")
    
    # 可在 config 里加权重策略等参数，暂时设为恒定
    return fit_model(
        wavelengths[clipped_mask],
        flux_obs[clipped_mask],
        error_obs[clipped_mask],
        base_matrix[clipped_mask, :],
        YAV_flags,
        config=config
    )

"""
用法：
# 第一步：初始拟合 + clip 重拟合
fit_result = clip_and_refit(wavelengths, flux_obs, error_obs, base_matrix, YAV_flags)

# 第二步：Burn-In 精修（T=1）
final_result = burn_in_refit(
    wavelengths,
    flux_obs,
    error_obs,
    base_matrix,
    YAV_flags,
    clipped_mask=fit_result['clipped_mask']
)
"""