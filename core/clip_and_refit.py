# core.py 中追加

# 实现 STARLIGHT 中的 NSIGMA裁剪机制
"""
目标功能：
	•	找出 |Oλ - Mλ| > N × eλ 的像素点
	•	将这些点的 wλ = 0，即不参与拟合（设置 error = ∞）
	•	然后再次执行一次 fit_model()，得到最终结果

"""

import numpy as np
from core import fit_model


def clip_and_refit(wavelengths, flux_obs, error_obs, base_matrix,
                   YAV_flags=None, config=None, clip_method='NSIGMA', sig_clip_thresh=3.0):
    """
    Clip & Refit 阶段：对偏差过大的像素进行裁剪（将 error 设置为 ∞），再重新拟合

    参数：
        clip_method: 'NSIGMA'|'RELRES'|'ABSRES'|'NOCLIP'
        sig_clip_thresh: threshold，通常 = 3.0
    返回：
        与 fit_model 相同，但带有裁剪信息
    """
    # 第一次拟合
    result_1 = fit_model(wavelengths, flux_obs, error_obs, base_matrix, YAV_flags, config)
    M_lambda = result_1['M_lambda']

    # 计算残差
    residual = flux_obs - M_lambda
    mask_clip = np.zeros_like(error_obs, dtype=bool)

    if clip_method == 'NSIGMA':
        mask_clip = np.abs(residual) > sig_clip_thresh * error_obs
    elif clip_method == 'RELRES':
        rel_res = np.abs(residual / error_obs)
        clip_threshold = sig_clip_thresh * np.std(rel_res)
        mask_clip = rel_res > clip_threshold
    elif clip_method == 'ABSRES':
        clip_threshold = sig_clip_thresh * np.std(residual)
        mask_clip = np.abs(residual) > clip_threshold
    elif clip_method == 'NOCLIP':
        return result_1  # 不做裁剪
    else:
        raise ValueError(f"Unknown clip method: {clip_method}")

    # 构造新的误差谱（裁剪点置为无穷大 → 权重为 0）
    error_obs_clipped = error_obs.copy()
    error_obs_clipped[mask_clip] = np.inf

    # 第二次拟合
    result_2 = fit_model(wavelengths, flux_obs, error_obs_clipped, base_matrix, YAV_flags, config)

    result_2['clip_mask'] = mask_clip
    result_2['n_clipped'] = np.sum(mask_clip)
    result_2['M_lambda_1st'] = M_lambda
    return result_2

"""
使用方式:
from core import clip_and_refit

result = clip_and_refit(
    wavelengths=wl,
    flux_obs=flux,
    error_obs=err,
    base_matrix=base_matrix,
    YAV_flags=YAV_flags,
    config={'ext_law': 'CCM'},
    clip_method='NSIGMA',
    sig_clip_thresh=3.0
)

print("拟合成功:", result['success'])
print("A_V =", result['A_V'])
print("裁剪点数量:", result['n_clipped'])
"""

"""
裁剪后的拟合效果可以对比可视化:
import matplotlib.pyplot as plt

plt.figure(figsize=(10,4))
plt.plot(wl, flux, label="Observed", color='black')
plt.plot(wl, result['M_lambda_1st'], label="Before Clipping", alpha=0.6)
plt.plot(wl, result['M_lambda'], label="After Clipping", alpha=0.8)
plt.scatter(wl[result['clip_mask']], flux[result['clip_mask']], color='red', label="Clipped", s=10)
plt.legend()
plt.xlabel("Wavelength (Å)")
plt.ylabel("Flux")
plt.title("STARLIGHT 模型拟合效果")
plt.tight_layout()
plt.show()
"""