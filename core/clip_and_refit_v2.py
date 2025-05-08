import numpy as np

from core import fit_model

def clip_and_refit(wavelengths, flux_obs, error_obs, base_matrix, YAV_flags, config=None,
                   method='NSIGMA', threshold=3.0, min_clipped=1):
    """
    残差剪裁并重新拟合（支持 AYV）
    """

    # 第一次拟合
    result1 = fit_model(wavelengths, flux_obs, error_obs, base_matrix, YAV_flags, config)
    model1 = result1['M_lambda']
    residuals = flux_obs - model1

    # 初始掩码（True = 保留）
    mask = np.ones_like(wavelengths, dtype=bool)

    if method == 'NSIGMA':
        sig = error_obs
        mask = np.abs(residuals) <= threshold * sig
    elif method == 'RELRES':
        rel_res = residuals / error_obs
        std_rel = np.std(rel_res)
        mask = np.abs(rel_res) <= threshold * std_rel
    elif method == 'ABSRES':
        std_abs = np.std(residuals)
        mask = np.abs(residuals) <= threshold * std_abs
    elif method == 'NOCLIP':
        mask[:] = True
    else:
        raise ValueError(f"未知的剪裁方法: {method}")

    # 检查是否需要再次拟合
    if np.sum(~mask) >= min_clipped:
        print(f"[Clip&Refit] Clipped {np.sum(~mask)} pixels using {method} (threshold={threshold})")
        result2 = fit_model(wavelengths[mask], flux_obs[mask], error_obs[mask],
                            base_matrix[mask, :], YAV_flags, config)
        result2['clipped_mask'] = mask
        return result2
    else:
        print("[Clip&Refit] Clipping too small, skipping re-fit")
        result1['clipped_mask'] = np.ones_like(wavelengths, dtype=bool)
        return result1