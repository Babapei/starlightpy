import numpy as np

from core import fit_model


"""
在 Burn-In 阶段收敛之后，我们可以“丢掉那些几乎没啥贡献的分量”，比如：
xj < 0.02 的 base 分量，然后再拟合一次。


STARLIGHT 原始设计逻辑：
	•	默认方式是：丢掉所有 xj < threshold 的成分（比如 2%）
	•	或者：丢掉最小的若干个 xj，但总和不超过 threshold（CUMUL 模式）
	•	最后不能低于 3 个分量（或自定义最小成分数）

"""

def ex0_condensed_refit(wavelengths, flux_obs, error_obs, base_matrix, YAV_flags,
                        clipped_mask, previous_x, method='CUMUL', threshold=0.02, min_size=3,
                        config=None):
    """
    精简 base 分量后再拟合（EX0s 模拟）
    参数        类型        默认              含义
    method      str         'CUMUL'         剪裁方式（CUMUL / SMALL）
    threshold   float       0.02            最小保留贡献（<2% 会被删）
    min_size    int          3              最少保留基底数量
    pop_vector  np.array    可选            是否使用自定义参考×
    """
    x = np.array(previous_x)

    if method == 'CUMUL':
        sorted_idx = np.argsort(x)
        cum_sum = np.cumsum(x[sorted_idx])
        keep_idx = sorted_idx[cum_sum > threshold]
    elif method == 'SMALL':
        keep_idx = np.where(x >= threshold)[0]
    else:
        raise ValueError("method 必须是 'CUMUL' 或 'SMALL'")

    if len(keep_idx) < min_size:
        keep_idx = np.argsort(x)[-min_size:]

    print(f"[EX0] Condensed base: 保留 {len(keep_idx)} 个分量")

    # 筛选保留的 base 和 YAV
    base_reduced = base_matrix[:, keep_idx]
    YAV_reduced = [YAV_flags[i] for i in keep_idx]

    # 使用 clip 后的 mask 重新拟合
    return fit_model(
        wavelengths[clipped_mask],
        flux_obs[clipped_mask],
        error_obs[clipped_mask],
        base_reduced[clipped_mask, :],
        YAV_reduced,
        config=config
    )



"""
使用：
# 3. EX0 Condensed Base 精简
final_result = ex0_condensed_refit(
    wavelengths,
    flux_obs,
    error_obs,
    base_matrix,
    YAV_flags,
    clipped_mask=fit_result['clipped_mask'],
    previous_x=burned_result['x'],
    method='CUMUL',  # 或 SMALL
    threshold=0.02
)
"""