import sys
sys.path.append("/Users/xuwangweifan/VSC/starlight")  # 替换为你的项目根目录
from core import burn_in_refit, clip_and_refit, ex0_condensed_refit, fit_model


"""
输入预处理好的光谱 + base，自动完成以下四步拟合：

	1.	First Fit（初拟合）
	2.	Clip & Refit（去掉离谱点后再拟合）
	3.	Burn-In（精细收敛）
	4.	EX0s（精简基底后再次拟合）

"""

def fit_spectrum_pipeline(wavelengths, flux_obs, error_obs,
                          base_matrix, YAV_flags,
                          config=None,
                          mask=None,
                          method='CUMUL',
                          ex0_threshold=0.02,
                          ex0_min_size=3,
                          verbose=True):
    """
    主拟合流程控制函数：First Fit -> Clip & Refit -> Burn-In -> EX0s
    """

    if verbose: print("[Step 1] First fit...")

    first_fit_result = fit_model(
        wavelengths, flux_obs, error_obs, base_matrix, YAV_flags,
        config=config
    )

    if verbose: print("[Step 2] Clipping residuals...")

    clip_result = clip_and_refit(
        wavelengths, flux_obs, error_obs,
        base_matrix, YAV_flags,
        init_result=first_fit_result,
        config=config
    )

    if verbose: print("[Step 3] Burn-in refining...")

    burn_result = burn_in_refit(
        wavelengths, flux_obs, error_obs,
        base_matrix, YAV_flags,
        clipped_mask=clip_result['clipped_mask'],
        init_result=clip_result,
        config=config
    )

    if verbose: print("[Step 4] EX0s condensed base fitting...")

    final_result = ex0_condensed_refit(
        wavelengths, flux_obs, error_obs,
        base_matrix, YAV_flags,
        clipped_mask=clip_result['clipped_mask'],
        previous_x=burn_result['x'],
        method=method,
        threshold=ex0_threshold,
        min_size=ex0_min_size,
        config=config
    )

    if verbose:
        print("[Done] Final fit completed.")
        print(f"  χ² = {final_result['chi2']:.3f}")
        print(f"  Condensed base size = {len(final_result['x'])}")

    return final_result