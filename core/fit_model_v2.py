# core.py（替换原来的 fit_model）
# 这个版本拓展支持对ayv参数的拟合（第二种消光分量 AY_V（用于年龄较小成分，YAV_flag=1 时额外叠加）
# 先前版本目前 fit_model() 只支持拟合：•成分比例 x_j，•全局消光 A_V

"""
	•	引入第三个拟合参数 AY_V
	•	根据 YAV_flag 对每个成分设置：
	•	AV_j = AV （若 flag=0）
	•	AV_j = AV + AYV（若 flag=1）
	•	对 AY_V 添加合理的 bounds，比如 [0, 2.0]
"""



import numpy as np
from extinction.extinction import get_extinction_curve
from models import build_model
from scipy.optimize import minimize


def fit_model(wavelengths, flux_obs, error_obs, base_matrix, YAV_flags=None, config=None):
    """
    更新版：拟合 x_j, AV, AYV
    """
    W, N = base_matrix.shape
    if config is None:
        config = {}

    R_V = config.get('R_V', 3.1)
    ext_law = config.get('ext_law', 'CCM')
    AV_bounds = config.get('AV_bounds', (0.0, 3.0))
    AYV_bounds = config.get('AYV_bounds', (0.0, 2.0))
    norm_range = config.get('norm_range', (4010, 4060))

    q_lambda = get_extinction_curve(wavelengths, law=ext_law, R_V=R_V)
    λ_norm_mask = (wavelengths >= norm_range[0]) & (wavelengths <= norm_range[1])
    q_lambda0 = np.median(q_lambda[λ_norm_mask])

    norm_factor = np.median(flux_obs[λ_norm_mask])
    flux_obs_norm = flux_obs / norm_factor
    error_obs_norm = error_obs / norm_factor
    base_matrix_norm = base_matrix / np.median(base_matrix[λ_norm_mask])

    x0 = np.ones(N) / N
    A0 = 0.5
    AY0 = 0.1
    p0 = np.concatenate([x0, [A0, AY0]])

    bounds = [(0.0, 1.0)] * N + [AV_bounds, AYV_bounds]
    constraints = [{
        'type': 'eq',
        'fun': lambda p: np.sum(p[:N]) - 1.0
    }]

    def loss(p):
        x, AV, AYV = p[:N], p[N], p[N+1]
        M_lambda = build_model(
            x_j=x,
            base_matrix=base_matrix_norm,
            extinction_curve=q_lambda,
            A_V=AV,
            AY_V=AYV,
            YAV_flags=YAV_flags,
            q_lambda0=q_lambda0
        )
        chi2 = np.sum(((flux_obs_norm - M_lambda) / error_obs_norm) ** 2)
        return chi2

    opt = minimize(loss, p0, method='SLSQP', bounds=bounds, constraints=constraints)

    x_opt, AV_opt, AYV_opt = opt.x[:N], opt.x[N], opt.x[N+1]
    M_lambda_opt = build_model(x_opt, base_matrix_norm, q_lambda, AV_opt, AYV_opt, YAV_flags, q_lambda0)

    return {
        'x': x_opt,
        'A_V': AV_opt,
        'AY_V': AYV_opt,
        'M_lambda': M_lambda_opt * norm_factor,
        'chi2': opt.fun,
        'success': opt.success,
        'message': opt.message
    }

"""
使用示例：
result = fit_model(...)
print("A_V =", result['A_V'])
print("AY_V =", result['AY_V'])
"""