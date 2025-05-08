# core.py

import sys
sys.path.append("/Users/xuwangweifan/VSC/starlight")  # 替换为你的项目根目录


import numpy as np
from scipy.optimize import minimize
from extinction import get_extinction_curve
from models import build_model

def fit_model(wavelengths, flux_obs, error_obs, base_matrix, YAV_flags=None, config=None):
    """
    拟合模型：根据 base 和观测谱，拟合 (x_j, A_V)，输出 M_lambda 和 χ²

    参数：
        wavelengths: (W,) 波长数组
        flux_obs:    (W,) 观测光谱 Oλ
        error_obs:   (W,) 对应误差 eλ
        base_matrix: (W, N) base 光谱
        YAV_flags:   (N,) base 中每个分量是否使用 AYV（暂不使用）
        config: dict，配置项：
            'R_V', 'ext_law', 'AV_bounds', 'norm_range' 等

    返回：
        result = {
            'x': x_j,
            'A_V': A_V,
            'M_lambda': Mλ,
            'chi2': χ²,
            'success': True/False,
            ...
        }
    """
    W, N = base_matrix.shape
    if config is None:
        config = {}

    R_V = config.get('R_V', 3.1)
    ext_law = config.get('ext_law', 'CCM')
    AV_bounds = config.get('AV_bounds', (0.0, 3.0))
    norm_range = config.get('norm_range', (4010, 4060))  # λ范围做归一化

    # 生成 extinction curve q_lambda
    q_lambda = get_extinction_curve(wavelengths, law=ext_law, R_V=R_V)
    λ_norm_mask = (wavelengths >= norm_range[0]) & (wavelengths <= norm_range[1])
    q_lambda0 = np.median(q_lambda[λ_norm_mask])  # 用于构建 Mλ

    # 归一化 Oλ 和 base，避免幅度问题
    norm_factor = np.median(flux_obs[λ_norm_mask])
    flux_obs_norm = flux_obs / norm_factor
    error_obs_norm = error_obs / norm_factor
    base_matrix_norm = base_matrix / np.median(base_matrix[λ_norm_mask])

    # 初始值
    x0 = np.ones(N) / N
    A0 = 0.5
    p0 = np.concatenate([x0, [A0]])

    # 约束：x_j >= 0，∑x_j ≈ 1，AV ∈ [min, max]
    bounds = [(0.0, 1.0)] * N + [AV_bounds]
    constraints = [{
        'type': 'eq',
        'fun': lambda p: np.sum(p[:N]) - 1.0
    }]

    def loss(p):
        x, AV = p[:N], p[-1]
        M_lambda = build_model(
            x_j=x,
            base_matrix=base_matrix_norm,
            extinction_curve=q_lambda,
            A_V=AV,
            YAV_flags=YAV_flags,
            q_lambda0=q_lambda0
        )
        chi2 = np.sum(((flux_obs_norm - M_lambda) / error_obs_norm) ** 2)
        return chi2

    opt = minimize(loss, p0, method='SLSQP', bounds=bounds, constraints=constraints)

    x_opt = opt.x[:N]
    AV_opt = opt.x[-1]
    M_lambda_opt = build_model(x_opt, base_matrix_norm, q_lambda, AV_opt, YAV_flags, q_lambda0)

    return {
        'x': x_opt,
        'A_V': AV_opt,
        'M_lambda': M_lambda_opt * norm_factor,  # 还原单位
        'chi2': opt.fun,
        'success': opt.success,
        'message': opt.message
    }


"""
使用方式
from core import fit_model
from io import load_spectrum, load_base_master, load_base_spectra

# 加载数据
wl, flux, err, flag, _, _ = load_spectrum("your_spectrum.txt")
base_list = load_base_master("Base.BC03.N")
wl_base, base_matrix = load_base_spectra("BasesDir", base_list)
YAV_flags = np.array([b['YAV'] for b in base_list])

# 拟合
result = fit_model(wl, flux, err, base_matrix, YAV_flags)

print("拟合成功：", result['success'])
print("AV =", result['A_V'])
print("光度分量：", result['x'])

"""