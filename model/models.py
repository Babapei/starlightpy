# models.py

import numpy as np

def build_model(x_j, base_matrix, extinction_curve, A_V, AY_V=0.0, YAV_flags=None, q_lambda=None, q_lambda0=None):
    """
    构建模型光谱 Mλ

    参数:
        x_j:           (N,) 光度分量（每个成分的归一化权重）
        base_matrix:   (W, N) base 光谱矩阵，每列一个成分，每行为一个波长
        extinction_curve: (W,) qλ = Aλ / AV
        A_V:           全局消光
        AY_V:          额外 selective extinction（可选）
        YAV_flags:     (N,) 每个分量是否使用 AYV
        q_lambda0:     基准波长点的 qλ（默认取 qλ 中央值）

    返回:
        M_lambda:      (W,) 拟合模型光谱
    """
    n_components = base_matrix.shape[1]
    n_pixels = base_matrix.shape[0]

    if YAV_flags is None:
        YAV_flags = np.zeros(n_components, dtype=int)

    if q_lambda0 is None:
        q_lambda0 = np.median(extinction_curve)

    A_Vj = A_V + AY_V * (YAV_flags > 0).astype(float)
    delta_q = extinction_curve[:, np.newaxis] - q_lambda0
    extinction_factor = 10 ** (-0.4 * delta_q * A_Vj[np.newaxis, :])

    # 模拟光谱：按列组合 base * extinction * x_j
    M_lambda = np.sum(base_matrix * extinction_factor * x_j[np.newaxis, :], axis=1)
    return M_lambda