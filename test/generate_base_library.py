import numpy as np

def generate_base_library(wavelengths, n_base=5, seed=42):
    """
    生成一个简化的 base 矩阵：n_base 个成分，每个是一个不同形状的黑体-like 谱
    """
    np.random.seed(seed)
    base_matrix = []

    for i in range(n_base):
        T = np.random.uniform(3000, 10000)  # 类似恒星的温度
        spectrum = planck_like(wavelengths, T)
        spectrum /= np.median(spectrum)  # 归一化
        base_matrix.append(spectrum)

    base_matrix = np.array(base_matrix).T  # shape: [n_wavelengths, n_base]
    YAV_flags = np.zeros(n_base, dtype=int)  # 全部设为 0（无 AYV）

    return base_matrix, YAV_flags

def planck_like(wavelengths, T):
    """黑体近似谱，用来模拟 base 分量"""
    # λ(Å) 转换为 cm
    wl_cm = wavelengths * 1e-8
    # 简化的 Wien 近似 (忽略常数因子)：flux ∝ 1/λ^5 * exp(-c2 / (λT))
    c2 = 1.4388e-1  # cm·K
    flux = 1.0 / wl_cm**5 * np.exp(-c2 / (wl_cm * T))
    return flux