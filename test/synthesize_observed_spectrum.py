import numpy as np
import scipy.ndimage

def synthesize_observed_spectrum(wavelengths, base_matrix, x,
                                  AV=0.2, AYV=0.0,
                                  sigma=150, v_shift=0.0,
                                  SNR=20):
    """
    合成 Oλ 和 eλ：线性组合 base -> 加 reddening -> 加速度模糊 -> 加噪声
    """
    # 基础模型谱
    M0 = np.dot(base_matrix, x)

    # 加入 reddening
    q_lambda = compute_q_lambda(wavelengths, law='CCM')  # 默认 CCM
    A_lambda = AV * q_lambda
    M_red = M0 * 10**(-0.4 * A_lambda)

    # 加速度模糊（σ 为 km/s）
    delta_lambda = np.mean(np.diff(wavelengths))
    sigma_pix = sigma / 3e5 * (wavelengths.mean() / delta_lambda)  # 像素数
    M_conv = scipy.ndimage.gaussian_filter1d(M_red, sigma_pix)

    # 加随机噪声
    signal_level = np.median(M_conv)
    noise_std = signal_level / SNR
    noise = np.random.normal(0, noise_std, size=M_conv.shape)

    O_lambda = M_conv + noise
    e_lambda = np.full_like(O_lambda, noise_std)

    return O_lambda, e_lambda

def compute_q_lambda(wavelengths, law='CCM', Rv=3.1):
    """
    返回 q_lambda = A(λ)/A(V) 的数组（目前仅支持 CCM）
    """
    from scipy.interpolate import interp1d

    inv_lambda = 1e4 / wavelengths  # μm^-1

    # Cardelli et al. 1989 CCM law
    x = inv_lambda
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    mask1 = (x >= 0.3) & (x < 1.1)
    a[mask1] = 0.574 * x[mask1]**1.61
    b[mask1] = -0.527 * x[mask1]**1.61

    mask2 = (x >= 1.1) & (x <= 3.3)
    y = x[mask2] - 1.82
    a[mask2] = 1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4 + \
               0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
    b[mask2] = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4 - \
               0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7

    A_lambda_over_AV = a + b / Rv
    return A_lambda_over_AV