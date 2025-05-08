import numpy as np
import matplotlib.pyplot as plt


import sys
sys.path.append("/Users/xuwangweifan/VSC/starlight")  # 替换为你的项目根目录
from fit import fit_spectrum_pipeline
from test import generate_base_library, synthesize_observed_spectrum


# Step 1: 构造 λ 波长数组 & 伪基底
wavelengths = np.arange(3500, 7001, 1)  # 3500–7000 Å, step=1
base_matrix, YAV_flags = generate_base_library(wavelengths, n_base=5)  # 简单 base

# Step 2: 随机生成 xj 权重 + 模拟 Oλ（加 reddening、加噪）
true_x = np.array([0.4, 0.3, 0.2, 0.1, 0.0])
AV, AYV = 0.3, 0.0
v_shift, sigma_disp = 0.0, 150.0  # km/s
obs_flux, error_flux = synthesize_observed_spectrum(
    wavelengths, base_matrix, true_x,
    AV=AV, AYV=AYV, sigma=sigma_disp,
    SNR=20
)

# Step 3: 拟合
result = fit_spectrum_pipeline(
    wavelengths, obs_flux, error_flux,
    base_matrix, YAV_flags,
    config=None,  # 默认配置
    verbose=True
)

# Step 4: 画图
plt.figure(figsize=(10, 6))
plt.plot(wavelengths, obs_flux, label='Observed', lw=1)
plt.plot(wavelengths, result['M_lambda'], label='Model', lw=1.2)
plt.fill_between(wavelengths, obs_flux - error_flux, obs_flux + error_flux, color='gray', alpha=0.2, label='±1σ error')
plt.xlabel("Wavelength (Å)")
plt.ylabel("Flux")
plt.title("Spectral Fit")
plt.legend()
plt.tight_layout()
plt.show()

# Step 5: 输出 xj
plt.figure()
plt.bar(np.arange(len(result['x'])), result['x'])
plt.title("Population vector xj")
plt.xlabel("Component index")
plt.ylabel("Flux fraction")
plt.tight_layout()
plt.show()