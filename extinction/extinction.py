# extinction.py

# 红化定律函数（消光函数）
# 后续可以加 CAL, GD1（已补充）, HZ1... 等分支

import numpy as np
# import scipy.interpolate as si
from scipy.interpolate import interp1d


def extinction_ccm89(wavelengths, R_V=3.1):
    """
    Cardelli, Clayton & Mathis (1989) extinction law.
    输入波长单位：Å
    返回 A(λ)/A_V，即 q_lambda
    """
    wl_um = wavelengths * 1e-4  # 转换为 μm
    x = 1.0 / wl_um  # 单位：μm^-1

    a = np.zeros_like(x)
    b = np.zeros_like(x)

    # IR & optical range: 0.3 < x < 1.1  (λ > ~0.91μm)
    idx_ir = (x >= 0.3) & (x < 1.1)
    a[idx_ir] = 0.574 * x[idx_ir] ** 1.61
    b[idx_ir] = -0.527 * x[idx_ir] ** 1.61

    # Optical/NIR: 1.1 <= x < 3.3
    idx_opt = (x >= 1.1) & (x < 3.3)
    y = x[idx_opt] - 1.82
    a[idx_opt] = (
        1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 +
        0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
    )
    b[idx_opt] = (
        1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 -
        5.38434 * y**4 - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7
    )

    # UV: 3.3 <= x <= 8.0
    idx_uv = (x >= 3.3) & (x <= 8.0)
    a[idx_uv] = 1.752 - 0.316 * x[idx_uv] - 0.104 / (
        (x[idx_uv] - 4.67)**2 + 0.341
    )
    b[idx_uv] = -3.090 + 1.825 * x[idx_uv] + 1.206 / (
        (x[idx_uv] - 4.62)**2 + 0.263
    )

    # Far-UV: 8.0 < x <= 10
    idx_fuv = (x > 8.0) & (x <= 10)
    a[idx_fuv] = -1.073 - 0.628 * (x[idx_fuv] - 8) + 0.137 * (x[idx_fuv] - 8)**2 - 0.070 * (x[idx_fuv] - 8)**3
    b[idx_fuv] = 13.670 + 4.257 * (x[idx_fuv] - 8) - 0.420 * (x[idx_fuv] - 8)**2 + 0.374 * (x[idx_fuv] - 8)**3

    A_lambda_div_AV = a + b / R_V
    return A_lambda_div_AV


def get_extinction_curve(wavelengths, law='CCM', R_V=3.1, custom_law_path=None):
    """
    通用接口：根据波长返回 extinction curve q_lambda（Aλ/A_V）
    """
    law = law.upper()

    if law == 'CCM':
        return extinction_ccm89(wavelengths, R_V)
    
    # 后续可以加 CAL, GD1, HZ1... 等分支
    else:
        raise NotImplementedError(f"Reddening law {law} not implemented yet.")
    

"""
添加其余红花定律前
使用方式：
from extinction import get_extinction_curve

wavelengths = np.arange(3500, 7501, 1)  # 以 Å 为单位
q_lambda = get_extinction_curve(wavelengths, law='CCM', R_V=3.1)
"""   

# 加 CAL, GD1
"""
注意：
先添加的这几个
	CAL：Calzetti et al. 2000
	GD1/GD2/GD3：使用手动表格或硬编码插值（简化）
	HZ1~HZ5：支持结构但先不写具体表达式（你有 HyperZ 表格时再加）
"""

def extinction_calzetti(wavelengths, R_V=4.05):
    """
    Calzetti et al. (2000) attenuation law
    波长单位：Å
    """
    wl = wavelengths / 1e4  # 转成 μm
    k = np.zeros_like(wl)

    # 0.12–0.63 μm
    mask1 = (wl >= 0.12) & (wl <= 0.63)
    k[mask1] = (
        2.659 * (-2.156 + 1.509/wl[mask1] - 0.198/(wl[mask1]**2) + 0.011/(wl[mask1]**3)) + R_V
    )
    
    # 0.63–2.2 μm
    mask2 = (wl > 0.63) & (wl <= 2.2)
    k[mask2] = (
        2.659 * (-1.857 + 1.040/wl[mask2]) + R_V
    )

    return k / R_V

def extinction_gordon(wavelengths, table='GD1'):
    """
    Gordon et al. (2003) extinction laws: SMC/LMC
    输入 table='GD1'|'GD2'|'GD3'
    """
    # Gordon+2003 Table 4: λ [μm], A(λ)/A_V
    law_tables = {
        'GD1': [  # SMC Bar
            (0.125, 7.75),
            (0.150, 5.91),
            (0.175, 4.89),
            (0.200, 4.16),
            (0.250, 3.20),
            (0.300, 2.62),
            (0.440, 1.76),
            (0.550, 1.00),
            (0.700, 0.71),
            (0.900, 0.43),
        ],
        'GD2': [  # LMC2 supershell
            (0.125, 6.95),
            (0.150, 5.58),
            (0.175, 4.67),
            (0.200, 4.01),
            (0.250, 3.12),
            (0.300, 2.55),
            (0.440, 1.75),
            (0.550, 1.00),
            (0.700, 0.70),
            (0.900, 0.46),
        ],
        'GD3': [  # LMC average
            (0.125, 5.61),
            (0.150, 4.54),
            (0.175, 3.80),
            (0.200, 3.29),
            (0.250, 2.53),
            (0.300, 2.05),
            (0.440, 1.58),
            (0.550, 1.00),
            (0.700, 0.72),
            (0.900, 0.51),
        ],
    }

    data = law_tables.get(table.upper())
    if data is None:
        raise ValueError(f"Gordon table {table} not recognized.")

    λ_tab, A_Av_tab = zip(*data)
    f = interp1d(λ_tab, A_Av_tab, kind='linear', bounds_error=False, fill_value='extrapolate')

    wl_um = wavelengths * 1e-4
    return f(wl_um)

def get_extinction_curve(wavelengths, law='CCM', R_V=3.1, custom_law_path=None):
    law = law.upper()

    if law == 'CCM':
        return extinction_ccm89(wavelengths, R_V)
    elif law == 'CAL':
        return extinction_calzetti(wavelengths, R_V)
    elif law in {'GD1', 'GD2', 'GD3'}:
        return extinction_gordon(wavelengths, law)
    elif law.startswith('HZ'):
        raise NotImplementedError(f"HyperZ law {law} not implemented yet.")
    elif custom_law_path is not None:
        # 用户自定义表格
        data = np.loadtxt(custom_law_path)
        wl_tab, A_Av_tab = data[:, 0], data[:, 1]
        f = interp1d(wl_tab, A_Av_tab, kind='linear', bounds_error=False, fill_value='extrapolate')
        return f(wavelengths)
    else:
        raise ValueError(f"Reddening law '{law}' not recognized.")
    

"""
使用示例：
from extinction import get_extinction_curve

wl = np.arange(3000, 8001)  # 单位：Å

q_ccm = get_extinction_curve(wl, law='CCM')
q_cal = get_extinction_curve(wl, law='CAL')
q_gd1 = get_extinction_curve(wl, law='GD1')
"""