# io.py
# 光谱文件读取函数
# 此函数是读取数据输入文件(cxt文件)

import numpy as np
from scipy.ndimage import uniform_filter1d


def load_spectrum(filename: str,
                  is_header=True,
                  smooth=False,
                  smooth_window=5,
                  err_min=1e-5,
                  err_max=1e5):
    """
    加载 STARLIGHT 光谱文件，支持自动格式识别和异常值保护。

    参数:
        filename: 光谱文件路径
        is_header: 是否跳过第一行
        smooth: 是否对 flux 做滑动平均平滑
        smooth_window: 平滑窗口大小
        err_min: 最小允许 error（防止数值发散）
        err_max: 最大允许 error（防止异常影响）
    
    返回:
        wavelength: np.ndarray
        flux: np.ndarray
        error: np.ndarray
        flags: np.ndarray
        IsErrSpecAvailable: bool
        IsFlagSpecAvailable: bool
    """
    wavelengths = []
    fluxes = []
    errors = []
    flags = []

    with open(filename, 'r') as f:
        lines = f.readlines()
        if is_header:
            lines = lines[1:]

        for line in lines:
            if line.strip() == '':
                continue
            parts = line.strip().split()

            if len(parts) == 4:
                wl, fl, err, flag = parts
                wavelengths.append(float(wl))
                fluxes.append(float(fl))
                errors.append(float(err))
                flags.append(int(flag))
            elif len(parts) == 3:
                wl, fl, err = parts
                wavelengths.append(float(wl))
                fluxes.append(float(fl))
                errors.append(float(err))
                flags.append(0)
            elif len(parts) == 2:
                wl, fl = parts
                wavelengths.append(float(wl))
                fluxes.append(float(fl))
                errors.append(1.0)  # 默认 eλ
                flags.append(0)
            else:
                raise ValueError(f"Invalid spectrum line: {line.strip()}")

    wl = np.array(wavelengths)
    flux = np.array(fluxes)
    err = np.array(errors)
    flag = np.array(flags)

    # 对 error 做限制（避免过小或过大）
    err = np.clip(err, err_min, err_max)

    # 简单平滑（仅作用于 flux）
    if smooth:
        flux = uniform_filter1d(flux, size=smooth_window, mode='nearest')

    # 自动判断格式
    is_err_spec = not np.all(err == 1.0)
    is_flag_spec = not np.all(flag == 0)

    return wl, flux, err, flag, is_err_spec, is_flag_spec


"""
调试：
from io import load_spectrum

wl, flux, err, flag, has_err, has_flag = load_spectrum("0414.51901.393.cxt", smooth=True)

print("波长范围：", wl[0], "→", wl[-1])
print("是否有误差谱？", has_err)
print("是否有 flag？", has_flag)
"""