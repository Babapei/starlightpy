# io.py（继续添加）
# 这两个函数是读取base文件以及加载模版光谱

def load_base_master(master_path: str):
    """
    读取 Base Master 文件（如 Base.BC03.N）
    
    返回：
        - base_list: List of dict, 每个 base 包含文件名、age、Z、nickname、fstar、YAV、alphaFe
    """
    base_list = []

    with open(master_path, 'r') as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("Base master file is empty.")
    
    try:
        n_base = int(lines[0].strip())
    except ValueError:
        raise ValueError("First line of base master must be an integer (number of bases)")

    for line in lines[1:n_base+1]:
        parts = line.strip().split()
        if len(parts) < 7:
            raise ValueError(f"Base line malformed: {line}")
        
        base_list.append({
            "filename": parts[0],
            "age": float(parts[1]),
            "Z": float(parts[2]),
            "nickname": parts[3],
            "f_star": float(parts[4]),
            "YAV": int(parts[5]),
            "alpha_Fe": float(parts[6])
        })

    return base_list



import os

import numpy as np

def load_base_spectra(base_dir: str, base_list: list, verbose=False):
    """
    加载 base_list 中每个 base 的光谱文件，组合为 λ × N 的 base 矩阵。
    
    参数：
        base_dir: 所有 base 光谱文件所在文件夹
        base_list: load_base_master() 的输出
    返回：
        wl_base: 所有 base 共用的波长数组
        base_matrix: shape=(N_wavelengths, N_components)
    """
    all_spectra = []
    wl_ref = None

    for base in base_list:
        path = os.path.join(base_dir, base["filename"])

        wavelengths = []
        fluxes = []
        with open(path, 'r') as f:
            for line in f:
                if line.startswith("#") or line.strip() == "":
                    continue
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                wl, flux = float(parts[0]), float(parts[1])
                wavelengths.append(wl)
                fluxes.append(flux)

        wavelengths = np.array(wavelengths)
        fluxes = np.array(fluxes)

        # 第一次加载时记录波长轴
        if wl_ref is None:
            wl_ref = wavelengths
        else:
            if not np.allclose(wl_ref, wavelengths):
                raise ValueError(f"Base file {path} has different wavelength sampling.")

        all_spectra.append(fluxes)

        if verbose:
            print(f"Loaded base: {base['filename']}, shape={fluxes.shape}")

    base_matrix = np.column_stack(all_spectra)
    return wl_ref, base_matrix


"""
base_list = load_base_master("Base.BC03.N")
wl_base, base_matrix = load_base_spectra("BasesDir", base_list)

print("Base 分量数：", len(base_list))
print("Base 矩阵大小：", base_matrix.shape)  # (波长数, 成分数)

"""

"""
后续建议补充改进：
支持•gz压缩光谱
    对BC03 的压缩版支持更好
自动缓存 base 矩阵
    提升加载速度（如 pickle + hash）
对f_star 进行预处理
    便于 mass/light 之间转换
"""