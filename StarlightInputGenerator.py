import numpy as np
import os
from pathlib import Path
import warnings

class StarlightInputGenerator:
    def __init__(self, work_dir='./starlight_run'):
        """初始化工作目录"""
        self.work_dir = Path(work_dir).absolute()
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # 初始化所有输入组件
        self._init_spectrum()
        self._init_masks()
        self._init_base()
        self._init_config()
        self._init_grid()

    def _init_spectrum(self):
        """观测光谱参数"""
        self.obs_wave = None      # 波长数组 (Å)
        self.obs_flux = None      # 流量数组
        self.obs_error = None    # 误差数组
        self.obs_flags = None     # 标记数组
        self.has_header = False   # 是否包含头部行
        
    def _init_masks(self):
        """掩膜参数"""
        self.mask_regions = []    # 格式: [(start1, end1, weight1), ...]
        
    def _init_base(self):
        """基组参数"""
        self.base_components = []   # 基组元素列表，每个元素为字典
        self.base_dir = None        # 基组文件存储目录
        
    def _init_config(self):
        """配置参数 (仅列举关键参数，完整列表见手册)"""
        self.config_params = {
            'N_chains': 5,
            'dir_base': './BaseDir/',
            'IsClipActive': 1,
            'AV_min': 0.0,
            'AV_max': 3.0,
            'Is1stLineHeader': 0,
            'wei_nsig_threshold': 3.0,
        }
        
    def _init_grid(self):
        """网格文件参数"""
        self.grid_params = {
            'IsErrSpecAvailable': 1,
            'IsFlagSpecAvailable': 1,
            'syn_ini': 3500,
            'syn_fin': 8000,
            'dlsyn': 1,
        }

    def set_observation(self, wave, flux, error=None, flags=None, has_header=False):
        """设置观测光谱 (自动处理不同格式)"""
        # 单位检查
        if not all(np.diff(wave) > 0):
            warnings.warn("波长数组未排序或存在重复值！将自动排序。")
            sort_idx = np.argsort(wave)
            wave = wave[sort_idx]
            flux = flux[sort_idx]
            
        # 归一化检查
        if np.median(flux) < 1e-10:
            warnings.warn("流量值过小，请确认是否已归一化！")
            
        # 输入格式处理
        self.obs_wave = np.asarray(wave, dtype=float)
        self.obs_flux = np.asarray(flux, dtype=float)
        self.has_header = has_header
        
        # 处理误差和标记
        self.obs_error = np.full_like(flux, 0.05) if error is None else np.asarray(error)
        self.obs_flags = np.zeros_like(flux, dtype=int) if flags is None else np.asarray(flags)
        
        # 自动检测输入格式
        self.grid_params['IsErrSpecAvailable'] = int(error is not None)
        self.grid_params['IsFlagSpecAvailable'] = int(flags is not None)
        self.config_params['Is1stLineHeader'] = int(has_header)

    def add_mask_region(self, start, end, weight=0.0):
        """添加屏蔽区域 (支持多次调用)"""
        if start >= end:
            raise ValueError(f"无效波长范围: {start}-{end} Å")
        self.mask_regions.append((float(start), float(end), float(weight)))
        
    def _write_mask_file(self):
        """生成掩膜文件"""
        path = self.work_dir / 'Mask.custom'
        with open(path, 'w') as f:
            f.write(f"{len(self.mask_regions)}\n")
            for start, end, weight in self.mask_regions:
                f.write(f"{start:.1f} {end:.1f} {weight:.1f}\n")

    def set_base(self, base_dir, components):
        """
        设置基组参数
        :param base_dir: 基组光谱文件存储目录
        :param components: 列表，每个元素为字典格式：
            {
                'spec_file': 'bc2003_hr_m42_chab_ssp_020.spec',
                'age_yr': 1e9,
                'metallicity': 0.004,
                'nickname': 'age020_m42',
                'mass_frac': 1.0,
                'yav_flag': 0,
                'alpha_fe': 0.0
            }
        """
        self.base_dir = Path(base_dir).absolute()
        if not self.base_dir.exists():
            raise FileNotFoundError(f"基组目录不存在: {self.base_dir}")
            
        self.base_components = components
        
    def _write_base_files(self):
        """生成基组主文件"""
        # 生成主文件
        base_master_path = self.work_dir / 'Base.custom'
        with open(base_master_path, 'w') as f:
            f.write(f"{len(self.base_components)}\n")
            for comp in self.base_components:
                line = (
                    f"{comp['spec_file']} "
                    f"{comp['age_yr']:.5e} "
                    f"{comp['metallicity']:.5f} "
                    f"{comp['nickname']} "
                    f"{comp['mass_frac']:.4f} "
                    f"{comp['yav_flag']} "
                    f"{comp['alpha_fe']:.4f}\n"
                )
                f.write(line)
                
        # 验证基组光谱文件存在
        for comp in self.base_components:
            spec_path = self.base_dir / comp['spec_file']
            if not spec_path.exists():
                raise FileNotFoundError(f"基组光谱文件缺失: {spec_path}")
            
    def generate_all_inputs(self):
        """一键生成全部输入文件"""
        # 生成光谱文件
        self._write_spectrum_file()
        
        # 生成其他文件
        self._write_mask_file()
        self._write_base_files()
        self._write_config_file()
        self._write_grid_file()
        
        # 创建输出目录
        (self.work_dir / 'output').mkdir(exist_ok=True)
        
        print(f"所有输入文件已生成至: {self.work_dir}")

    def _write_spectrum_file(self):
        """生成观测光谱文件 (自动处理不同格式)"""
        path = self.work_dir / 'spectrum.cxt'
        cols = [self.obs_wave, self.obs_flux]
        
        # 根据可用性添加列
        if self.grid_params['IsErrSpecAvailable']:
            cols.append(self.obs_error)
        if self.grid_params['IsFlagSpecAvailable']:
            cols.append(self.obs_flags)
            
        data = np.column_stack(cols)
        header = "# STARLIGHT Input Spectrum" if self.has_header else ""
        np.savetxt(path, data, fmt=self._get_spectrum_fmt(), header=header)

    def _get_spectrum_fmt(self):
        """根据数据格式确定保存格式"""
        ncols = 2 + self.grid_params['IsErrSpecAvailable'] + self.grid_params['IsFlagSpecAvailable']
        fmt_dict = {
            2: ['%.2f', '%.6e'],
            3: ['%.2f', '%.6e', '%.6e'],
            4: ['%.2f', '%.6e', '%.6e', '%d']
        }
        return '  '.join(fmt_dict[ncols])