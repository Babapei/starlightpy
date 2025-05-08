"""
STARLIGHT 输入生成库
包含观测数据、配置、网格文件的生成与管理
"""

import os
import logging
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Union
from dataclasses import dataclass

# ----------------------
# 自定义异常类
# ----------------------
class StarlightError(Exception):
    """库级基础异常"""
    pass

class InvalidSpectralDataError(StarlightError):
    """光谱数据不合法异常"""
    pass

class FileGenerationError(StarlightError):
    """文件生成失败异常"""
    pass

# ----------------------
# 数据结构定义
# ----------------------
@dataclass
class SpectralData:
    """标准化光谱数据容器"""
    wavelength: np.ndarray  # 波长数组 (Å)
    flux: np.ndarray        # 流量数组 (erg/s/cm²/Å)
    error: Optional[np.ndarray] = None  # 误差数组
    flags: Optional[np.ndarray] = None  # 标志数组 (0=正常, ≥1=掩码)

# ----------------------
# 输入文件生成核心类
# ----------------------
class StarlightInputGenerator:
    """生成观测数据文件(.cxt)和掩膜文件"""
    
    def __init__(self, work_dir: Union[str, Path]):
        self.work_dir = Path(work_dir).resolve()
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(f"{__name__}.InputGenerator")
        
    def _validate_spectral_data(self, data: SpectralData) -> None:
        """执行严格的光谱数据校验"""
        # 波长校验
        if not np.allclose(np.diff(data.wavelength), 1.0):
            raise InvalidSpectralDataError("波长必须为1Å等间距采样")
            
        # 流量非负校验
        if np.any(data.flux < 0):
            self.logger.warning("流量包含负值，自动设置对应flag≥2")
            if data.flags is None:
                data.flags = np.zeros_like(data.flux, dtype=int)
            data.flags[data.flux < 0] = 2
        
        # 维度一致性校验
        if data.error is not None and len(data.error) != len(data.wavelength):
            raise InvalidSpectralDataError("误差数组维度与波长不一致")
        if data.flags is not None and len(data.flags) != len(data.wavelength):
            raise InvalidSpectralDataError("标志数组维度与波长不一致")
    
    def generate_observation_file(
        self,
        data: SpectralData,
        filename: str
    ) -> Path:
        """
        生成STARLIGHT观测文件(.cxt)
        
        参数:
            data: 标准化光谱数据
            filename: 输出文件名(不带后缀)
            
        返回:
            生成文件的Path对象
            
        异常:
            FileGenerationError: 文件生成失败时抛出
        """
        try:
            self._validate_spectral_data(data)
            
            output_path = self.work_dir / f"{filename}.cxt"
            header = (
                "# lambda flux [error] [flag]\n"
                "# Created by StarlightInputGenerator\n"
            )
            
            # 构建数据列
            columns = [data.wavelength, data.flux]
            if data.error is not None:
                columns.append(data.error)
            if data.flags is not None:
                columns.append(data.flags.astype(int))
                
            # 写入文件
            np.savetxt(
                output_path,
                np.column_stack(columns),
                header=header,
                fmt="%.3f %.6e" + (" %.6e" if data.error else "") + (" %d" if data.flags else "")
            )
            
            self.logger.info(f"成功生成观测文件: {output_path}")
            return output_path
            
        except Exception as e:
            self.logger.error(f"观测文件生成失败: {str(e)}")
            raise FileGenerationError(f"无法生成 {filename}.cxt") from e

# ----------------------
# 配置生成模块
# ----------------------
class StarlightConfigGenerator:
    """生成STARLIGHT配置文件(.config)"""
    
    _DEFAULT_CONFIG = {
        # 基础参数
        'N_chains': 5,
        'IsClipActive': 1,
        'clip_nsig': 3.0,
        
        # 消光参数
        'AV_min': 0.0,
        'AV_max': 3.0,
        'RV': 3.1,
        
        # 运动学参数
        'v0_low': -500.0,
        'v0_upp': 500.0,
        'vd_low': 50.0,
        'vd_upp': 300.0,
        
        # 物理约束
        'X_min': 0.0,
        'X_max': 1.0,
        'Y_min': 0.0,
        'Y_max': 1.0,
    }
    
    _PARAM_COMMENTS = {
        'N_chains': "马尔可夫链数量 (推荐5-10)",
        'IsClipActive': "异常值剔除 (0=关闭, 1=弱剔除, 2=强剔除)",
        'AV_max': "最大消光值 (根据目标类型调整)"
    }
    
    def __init__(self):
        self.params = self._DEFAULT_CONFIG.copy()
        self.logger = logging.getLogger(f"{__name__}.ConfigGenerator")
    
    def set_parameter(self, name: str, value: Union[float, int]) -> None:
        """安全设置配置参数"""
        if name not in self._DEFAULT_CONFIG:
            raise KeyError(f"无效参数名: {name}. 有效参数: {list(self._DEFAULT_CONFIG.keys())}")
        self.params[name] = value
    
    def generate_config_file(self, output_path: Union[str, Path]) -> Path:
        """生成配置文件"""
        try:
            output_path = Path(output_path).resolve()
            lines = []
            
            for param, value in self.params.items():
                comment = self._PARAM_COMMENTS.get(param, "详见STARLIGHT手册4.4节")
                lines.append(f"{param} = {value}  # {comment}")
            
            with open(output_path, 'w') as f:
                f.write("# STARLIGHT 配置文件\n")
                f.write("# 自动生成 - 请勿手动修改\n")
                f.write("\n".join(lines))
            
            self.logger.info(f"成功生成配置文件: {output_path}")
            return output_path
            
        except Exception as e:
            self.logger.error(f"配置文件生成失败: {str(e)}")
            raise FileGenerationError(f"无法生成 {output_path.name}") from e

# ----------------------
# 网格文件生成模块
# ----------------------
@dataclass
class FitTask:
    """单个拟合任务参数容器"""
    obs_file: str          # 观测文件名
    config_file: str       # 配置文件
    base_file: str         # 基组文件
    mask_file: str         # 掩膜文件
    reddening_law: str     # 消光曲线 (e.g., 'CCM')
    v0_start: float = 0.0  # 初始速度偏移 (km/s)
    vd_start: float = 150.0  # 初始速度弥散 (km/s)

class StarlightGridGenerator:
    """生成STARLIGHT网格文件(grid.in)"""
    
    def __init__(self, work_dir: Union[str, Path]):
        self.work_dir = Path(work_dir).resolve()
        self.tasks: List[FitTask] = []
        self.global_params = {
            'seed': 123456,
            'llow_SN': 4730.0,
            'lupp_SN': 4780.0,
            'syn_ini': 3400.0,
            'syn_fin': 8900.0,
            'dlsyn': 1.0,
            'fit_option': 'FIT'
        }
        self.logger = logging.getLogger(f"{__name__}.GridGenerator")
    
    def add_task(self, task: FitTask) -> None:
        """添加拟合任务"""
        # 验证消光曲线
        valid_laws = {'CCM', 'CAL', 'GD1', 'GD2', 'GD3', 'HZ1', 'HZ2', 'HZ3', 'HZ4', 'HZ5'}
        if task.reddening_law not in valid_laws:
            raise ValueError(f"无效消光曲线: {task.reddening_law}. 有效选项: {valid_laws}")
        
        self.tasks.append(task)
        self.logger.debug(f"添加任务: {task.obs_file}")
    
    def generate_grid_file(self, filename: str = "grid.in") -> Path:
        """生成完整网格文件"""
        try:
            output_path = self.work_dir / filename
            
            # 生成头部
            header = [
                f"{len(self.tasks)}  # 拟合任务数",
                f"{self.work_dir / 'bases/'}  # 基组目录",
                f"{self.work_dir / 'obs/'}  # 观测文件目录",
                f"{self.work_dir / 'masks/'}  # 掩膜目录",
                f"{self.work_dir / 'output/'}  # 输出目录",
                f"{self.global_params['seed']}  # 随机种子",
                f"{self.global_params['llow_SN']}  # SNR计算下限波长",
                f"{self.global_params['lupp_SN']}  # SNR计算上限波长",
                f"{self.global_params['syn_ini']}  # 拟合起始波长",
                f"{self.global_params['syn_fin']}  # 拟合结束波长",
                f"{self.global_params['dlsyn']}  # 波长采样间隔",
                "1.0  # chi2缩放因子",
                f"{self.global_params['fit_option']}  # 运动学拟合选项(FIT/FXK)",
                "1  # 存在误差谱",
                "1  # 存在标志谱"
            ]
            
            # 生成任务行
            task_lines = []
            for task in self.tasks:
                line = " ".join([
                    task.obs_file,
                    task.config_file,
                    task.base_file,
                    task.mask_file,
                    task.reddening_law,
                    f"{task.v0_start:.1f}",
                    f"{task.vd_start:.1f}",
                    f"output_{Path(task.obs_file).stem}.out"
                ])
                task_lines.append(line)
            
            # 写入文件
            with open(output_path, 'w') as f:
                f.write("\n".join(header) + "\n")
                f.write("\n".join(task_lines))
            
            self.logger.info(f"成功生成网格文件: {output_path}")
            return output_path
            
        except Exception as e:
            self.logger.error(f"网格文件生成失败: {str(e)}")
            raise FileGenerationError(f"无法生成 {filename}") from e

# ----------------------
# 使用示例
# ----------------------
if __name__ == "__main__":
    # 配置日志
    logging.basicConfig(level=logging.INFO)
    
    # 初始化工作目录
    work_dir = Path("starlight_run")
    work_dir.mkdir(exist_ok=True)
    
    try:
        # 生成观测文件
        input_gen = StarlightInputGenerator(work_dir / "obs")
        data = SpectralData(
            wavelength=np.arange(3400, 8901),
            flux=np.random.randn(5501).clip(0),
            error=0.1 * np.ones(5501),
            flags=np.zeros(5501, dtype=int)
        )
        obs_file = input_gen.generate_observation_file(data, "sample_galaxy")
        
        # 生成配置文件
        config_gen = StarlightConfigGenerator()
        config_gen.set_parameter('AV_max', 2.5)
        config_file = config_gen.generate_config_file(work_dir / "configs/sample.config")
        
        # 配置网格文件
        grid_gen = StarlightGridGenerator(work_dir)
        grid_gen.add_task(FitTask(
            obs_file=obs_file.name,
            config_file=config_file.name,
            base_file="Base.BC03.N",
            mask_file="Mask.general",
            reddening_law='CCM'
        ))
        grid_file = grid_gen.generate_grid_file()
        
        print(f"成功生成所有输入文件于: {work_dir}")
        
    except StarlightError as e:
        print(f"STARLIGHT 错误发生: {str(e)}")
        exit(1)