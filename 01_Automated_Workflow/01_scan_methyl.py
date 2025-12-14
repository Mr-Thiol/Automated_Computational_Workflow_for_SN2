# run_batch_safe.py
# 带有【断点续传】功能的批量计算脚本

import os
import shutil
import time
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms
import numpy as np

# 1. 导入工具箱
try:
    from chem_utils import *
except ImportError:
    print("❌ 错误：找不到 chem_utils.py。")
    exit()

# 硬件配置
N_PROC = 6
MEM_SIZE = "4GB"

def get_rotation_matrix_to_align_vector(vec_start, vec_end, target_axis=np.array([0, 0, 1])):
    """计算旋转矩阵：将向量旋转到 Z 轴"""
    v_current = np.array(vec_end) - np.array(vec_start)
    norm = np.linalg.norm(v_current)
    if norm < 1e-6: return np.eye(4) # 避免除零
    v_current = v_current / norm
    
    if np.allclose(v_current, target_axis): return np.eye(4)
    if np.allclose(v_current, -target_axis): return np.array([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,1]])
    
    axis = np.cross(v_current, target_axis)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-6: return np.eye(4) # 应该不会发生
    axis = axis / axis_norm
    
    angle = np.arccos(np.dot(v_current, target_axis))
    c, s = np.cos(angle), np.sin(angle)
    t = 1 - c
    x, y, z = axis
    
    return np.array([
        [t*x*x + c,   t*x*y - z*s, t*x*z + y*s, 0],
        [t*x*y + z*s, t*y*y + c,   t*y*z - x*s, 0],
        [t*x*z - y*s, t*y*z + x*s, t*z*z + c,   0],
        [0,           0,           0,           1]
    ])


# --- 核心保底函数 ---
def check_if_job_finished(log_filename):
    """
    检查高斯任务是否已经【正常结束】。
    如果不正常（没跑完、报错），返回 False，脚本就会重跑它。
    """
    if not os.path.exists(log_filename):
        return False
    
    try:
        # 只读取文件最后几行，避免读取大文件占用内存
        with open(log_filename, 'rb') as f:
            try:
                f.seek(-2048, os.SEEK_END) # 倒退 2KB 读取
            except OSError:
                f.seek(0) # 文件太小就从头读
            
            last_content = f.read().decode('utf-8', errors='ignore')
            
            # Gaussian 的成功标志
            if "Normal termination" in last_content:
                return True
    except Exception:
        return False
    
    return False

# --- 构建函数 (保持不变) ---
def build_aligned_system(amine_smiles):
    # A. 构建胺
    mol_amine = Chem.MolFromSmiles(amine_smiles)
    mol_amine = Chem.AddHs(mol_amine)
    AllChem.EmbedMolecule(mol_amine, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_amine)
    
    conf = mol_amine.GetConformer()
    # 找 N 原子
    n_atom = [a for a in mol_amine.GetAtoms() if a.GetSymbol() == 'N'][0]
    n_idx = n_atom.GetIdx()
    
    # 计算几何中心用于对齐
    neighbors = n_atom.GetNeighbors()
    if len(neighbors) > 0:
        centroid = np.mean([np.array(conf.GetAtomPosition(nbr.GetIdx())) for nbr in neighbors], axis=0)
    else:
        centroid = np.array([0.0, 0.0, -1.0]) # NH3 case

    # 旋转对齐
    rot_mat = get_rotation_matrix_to_align_vector(centroid, np.array(conf.GetAtomPosition(n_idx)), np.array([0, 0, 1]))
    rdMolTransforms.TransformConformer(conf, rot_mat)
    
    # 平移归零
    pos = conf.GetAtomPosition(n_idx)
    shift = Geometry.Point3D(-pos.x, -pos.y, -pos.z)
    for i in range(mol_amine.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Geometry.Point3D(p.x+shift.x, p.y+shift.y, p.z+shift.z))
        
    # B. 构建 MeI
    mol_mei = Chem.MolFromSmiles("CI")
    mol_mei = Chem.AddHs(mol_mei)
    AllChem.EmbedMolecule(mol_mei, randomSeed=42)
    conf_mei = mol_mei.GetConformer()
    c_idx = [a.GetIdx() for a in mol_mei.GetAtoms() if a.GetSymbol() == 'C'][0]
    i_idx = [a.GetIdx() for a in mol_mei.GetAtoms() if a.GetSymbol() == 'I'][0]
    h_idxs = [a.GetIdx() for a in mol_mei.GetAtoms() if a.GetSymbol() == 'H']
    
    conf_mei.SetAtomPosition(c_idx, Geometry.Point3D(0.0, 0, 3.0)) 
    conf_mei.SetAtomPosition(i_idx, Geometry.Point3D(0.0, 0, 5.2))
    offsets = [(0.0, 1.03, 2.65), (0.89, -0.51, 2.65), (-0.89, -0.51, 2.65)]
    for idx, off in zip(h_idxs, offsets):
        conf_mei.SetAtomPosition(idx, Geometry.Point3D(off[0], off[1], off[2]))
        
    return Chem.CombineMols(mol_amine, mol_mei)

# --- 主程序 ---
if __name__ == "__main__":
    
    if not shutil.which("g16"):
        print("❌ 找不到 g16，请检查环境！")
        exit()

    # 任务列表 (建议先用这三个跑通)
    target_amines = [
        ("MeNH2", "CN"),
        ("Me2NH", "CNC"),
        ("Me3N",  "CN(C)C"),
        ("EtNH2", "CCN"),
        # ("iPrNH2", "CC(C)N"), # 想跑更多可以随时加
    ]

    # 打开 CSV 文件准备记录 (使用 append 模式 'a'，防止覆盖之前的记录)
    csv_file = "batch_results.csv"
    if not os.path.exists(csv_file):
        with open(csv_file, "w") as f:
            f.write("Name,Barrier_kcal_mol,Note\n")

    print(f">>> 启动安全版批量脚本 (支持断点续传)...")
    
    for name, smiles in target_amines:
        print(f"\n{'-'*20} {name} {'-'*20}")
        
        scan_file = f"{name}_scan.gjf"
        log_file = f"{name}_scan.log"
        
        # === 🛡️ 保底检查 1: Scan 算完没？ ===
        if check_if_job_finished(log_file):
            print(f"⏩ {name} 的扫描任务已完成 (Normal termination)，跳过计算。")
            # 即使跳过，也要尝试读取一下结果，为了算后面的 TS
            # (这里省略读取逻辑，默认如果 Scan 完了，我们假设它没问题)
            # 如果你想做得更细，可以在这里读一下 points
        else:
            print(f"🚀 开始计算 {name} (Scan)...")
            try:
                # 1. 构建
                mol_sys = build_aligned_system(smiles)
                # 自动找序号
                conf = mol_sys.GetConformer()
                n_idx, c_idx = -1, -1
                for atom in mol_sys.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    if atom.GetSymbol() == 'N' and abs(pos.z) < 0.2: n_idx = atom.GetIdx() + 1
                    if atom.GetSymbol() == 'C' and abs(pos.z - 3.0) < 0.5: c_idx = atom.GetIdx() + 1
                
                # 2. 运行 Scan
                scan_lines = [f"{n_idx} {c_idx} S 15 -0.10"]
                create_advanced_scan_gaussian_input(
                    mol_sys, f"{name}_Scan", scan_lines, charge=0,
                    nproc=N_PROC, mem=MEM_SIZE,
                    method="#PM7 opt(modredundant)", 
                    filename=scan_file
                )
                run_gaussian_job(scan_file) # 这是一个阻塞操作，算完才会往下走
                
            except Exception as e:
                print(f"❌ 构建/提交阶段出错: {e}")
                continue # 出错就跳过这个分子，别卡死

        # === 🛡️ 保底检查 2: 结果记录 ===
        # 无论刚才是在算，还是刚跳过，我们都尝试读取 Log 来记录结果
        # 这样即使你删了 CSV，只要 Log 在，重跑一遍脚本就能自动重建 CSV
        try:
            points = read_scan_output(log_file)
            if points:
                max_point = max(points, key=lambda x: x['energy'])
                barrier = (max_point['energy'] - points[0]['energy']) * 627.5
                print(f"📊 {name} 能垒: {barrier:.2f} kcal/mol")
                
                # 实时写入 CSV (防止脚本最后崩溃导致数据没存)
                with open(csv_file, "a") as f:
                    f.write(f"{name},{barrier:.4f},Scan_Done\n")
            else:
                print(f"⚠️ {log_file} 读取失败或为空 (可能计算崩了)。")
        except Exception as e:
            print(f"⚠️ 结果解析出错: {e}")

    print("\n>>> 所有任务检查/计算完毕！")