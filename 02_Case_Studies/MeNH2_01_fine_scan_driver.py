# retry_MeNH2_fine.py
# 专门用于修复 MeNH2 的精细扫描脚本

import os
import shutil
import sys
import numpy as np
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms

# --- 导入工具箱 ---
try:
    from chem_utils import *
except ImportError:
    try:
        from chem_utils_new import *
    except ImportError:
        print("❌ 找不到 chem_utils.py")
        sys.exit()

# 硬件配置 (MeNH2 很小，单核跑都很快，用 6 核飞快)
N_PROC = 6
MEM_SIZE = "4GB"

# --- 核心对齐函数 (保持一致性) ---
def get_rotation_matrix_to_align_vector(vec_start, vec_end, target_axis=np.array([0, 0, 1])):
    v_current = np.array(vec_end) - np.array(vec_start)
    norm = np.linalg.norm(v_current)
    if norm < 1e-6: return np.eye(4)
    v_current = v_current / norm
    if np.allclose(v_current, target_axis): return np.eye(4)
    if np.allclose(v_current, -target_axis): return np.array([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,1]])
    axis = np.cross(v_current, target_axis)
    axis = axis / np.linalg.norm(axis)
    angle = np.arccos(np.dot(v_current, target_axis))
    c, s = np.cos(angle), np.sin(angle)
    t = 1 - c
    x, y, z = axis
    return np.array([[t*x*x+c, t*x*y-z*s, t*x*z+y*s, 0],
                     [t*x*y+z*s, t*y*y+c, t*y*z-x*s, 0],
                     [t*x*z-y*s, t*y*z+x*s, t*z*z+c, 0],
                     [0, 0, 0, 1]])

def build_menh2_system():
    """专门构建 MeNH2 + MeI"""
    print(">>> 正在构建 MeNH2 + MeI 体系...")
    mol_amine = Chem.MolFromSmiles("CN") # MeNH2
    mol_amine = Chem.AddHs(mol_amine)
    AllChem.EmbedMolecule(mol_amine, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_amine)
    
    conf = mol_amine.GetConformer()
    n_atom = [a for a in mol_amine.GetAtoms() if a.GetSymbol() == 'N'][0]
    n_idx = n_atom.GetIdx()
    
    # 找几何中心
    neighbors = n_atom.GetNeighbors()
    centroid = np.mean([np.array(conf.GetAtomPosition(nbr.GetIdx())) for nbr in neighbors], axis=0)
    
    # 对齐
    rot_mat = get_rotation_matrix_to_align_vector(centroid, np.array(conf.GetAtomPosition(n_idx)))
    rdMolTransforms.TransformConformer(conf, rot_mat)
    
    # 平移 N 到原点
    pos = conf.GetAtomPosition(n_idx)
    shift = Geometry.Point3D(-pos.x, -pos.y, -pos.z)
    for i in range(mol_amine.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Geometry.Point3D(p.x+shift.x, p.y+shift.y, p.z+shift.z))
    
    # MeI
    mol_mei = Chem.MolFromSmiles("CI")
    mol_mei = Chem.AddHs(mol_mei)
    AllChem.EmbedMolecule(mol_mei, randomSeed=42)
    conf_mei = mol_mei.GetConformer()
    c_idx = [a.GetIdx() for a in mol_mei.GetAtoms() if a.GetSymbol() == 'C'][0]
    i_idx = [a.GetIdx() for a in mol_mei.GetAtoms() if a.GetSymbol() == 'I'][0]
    h_idxs = [a.GetIdx() for a in mol_mei.GetAtoms() if a.GetSymbol() == 'H']
    
    # 【关键修改】起始位置设为 2.4 埃，而不是 3.0
    # 这样我们可以扫得更细，而不浪费时间在远处
    start_dist = 2.4
    conf_mei.SetAtomPosition(c_idx, Geometry.Point3D(0.0, 0, start_dist)) 
    conf_mei.SetAtomPosition(i_idx, Geometry.Point3D(0.0, 0, start_dist + 2.2))
    
    # 摆放 H
    offsets = [(0.0, 1.03, start_dist-0.35), (0.89, -0.51, start_dist-0.35), (-0.89, -0.51, start_dist-0.35)]
    for idx, off in zip(h_idxs, offsets):
        conf_mei.SetAtomPosition(idx, Geometry.Point3D(off[0], off[1], off[2]))
        
    return Chem.CombineMols(mol_amine, mol_mei)

# --- 主程序 ---
if __name__ == "__main__":
    name = "MeNH2_Fine" # 改个名字，防止覆盖原来的
    
    # 1. 构建
    mol_sys = build_menh2_system()
    
    # 找序号
    conf = mol_sys.GetConformer()
    n_idx, c_idx = -1, -1
    for atom in mol_sys.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        if atom.GetSymbol() == 'N' and abs(pos.z) < 0.2: n_idx = atom.GetIdx() + 1
        if atom.GetSymbol() == 'C' and pos.z > 1.0: c_idx = atom.GetIdx() + 1
    
    # 2. 精细扫描 (Fine Scan)
    # 从 2.4 扫到 1.6，共 0.8 距离
    # 步长 0.04 (比之前细 2.5倍)，共 20 步
    # 这样肯定能抓到最高点
    print(f"🚀 启动精细扫描 (Start=2.4Å, Step=-0.04Å)...")
    
    scan_file = f"{name}_scan.gjf"
    scan_lines = [f"{n_idx} {c_idx} S 20 -0.04"] 
    
    create_advanced_scan_gaussian_input(
        mol=mol_sys, mol_name=f"{name}_Scan", 
        scan_lines=scan_lines, charge=0,
        nproc=N_PROC, mem=MEM_SIZE,
        method="#PM7 opt(modredundant)", 
        filename=scan_file
    )
    
    if run_gaussian_job(scan_file):
        # 3. 提取最高点
        points = read_scan_output(f"{name}_scan.log")
        if points:
            max_point = max(points, key=lambda x: x['energy'])
            print(f"📊 锁定精细最高点: Step {max_point['step']}, Energy={max_point['energy']:.5f}")
            
            # 4. 再次尝试 TS 优化
            print(f"🚀 启动 TS 优化 (PM7 calcall)...")
            ts_guess, _ = extract_scan_point_to_mol(points, max_point['step'], mol_sys, show_3d=False)
            
            ts_file = f"{name}_ts.gjf"
            ts_log = f"{name}_ts.log"
            
            # 使用 calcall 强力计算力常数
            create_gaussian_input_advanced(
                mol=ts_guess, mol_name=f"{name}_TS", 
                scan_lines=[], charge=0,
                nproc=N_PROC, mem=MEM_SIZE,
                method="#PM7 opt(ts,calcall,noeig) freq", 
                filename=ts_file
            )
            
            run_gaussian_job(ts_file)
            
            # 5. 验证
            res = read_gaussian16_output_opt(ts_log)
            freqs = res.get('frequencies', [])
            if freqs and freqs[0] < 0:
                print(f"✅✅✅ 成功救回 MeNH2! 虚频: {freqs[0]:.2f} cm^-1")
            else:
                print(f"❌ 依然失败。放弃 MeNH2，专心做 Me2NH 的 Pre 吧。")
        else:
            print("❌ 扫描读取失败")