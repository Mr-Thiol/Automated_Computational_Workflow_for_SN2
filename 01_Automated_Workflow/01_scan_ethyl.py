# batch_Et_series.py
# 专门用于修复 EtNH2 并计算 Et2NH, NEt3

import os
import shutil
import sys
import numpy as np
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms

# 导入工具箱
try:
    from chem_utils import *
except ImportError:
    try:
        from chem_utils_new import *
    except ImportError:
        print("❌ 找不到 chem_utils.py")
        sys.exit()

# 硬件配置
N_PROC = 6
MEM_SIZE = "4GB"

# --- 核心对齐函数 (保持一致) ---
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

def build_aligned_system_fine(amine_smiles):
    """构建 N...C-I 体系，专门针对精细扫描优化初始距离"""
    # A. 胺
    mol_amine = Chem.MolFromSmiles(amine_smiles)
    mol_amine = Chem.AddHs(mol_amine)
    AllChem.EmbedMolecule(mol_amine, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_amine)
    
    conf = mol_amine.GetConformer()
    # 找 N 原子
    n_atom = [a for a in mol_amine.GetAtoms() if a.GetSymbol() == 'N'][0]
    n_idx = n_atom.GetIdx()
    
    # 计算几何中心 (用于确定 N 的朝向)
    neighbors = n_atom.GetNeighbors()
    if len(neighbors) > 0:
        centroid = np.mean([np.array(conf.GetAtomPosition(nbr.GetIdx())) for nbr in neighbors], axis=0)
    else:
        centroid = np.array([0.0, 0.0, -1.0]) 

    # 旋转对齐到 Z 轴
    rot_mat = get_rotation_matrix_to_align_vector(centroid, np.array(conf.GetAtomPosition(n_idx)), np.array([0, 0, 1]))
    rdMolTransforms.TransformConformer(conf, rot_mat)
    
    # 平移 N 到原点
    pos = conf.GetAtomPosition(n_idx)
    shift = Geometry.Point3D(-pos.x, -pos.y, -pos.z)
    for i in range(mol_amine.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Geometry.Point3D(p.x+shift.x, p.y+shift.y, p.z+shift.z))
        
    # B. 碘甲烷 (放置在 2.4 埃处，开始精细扫描)
    mol_mei = Chem.MolFromSmiles("CI")
    mol_mei = Chem.AddHs(mol_mei)
    AllChem.EmbedMolecule(mol_mei, randomSeed=42)
    conf_mei = mol_mei.GetConformer()
    c_idx = [a.GetIdx() for a in mol_mei.GetAtoms() if a.GetSymbol() == 'C'][0]
    i_idx = [a.GetIdx() for a in mol_mei.GetAtoms() if a.GetSymbol() == 'I'][0]
    h_idxs = [a.GetIdx() for a in mol_mei.GetAtoms() if a.GetSymbol() == 'H']
    
    # 初始距离设为 2.4 (比之前的 3.0 近，节省计算量且更准)
    start_dist = 2.4
    conf_mei.SetAtomPosition(c_idx, Geometry.Point3D(0.0, 0, start_dist)) 
    conf_mei.SetAtomPosition(i_idx, Geometry.Point3D(0.0, 0, start_dist + 2.14)) # C-I 键长约 2.14
    
    # 调整 MeI 的 H 方向，避免与胺撞车
    offsets = [(0.0, 1.03, start_dist-0.35), (0.89, -0.51, start_dist-0.35), (-0.89, -0.51, start_dist-0.35)]
    for idx, off in zip(h_idxs, offsets):
        conf_mei.SetAtomPosition(idx, Geometry.Point3D(off[0], off[1], off[2]))
        
    return Chem.CombineMols(mol_amine, mol_mei)

# --- 主程序 ---
if __name__ == "__main__":
    
    if not shutil.which("g16"):
        print("❌ 找不到 g16")
        exit()

    # 定义要跑的任务：名字, SMILES
    tasks = [
        ("EtNH2", "CCN"),       # 重跑
        ("Et2NH", "CCNCC"),     # 新增
        ("NEt3",  "CCN(CC)CC")  # 新增 (重头戏)
    ]
    
    print(f">>> 启动 Et 系列补完计划 (精细扫描模式)...")

    for name, smiles in tasks:
        print(f"\n{'='*40}")
        print(f"正在处理: {name}")
        print(f"{'='*40}")
        
        scan_file = f"{name}_scan.gjf"
        ts_file = f"{name}_ts.gjf"
        ts_log = f"{name}_ts.log"

        try:
            # 1. 构建精细对齐体系
            mol_sys = build_aligned_system_fine(smiles)
            
            # 找原子序号
            conf = mol_sys.GetConformer()
            n_idx, c_idx = -1, -1
            for atom in mol_sys.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                if atom.GetSymbol() == 'N' and abs(pos.z) < 0.2: n_idx = atom.GetIdx() + 1
                # 这里的 C 是 MeI 的 C，Z 坐标应该在 2.4 左右
                if atom.GetSymbol() == 'C' and pos.z > 1.5: c_idx = atom.GetIdx() + 1
            
            # 2. 精细扫描 (从 2.4 扫到 1.6，步长 0.04)
            print("🚀 启动 Scan...")
            scan_lines = [f"{n_idx} {c_idx} S 20 -0.04"]
            
            create_advanced_scan_gaussian_input(
                mol=mol_sys, mol_name=f"{name}_Scan", 
                scan_lines=scan_lines, charge=0,
                nproc=N_PROC, mem=MEM_SIZE,
                method="#PM7 opt(modredundant)", 
                filename=scan_file
            )
            
            run_gaussian_job(scan_file)
            
            # 3. 提取最高点
            points = read_scan_output(f"{name}_scan.log")
            if points:
                max_point = max(points, key=lambda x: x['energy'])
                print(f"📊 锁定最高点: Step {max_point['step']}, E={max_point['energy']:.5f}")
                
                # 4. 提交 TS 优化
                print("🚀 启动 TS Optimization...")
                ts_guess, _ = extract_scan_point_to_mol(points, max_point['step'], None, show_3d=False)
                
                create_gaussian_input_advanced(
                    mol=ts_guess, mol_name=f"{name}_TS", 
                    # scan_lines=[], # 记得这里不需要 scan_lines
                    charge=0,
                    nproc=N_PROC, mem=MEM_SIZE,
                    method="#PM7 opt(ts,calcall,noeig) freq", 
                    filename=ts_file
                )
                
                run_gaussian_job(ts_file)
                
                # 5. 验证
                res = read_gaussian16_output_opt(ts_log)
                freqs = res.get('frequencies', [])
                if freqs and freqs[0] < 0:
                     print(f"✅ 成功! 虚频: {freqs[0]:.2f} cm^-1")
                else:
                     print(f"❌ 无虚频。可能原因：位阻太大导致无法形成有效 TS (这本身也是 NEt3 的结论!)")
            else:
                print("❌ 扫描失败，无数据。")
                
        except Exception as e:
            print(f"❌ 发生错误: {e}")
            import traceback
            traceback.print_exc()

    print("\n>>> 所有 Et 系列任务完成！")