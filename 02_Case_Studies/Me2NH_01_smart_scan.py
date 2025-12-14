import os
import sys
import numpy as np
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms

# 导入工具箱
try:
    from chem_utils import *
except:
    from chem_utils_new import *

N_PROC = 6
MEM = "4GB"
NAME = "Me2NH_Smart"

def get_lone_pair_vector(mol, n_idx):
    """计算 N 原子的孤对电子方向 (基于邻居原子的反向平均)"""
    conf = mol.GetConformer()
    n_pos = np.array(conf.GetAtomPosition(n_idx))
    
    atom = mol.GetAtomWithIdx(n_idx)
    neighbors = atom.GetNeighbors()
    
    vec_sum = np.array([0.0, 0.0, 0.0])
    for nbr in neighbors:
        p = np.array(conf.GetAtomPosition(nbr.GetIdx()))
        v = p - n_pos
        v = v / np.linalg.norm(v)
        vec_sum += v
        
    # 孤对电子方向大致是键向量之和的反方向
    lp_vec = -vec_sum
    lp_vec = lp_vec / np.linalg.norm(lp_vec)
    return lp_vec

def build_smart_system():
    print(">>> 1. 构建智能对齐的 Me2NH + MeI 体系...")
    
    # --- A. 准备 Me2NH ---
    amine = Chem.MolFromSmiles("CNC")
    amine = Chem.AddHs(amine)
    AllChem.EmbedMolecule(amine, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(amine)
    conf_amine = amine.GetConformer()
    
    # 找 N
    n_atom = [a for a in amine.GetAtoms() if a.GetSymbol() == 'N'][0]
    n_idx = n_atom.GetIdx()
    
    # 算孤对电子方向
    lp_vec = get_lone_pair_vector(amine, n_idx)
    
    # 旋转矩阵：把孤对电子转到 Z 轴正方向 (0,0,1)
    # 这样 N 指向 Z+，我们要把 MeI 放在 Z+ 处
    rot_mat = get_rotation_matrix_to_align_vector(
        np.array([0,0,0]), lp_vec, target_axis=np.array([0,0,1])
    )
    rdMolTransforms.TransformConformer(conf_amine, rot_mat)
    
    # 平移 N 到原点
    n_pos = conf_amine.GetAtomPosition(n_idx)
    shift = Geometry.Point3D(-n_pos.x, -n_pos.y, -n_pos.z)
    for i in range(amine.GetNumAtoms()):
        p = conf_amine.GetAtomPosition(i)
        conf_amine.SetAtomPosition(i, Geometry.Point3D(p.x+shift.x, p.y+shift.y, p.z+shift.z))
        
    # --- B. 准备 MeI ---
    mei = Chem.MolFromSmiles("CI")
    mei = Chem.AddHs(mei)
    AllChem.EmbedMolecule(mei, randomSeed=42)
    conf_mei = mei.GetConformer()
    
    c_idx = [a.GetIdx() for a in mei.GetAtoms() if a.GetSymbol() == 'C'][0]
    i_idx = [a.GetIdx() for a in mei.GetAtoms() if a.GetSymbol() == 'I'][0]
    
    # --- C. 摆放 (关键微扰) ---
    # 我们把 MeI 的 C 放在 Z 轴上，距离 2.3
    # 但是为了给一点微扰（避免完美直线），我们在 X 轴偏移一点点 (0.05)
    # 这样 N-C-I 角度大概是 178 度左右
    
    start_dist = 2.3
    perturbation = 0.08 # 偏移量
    
    # 摆放 C
    conf_mei.SetAtomPosition(c_idx, Geometry.Point3D(perturbation, 0, start_dist))
    
    # 摆放 I (C-I 键长 2.14)
    # I 也要顺着微扰方向摆，或者摆在轴上制造更明显的角度
    # 这里我们把 I 放在更远处，保持 C-I 也是竖直的，这样 N...C-I 就自然弯曲了
    conf_mei.SetAtomPosition(i_idx, Geometry.Point3D(perturbation, 0, start_dist + 2.14))
    
    # 摆放 H (调整方向避免撞车)
    h_idxs = [a.GetIdx() for a in mei.GetAtoms() if a.GetSymbol() == 'H']
    offsets = [(0, 1.0, -0.4), (0.86, -0.5, -0.4), (-0.86, -0.5, -0.4)]
    for k, h_id in enumerate(h_idxs):
        dx, dy, dz = offsets[k]
        conf_mei.SetAtomPosition(h_id, Geometry.Point3D(perturbation+dx, dy, start_dist+dz))

    # 组合
    combo = Chem.CombineMols(amine, mei)
    return combo

# --- 对齐向量辅助函数 (防止报错) ---
def get_rotation_matrix_to_align_vector(vec_start, vec_end, target_axis=np.array([0, 0, 1])):
    # 简单的实现，只用 vec_end 即可，因为 vec_start 通常是原点
    if isinstance(vec_end, Geometry.Point3D): vec_end = np.array([vec_end.x, vec_end.y, vec_end.z])
    
    v_current = vec_end
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

# --- 主程序 ---
if __name__ == "__main__":
    # 1. 构建
    mol = build_smart_system()
    
    # 找原子序号
    atoms = mol.GetAtoms()
    # 假设前部分是胺，后部分是 MeI (CombineMols 的顺序)
    # N 在原点附近
    conf = mol.GetConformer()
    n_idx = -1
    c_idx = -1
    
    for a in atoms:
        pos = conf.GetAtomPosition(a.GetIdx())
        if a.GetSymbol() == 'N' and abs(pos.z) < 0.5: n_idx = a.GetIdx() + 1
        if a.GetSymbol() == 'C' and pos.z > 1.0: c_idx = a.GetIdx() + 1
        
    print(f"    锁定原子: N={n_idx}, C(MeI)={c_idx}")
    
    # 2. 设置扫描
    # Start=2.3 (当前位置), End=1.6
    # Diff = 0.7
    # Step = 0.05
    # Steps = 14
    
    print(f">>> 2. 启动智能扫描 (2.3 -> 1.6 Å)...")
    scan_file = f"{NAME}_scan.gjf"
    
    # 关键：扫描 N-C 键
    scan_lines = [f"{n_idx} {c_idx} S 14 -0.05"]
    
    create_advanced_scan_gaussian_input(
        mol=mol, mol_name=f"{NAME}_Scan",
        scan_lines=scan_lines, charge=0,
        nproc=N_PROC, mem=MEM,
        method="#PM7 opt(modredundant)",
        filename=scan_file
    )
    
    # 运行
    if run_gaussian_job(scan_file):
        # 3. 提取与优化
        points = read_scan_output(f"{NAME}_scan.log")
        if points:
            max_point = max(points, key=lambda x: x['energy'])
            print(f"📊 锁定最高点: Step {max_point['step']}, Energy={max_point['energy']:.5f}")
            
            # 4. 提交最终 TS 优化
            print(f">>> 3. 提交 TS 优化 (CalcAll)...")
            ts_guess, _ = extract_scan_point_to_mol(points, max_point['step'], mol, show_3d=False)
            
            ts_file = f"{NAME}_ts.gjf"
            ts_log = f"{NAME}_ts.log"
            
            create_gaussian_input_advanced(
                mol=ts_guess, mol_name=f"{NAME}_TS",
                scan_lines=[], charge=0,
                nproc=N_PROC, mem=MEM,
                method="#PM7 opt(ts,calcall,noeig) freq",
                filename=ts_file
            )
            
            run_gaussian_job(ts_file)
            
            # 5. 验证
            res = read_gaussian16_output_opt(ts_log)
            freqs = res.get('frequencies', [])
            if freqs and freqs[0] < -200:
                print(f"✅✅✅ 最终胜利！Me2NH 虚频: {freqs[0]:.2f} cm^-1")
            else:
                print(f"❌ 依然困难。Freqs: {freqs}")