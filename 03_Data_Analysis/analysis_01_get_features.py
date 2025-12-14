import pandas as pd
import numpy as np
import os
import math
from rdkit import Chem, Geometry

# ================= 配置 =================
# 文件名映射
filename_map = {
    "Et3N": "NEt3", 
}
# =======================================

def get_angle(c0, c1, c2):
    """计算三个点 c1-c0-c2 的夹角 (c0是中心)"""
    v1 = np.array(c1) - np.array(c0)
    v2 = np.array(c2) - np.array(c0)
    
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    
    if norm1 < 1e-6 or norm2 < 1e-6: return 0.0
    
    dot = np.dot(v1, v2)
    cos_theta = dot / (norm1 * norm2)
    # 防止数值误差导致超出范围
    cos_theta = max(min(cos_theta, 1.0), -1.0)
    
    return math.degrees(math.acos(cos_theta))

def extract_reactant_geometry(logfile):
    """
    从 Scan Log 提取第一步(反应物)的几何结构
    返回一个原子列表 [{'symbol': 'N', 'coords': [x,y,z]}, ...]
    """
    if not os.path.exists(logfile): return None
    
    atoms = []
    reading_geom = False
    
    with open(logfile, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    # Gaussian 16 提取 Input Orientation 或 Standard Orientation
    # 我们找第一步的 Standard Orientation
    found_first = False
    
    start_idx = -1
    for i, line in enumerate(lines):
        if "Standard orientation:" in line:
            start_idx = i + 5 # 跳过表头
            found_first = True
            break
            
    if not found_first: return None
    
    # 读取坐标
    for i in range(start_idx, len(lines)):
        line = lines[i]
        if "----------------" in line: break # 结束
        parts = line.split()
        if len(parts) >= 6:
            # Center AtomicNumber ... X Y Z
            # AtomicNumber 转 Symbol 需要个表，或者我们直接通过原子序数判断
            # 6=C, 1=H, 7=N, 53=I
            atom_num = int(parts[1])
            x = float(parts[3])
            y = float(parts[4])
            z = float(parts[5])
            
            symbol = "X"
            if atom_num == 1: symbol = 'H'
            elif atom_num == 6: symbol = 'C'
            elif atom_num == 7: symbol = 'N'
            elif atom_num == 53: symbol = 'I'
            
            atoms.append({'symbol': symbol, 'coords': [x, y, z]})
            
    return atoms

def calculate_n_geometry(atoms):
    """计算 N 原子的几何参数"""
    # 1. 找 N 原子
    n_atom = None
    n_idx = -1
    for i, a in enumerate(atoms):
        if a['symbol'] == 'N':
            n_atom = a
            n_idx = i
            break
            
    if n_atom is None: return None
    
    # 2. 找 N 的邻居 (在 1.6 埃以内的算成键)
    neighbors = []
    n_coords = n_atom['coords']
    
    for i, a in enumerate(atoms):
        if i == n_idx: continue
        # 不计算 I (离得很远)
        if a['symbol'] == 'I': continue
        
        dist = np.linalg.norm(np.array(a['coords']) - np.array(n_coords))
        if dist < 1.65: # C-N 约 1.47, N-H 约 1.01
            neighbors.append(a)
            
    # 3. 计算参数
    if len(neighbors) < 2: 
        return None # 至少要 2 个邻居才能算角度 (虽然伯胺有2个H)
    
    angles = []
    bond_lengths = []
    
    for k in range(len(neighbors)):
        # 键长
        d = np.linalg.norm(np.array(neighbors[k]['coords']) - np.array(n_coords))
        bond_lengths.append(d)
        
        # 键角
        for m in range(k+1, len(neighbors)):
            ang = get_angle(n_coords, neighbors[k]['coords'], neighbors[m]['coords'])
            angles.append(ang)
            
    if not angles: return None
    
    return {
        "Sum_Angles": sum(angles),
        "Avg_Angle": sum(angles) / len(angles),
        "Max_Angle": max(angles),
        "Avg_Bond_Len": sum(bond_lengths) / len(bond_lengths),
        "N_Coordination": len(neighbors) # 配位数
    }

def find_correct_log(base_name):
    """找 scan log"""
    target_name = filename_map.get(base_name, base_name)
    candidates = [
        f"{target_name}_Fine_scan.log",
        f"{target_name}_scan.log",
        f"{base_name}_Fine_scan.log",
        f"{base_name}_scan.log",
        f"{target_name}_Smart_scan.log",
        f"{base_name}_Smart_scan.log"
    ]
    for f in candidates:
        if os.path.exists(f): return f
    return None

# --- 主程序 ---
try:
    df_ea = pd.read_csv("Final_Report_Automated.csv")
except:
    print("❌ 找不到 CSV")
    exit()

data = []
print(">>> 提取 N 中心几何特征 (Geometric Features)...")

for index, row in df_ea.iterrows():
    name = row['Name']
    
    # 兼容 Ea
    if 'Ea_kcal' in row: ea = row['Ea_kcal']
    elif 'Ea_kJ' in row: ea = row['Ea_kJ'] / 4.184
    else: ea = 0
    
    found_log = find_correct_log(name)
    
    if found_log:
        atoms = extract_reactant_geometry(found_log)
        if atoms:
            geo_feats = calculate_n_geometry(atoms)
            if geo_feats:
                print(f"  ✅ {name:<6} | Sum_Angles: {geo_feats['Sum_Angles']:.1f}° | Avg_Len: {geo_feats['Avg_Bond_Len']:.3f} Å")
                
                # 读取之前的 HOMO 等数据 (从 SISSO_Input_Advanced.csv 合并更好，这里简单起见只生成几何的)
                # 我们建议你把这个 merge 进去
                
                row_data = {
                    "Name": name,
                    "Ea": ea,
                }
                row_data.update(geo_feats)
                data.append(row_data)
            else:
                print(f"  ⚠️ {name}: 几何分析失败 (邻居太少?)")
        else:
             print(f"  ⚠️ {name}: 无法从 {found_log} 提取结构")
    else:
        print(f"  ❌ {name}: 找不到 log")

df_geo = pd.DataFrame(data)

# 如果存在之前的 Advanced Input，合并之
if os.path.exists("SISSO_Input_Advanced.csv"):
    df_old = pd.read_csv("SISSO_Input_Advanced.csv")
    # Merge
    df_final = pd.merge(df_old, df_geo, on=["Name", "Ea"], how="inner")
else:
    df_final = df_geo

df_final.to_csv("SISSO_Input_Geometric.csv", index=False)
print(f"\n>>> 结果已保存至 SISSO_Input_Geometric.csv")
print("    现在你可以用这些 '硬核' 的几何参数去跑 LOO 了！")