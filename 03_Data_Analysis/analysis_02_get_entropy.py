import pandas as pd
import os
import re
from rdkit import Chem
from rdkit.Chem import Descriptors, GraphDescriptors, Lipinski

# ================= 配置 =================
filename_map = {
    "Et3N": "NEt3", 
}
molecules_smiles = {
    "NH3": "N",
    "MeNH2": "CN",
    "Me2NH": "CNC",
    "Me3N": "CN(C)C",
    "EtNH2": "CCN",
    "Et2NH": "CCNCC",
    "Et3N": "CCN(CC)CC", 
    "NEt3": "CCN(CC)CC"  
}

def get_entropy_from_log(logfile):
    """从 Log 文件提取总熵 (Total Entropy) - 柔性的硬核指标"""
    if not os.path.exists(logfile): return None
    
    with open(logfile, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
        
    # Gaussian 格式: 
    # Total                      86.543    (cal/mol.K)
    # 或者寻找 "Cal/Mol-Kelvin" 下面的 Total
    
    for i, line in enumerate(lines):
        if "Total" in line and "Cal/Mol-Kelvin" not in line:
            # 检查上一行或上几行是否有 Entropy 标识
            # 简单粗暴法：找包含 Total 且后面跟着数字的行，通常在 thermochemistry 部分
            # 更稳健的方法：
            pass
            
    # 从后往前找 "Total                      " 
    # 位于 "Molar capacity at constant volume" 附近
    
    entropy = None
    # 倒序搜索更靠近结尾的 Freq 分析
    for line in reversed(lines):
        if "Total" in line and len(line.split()) >= 2:
            try:
                # 尝试提取行尾的数值
                parts = line.split()
                # 典型格式: Total                      75.234
                # 或者是 CV, S 的表格
                val = float(parts[-1])
                # 熵通常在 40 - 150 之间
                if 30 < val < 200:
                    entropy = val
                    break
            except:
                pass
    
    # 如果没找到，尝试正则匹配更严格的格式
    if entropy is None:
        content = "".join(lines)
        # 匹配像 "Total                      86.543" 这样的行，前面可能有 S 或 Entropy 字眼
        # Gaussian 输出里通常是三列: E (Thermal), CV, S
        # 我们找 S 那一列 (最后一列)
        match = re.search(r'Total\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)', content)
        if match:
            entropy = float(match.group(1))
            
    return entropy

def find_correct_log(name):
    target = filename_map.get(name, name)
    candidates = [
        f"{target}_Ultimate.log",
        f"{target}_Smart_ts.log",
        f"{target}_ts.log",
        f"{target}_Fine_ts.log",
        f"{name}_ts.log"
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
print(">>> 提取柔性特征 (Flexibility: Entropy & Topology)...")

for index, row in df_ea.iterrows():
    name = row['Name']
    
    # 兼容 Ea
    if 'Ea_kcal' in row: ea = row['Ea_kcal']
    elif 'Ea_kJ' in row: ea = row['Ea_kJ'] / 4.184
    else: ea = 0
    
    smiles_key = filename_map.get(name, name)
    smiles = molecules_smiles.get(smiles_key)
    
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        # 1. 拓扑柔性
        rot_bonds = Lipinski.NumRotatableBonds(mol)
        kappa3 = GraphDescriptors.Kappa3(mol) # 形状指数
        hall_kier = GraphDescriptors.HallKierAlpha(mol) # 长链指数
        
        # 2. 物理柔性 (熵)
        found_log = find_correct_log(name)
        entropy = get_entropy_from_log(found_log) if found_log else None
        
        # 3. 之前的电子/位阻特征 (为了 SISSO 方便，重新算一下或者Merge)
        # 这里简单算几个核心的
        chi1v = GraphDescriptors.Chi1v(mol)
        
        if entropy is not None:
            print(f"  ✅ {name:<6} | S(Entropy): {entropy:.2f} | RotBonds: {rot_bonds} | Kappa3: {kappa3:.2f}")
            
            data.append({
                "Name": name,
                "Ea": ea,
                "S_Total": entropy,       # 核心柔性指标
                "Num_Rotors": rot_bonds,  # 简单柔性指标
                "HallKier": hall_kier,    # 拓扑形状
                "Kappa3": kappa3,         # 形状(球vs线)
                "Chi1v": chi1v            # 之前的位阻
            })
        else:
             print(f"  ⚠️ {name}: 没找到熵数据 (检查 {found_log})")

df_flex = pd.DataFrame(data)

# 尝试合并之前的高级特征 (HOMO, Q_N 等)
if os.path.exists("SISSO_Input_Advanced.csv"):
    df_old = pd.read_csv("SISSO_Input_Advanced.csv")
    # 只取 HOMO, Q_N, Num_H, TPSA 这些
    cols_to_keep = [c for c in df_old.columns if c not in df_flex.columns and c != "Name" and c != "Ea"]
    df_merged = pd.concat([df_flex, df_old[cols_to_keep]], axis=1) # 简单假设顺序一致(通常一致)
    # 更安全的 Merge:
    df_merged = pd.merge(df_flex, df_old[["Name", "HOMO", "Q_N", "Num_H", "TPSA"]], on="Name", how="left")
else:
    df_merged = df_flex

df_merged.to_csv("SISSO_Input_Flexibility.csv", index=False)
print(f"\n>>> 结果已保存至 SISSO_Input_Flexibility.csv")