import os
import glob
import re
import pandas as pd

# ================= 配置 =================
# 输出文件名
OUTPUT_CSV = "Final_Report_Automated.csv"
# 能量转换因子 (Hartree -> kJ/mol)
H_TO_KJ = 2625.4996
# =======================================

def extract_energy_and_freq(logfile):
    """
    从 TS log 中提取：
    1. 最终能量 (SCF Done)
    2. 虚频 (Frequencies)
    """
    if not os.path.exists(logfile): return None, None
    
    energy = None
    freqs = []
    
    try:
        with open(logfile, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
            
        # 倒序查找最后的 SCF Energy
        for line in reversed(lines):
            if "SCF Done:" in line:
                # 匹配浮点数
                match = re.search(r'SCF Done:.*=\s*([-\d\.]+(?:[eE][-+]?\d+)?)', line)
                if match:
                    energy = float(match.group(1))
                    break
        
        # 正序查找频率 (Frequencies --)
        # 只要有一个负的就行，取第一个
        for line in lines:
            if "Frequencies --" in line:
                # 提取行内所有数字
                nums = re.findall(r'[-\d\.]+', line)
                # 过滤掉非数字字符并转换
                current_freqs = [float(x) for x in nums if x != '-' and x != '--']
                freqs.extend(current_freqs)
                # 如果找到了虚频（负数），通常在第一行
                if freqs and freqs[0] < 0:
                    break
                    
    except Exception as e:
        print(f"  ⚠️ 解析出错 {logfile}: {e}")
        
    return energy, freqs

def extract_scan_endpoints(logfile):
    """
    从 Scan log 中提取：
    1. 起点能量 (Step 1) -> 作为 Reactant
    2. 终点能量 (Last Step) -> 作为 Product
    """
    if not os.path.exists(logfile): return None, None
    
    energies = []
    try:
        with open(logfile, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if "SCF Done:" in line:
                    match = re.search(r'SCF Done:.*=\s*([-\d\.]+(?:[eE][-+]?\d+)?)', line)
                    if match:
                        energies.append(float(match.group(1)))
    except:
        pass
        
    if not energies: return None, None
    
    return energies[0], energies[-1]

# ================= 主程序 =================

print(f">>> 启动全自动数据整理脚本...")

# 1. 找到所有 TS 文件 (以此为核心)
ts_files = glob.glob("*_ts.log")
print(f"发现 {len(ts_files)} 个过渡态文件: {ts_files}")

data_rows = []

for ts_log in ts_files:
    # 提取基础名字
    # 例如: "MeNH2_Fine_ts.log" -> base_name="MeNH2_Fine"
    # 例如: "NH3_ts.log" -> base_name="NH3"
    base_name = ts_log.replace("_ts.log", "")
    print(f"\n处理分子: {base_name}")
    
    # A. 提取 TS 数据
    ts_energy, freqs = extract_energy_and_freq(ts_log)
    
    if ts_energy is None:
        print("  ❌ TS 能量读取失败")
        continue
        
    imag_freq = freqs[0] if (freqs and freqs[0] < 0) else 0.0
    print(f"  ✅ TS能量: {ts_energy:.5f} | 虚频: {imag_freq:.2f}")
    
    # B. 寻找对应的 Scan 文件以获取 Reactant/Product 能量
    # 策略：优先找同名的 _scan.log
    scan_log = f"{base_name}_scan.log"
    
    # 如果没找到，且名字里有 Fine，尝试找去掉 Fine 的粗扫文件 (为了更准的反应物基准)
    # 比如 MeNH2_Fine_ts.log，优先找 MeNH2_Fine_scan.log
    # 如果没找到，再试 MeNH2_scan.log (看你具体情况，这里先默认找同名)
    
    reactant_e = None
    product_e = None
    scan_source = "None"
    
    if os.path.exists(scan_log):
        reactant_e, product_e = extract_scan_endpoints(scan_log)
        scan_source = scan_log
    else:
        # 尝试去掉 _Fine 的情况 (针对 MeNH2 这种修正过的)
        if "_Fine" in base_name:
            simple_name = base_name.replace("_Fine", "")
            simple_scan = f"{simple_name}_scan.log"
            if os.path.exists(simple_scan):
                # 注意：如果用粗扫文件，可能需要确认是否兼容
                # 这里为了自动化，先尝试读取
                r_e, p_e = extract_scan_endpoints(simple_scan)
                if reactant_e is None: reactant_e = r_e # 优先取反应物
                # Product 还是建议用 Fine scan 的结尾，如果 Fine scan 没找到，那就没办法了
                scan_source = f"{simple_scan} (Reactant only)"
    
    # C. 计算 Ea 和 dE
    row = {
        "Name": base_name,
        "TS_Energy": ts_energy,
        "Imag_Freq": imag_freq,
        "Reactant_E": reactant_e,
        "Product_E": product_e,
        "Ea_kcal": None,
        "dE_kcal": None,
        "Note": ""
    }
    
    if reactant_e is not None:
        ea = (ts_energy - reactant_e) * H_TO_KJ
        row["Ea_kcal"] = round(ea, 2)
        print(f"  📊 活化能 Ea: {ea:.2f} kcal/mol")
    else:
        row["Note"] += "No_Scan_Found "
        print("  ⚠️ 未找到扫描文件，无法计算 Ea")
        
    if reactant_e is not None and product_e is not None:
        de = (product_e - reactant_e) * H_TO_KJ
        row["dE_kcal"] = round(de, 2)
    
    data_rows.append(row)

# 3. 汇总保存
df = pd.DataFrame(data_rows)

# 调整列顺序
cols = ["Name", "Ea_kcal", "dE_kcal", "Imag_Freq", "TS_Energy", "Reactant_E", "Product_E", "Note"]
df = df[cols]

# 排序 (可选)
df = df.sort_values(by="Name")

print("\n" + "="*40)
print("最终报表预览")
print("="*40)
print(df)

df.to_csv(OUTPUT_CSV, index=False)
print(f"\n>>> ✅ 结果已保存至: {OUTPUT_CSV}")