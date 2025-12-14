# batch_phase2_ts.py
# Workflow 第二阶段：读取扫描结果，进行 TS 精修和验证

import os
import shutil
import glob
from rdkit import Chem

# 导入工具箱 (确保目录下有 chem_utils.py 或 chem_utils_new.py)
try:
    from chem_utils import * 
except ImportError:
    try:
        from chem_utils_new import *
    except:
        print("❌ 找不到 chem_utils.py")
        exit()

# 硬件配置
N_PROC = 6
MEM_SIZE = "4GB"

# --- 核心辅助函数 ---
def check_normal_termination(logfile):
    """检查 Log 是否正常结束"""
    if not os.path.exists(logfile): return False
    try:
        with open(logfile, 'rb') as f:
            f.seek(-2048, os.SEEK_END)
            if "Normal termination" in f.read().decode('utf-8', errors='ignore'):
                return True
    except: pass
    return False

def get_scan_max_structure(scan_log, original_mol):
    """从扫描 Log 中提取能量最高的结构"""
    points = read_scan_output(scan_log)
    if not points: return None
    
    # 找最高点
    max_point = max(points, key=lambda x: x['energy'])
    print(f"    -> 锁定扫描最高点: Step {max_point['step']}, Energy={max_point['energy']:.5f}")
    
    # 提取结构 (利用 chem_utils 里的函数)
    # 注意：这里我们不需要显示 3D 图，show_3d=False
    ts_guess_mol, _ = extract_scan_point_to_mol(points, max_point['step'], original_mol, show_3d=False)
    return ts_guess_mol

# --- 主程序 ---
if __name__ == "__main__":
    if not shutil.which("g16"):
        print("❌ 找不到 g16")
        exit()

    print(">>> 启动 Phase 2: TS 精修与验证...")
    
    # 1. 自动搜索当前文件夹下所有的 "_scan.log" 文件
    # 这样你不需要手动改列表，Phase 1 跑成几个，这里就接着算几个
    scan_logs = glob.glob("*_scan.log")
    
    if not scan_logs:
        print("❌ 未找到任何扫描日志 (*_scan.log)。请先运行 Phase 1。")
        exit()

    results_csv = "ts_results.csv"
    if not os.path.exists(results_csv):
        with open(results_csv, "w") as f:
            f.write("Name,TS_Energy,Imag_Freq,Barrier_kcal_mol,Status\n")

    for scan_log in scan_logs:
        # 从文件名提取分子名 (例如 "MeNH2_scan.log" -> "MeNH2")
        name = scan_log.replace("_scan.log", "")
        print(f"\n{'-'*15} 处理分子: {name} {'-'*15}")
        
        ts_file = f"{name}_ts.gjf"
        ts_log  = f"{name}_ts.log"
        
        # --- 步骤 A: 检查是否已经算完 ---
        if check_normal_termination(ts_log):
            print(f"⏩ {name} TS 优化已完成，跳过。")
            continue
            
        # --- 步骤 B: 提取初猜结构 ---
        # 我们需要重新构建一下分子对象来作为模板 (为了保留原子序数等信息)
        # 简单粗暴的方法：读取 scan log 里的 SMILES 或者直接读取 scan log 的第一帧
        # 这里为了稳健，我们尝试从 scan log 提取
        
        try:
            # 这里的 original_mol 设为 None，让提取函数自己根据坐标创建新分子
            # 只要原子对上了就行
            ts_guess = get_scan_max_structure(scan_log, None)
            
            if not ts_guess:
                print(f"❌ 无法从 {scan_log} 提取结构，跳过。")
                continue
                
            # --- 步骤 C: 生成 TS 任务 ---
            print(f"🚀 提交 TS 优化任务...")
            
            # 关键：TS 优化必须用 calcall (算力常数) 保证准确，用 noeig 减少输出
            create_gaussian_input_advanced(
                ts_guess, 
                f"{name}_TS", 
                charge=0,
                nproc=N_PROC, 
                mem=MEM_SIZE,
                method="#PM7 opt(ts,calcall,noeig) freq", # 加 freq 确认虚频
                filename=ts_file
            )
            
            # --- 步骤 D: 运行 ---
            run_gaussian_job(ts_file)
            
            # --- 步骤 E: 验证虚频 (Quality Control) ---
            # 读取结果
            res = read_gaussian16_output_opt(ts_log)
            freqs = res.get('frequencies', [])
            energy = res.get('energy', 0)
            
            if freqs and freqs[0] < 0:
                img_freq = freqs[0]
                status = "Success"
                if len(freqs) > 1 and freqs[1] < 0:
                    status = "Warning_Multiple_Imag" # 多个虚频，可能不对
                
                print(f"✅ TS 验证通过! 虚频: {img_freq:.2f} cm^-1")
                
                # 记录到 CSV
                with open(results_csv, "a") as f:
                    f.write(f"{name},{energy:.6f},{img_freq:.2f},CALC_LATER,{status}\n")
            else:
                print(f"❌ TS 验证失败: 没有虚频或优化未收敛。Freqs: {freqs}")
                
        except Exception as e:
            print(f"❌ 处理出错: {e}")

    print("\n>>> Phase 2 全部结束！Workflow 闭环！")