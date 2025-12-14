# batch_phase3_irc_fix.py
# Workflow 最终章：批量运行 IRC，验证反应路径 (修复版)

import os
import shutil
import glob
from rdkit import Chem

# 导入工具箱
try:
    from chem_utils import *
except ImportError:
    try:
        from chem_utils_new import *
    except ImportError:
        print("❌ 找不到 chem_utils.py")
        exit()

# 硬件配置
N_PROC = 6
MEM_SIZE = "4GB"

# 要跑 IRC 的分子列表
TARGETS = [
    "MeNH2", "Me2NH", "Me3N",
    "EtNH2", "Et2NH", "NEt3"
    # "NH3" 
]

# --- 核心辅助函数 ---
def find_best_ts_log(mol_name):
    """
    智能寻找最好的 TS log 文件
    优先找 _Fine_ts.log，其次找 _ts.log
    """
    # 1. 尝试 Fine 版本 (修复版)
    fine_log = f"{mol_name}_Fine_ts.log"
    if os.path.exists(fine_log):
        return fine_log
    
    # 2. 尝试普通版本
    normal_log = f"{mol_name}_ts.log"
    if os.path.exists(normal_log):
        return normal_log
        
    return None

def check_irc_done(irc_log):
    """检查 IRC 是否已经跑完"""
    if not os.path.exists(irc_log): return False
    try:
        with open(irc_log, 'rb') as f:
            content = f.read().decode('utf-8', errors='ignore')
            if "Normal termination" in content:
                return True
    except: pass
    return False

# --- 主程序 ---
if __name__ == "__main__":
    if not shutil.which("g16"):
        print("❌ 找不到 g16")
        exit()

    print(">>> 启动 Phase 3: 批量 IRC 验证 (修复版)...")
    
    for name in TARGETS:
        print(f"\n{'-'*15} 处理分子: {name} {'-'*15}")
        
        # 1. 寻找 TS 源文件
        ts_log = find_best_ts_log(name)
        if not ts_log:
            print(f"⚠️ 跳过: 找不到 {name} 的 TS 结果文件。")
            continue
            
        print(f"  📖 读取 TS 结构: {ts_log}")
        
        # 2. 检查 IRC 是否已完成
        base_name = ts_log.replace("_ts.log", "")
        irc_file = f"{base_name}_irc.gjf"
        irc_log = f"{base_name}_irc.log"
        
        if check_irc_done(irc_log):
            print(f"  ⏩ {irc_log} 已完成，跳过。")
            continue
            
        # 3. 提取 TS 结构
        try:
            res = read_gaussian16_output_opt(ts_log)
            if not res or not res.get('geometry'):
                print("  ❌ 无法从 log 提取几何结构")
                continue
                
            ts_mol = create_mol_from_geometry(res['geometry'])
            if not ts_mol:
                print("  ❌ 结构转换失败")
                continue
            
            # 4. 生成 IRC 输入文件
            print(f"  🚀 提交 IRC 计算...")
            
            # 【修复点】这里删除了 scan_lines=[]
            create_gaussian_input_advanced(
                mol=ts_mol, 
                mol_name=f"{base_name}_IRC", 
                # scan_lines=[],  <-- 罪魁祸首已删除
                charge=0,
                nproc=N_PROC, 
                mem=MEM_SIZE,
                method="#PM7 IRC(CalcFC, MaxPoints=20, StepSize=5)", 
                filename=irc_file
            )
            
            # 5. 运行
            run_gaussian_job(irc_file)
            
            if check_irc_done(irc_log):
                print(f"  ✅ IRC 成功完成！")
            else:
                print(f"  ⚠️ IRC 可能未完全跑完 (请检查 log)")
                
        except Exception as e:
            print(f"  ❌ 出错: {e}")
            import traceback
            traceback.print_exc()

    print("\n>>> Phase 3 全部结束！所有证据链已闭环！")