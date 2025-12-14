import os
import shutil
try:
    from chem_utils import *
except:
    try:
        from chem_utils_new import *
    except:
        print("❌ 找不到 chem_utils")
        exit()

N_PROC = 6
MEM = "4GB"
# 注意：这里要跟刚才 salvaging 脚本里生成的 TS 名字对应
# 刚才生成的是 Me2NH_Smart_ts.gjf -> Me2NH_Smart_ts.log
BASE_NAME = "Me2NH_Smart" 

print(f">>> 启动 Me2NH (Smart) IRC 计算...")

ts_log = f"{BASE_NAME}_ts.log"
if not os.path.exists(ts_log):
    print(f"❌ 找不到 {ts_log}，请确认文件名")
    exit()

print(f"  📖 读取 TS 结构: {ts_log}")

try:
    # 提取结构
    res = read_gaussian16_output_opt(ts_log)
    if not res or not res.get('geometry'):
        print("  ❌ 无法提取几何结构")
        exit()
        
    ts_mol = create_mol_from_geometry(res['geometry'])
    
    # 生成 IRC 输入文件
    irc_file = f"{BASE_NAME}_irc.gjf"
    
    print(f"  🚀 提交 IRC 计算...")
    create_gaussian_input_advanced(
        mol=ts_mol, 
        mol_name=f"{BASE_NAME}_IRC", 
        # scan_lines=[], # 确保没有 scan 参数
        charge=0,
        nproc=N_PROC, 
        mem=MEM,
        method="#PM7 IRC(CalcFC, MaxPoints=20, StepSize=5)", 
        filename=irc_file
    )
    
    # 运行
    run_gaussian_job(irc_file)
    print(f"  ✅ IRC 计算提交完成！请等待 Me2NH_Smart_irc.log 生成。")

except Exception as e:
    print(f"❌ 出错: {e}")