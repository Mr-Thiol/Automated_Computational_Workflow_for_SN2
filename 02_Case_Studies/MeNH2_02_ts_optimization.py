# finish_MeNH2_ts.py
# 专门用于完成 MeNH2 的最后一步 TS 优化

import os
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
name = "MeNH2_Fine"

print(f">>> 继续完成 MeNH2 的 TS 优化...")

# 1. 读取刚刚跑完的精细扫描日志
scan_log = f"{name}_scan.log"
if not os.path.exists(scan_log):
    print(f"❌ 找不到 {scan_log}，请确认上一波扫描是否生成了日志。")
    exit()

points = read_scan_output(scan_log)

if points:
    # 2. 再次锁定最高点 (Step 15)
    max_point = max(points, key=lambda x: x['energy'])
    print(f"📊 重新锁定最高点: Step {max_point['step']}, Energy={max_point['energy']:.5f}")
    
    # 3. 提取结构
    # 注意：这里 original_mol 传 None，让函数自动创建
    ts_guess, _ = extract_scan_point_to_mol(points, max_point['step'], None, show_3d=False)
    
    # 4. 提交 TS 优化 (修复了报错的地方！)
    ts_file = f"{name}_ts.gjf"
    ts_log = f"{name}_ts.log"
    
    print(f"🚀 启动 TS 优化 (PM7 calcall)...")
    
    # 【关键修复】删除了 scan_lines=[] 参数
    create_gaussian_input_advanced(
        mol=ts_guess, 
        mol_name=f"{name}_TS", 
        # scan_lines=[],  <-- 之前就是这里多了这一行，删掉就好了
        charge=0,
        nproc=N_PROC, 
        mem=MEM_SIZE,
        method="#PM7 opt(ts,calcall,noeig) freq", 
        filename=ts_file
    )
    
    # 5. 运行计算
    run_gaussian_job(ts_file)
    
    # 6. 验证结果
    res = read_gaussian16_output_opt(ts_log)
    freqs = res.get('frequencies', [])
    
    if freqs and freqs[0] < 0:
        print(f"✅✅✅ 最终胜利！MeNH2 救回来了！虚频: {freqs[0]:.2f} cm^-1")
    else:
        print(f"❌ 依然没有虚频。Freqs: {freqs}")
        print("建议：不用纠结了，直接用 Me2NH 做 Pre。")

else:
    print("❌ 读取扫描日志失败")