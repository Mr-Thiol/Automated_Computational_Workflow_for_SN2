import os
import sys

# 导入工具箱
try:
    from chem_utils import *
except:
    from chem_utils_new import *

N_PROC = 6
MEM = "4GB"
NAME = "Me2NH_Smart"
SCAN_LOG = f"{NAME}_scan.log"

print(f">>> 启动 Me2NH 抢救计划...")

# 1. 尝试读取“失败”的日志
if not os.path.exists(SCAN_LOG):
    print(f"❌ 找不到 {SCAN_LOG}，无法抢救。")
    sys.exit()

print(f"    正在解析 {SCAN_LOG}...")
# read_scan_output 通常比较健壮，能读多少读多少
points = read_scan_output(SCAN_LOG)

if not points:
    print("❌ 读取失败：日志里似乎没有包含任何完整的优化步骤。")
    print("   原因推测：可能第一步就崩了 (SCF不收敛)。")
    print("   建议：直接放弃 Me2NH，使用 SISSO 剔除方案。")
    sys.exit()

# 2. 寻找最高点
num_points = len(points)
max_point = max(points, key=lambda x: x['energy'])
max_step = max_point['step']

print(f"✅ 成功提取到 {num_points} 个扫描点！")
print(f"📊 锁定最高点: Step {max_step}, Energy={max_point['energy']:.5f}")

# 简单的判断：如果最高点不是第一点也不是最后一点，说明确实翻山了
if max_step > 1 and max_step < num_points:
    print("✨ 好消息！扫描成功翻越了能垒！TS 结构就在这里！")
elif max_step == 1:
    print("⚠️ 警告：最高点在起点。可能是还在爬坡就崩了，或者本来就是下坡。")
elif max_step == num_points:
    print("⚠️ 警告：最高点在终点。可能还没爬到顶就崩了。")

# 3. 提取结构并提交 TS 优化
print(f">>> 提取 Step {max_step} 结构，提交 TS 优化 (Calcall)...")

# 注意：extract_scan_point_to_mol 需要 original_mol 来推断连接，这里传 None 让它自动感知
ts_guess, _ = extract_scan_point_to_mol(points, max_step, None, show_3d=False)

ts_file = f"{NAME}_ts.gjf"
ts_log = f"{NAME}_ts.log"

create_gaussian_input_advanced(
    mol=ts_guess, 
    mol_name=f"{NAME}_TS",
    # scan_lines=[], 
    charge=0,
    nproc=N_PROC, 
    mem=MEM,
    # 关键：用 calcall 算准 Hessian，noeig 防止它乱跑
    method="#PM7 opt(ts,calcall,noeig) freq",
    filename=ts_file
)

# 4. 运行
run_gaussian_job(ts_file)

# 5. 最终开奖
res = read_gaussian16_output_opt(ts_log)
freqs = res.get('frequencies', [])
energy = res.get('energy', 0)

print("\n" + "="*40)
print("最终抢救结果")
print("="*40)

if freqs and freqs[0] < -100:
    print(f"🎉🎉🎉 奇迹发生了！Me2NH 救回来了！")
    print(f"虚频: {freqs[0]:.2f} cm^-1")
    print(f"能量: {energy:.5f} Hartree")
    print("快把这个数据填进 CSV，不需要剔除离群点了！")
else:
    print(f"❌ 还是不行。虚频: {freqs}")
    print("Plan B 启动：请在报告中剔除 Me2NH (SISSO R2=0.68)。")