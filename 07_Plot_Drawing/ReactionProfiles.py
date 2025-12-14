import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

# ==============================================================================
# 🎛️ [微调参数区] - 请在这里修改你要画的分子和标签位置
# ==============================================================================

# 1. 选妃：你想画哪个分子？(必须是 CSV 第一列里有的名字)
# 选项: "NH3", "MeNH2", "Me2NH", "Me3N", "EtNH2", "Et2NH", "NEt3"
TARGET_MOLECULE = "NEt3"  # <--- 在这里修改名字

# 2. 标签位置微调 (Offset)
# 格式: (x_shift, y_shift)
# x_shift: 正数向右，负数向左
# y_shift: 正数向上，负数向下 (单位是能量 kJ/mol)

# "Reactants" 文字位置
OFFSET_REACTANT_TEXT = (0.0, -5.0) 

# "TS" 文字位置
OFFSET_TS_TEXT = (0.0, 5.0)       

# "Products" 文字位置
OFFSET_PRODUCT_TEXT = (0.0, -5.0) 

# Ea (红色) 相关
X_POS_EA_ARROW = 1.45             # 红色箭头的 X 轴位置 (TS峰左侧一点点)
OFFSET_EA_LABEL = (0.05, -0.5)     # Ea 数值文字相对于箭头的偏移

# dE (绿色) 相关
X_POS_DE_ARROW = 2.8              # 绿色箭头的 X 轴位置 (图的最右侧)
OFFSET_DE_LABEL = (-0.03, 0.05)     # dE 数值文字相对于箭头的偏移 (负数表示文字在箭头左边)

# ==============================================================================
# 下面的代码通常不需要动，除非你想改配色或线宽
# ==============================================================================

INPUT_CSV = "Final_Report_Automated.csv"
COLOR_CURVE = "#000000"   
COLOR_FILL = "#E8F4F8"    
COLOR_LEVEL = "#333333"   
COLOR_EA = "#D62728"      
COLOR_DE = "#2CA02C"      
LINE_WIDTH = 2.0          

def plot_single_perfect(name, ea_kcal, de_kcal):
    # 单位转换
    Ea = ea_kcal * 4.184
    dE = de_kcal * 4.184
    
    print(f"Plotting {name}: Ea={Ea:.1f} kJ, dE={dE:.1f} kJ")

    # --- 1. 构建 PCHIP 曲线 ---
    # 针对吸热/放热反应自动调整辅助点，确保曲线平滑
    if dE > Ea * 0.8: # 如果产物能量很高（接近TS），需要调整右侧下降弧度
        x_anchors = np.array([0.0, 0.6, 1.4, 1.5, 1.6, 2.4, 3.0])
        y_anchors = np.array([0.0, 0.0, Ea,  Ea,  Ea,  dE,  dE])
    else:
        # 标准情况
        x_anchors = np.array([0.0, 0.6, 1.45, 1.5, 1.55, 2.4, 3.0])
        y_anchors = np.array([0.0, 0.0, Ea,   Ea,  Ea,   dE,  dE])
    
    interpolator = PchipInterpolator(x_anchors, y_anchors)
    x_smooth = np.linspace(0, 3.0, 500)
    y_smooth = interpolator(x_smooth)

    # --- 2. 绘图 ---
    fig, ax = plt.subplots(figsize=(5, 4), dpi=150) # 预览时 dpi 低一点快一点
    
    # 曲线与填充
    ax.plot(x_smooth, y_smooth, color=COLOR_CURVE, linewidth=LINE_WIDTH, zorder=3)
    ax.fill_between(x_smooth, 0, y_smooth, color=COLOR_FILL, alpha=0.6, zorder=1)
    
    # 能级平台
    ax.hlines(0, 0.0, 0.7, colors=COLOR_LEVEL, linewidth=1.5, zorder=4)
    ax.hlines(Ea, 1.3, 1.7, colors=COLOR_LEVEL, linewidth=1.5, zorder=4)
    ax.hlines(dE, 2.3, 3.0, colors=COLOR_LEVEL, linewidth=1.5, zorder=4)
    
    # --- 3. 标签 (应用微调参数) ---
    
    # Reactants
    ax.text(0.35 + OFFSET_REACTANT_TEXT[0], 
            10 + OFFSET_REACTANT_TEXT[1], 
            "Reactants", ha='center', va='top', fontsize=10, family='sans-serif')
    
    # TS
    ax.text(1.5 + OFFSET_TS_TEXT[0], 
            Ea - 5 + OFFSET_TS_TEXT[1], 
            "TS", ha='center', va='bottom', fontsize=10) # , fontweight='bold')

    # Products
    ax.text(2.65 + OFFSET_PRODUCT_TEXT[0], 
            dE + 3 + OFFSET_PRODUCT_TEXT[1], 
            "Products", ha='center', va='top', fontsize=10, family='sans-serif')

    # --- 4. 箭头与数值 ---
    
    # Ea 箭头
    ax.annotate('', xy=(X_POS_EA_ARROW, Ea), xytext=(X_POS_EA_ARROW, 0),
                arrowprops=dict(arrowstyle='<->', color=COLOR_EA, lw=1.2))
    # Ea 数值
    ax.text(X_POS_EA_ARROW + OFFSET_EA_LABEL[0], 
            Ea/2 + OFFSET_EA_LABEL[1], 
            f"$E_a={Ea:.1f}$", color=COLOR_EA, ha='left', va='center', fontsize=9,
            bbox=dict(boxstyle="square,pad=0.0", fc="white", ec="none", alpha=0.7))

    # dE 箭头
    ax.annotate('', xy=(X_POS_DE_ARROW, dE), xytext=(X_POS_DE_ARROW, 0),
                arrowprops=dict(arrowstyle='<->', color=COLOR_DE, lw=1.2))
    # dE 数值
    # 智能判断：如果文字偏移是正数，放右边；负数放左边
    align_h = 'left' if OFFSET_DE_LABEL[0] > 0 else 'right'
    ax.text(X_POS_DE_ARROW + OFFSET_DE_LABEL[0], 
            dE/2 + OFFSET_DE_LABEL[1], 
            f"$\\Delta E={dE:.1f}$", color=COLOR_DE, ha=align_h, va='center', fontsize=9,
            bbox=dict(boxstyle="square,pad=0.0", fc="white", ec="none", alpha=0.7))

    # 装饰
    ax.set_ylabel("Energy (kJ/mol)", fontsize=11)
    ax.set_xticks([])
    ax.set_xlabel("Reaction Coordinate", fontsize=11)
    ax.set_title(r"Reaction Profile: Et$_3$N + MeI", fontsize=12) #, fontweight='bold')
    
    # 自动缩放 Y 轴
    y_all = np.concatenate([y_anchors, [0, Ea, dE]])
    margin = (y_all.max() - y_all.min()) * 0.25
    ax.set_ylim(y_all.min() - margin, y_all.max() + margin)
    ax.set_xlim(0, 3.2)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    
    # ==========================================================================
    # 💾 [保存区] - 只有当你对 plt.show() 的结果满意时，取消下面两行的注释
    # ==========================================================================
    
    save_name = f"Final_Profile_{name}.png"
    plt.savefig(save_name, dpi=300)   # <--- 取消这个注释来保存
    print(f"✅ Saved: {save_name}")

    plt.show()  # <--- 先看图，满意再保存

# --- 主程序 ---
try:
    df = pd.read_csv(INPUT_CSV)
    row = df[df['Name'] == TARGET_MOLECULE]
    if not row.empty:
        # 兼容列名
        ea_val = row['Ea_kcal'].values[0] if 'Ea_kcal' in df.columns else row['Ea_kJ'].values[0]/4.184
        de_val = row['dE_kcal'].values[0] if 'dE_kcal' in df.columns else row['dE_kJ'].values[0]/4.184
        
        plot_single_perfect(TARGET_MOLECULE, ea_val, de_val)
    else:
        print(f"❌ 找不到分子: {TARGET_MOLECULE}，请检查名字是否写对。")
        print("可用名字:", df['Name'].values)
except Exception as e:
    print(f"❌ 读取 CSV 失败: {e}")