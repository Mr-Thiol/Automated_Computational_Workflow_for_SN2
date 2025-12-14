import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_absolute_error
import re  # 引入正则表达式库，用于处理化学式下标

# ================= 配置 =================
INPUT_FILE = "SISSO_Input_Flexibility.csv"
FEATURES = ['Q_N', 'Num_H', 'Kappa3']
TARGET = 'Ea'

# --- 标签位置微调配置 ---
# 依然使用原始名字作为 Key (例如 'Me3N')
# 格式：'原始名': (x偏移, y偏移)
LABEL_OFFSETS = {
    'Me3N':   (0, -15),
    'NEt3':   (10, -10),
    'NH3':    (10, 0),
    'Me2NH':  (-30, 10),
    'Et2NH': (0, -15),
    # 在这里根据生成的图片继续添加需要调整的分子...
}
DEFAULT_OFFSET = (5, 5) 
# =======================================

def format_chemical_formula(text):
    """
    将化学式字符串转换为 LaTeX 格式
    输入: "Me3N"
    输出: "$\mathrm{Me_{3}N}$"
    """
    # 1. 查找所有数字，替换为 _{数字}
    latex_text = re.sub(r'(\d+)', r'_{\1}', text)
    # 2. 包裹在 mathrm 中保证是正体，外层包裹 $ 触发 LaTeX
    return f"$\\mathrm{{{latex_text}}}$"

def plot_parity():
    print(f">>> 正在绘制最终拟合图 (Parity Plot)...")
    
    try:
        df = pd.read_csv(INPUT_FILE)
    except:
        print(f"❌ 找不到 {INPUT_FILE}，正在生成测试数据演示效果...")
        # 生成假数据演示
        data = {
            'Name': ['Me3N', 'Me2NH', 'MeNH2', 'Et3N', 'Et2NH', 'EtNH2', 'NH3', 'MeEtNH'],
            'Q_N': np.random.rand(8), 'Num_H': np.random.rand(8), 'Kappa3': np.random.rand(8),
            'Ea': np.random.uniform(10, 30, 8)
        }
        df = pd.DataFrame(data)

    # --- 1. 数据准备 ---
    df_clean = df.dropna(subset=FEATURES + [TARGET]).copy()
    X = df_clean[FEATURES].values
    y_true = df_clean[TARGET].values
    names = df_clean['Name'].values
    
    # --- 2. 训练模型 ---
    model = LinearRegression()
    model.fit(X, y_true)
    y_pred = model.predict(X)
    r2 = r2_score(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    
    # --- 3. 绘图设置 ---
    plt.figure(figsize=(7, 7), dpi=300)
    
    # 坐标范围
    all_vals = np.concatenate([y_true, y_pred])
    min_val = np.min(all_vals) - 3 
    max_val = np.max(all_vals) + 3
    plt.plot([min_val, max_val], [min_val, max_val], color='gray', linestyle='--', alpha=0.5, label='Perfect Fit')
    
    # --- 4. 分组绘制 (图例) ---
    def get_category(name):
        if "Me" in name and "N" in name: return "Methyl Series", '#1f77b4'
        if "Et" in name: return "Ethyl Series", '#d62728'
        return "Others", 'black'
    
    groups = {}
    for i, nm in enumerate(names):
        cat, color = get_category(nm)
        if cat not in groups:
            groups[cat] = {'x': [], 'y': [], 'c': color}
        groups[cat]['x'].append(y_true[i])
        groups[cat]['y'].append(y_pred[i])
        
    for cat, data in groups.items():
        plt.scatter(data['x'], data['y'], color=data['c'], label=cat, 
                    s=150, alpha=0.8, edgecolors='k', zorder=5)

    # --- 5. 标注文字 (含化学式格式化) ---
    for i, original_name in enumerate(names):
        # 转换名字格式：Me3N -> Me_3N
        display_text = format_chemical_formula(original_name)
        
        # 查找位置偏移 (使用原始名字查找 Key)
        if original_name in LABEL_OFFSETS:
            offset = LABEL_OFFSETS[original_name]
            arrow_props = dict(arrowstyle='-', color='gray', lw=0.5)
        else:
            offset = DEFAULT_OFFSET
            arrow_props = None

        plt.annotate(display_text, 
                     xy=(y_true[i], y_pred[i]), 
                     xytext=offset, 
                     textcoords='offset points', 
                     fontsize=11, # 稍微加大一点字体，因为LaTeX显示会显小
                     arrowprops=arrow_props,
                     zorder=10)

    # --- 6. 装饰 ---
    plt.legend(loc='lower right', frameon=True, fancybox=True, framealpha=0.9)
    
    # 标题去掉了 fontweight='bold'
    plt.title(f"Machine Learning Prediction of $E_a$ \n($R^2={r2:.3f}$, MAE={mae:.2f} kcal/mol)", fontsize=14)
    
    plt.xlabel("Calculated $E_a$ kcal/mol", fontsize=12)
    plt.ylabel("Predicted $E_a$ kcal/mol", fontsize=12)
    plt.xlim(min_val, max_val)
    plt.ylim(min_val, max_val)
    plt.grid(True, linestyle=':', alpha=0.6)
    
    # 水印
    formula_text = f"$E_a = {model.intercept_:.0f} + {model.coef_[0]:.1f} Q_N + {model.coef_[1]:.1f} N_H + {model.coef_[2]:.2f} \\kappa_3$"
    plt.text(0.05, 0.95, formula_text, transform=plt.gca().transAxes, fontsize=10, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    plt.tight_layout()
    
    save_name = "Final_ML_Parity_Plot_Chem.png"
    plt.savefig(save_name)
    print(f"✅ 图片已自动保存为: {save_name}")
    plt.show()

if __name__ == "__main__":
    plot_parity()