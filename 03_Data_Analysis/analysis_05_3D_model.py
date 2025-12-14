import pandas as pd
import numpy as np
import itertools
from sklearn.linear_model import LinearRegression

# ================= 配置 =================
INPUT_FILE = "SISSO_Input_Flexibility.csv"
TARGET = "Ea"
# 候选特征池：把最好的都放进来
FEATURES = [
    "HOMO", "Q_N",       # 电子
    "Chi1v", "Num_H",    # 位阻/几何
    "S_Total", "Kappa3"  # 柔性
]
# =======================================

def run_3d_sisso():
    print(f">>> 启动 3D-SISSO (寻找 3 个描述符的线性组合)...")
    print(f"    警告：样本数(7)较少，3D 模型存在过拟合风险。")
    print(f"    目标：寻找 Ea = c1*Electronic + c2*Steric + c3*Flexibility")
    
    try:
        df = pd.read_csv(INPUT_FILE)
    except:
        print("❌ 找不到 input 文件")
        return
        
    # 剔除缺失值
    df = df.dropna(subset=FEATURES)
    y = df[TARGET].values
    X_dict = {f: df[f].values for f in FEATURES}
    names = df['Name'].values
    
    # 构建特征组合 (C(n, 3))
    combos = list(itertools.combinations(FEATURES, 3))
    print(f"    正在评估 {len(combos)} 种 3D 组合...")
    
    results = []
    
    for c in combos:
        f1, f2, f3 = c
        # 构建矩阵
        X_mat = np.column_stack([X_dict[f1], X_dict[f2], X_dict[f3]])
        
        reg = LinearRegression().fit(X_mat, y)
        r2 = reg.score(X_mat, y)
        
        # 记录结果
        formula = f"{reg.intercept_:.1f} + {reg.coef_[0]:.2f}*{f1} + {reg.coef_[1]:.2f}*{f2} + {reg.coef_[2]:.2f}*{f3}"
        results.append((r2, formula, c))

    # 排序
    results.sort(key=lambda x: x[0], reverse=True)
    
    print("\n" + "="*80)
    print(f"3D-SISSO Top 5 Models")
    print("="*80)
    
    for k in range(min(5, len(results))):
        r2, form, feats = results[k]
        print(f"Rank {k+1} | R2 = {r2:.4f}")
        print(f"         Formula: {form}")
        print(f"         Features: {feats}")
        print("-" * 60)
        
    # 物理意义检查
    top_feats = results[0][2]
    print("\n[物理意义自检]")
    has_elec = any(f in top_feats for f in ["HOMO", "Q_N"])
    has_steric = any(f in top_feats for f in ["Chi1v", "Num_H"])
    has_flex = any(f in top_feats for f in ["S_Total", "Kappa3", "HallKier"])
    
    if has_elec and has_steric and has_flex:
        print("🎉 完美！该模型同时包含了 电子、位阻 和 柔性 三要素！")
        print("   这就是传说中的“大统一模型”。")
    else:
        print("⚠️ 模型偏科了。可能某些效应之间存在强共线性（比如位阻大通常柔性也差）。")

run_3d_sisso()