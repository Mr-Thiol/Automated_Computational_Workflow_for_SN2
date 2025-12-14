import pandas as pd
import numpy as np
import itertools
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# ================= 配置 =================
INPUT_FILE = "SISSO_Input.csv"
TARGET = "Ea"
FEATURES = ["HOMO", "MW", "Chi1v", "HeavyAtoms"]
OPS = ['+', '-', '*', '/'] # 允许的运算符
# =======================================

def run_mini_sisso():
    print(f">>> 启动 Mini-SISSO 筛选...")
    
    try:
        df = pd.read_csv(INPUT_FILE)
    except:
        print("❌ 找不到 input 文件")
        return

    # 归一化/标准化数据 (可选，为了数值稳定，这里先不做，保持物理意义)
    X_raw = df[FEATURES].values
    y = df[TARGET].values
    names = df['Name'].values
    
    best_r2 = -100
    best_desc = ""
    best_model = None
    
    results = []

    print(f"  正在构建特征空间 (原始特征 {len(FEATURES)} 个)...")
    
    # 1. 遍历 1D 描述符 (原始特征)
    for i, f in enumerate(FEATURES):
        x_curr = X_raw[:, i].reshape(-1, 1)
        reg = LinearRegression().fit(x_curr, y)
        r2 = reg.score(x_curr, y)
        results.append((r2, f, reg.coef_[0], reg.intercept_))

    # 2. 遍历 2D 组合特征 (SISSO Level 1: A+B, A/B, A*B...)
    # 复杂度 O(N^2)，这里特征很少，瞬间跑完
    for i in range(len(FEATURES)):
        for j in range(len(FEATURES)):
            feat_a = FEATURES[i]
            feat_b = FEATURES[j]
            val_a = X_raw[:, i]
            val_b = X_raw[:, j]
            
            # 生成组合
            combos = [
                (val_a + val_b, f"({feat_a} + {feat_b})"),
                (val_a * val_b, f"({feat_a} * {feat_b})"),
                (val_a - val_b, f"({feat_a} - {feat_b})"),
                (np.abs(val_a - val_b), f"|{feat_a} - {feat_b}|")
            ]
            
            # 除法要小心分母为0
            if not np.any(np.abs(val_b) < 1e-6):
                combos.append((val_a / val_b, f"({feat_a} / {feat_b})"))
            
            for x_new, desc in combos:
                x_new = x_new.reshape(-1, 1)
                reg = LinearRegression().fit(x_new, y)
                r2 = reg.score(x_new, y)
                results.append((r2, desc, reg.coef_[0], reg.intercept_))

    # 3. 排序并输出 Top 5
    results.sort(key=lambda x: x[0], reverse=True)
    
    print("\n" + "="*50)
    print(f"SISSO 筛选结果 (Target: {TARGET})")
    print("="*50)
    print(f"{'Rank':<5} | {'R2 Score':<10} | {'Descriptor (Formula)':<30} | {'Coeff'}")
    print("-" * 70)
    
    for k in range(min(10, len(results))):
        r2, desc, coef, intercept = results[k]
        print(f"{k+1:<5} | {r2:.4f}     | {desc:<30} | {coef:.2e}")

    # 4. 输出最佳公式的解释
    top_r2, top_desc, top_coef, top_int = results[0]
    print("\n>>> 最佳物理模型发现:")
    print(f"Ea = {top_coef:.2f} * {top_desc} + {top_int:.2f}")
    print(f"拟合优度 R2 = {top_r2:.4f}")
    
    # 5. 简单的物理意义推测
    print("\n[物理意义解读]")
    if "HOMO" in top_desc and "Chi1v" in top_desc:
        if "/" in top_desc:
            print("公式包含 (位阻 / 电子) 或反之，这完美体现了两个效应的竞争关系！")
        elif "*" in top_desc:
            print("公式体现了位阻和电子效应的耦合作用。")

run_mini_sisso()