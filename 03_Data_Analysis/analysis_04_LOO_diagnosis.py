import pandas as pd
import numpy as np
import warnings
from sklearn.linear_model import LinearRegression

# 忽略 log(0) 或除以 0 的警告
warnings.filterwarnings('ignore')

# ================= 配置 =================
INPUT_FILE = "SISSO_Input_Geometric.csv"
TARGET = "Ea"
# 基础特征
BASE_FEATURES = ["HOMO", "Chi1v", "Sum_Angles", "Avg_Bond_Len", "TPSA"]
# =======================================

def generate_feature_space(df, features):
    """
    构建扩展特征空间 (包含一元和二元变换)
    """
    X_dict = {}
    names = []
    
    # 1. 原始特征
    for f in features:
        val = df[f].values
        X_dict[f] = val
        names.append(f)
        
        # 2. 一元变换 (Unary Operators)
        # 倒数 1/x
        if not np.any(np.abs(val) < 1e-6):
            desc = f"(1/{f})"
            X_dict[desc] = 1 / val
            names.append(desc)
        
        # 平方 x^2
        desc = f"({f}^2)"
        X_dict[desc] = val**2
        names.append(desc)
        
        # 指数 exp(x) (注意数值爆炸，先归一化处理一下再做exp意义不大，这里直接做可能溢出，先跳过或做简单的)
        # 对数 log(|x|)
        if np.all(np.abs(val) > 1e-6):
            desc = f"log(|{f}|)"
            X_dict[desc] = np.log(np.abs(val))
            names.append(desc)

    # 3. 二元变换 (Binary Operators)
    # 转换为列表方便遍历
    base_keys = list(X_dict.keys()) # 包含了一元变换后的特征
    
    # 为了防止特征爆炸，我们只对原始特征做二元组合，或者对一元后的做？
    # 策略：只用原始特征做二元组合，加上一元变换的特征，作为总池子
    # 否则 (1/A) * (B^2) 这种组合太多了
    
    combined_X = {}
    combined_names = []
    
    # 先把已有的加进去
    for k in base_keys:
        combined_X[k] = X_dict[k]
        combined_names.append(k)
    
    n_base = len(features)
    for i in range(n_base):
        for j in range(i, n_base): # 允许自组合
            f_a = features[i]
            f_b = features[j]
            v_a = df[f_a].values
            v_b = df[f_b].values
            
            # A + B
            d = f"({f_a}+{f_b})"
            combined_X[d] = v_a + v_b
            combined_names.append(d)
            
            # A - B
            if i != j:
                d = f"|{f_a}-{f_b}|"
                combined_X[d] = np.abs(v_a - v_b)
                combined_names.append(d)
            
            # A * B
            d = f"({f_a}*{f_b})"
            combined_X[d] = v_a * v_b
            combined_names.append(d)
            
            # A / B
            if not np.any(np.abs(v_b) < 1e-6):
                d = f"({f_a}/{f_b})"
                combined_X[d] = v_a / v_b
                combined_names.append(d)
            
            # B / A
            if i != j and not np.any(np.abs(v_a) < 1e-6):
                d = f"({f_b}/{f_a})"
                combined_X[d] = v_b / v_a
                combined_names.append(d)

    return combined_X, combined_names

def run_loo_sisso():
    print(f">>> 启动留一法 (LOO) SISSO 诊断...")
    print(f"    运算符池: +, -, *, /, ^2, 1/x, log")
    
    try:
        df_full = pd.read_csv(INPUT_FILE)
    except:
        print("❌ 找不到 input 文件")
        return

    # 生成全量特征空间
    X_dict, feature_names = generate_feature_space(df_full, BASE_FEATURES)
    print(f"    特征空间大小: {len(feature_names)} 个候选公式")
    
    y_full = df_full[TARGET].values
    mol_names = df_full['Name'].values
    
    # 定义留一循环
    # 包含 "None" 表示不剔除任何点
    loo_targets = [None] + list(range(len(df_full)))
    
    print("\n" + "="*80)
    print(f"{'Excluded Molecule':<20} | {'Best R2':<8} | {'Best Formula (Model)':<40}")
    print("="*80)
    
    results_summary = []

    for idx in loo_targets:
        if idx is None:
            exclude_name = "None (All Used)"
            mask = np.ones(len(df_full), dtype=bool) # 全选
        else:
            exclude_name = mol_names[idx]
            mask = np.ones(len(df_full), dtype=bool)
            mask[idx] = False # 剔除第 idx 个
            
        y_train = y_full[mask]
        
        best_r2 = -100
        best_desc = ""
        best_coef = 0
        
        # 遍历所有特征找最佳
        for fname in feature_names:
            val = X_dict[fname][mask]
            
            # 简单的清洗：去inf/nan
            if np.isinf(val).any() or np.isnan(val).any(): continue
            if np.std(val) < 1e-9: continue # 方差为0无意义
            
            x_in = val.reshape(-1, 1)
            reg = LinearRegression().fit(x_in, y_train)
            r2 = reg.score(x_in, y_train)
            
            if r2 > best_r2:
                best_r2 = r2
                best_desc = fname
                best_coef = reg.coef_[0]
        
        # 记录
        print(f"{exclude_name:<20} | {best_r2:.4f}   | Ea = {best_coef:+.2e} * {best_desc}")
        results_summary.append((exclude_name, best_r2, best_desc))

    # ----------------------------------------------------
    # 自动总结
    # ----------------------------------------------------
    print("\n" + "="*80)
    print(">>> 诊断报告 (Diagnosis Report)")
    
    # 1. 找基准 (None)
    base_res = [x for x in results_summary if "None" in x[0]][0]
    base_r2 = base_res[1]
    
    # 2. 找提升最大的
    # 排除 None
    loo_res = [x for x in results_summary if "None" not in x[0]]
    # 按 R2 排序
    loo_res.sort(key=lambda x: x[1], reverse=True)
    
    best_loo = loo_res[0]
    improvement = best_loo[1] - base_r2
    
    print(f"1. 原始模型 (全样本) R2: {base_r2:.4f}")
    print(f"2. 最佳留一模型 (剔除 {best_loo[0]}) R2: {best_loo[1]:.4f}")
    print(f"3. 性能提升: {improvement:+.4f}")
    
    if improvement > 0.3:
        print(f"\n[结论] 🚨 发现严重离群点: {best_loo[0]}！")
        print(f"       剔除它之后，物理规律变得非常清晰 ({best_loo[2]})。")
        print(f"       建议：在最终报告中剔除 {best_loo[0]} 并单独讨论其异常原因。")
    elif improvement > 0.1:
        print(f"\n[结论] ⚠️ 疑似离群点: {best_loo[0]}。")
        print(f"       剔除它能显著改善模型，建议深入检查其计算数据。")
    else:
        print(f"\n[结论] ✅ 数据集比较稳健，没有单一离群点能完全解释低 R2。")
        print(f"       可能是物理规律本身就很复杂，或者需要更多特征。")

run_loo_sisso()