import os
import sys
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms
import numpy as np

# 导入工具箱
try:
    from chem_utils import *
except ImportError:
    try:
        from chem_utils_new import *
    except ImportError:
        print("❌ 找不到 chem_utils.py")
        sys.exit()

N_PROC = 6
MEM_SIZE = "4GB"
NAME = "Me2NH_Ultimate"

def build_perfect_ts_guess():
    print(">>> 正在手动组装 Me2NH 完美过渡态初猜...")
    
    # 1. 构建二甲胺 (Me2NH)
    amine = Chem.MolFromSmiles("CNC")
    amine = Chem.AddHs(amine)
    AllChem.EmbedMolecule(amine, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(amine)
    conf = amine.GetConformer()
    
    # 2. 找到 N 和它的邻居
    atoms = amine.GetAtoms()
    n_idx = [a.GetIdx() for a in atoms if a.GetSymbol() == 'N'][0]
    
    # 3. 对齐：让 N 位于原点，且孤对电子方向指向 Z 轴正方向
    # 简单做法：算出 N 的几何中心和 N 的向量，反向延长
    # 或者简单点：让 N 在原点，两个 C 尽量往 Z 轴负方向压
    # 这里我们用最简单的方法：先平移 N 到原点
    pos_n = conf.GetAtomPosition(n_idx)
    shift = Geometry.Point3D(-pos_n.x, -pos_n.y, -pos_n.z)
    
    for i in range(amine.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Geometry.Point3D(p.x+shift.x, p.y+shift.y, p.z+shift.z))
    
    # 4. 构建 碘甲烷 (MeI)
    mei = Chem.MolFromSmiles("CI")
    mei = Chem.AddHs(mei)
    AllChem.EmbedMolecule(mei, randomSeed=42)
    conf_mei = mei.GetConformer()
    
    c_idx = [a.GetIdx() for a in mei.GetAtoms() if a.GetSymbol() == 'C'][0]
    i_idx = [a.GetIdx() for a in mei.GetAtoms() if a.GetSymbol() == 'I'][0]
    
    # 5. 摆放 MeI
    # S_N2 关键几何参数：
    # N...C 距离约 2.2 埃
    # C...I 距离约 2.6 埃
    # N-C-I 成 180 度直线
    
    # 把 MeI 的 C 放在 (0, 0, 2.2)
    conf_mei.SetAtomPosition(c_idx, Geometry.Point3D(0, 0, 2.2))
    
    # 把 MeI 的 I 放在 (0, 0, 4.8) (2.2 + 2.6)
    conf_mei.SetAtomPosition(i_idx, Geometry.Point3D(0, 0, 4.8))
    
    # 把 MeI 的 H 放在 C 附近平面上
    # 这一步稍微粗糙点没关系，优化会修好的
    h_idxs = [a.GetIdx() for a in mei.GetAtoms() if a.GetSymbol() == 'H']
    # 将 H 稍微往后（Z轴正向）推一点，呈现平面状
    offsets = [
        (0.0, 1.0, 2.3), 
        (0.86, -0.5, 2.3), 
        (-0.86, -0.5, 2.3)
    ]
    for k, h_id in enumerate(h_idxs):
        if k < 3:
            conf_mei.SetAtomPosition(h_id, Geometry.Point3D(offsets[k][0], offsets[k][1], offsets[k][2]))
            
    # 6. 组合分子
    combo = Chem.CombineMols(amine, mei)
    
    return combo

# --- 主流程 ---
if __name__ == "__main__":
    # 1. 生成结构
    mol = build_perfect_ts_guess()
    
    # 2. 生成输入文件
    # 必须用 calcall，因为我们是手摆的结构，Hessian 矩阵必须算准才能引导方向
    print(f"🚀 提交 Me2NH 强力修复任务 (不扫描，直接 Opt TS)...")
    
    create_gaussian_input_advanced(
        mol=mol, 
        mol_name=NAME, 
        # scan_lines=[], 
        charge=0, # 中性体系
        nproc=N_PROC, 
        mem=MEM_SIZE,
        method="#PM7 opt(ts,calcall,noeig) freq", 
        filename=f"{NAME}.gjf"
    )
    
    # 3. 运行
    run_gaussian_job(f"{NAME}.gjf")
    
    # 4. 验证
    log_file = f"{NAME}.log"
    res = read_gaussian16_output_opt(log_file)
    freqs = res.get('frequencies', [])
    energy = res.get('energy', 0)
    
    print("\n" + "="*40)
    print("修复结果验证")
    print("="*40)
    
    if freqs and freqs[0] < -400:
        print(f"✅✅✅ 成功！虚频: {freqs[0]:.2f} cm^-1")
        print("这才是真正的 S_N2 过渡态！")
        print(f"TS 能量: {energy:.5f}")
        
        # 估算一下 Ea (假设 Reactant 能量大概是 -0.002 左右)
        # 之前的 MeNH2 TS 是 0.04，EtNH2 是 0.02
        # 这个 Me2NH 应该也在 0.02 - 0.08 之间
        if energy < 0.1:
            print("⚡ 能量值正常！不再是 229 kJ/mol 那个怪物了。")
        else:
            print("⚠️ 能量依然偏高，请仔细检查。")
            
    else:
        print(f"❌ 失败。虚频: {freqs}")
        print("建议：如果还是不行，就用 SISSO 剔除法处理。")