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
NAME = "Me3N_Ultimate"

def build_perfect_ts_guess():
    print(">>> 正在手动组装 Me3N 完美过渡态初猜 (优化甲基取向)...")
    
    # 1. 构建三甲胺 (Me3N)
    amine = Chem.MolFromSmiles("CN(C)C")
    amine = Chem.AddHs(amine)
    # 使用 ETKDG 算法生成更合理的初始构象
    AllChem.EmbedMolecule(amine, randomSeed=42, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    AllChem.MMFFOptimizeMolecule(amine)
    conf = amine.GetConformer()
    
    # 2. 找到 N 原子
    atoms = amine.GetAtoms()
    n_idx = [a.GetIdx() for a in atoms if a.GetSymbol() == 'N'][0]
    
    # 3. 对齐：平移 N 到原点
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
    
    # 5. 摆放 MeI (TS 典型距离)
    # 对于位阻大的体系，N...C 距离通常比 2.2 稍远一点点，设为 2.25
    conf_mei.SetAtomPosition(c_idx, Geometry.Point3D(0, 0, 2.25))
    conf_mei.SetAtomPosition(i_idx, Geometry.Point3D(0, 0, 4.90)) # 2.25 + 2.65
    
    # 6. 关键：摆放 MeI 的 H (平面构型)，并旋转以错开位阻
    # 我们让 MeI 的 H 处于 "Staggered" (交叉) 构象，而不是 "Eclipsed" (重叠)
    h_idxs = [a.GetIdx() for a in mei.GetAtoms() if a.GetSymbol() == 'H']
    
    # 这里我们简单地摆成平面三角，Gaussian 优化时会自动微调旋转
    offsets = [
        (0.0, 1.0, 2.35), 
        (0.86, -0.5, 2.35), 
        (-0.86, -0.5, 2.35)
    ]
    for k, h_id in enumerate(h_idxs):
        if k < 3:
            conf_mei.SetAtomPosition(h_id, Geometry.Point3D(offsets[k][0], offsets[k][1], offsets[k][2]))
            
    # 7. 组合
    combo = Chem.CombineMols(amine, mei)
    return combo

if __name__ == "__main__":
    # 1. 生成结构
    mol = build_perfect_ts_guess()
    
    # 2. 提交任务
    # 直接 Opt TS，跳过扫描
    print(f"🚀 提交 Me3N 强力修复任务...")
    
    create_gaussian_input_advanced(
        mol=mol, 
        mol_name=NAME, 
        # scan_lines=[], 
        charge=0,
        nproc=N_PROC, 
        mem=MEM_SIZE,
        # 用 calcall 算准二阶导数，引导它走向正确的鞍点
        method="#PM7 opt(ts,calcall,noeig) freq", 
        filename=f"{NAME}.gjf"
    )
    
    # 3. 运行
    run_gaussian_job(f"{NAME}.gjf")
    
    # 4. 验证结果
    log_file = f"{NAME}.log"
    res = read_gaussian16_output_opt(log_file)
    freqs = res.get('frequencies', [])
    energy = res.get('energy', 0)
    
    print("\n" + "="*40)
    print("Me3N 修复结果")
    print("="*40)
    
    if freqs and freqs[0] < -400:
        print(f"✅ TS 优化成功！")
        print(f"虚频: {freqs[0]:.2f} cm^-1")
        print(f"能量: {energy:.5f} Hartree")
        print(">>> 请一定要去 CSV 里重新算一下 Ea！")
        print("    Ea = (这个能量 - Me3N_scan.log的第一步能量) * 627.5")
    else:
        print(f"❌ 优化失败或未找到 TS。虚频: {freqs}")