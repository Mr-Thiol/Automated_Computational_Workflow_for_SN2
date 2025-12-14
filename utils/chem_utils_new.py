##库，通用参数
import os
import subprocess
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
import py3Dmol
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG, display
from rdkit.Geometry.rdGeometry import Point3D

# 设置RDKit显示选项
IPythonConsole.ipython_useSVG = True

#函数 mol3D
def generate_3d_structure_single(smiles, name):
    """生成分子的3D结构并优化"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"无法解析SMILES: {smiles}")
        return None
    
    # 添加氢原子
    mol = Chem.AddHs(mol)
    
    # 生成3D坐标
    AllChem.EmbedMolecule(mol, randomSeed=42)
    
    # 优化结构
    AllChem.MMFFOptimizeMolecule(mol)
    
    print(f"{name} SMILES: {smiles}")
    print(f"{name} 分子式: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
    
    return mol

# 如果需要可视化3D结构，可以使用以下代码
def view_3d_structure(mol, name="3D Structure", 
                                    font_color='black', font_size=16, offset=0.2):
    """在3D中查看分子结构"""
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
    if mmff_props is not None:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props)
        if ff is not None:
            energy = ff.CalcEnergy()
            print(f"分子能量: {energy:.2f} kcal/mol")
                    
    mb = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=600, height=600)
    view.addModel(mb, 'mol')
    view.setStyle({'stick': {'radius': 0.12}, 'sphere': {'radius': 0.4}})
    view.setBackgroundColor('0xeeeeee')
    # 添加原子序号标签 - 使用更大的字体模拟加粗
    conformer = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        # 获取原子坐标
        atom = mol.GetAtomWithIdx(i)
        atom_pos = conformer.GetAtomPosition(i)
        
        # 计算标签位置
        label_x = atom_pos.x + offset
        label_y = atom_pos.y + offset
        label_z = atom_pos.z + offset

       
        # 使用更大的字体模拟加粗效果
        # 修复：直接使用原子索引 i+1 作为标签
        view.addLabel(f"{i+1}",  # 使用原子索引作为标签，从1开始计数
                     {'position': {'x': label_x, 'y': label_y, 'z': label_z},
                      'fontColor': font_color,
                      'backgroundColor': 'transparent',
                      'backgroundOpacity': 0,
                      'fontSize': font_size,
                      'borderColor': 'transparent',
                      'borderWidth': 0,
                      'inFront': True})
    view.zoomTo()
    view.show()
    print(f"{name} 3D结构:")
    return view

def view_optimized_structure_elegant0(results, name="Optimized Structure", 
                                    font_color='black', font_size=16, offset=0.2):
    """使用完全透明的背景显示原子序号，通过增加字体大小模拟加粗效果"""
    geometry = results.get('geometry')
    if not geometry:
        print("未找到优化后的几何结构")
        return None
    
    # 创建XYZ格式的字符串
    xyz_content = f"{len(geometry)}\n\n"
    for atom in geometry:
        xyz_content += f"{atom['element']} {atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f}\n"
    
    # 使用py3Dmol显示
    view = py3Dmol.view(width=600, height=600)
    view.addModel(xyz_content, "xyz")
    
    # 设置分子样式
    view.setStyle({'stick': {'radius': 0.12}, 'sphere': {'radius': 0.4}})
    
    # 添加原子序号标签 - 使用更大的字体模拟加粗
    for i, atom in enumerate(geometry):
        # 计算标签位置（稍微偏移原子中心）
        label_x = atom['x'] + offset
        label_y = atom['y'] + offset
        label_z = atom['z'] + offset
        
        # 使用更大的字体模拟加粗效果
        view.addLabel(f"{atom['atom_num']}",
                     {'position': {'x': label_x, 'y': label_y, 'z': label_z},
                      'fontColor': font_color,
                      'backgroundColor': 'transparent',
                      'backgroundOpacity': 0,
                      'fontSize': font_size,
                      'borderColor': 'transparent',
                      'borderWidth': 0,
                      'inFront': True})
    
    view.zoomTo()
    
    print(f"\n显示 {name}:")
    return view    

def print_molecular_coordinates(mol, name):
    """打印分子坐标"""
    print(f"\n{name} 的3D坐标:")
    print("原子\tX\t\tY\t\tZ")
    conf = mol.GetConformer()
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        print(f"{symbol}\t{pos.x:.6f}\t{pos.y:.6f}\t{pos.z:.6f}")

def create_gaussian_input_advanced(mol, mol_name, charge=0, multiplicity=1, 
                                  method="#PM7 opt(calcall)", filename=None,
                                  nproc=2, mem="2GB"):
    """
    高级版本的Gaussian输入文件生成函数
    """
    if filename is None:
        # 移除特殊字符创建安全文件名
        safe_name = "".join(c for c in mol_name if c.isalnum() or c in ('_', '-')).rstrip()
        filename = f"{safe_name}.gjf"
    
    conf = mol.GetConformer()
    coordinates = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        atom = mol.GetAtomWithIdx(i)
        symbol = atom.GetSymbol()
        coordinates.append(f"{symbol:2s} {pos.x:14.8f} {pos.y:14.8f} {pos.z:14.8f}")
    
    with open(filename, 'w') as f:
        f.write(f"%nproc={nproc}\n")
        f.write(f"%mem={mem}\n")
        f.write("\n")
        f.write(f"{method}\n")
        f.write("\n")
        f.write(f"title: {mol_name}\n")
        f.write("\n")
        f.write(f"{charge},{multiplicity}\n")
        for coord in coordinates:
            f.write(coord + "\n")
        f.write("\n")
    
    print(f"Gaussian输入文件已保存为: {filename}")
    return filename

def run_gaussian_job(gjf_file, log_file=None):
    """
    运行Gaussian 16计算
    
    参数:
    gjf_file: 输入文件路径(.gjf)
    log_file: 输出文件路径(.log)，如果为None则自动生成
    """
    if log_file is None:
        log_file = gjf_file.replace('.gjf', '.log')
    
    try:
        # 方法1: 直接调用g16命令
        cmd = ['g16.exe ', gjf_file, log_file]
        
        # 或者方法2: 使用输入输出重定向
        # cmd = f"g16 < {gjf_file} > {log_file}"
        
        print(f"执行命令: {' '.join(cmd)}")
        #os.system(' '.join(cmd))
        # 执行计算
        result = subprocess.run(cmd, 
                              capture_output=True, 
                              text=True, 
                              check=True)
        
        print("计算完成!")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"计算失败: {e}")
        print(f"错误输出: {e.stderr}")
        return False
    except FileNotFoundError:
        print("错误: 未找到g16命令。请确保Gaussian 16已正确安装并配置环境变量。")
        return False


def print_optimized_geometry(results):
    """打印优化后的分子结构坐标"""
    geometry = results.get('geometry')
    if not geometry:
        print("未找到优化后的几何结构")
        return
    
    print("\n优化后的分子结构坐标:")
    print("原子\t元素\tX\t\tY\t\tZ")
    for atom in geometry:
        print(f"{atom['atom_num']}\t{atom['element']}\t{atom['x']:.6f}\t{atom['y']:.6f}\t{atom['z']:.6f}")

def create_mol_with_atom_numbers(mol):
    """创建带有原子序号的分子图像"""
    # 复制分子
    mol_copy = Chem.Mol(mol)
    
    # 设置原子标签为原子序号
    for atom in mol_copy.GetAtoms():
        atom.SetProp('atomNote', str(atom.GetIdx() + 1))
    
    return mol_copy

def view_optimized_structure_combined(results, original_mol, name="Optimized Structure"):
    """结合2D和3D显示优化后的分子结构"""
    geometry = results.get('geometry')
    if not geometry:
        print("未找到优化后的几何结构")
        return None
    
    # 更新分子坐标
    optimized_mol = update_mol_coordinates(original_mol, results)
    
    # 创建带有原子序号的2D图像
    mol_with_numbers = create_mol_with_atom_numbers(optimized_mol)
    
    # 生成2D图像
    drawer = rdMolDraw2D.MolDraw2DSVG(600, 600)
    drawer.DrawMolecule(mol_with_numbers)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    
    # 显示2D图像
    display(SVG(svg))
    
    # 创建XYZ格式的字符串
    xyz_content = f"{len(geometry)}\n\n"
    for atom in geometry:
        xyz_content += f"{atom['element']} {atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f}\n"
    
    # 使用py3Dmol显示3D结构
    view = py3Dmol.view(width=600, height=600)
    view.addModel(xyz_content, "xyz")
    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'radius': 0.5}})
    view.zoomTo()
    
    print(f"\n显示 {name}:")
    return view

def view_optimized_structure_elegant(results, name="Optimized Structure", 
                                    font_color='black', font_size=16, offset=0.2):
    """使用完全透明的背景显示原子序号，通过增加字体大小模拟加粗效果"""
    geometry = results.get('geometry')
    if not geometry:
        print("未找到优化后的几何结构")
        return None
    
    # 创建XYZ格式的字符串
    xyz_content = f"{len(geometry)}\n\n"
    for atom in geometry:
        xyz_content += f"{atom['element']} {atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f}\n"
    
    # 使用py3Dmol显示
    view = py3Dmol.view(width=600, height=600)
    view.addModel(xyz_content, "xyz")
    
    # 设置分子样式
    view.setStyle({'stick': {'radius': 0.12}, 'sphere': {'radius': 0.4}})
    
    # 添加原子序号标签 - 使用更大的字体模拟加粗
    for i, atom in enumerate(geometry):
        # 计算标签位置（稍微偏移原子中心）
        label_x = atom['x'] + offset
        label_y = atom['y'] + offset
        label_z = atom['z'] + offset
        
        # 使用更大的字体模拟加粗效果
        view.addLabel(f"{atom['atom_num']}",
                     {'position': {'x': label_x, 'y': label_y, 'z': label_z},
                      'fontColor': font_color,
                      'backgroundColor': 'transparent',
                      'backgroundOpacity': 0,
                      'fontSize': font_size,
                      'borderColor': 'transparent',
                      'borderWidth': 0,
                      'inFront': True})
    
    view.zoomTo()
    
    print(f"\n显示 {name}:")
    return view

def generate_3d_structure(smiles, name, separation=10.0):
    """生成分子的3D结构并优化，如果包含多个分子则分开它们"""
    # 检查SMILES是否包含多个分子（用点号分隔）
    if '.' in smiles:
        # 分割SMILES字符串
        smiles_list = smiles.split('.')
        mols = []
        
        # 为每个分子生成3D结构
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                print(f"无法解析SMILES: {smi}")
                continue
            
            # 添加氢原子
            mol = Chem.AddHs(mol)
            
            # 生成3D坐标
            AllChem.EmbedMolecule(mol, randomSeed=42)
            
            # 优化结构并获取能量
            result = AllChem.MMFFOptimizeMolecule(mol)
            
            # 获取力场并计算能量
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
            if mmff_props is not None:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props)
                if ff is not None:
                    energy = ff.CalcEnergy()
                    print(f"分子 {i+1} 优化后能量: {energy:.2f} kcal/mol")
            
            # 在Z轴上平移分子，使它们分开
            conf = mol.GetConformer()
            for j in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(j)
                conf.SetAtomPosition(j, (pos.x, pos.y, pos.z + i * separation))
            
            mols.append(mol)
            print(f"分子 {i+1} SMILES: {smi}")
            print(f"分子 {i+1} 分子式: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
        
        # 合并所有分子
        if mols:
            combined_mol = mols[0]
            for mol in mols[1:]:
                combined_mol = Chem.CombineMols(combined_mol, mol)
            
            print(f"{name} 包含 {len(mols)} 个分子，已在Z轴上分开")
            return combined_mol
        else:
            print(f"无法为 {name} 生成任何分子")
            return None
    else:
        # 单个分子的处理
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"无法解析SMILES: {smiles}")
            return None
        
        # 添加氢原子
        mol = Chem.AddHs(mol)
        
        # 生成3D坐标
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # 优化结构并获取能量
        result = AllChem.MMFFOptimizeMolecule(mol)
        
        # 获取力场并计算能量
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
        if mmff_props is not None:
            ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props)
            if ff is not None:
                energy = ff.CalcEnergy()
                print(f"{name} MMFF94优化后能量: {energy:.2f} kcal/mol")
        
        #print(f"{name} SMILES: {smiles}")
        #print(f"{name} 分子式: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
        
        return mol

        
def generate_3d_structure_advanced(smiles, name, positions=None):
    """生成分子的3D结构并优化，可以自定义每个分子的位置"""
    # 检查SMILES是否包含多个分子（用点号分隔）
    if '.' in smiles:
        # 分割SMILES字符串
        smiles_list = smiles.split('.')
        mols = []
        
        # 为每个分子生成3D结构
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                print(f"无法解析SMILES: {smi}")
                continue
            
            # 添加氢原子
            mol = Chem.AddHs(mol)
            
            # 生成3D坐标
            AllChem.EmbedMolecule(mol, randomSeed=42)
            
            # 优化结构
            AllChem.MMFFOptimizeMolecule(mol)
            
            # 平移分子到指定位置
            if positions and i < len(positions):
                dx, dy, dz = positions[i]
            else:
                # 默认在Z轴上分开
                dx, dy, dz = 0, 0, i * 10.0
            
            conf = mol.GetConformer()
            for j in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(j)
                conf.SetAtomPosition(j, (pos.x + dx, pos.y + dy, pos.z + dz))
            
            mols.append(mol)
            print(f"分子 {i+1} SMILES: {smi}")
            print(f"分子 {i+1} 分子式: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
            if positions and i < len(positions):
                print(f"分子 {i+1} 位置: ({dx}, {dy}, {dz})")
        
        # 合并所有分子
        if mols:
            combined_mol = mols[0]
            for mol in mols[1:]:
                combined_mol = Chem.CombineMols(combined_mol, mol)
            
            print(f"{name} 包含 {len(mols)} 个分子")
            return combined_mol
        else:
            print(f"无法为 {name} 生成任何分子")
            return None
    else:
        # 单个分子的处理
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"无法解析SMILES: {smiles}")
            return None
        
        # 添加氢原子
        mol = Chem.AddHs(mol)
        
        # 生成3D坐标
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # 优化结构
        AllChem.MMFFOptimizeMolecule(mol)
        
        print(f"{name} SMILES: {smiles}")
        print(f"{name} 分子式: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
        
        return mol
        
def update_mol_coordinates(original_mol, results):
    """更新原始分子的3D坐标为优化后的坐标"""
    geometry = results.get('geometry')
    if not geometry:
        print("未找到优化后的几何结构")
        return None
    
    # 复制原始分子
    optimized_mol = Chem.Mol(original_mol)
    
    # 获取构象
    if optimized_mol.GetNumConformers() == 0:
        # 如果没有构象，创建一个
        conf = Chem.Conformer(optimized_mol.GetNumAtoms())
        optimized_mol.AddConformer(conf)
    
    conf = optimized_mol.GetConformer(0)
    
    # 更新坐标
    for i, atom in enumerate(geometry):
        if i < optimized_mol.GetNumAtoms():
            conf.SetAtomPosition(i, (atom['x'], atom['y'], atom['z']))
        else:
            print(f"警告: 几何结构中的原子数 ({len(geometry)}) 与分子中的原子数 ({optimized_mol.GetNumAtoms()}) 不匹配")
            break
    
    return optimized_mol

def user_update_mol_coordinates(newzb, original_mol, results=None):
    """
    使用用户自定义坐标更新分子
    
    参数:
    newzb: 用户提供的坐标列表，格式可以是:
        ["1 C -4.652026 0.000000 -0.000002", ...] 或
        ["C -4.652026 0.000000 -0.000002", ...]
    original_mol: 原始分子对象
    results: 可选参数，保持与之前代码的兼容性
    """
    # 解析用户输入的坐标
    parsed_coords = []
    
    for line in newzb:
        # 移除可能的空白字符并分割
        parts = line.strip().split()
        
        if len(parts) >= 4:
            # 格式: "1 C -4.652026 0.000000 -0.000002" 或 "C -4.652026 0.000000 -0.000002"
            # 检查第一个部分是否是数字
            if parts[0].isdigit():
                # 格式: "1 C -4.652026 0.000000 -0.000002"
                element = parts[1]
                try:
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    parsed_coords.append((element, x, y, z))
                except (ValueError, IndexError):
                    print(f"无法解析坐标行: {line}")
                    continue
            else:
                # 格式: "C -4.652026 0.000000 -0.000002"
                element = parts[0]
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    parsed_coords.append((element, x, y, z))
                except (ValueError, IndexError):
                    print(f"无法解析坐标行: {line}")
                    continue
        else:
            print(f"无法识别的坐标格式: {line}")
    
    # 检查原子数量是否匹配
    if len(parsed_coords) != original_mol.GetNumAtoms():
        print(f"警告: 提供的坐标数量 ({len(parsed_coords)}) 与分子原子数 ({original_mol.GetNumAtoms()}) 不匹配")
        # 可以选择返回或继续
    
    # 更新分子坐标
    conf = original_mol.GetConformer()
    
    for i, (element, x, y, z) in enumerate(parsed_coords):
        if i < original_mol.GetNumAtoms():
            # 检查元素是否匹配（可选）
            atom = original_mol.GetAtomWithIdx(i)
            if atom.GetSymbol() != element:
                print(f"警告: 原子 {i} 元素不匹配，分子中为 {atom.GetSymbol()}，提供的为 {element}")
            
            # 更新坐标
            conf.SetAtomPosition(i, Point3D(x, y, z))
    
    print(f"成功更新 {len(parsed_coords)} 个原子的坐标")
    return original_mol

#函数
def read_gaussian16_output_opt(filename):
#def read_gaussian16_output(filename):
    """
    解析Gaussian16输出文件，提取关键计算结果
    """
    results = {
        'geometry': None,
        'energy': None,
        'max_force': None,
        'homo': None,
        'lumo': None,
        'mulliken_charges': None,
        'polarizability': None,
        'frequencies': None,
        'enthalpy': None,
        'free_energy': None,
        'termination_status': 'Normal termination of Gaussian 16',
        'error_info': None
    }
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # 1. 使用逐行扫描的方法提取最后一个Standard orientation部分
        geometry_sections = []
        current_section = []
        in_geometry_section = False
        found_header = False
        dash_count = 0
        
        for i, line in enumerate(lines):
            # 检查是否进入Standard orientation部分
            if 'Standard orientation:' in line:
                # 如果已经在一个几何结构部分中，先保存它
                if in_geometry_section and current_section:
                    geometry_sections.append(current_section)
                
                # 开始新的几何结构部分
                current_section = [line]
                in_geometry_section = True
                found_header = False
                dash_count = 0
            elif in_geometry_section:
                current_section.append(line)
                
                # 检查是否找到表头分隔线
                if not found_header and '----' in line:
                    dash_count += 1
                    if dash_count == 2:  # 找到第二个分隔线，开始提取坐标
                        found_header = True
                
                # 检查是否到达几何结构部分的结束（找到第三个分隔线）
                if found_header and '----' in line and dash_count >= 3:
                    geometry_sections.append(current_section)
                    current_section = []
                    in_geometry_section = False
                    found_header = False
                    dash_count = 0
        
        # 添加最后一个部分（如果有）
        if in_geometry_section and current_section:
            geometry_sections.append(current_section)
        
        print(f"找到 {len(geometry_sections)} 个几何结构部分")
        
        # 提取最后一个几何结构部分的原子坐标
        if geometry_sections:
            last_geometry = geometry_sections[-1]
            atoms = []
            
            # 跳过表头，找到原子坐标行
            in_coordinates = False
            dash_encountered = 0
            
            for line in last_geometry:
                # 计算遇到的分隔线数量
                if '----' in line:
                    dash_encountered += 1
                    # 第二个分隔线后开始坐标，第三个分隔线结束坐标
                    if dash_encountered == 2:
                        in_coordinates = True
                        continue
                    elif dash_encountered == 3:
                        in_coordinates = False
                        break
                
                # 提取坐标行
                if in_coordinates and line.strip():
                    parts = line.split()
                    # 检查是否是坐标行：应该有6个字段，前三个是整数，后三个是浮点数
                    if len(parts) >= 6:
                        try:
                            atom_num = int(parts[0])
                            element_num = int(parts[1])
                            x, y, z = map(float, parts[3:6])
                            
                            # 将原子序数转换为元素符号
                            element_map = {
                                1: 'H', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 
                                9: 'F',14: 'Si',  15: 'P', 16: 'S', 17: 'Cl',
                                35: 'Br', 53: 'I'
                            }
                            element = element_map.get(element_num, f"X{element_num}")
                            
                            atoms.append({
                                'atom_num': atom_num,
                                'element': element,
                                'x': x, 
                                'y': y, 
                                'z': z
                            })
                        except (ValueError, IndexError) as e:
                            # 忽略解析错误的行
                            continue
            
            if atoms:
                results['geometry'] = atoms
                print(f"成功提取 {len(atoms)} 个原子的坐标")
            else:
                print("在最后一个几何结构部分中未找到有效的原子坐标")
                # 打印部分内容用于调试
                print("最后一个几何结构部分的内容:")
                for i, line in enumerate(last_geometry):
                    print(f"{i}: {line}")
        else:
            print("未找到任何几何结构部分")
        
        # 将内容合并为字符串用于其他提取
        content = ''.join(lines)
        
        # 2. 读取最后一次能量 (SCF Done)
        scf_matches = re.findall(r'SCF Done:.*?=\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)', content)
        #scf_matches = re.findall(r'SCF Done:\s*E\([^)]+\)\s*=\s*([-\d\.]+[Ee][+-]?\d+)', content)
        #scf_matches = re.findall(r'SCF Done:\s*E\([^)]+\)\s*=\s*([-\d\.]+)', content)
        if scf_matches:
            results['energy'] = float(scf_matches[-1])
        
        # 3. 读取最后一次最大力 (Cartesian Forces: Max)
        force_sections = re.findall(r'Cartesian Forces:.*?Max\s*([-\d\.]+)', content, re.DOTALL)
        if force_sections:
            results['max_force'] = float(force_sections[-1])
        
        # 4. 读取HOMO和LUMO - 精确匹配格式
        # 找到所有轨道能级部分
        orbital_sections = re.findall(
            r'Alpha  occ\. eigenvalues --\s*([-\d\.\s-]+)\s*\n\s*Alpha virt\. eigenvalues --\s*([-\d\.\s-]+)',
            content
        )
        
        if orbital_sections:
            # 取最后一组轨道能级
            last_homo_line, last_lumo_line = orbital_sections[-1]
            
            # 提取HOMO (最后一个占据轨道)
            homo_numbers = re.findall(r'[-\d\.]+', last_homo_line.strip())
            if homo_numbers:
                results['homo'] = float(homo_numbers[-1])
            
            # 提取LUMO (第一个虚拟轨道)
            lumo_numbers = re.findall(r'[-\d\.]+', last_lumo_line.strip())
            if lumo_numbers:
                results['lumo'] = float(lumo_numbers[0])
        
        # 5. 读取最后一次Mulliken电荷 - 包含原子序号和元素符号
        mulliken_sections = re.findall(r'Mulliken charges:.*?(?=Sum of Mulliken charges|\n\n|\Z)', content, re.DOTALL)
        if mulliken_sections:
            last_mulliken = mulliken_sections[-1]
            mulliken_data = []  # 存储原子序号、元素符号和电荷
            
            for line in last_mulliken.split('\n'):
                # 匹配格式: "  1  C     0.266044" 或 "  1     0.266044"
                if re.match(r'^\s*\d+\s+[A-Za-z]+\s+[-\d\.]+', line):
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            atom_num = int(parts[0])
                            element = parts[1]
                            charge = float(parts[2])
                            mulliken_data.append({
                                'atom_num': atom_num,
                                'element': element,
                                'charge': charge
                            })
                        except (ValueError, IndexError):
                            continue
                # 匹配没有元素符号的格式: "  1     0.266044"
                elif re.match(r'^\s*\d+\s+[-\d\.]+', line):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            atom_num = int(parts[0])
                            charge = float(parts[1])
                            # 尝试从几何结构中获取元素符号
                            element = "?"
                            if results['geometry']:
                                for atom in results['geometry']:
                                    if atom['atom_num'] == atom_num:
                                        element = atom['element']
                                        break
                            
                            mulliken_data.append({
                                'atom_num': atom_num,
                                'element': element,
                                'charge': charge
                            })
                        except (ValueError, IndexError):
                            continue
            
            results['mulliken_charges'] = mulliken_data
        
        # 6. 读取完整的极化率张量 (6个分量)
        # 查找完整的极化率行
        polar_match = re.search(r'Exact polarizability:\s*([-\d\.\s]+)', content)
        if polar_match:
            polar_line = polar_match.group(1).strip()
            polar_numbers = re.findall(r'[-\d\.]+', polar_line)
            if len(polar_numbers) >= 6:
                results['polarizability'] = [float(x) for x in polar_numbers[:6]]
            elif len(polar_numbers) >= 3:
                # 如果只有3个分量，则复制为6个分量 (xx, yy, zz, xy, xz, yz)
                results['polarizability'] = [float(x) for x in polar_numbers] + [0.0] * (6 - len(polar_numbers))
        
        # 7. 读取第一次出现的三个频率
        freq_match = re.search(r'Frequencies --\s*([-\d\.\s-]+)', content)
        if freq_match:
            first_freq_line = freq_match.group(1).strip()
            freq_numbers = re.findall(r'[-\d\.]+', first_freq_line)
            if len(freq_numbers) >= 3:
                results['frequencies'] = [float(x) for x in freq_numbers[:3]]
        
        # 8. 读取焓
        enthalpy_matches = re.findall(r'Sum of electronic and thermal Enthalpies=\s*([-\d\.]+)', content)
        if enthalpy_matches:
            results['enthalpy'] = float(enthalpy_matches[-1])
        
        # 9. 读取自由能
        free_energy_matches = re.findall(r'Sum of electronic and thermal Free Energies=\s*([-\d\.]+)', content)
        if free_energy_matches:
            results['free_energy'] = float(free_energy_matches[-1])
        
        # 10. 检查计算是否正常结束
        if 'Normal termination of Gaussian 16' not in content:
            results['termination_status'] = 'Error termination'
            # 尝试提取错误信息
            error_match = re.search(r'l\d+\.exe.*?error', content, re.IGNORECASE)
            if error_match:
                results['error_info'] = error_match.group()
            else:
                # 查找其他错误模式
                error_match = re.search(r'Error termination.*', content)
                if error_match:
                    results['error_info'] = error_match.group()
        
        return results
        
    except FileNotFoundError:
        return {'error': f'文件 {filename} 未找到'}
    except Exception as e:
        return {'error': f'解析文件时出错: {str(e)}'}
        

def print_gaussian_results(results):
    """格式化打印Gaussian计算结果"""
    if 'error' in results:
        print("错误:", results['error'])
        return
    
    print("计算状态:", results.get('termination_status', 'Unknown'))
    print("info状态:", results.get('error_info', ''))
    print("能量 (Hartree):", results.get('energy'))
    print("最大力:", results.get('max_force'))
    print("HOMO:", results.get('homo'))
    print("LUMO:", results.get('lumo'))
    
    # 计算能隙
    if results.get('homo') is not None and results.get('lumo') is not None:
        gap = results['lumo'] - results['homo']
        print("能隙:", gap)
    else:
        print("能隙: N/A")
        
    print("精确极化率:", results.get('polarizability'))
    print("频率 (cm⁻¹):", results.get('frequencies'))
    print("焓 (Hartree):", results.get('enthalpy'))
    print("自由能 (Hartree):", results.get('free_energy'))
    
    # 打印格式化的Mulliken电荷
    mulliken_data = results.get('mulliken_charges')
    if mulliken_data:
        print("\nMulliken电荷:")
        for item in mulliken_data:
            print(f"{item['atom_num']:>5}  {item['element']:>2}  {item['charge']:>10.6f}")

def create_advanced_scan_gaussian_input(mol, mol_name, scan_lines,
                                      charge=0, multiplicity=1, filename=None,
                                      method="#PM7 opt(modredundant,calcall)",
                                      nproc=2, mem="2GB"):
    """
    高级版本的柔性扫描输入文件生成函数，接受完整的扫描行作为参数
    
    参数:
    mol: RDKit分子对象（必须包含3D坐标）
    mol_name: 分子名称，用于标题
    scan_lines: 扫描行列表，每个元素是一个完整的扫描指令字符串
    charge: 分子电荷，默认为0
    multiplicity: 自旋多重度，默认为1
    filename: 输出文件名，如果为None则自动生成
    method: 计算方法，默认为"#PM7 opt(modredundant,calcall)"
    nproc: 处理器数量
    mem: 内存大小
    """
    if filename is None:
        safe_name = "".join(c for c in mol_name if c.isalnum() or c in ('_', '-')).rstrip()
        filename = f"{safe_name}_scan.gjf"
    
    # 获取分子坐标
    conf = mol.GetConformer()
    coordinates = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        atom = mol.GetAtomWithIdx(i)
        symbol = atom.GetSymbol()
        coordinates.append(f"{symbol:2s} {pos.x:14.8f} {pos.y:14.8f} {pos.z:14.8f}")
    
    # 写入Gaussian输入文件
    with open(filename, 'w') as f:
        # 全局参数
        f.write(f"%nproc={nproc}\n")
        f.write(f"%mem={mem}\n")
        f.write("\n")
        
        # 计算级别
        f.write(f"{method}\n")
        f.write("\n")
        
        # 标题
        f.write(f"title: {mol_name} - Flexible Scan\n")
        f.write("\n")
        
        # 电荷和多重度
        f.write(f"{charge},{multiplicity}\n")
        
        # 3D坐标
        for coord in coordinates:
            f.write(coord + "\n")
        
        # 柔性扫描控制行
        f.write("\n")  # 空一行
        for scan_line in scan_lines:
            f.write(f"{scan_line}\n")
            print(f"扫描行: {scan_line}")
    
    print(f"柔性扫描Gaussian输入文件已保存为: {filename}")
    return filename


#函数
def read_scan_output(filename):
    """
    从Gaussian扫描输出文件中读取每个扫描点的最优结构和能量
    使用状态跟踪的方法
    """
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        scan_points = []  # 存储所有扫描点
        cur_geometry = None  # 当前几何结构
        cur_energy = None  # 当前能量
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # 检查是否是标准坐标部分的开头
            if 'Standard orientation:' in line:
                #print(line)
                # 提取几何结构
                cur_geometry = extract_geometry_from_position(lines, i)
                # 跳过已经处理的行
                i = find_end_of_geometry(lines, i)
                continue
            
            # 检查是否是SCF能量行
            elif 'SCF Done:' in line:
                #print(line)
                #match = re.search(r'SCF Done:\s*E\([^)]+\)\s*=\s*([-\d\.]+)', line)
                match = re.search(r'SCF Done:\s*E\([^)]+\)\s*=\s*([-\d\.]+(?:[eE][-+]?\d+)?)', line)
                if match:
                    cur_energy = float(match.group(1))
            
            # 检查是否是优化完成标记
            elif 'Optimization completed.' in line:
                #print(line,cur_energy)
                #print(cur_geometry)
                # 如果当前有几何结构和能量，则保存为一个扫描点
                if cur_geometry is not None and cur_energy is not None:
                    scan_points.append({
                        'step': len(scan_points) + 1,
                        'energy': cur_energy,
                        'geometry': cur_geometry
                    })
                    #print(f"找到扫描点 {len(scan_points)}: 能量 = {cur_energy:.8f}")
                
                # 重置当前状态（可选，取决于是否需要）
                # cur_geometry = None
                # cur_energy = None
            
            i += 1
        
        #print(f"成功提取 {len(scan_points)} 个扫描点的数据")
        return scan_points
        
    except FileNotFoundError:
        print(f"错误: 文件 {filename} 未找到")
        return []
    except Exception as e:
        print(f"解析扫描输出文件时出错: {e}")
        import traceback
        traceback.print_exc()
        return []
        
def find_end_of_geometry(lines, start_idx):
    """
    找到几何结构部分的结束位置
    """
    # 查找第二个分隔线（表示几何结构表格结束）
    dash_count = 0
    for i in range(start_idx, min(len(lines), start_idx + 100)):
        if '---' in lines[i]:
            dash_count += 1
            if dash_count == 2:  # 第二个分隔线表示表格结束
                #print(f"几何结构表格在第{i}行结束")
                return i + 1  # 返回下一行的索引
    
    # 如果没有找到第二个分隔线，返回一个合理的值
    #print(f"在从第{start_idx}行开始的100行内未找到第二个分隔线，返回默认位置")
    return min(len(lines), start_idx + 100)
    
def extract_geometry_from_position(lines, start_idx):
    """
    从指定行开始提取几何结构
    """
    atoms = []
    
    # 找到三个分隔线
    dash_count = 0
    data_start = None
    for i in range(start_idx, min(len(lines), start_idx + 100)):
        if '---' in lines[i]:
            dash_count += 1
            if dash_count == 2:
                data_start = i + 1  # 第二个分隔线之后是数据开始
            elif dash_count == 3:
                data_end = i
                break
    else:
        # 如果循环正常结束，说明没有找到三个分隔线
        print(f"在从第{start_idx}行开始的100行内未找到三个分隔线")
        return None

    if data_start is None:
        print("未找到数据开始位置")
        return None

    #print(f"几何结构数据从第{data_start}行开始，到第{data_end}行结束")
    
    # 提取原子坐标
    for i in range(data_start, data_end):
        line = lines[i].strip()
        # 匹配原子坐标行
        if re.match(r'^\s*\d+\s+\d+\s+\d+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+', line):
            parts = line.split()
            if len(parts) >= 6:
                try:
                    atom_num = int(parts[0])
                    element_num = int(parts[1])
                    x, y, z = map(float, parts[3:6])
                    
                    # 将原子序数转换为元素符号
                    element_map = {
                        1: 'H', 6: 'C', 7: 'N', 8: 'O', 
                        9: 'F', 15: 'P', 16: 'S', 17: 'Cl',
                        35: 'Br', 53: 'I'
                    }
                    element = element_map.get(element_num, f"X{element_num}")
                    
                    atoms.append({
                        'atom_num': atom_num,
                        'element': element,
                        'x': x, 
                        'y': y, 
                        'z': z
                    })
                except (ValueError, IndexError) as e:
                    print(f"解析原子坐标行时出错: {line}, 错误: {e}")
                    continue
    
    #print(f"提取到 {len(atoms)} 个原子")
    return atoms if atoms else None

    
def print_scan_results(scan_points):
    """
    打印扫描结果
    """
    if not scan_points:
        print("未找到扫描点数据")
        return
    
    print(f"找到 {len(scan_points)} 个扫描点:")
    print("步数\t能量 (Hartree)")
    for point in scan_points:
        print(f"{point['step']}\t{point['energy']:.8f}")
    
    # 打印第一个和最后一个扫描点的结构对比
    if len(scan_points) > 1:
        print(f"\n第一个扫描点的能量: {scan_points[0]['energy']:.8f} Hartree")
        print(f"最后一个扫描点的能量: {scan_points[-1]['energy']:.8f} Hartree")
        print(f"能量变化: {scan_points[-1]['energy'] - scan_points[0]['energy']:.8f} Hartree")

#函数
def extract_scan_point_to_mol(scan_points, step_number, original_mol=None, show_3d=True):
    """
    从扫描数据中提取特定扫描点的结构和能量，并转换为 RDKit 分子对象
    
    参数:
    scan_points: 扫描点数据列表
    step_number: 要提取的扫描点步数（从1开始）
    original_mol: 原始分子对象（可选，如果提供则用于更新坐标）
    show_3d: 是否显示3D结构
    
    返回:
    tuple: (分子对象, 能量值) 或 None（如果提取失败）
    """
    # 查找指定步数的扫描点
    target_point = None
    for point in scan_points:
        if point['step'] == step_number:
            target_point = point
            break
    
    if target_point is None:
        print(f"未找到扫描点 {step_number}")
        return None
    
    print(f"提取扫描点 {step_number}:")
    print(f"能量: {target_point['energy']:.8f} Hartree")
    
    # 如果有原始分子对象，使用它来更新坐标
    if original_mol is not None:
        # 创建一个字典格式的results，与update_mol_coordinates函数兼容
        results = {'geometry': target_point['geometry']}
        mol = update_mol_coordinates(original_mol, results)
        if mol is None:
            print("更新分子坐标失败")
            return None
    else:
        # 如果没有提供原始分子，创建一个新的分子对象
        mol = create_mol_from_geometry(target_point['geometry'])
        if mol is None:
            print("从几何结构创建分子失败")
            return None
    
    # 显示3D结构
    if show_3d:
        view = view_3d_structure(mol, f"扫描点 {step_number}")
        if view:
            view.show()
    
    return mol, target_point['energy']

def create_mol_from_geometry(geometry):
    """
    从几何结构创建RDKit分子对象
    
    参数:
    geometry: 几何结构数据
    
    返回:
    RDKit分子对象
    """
    try:
        # 创建空的分子对象
        mol = Chem.RWMol()
        
        # 添加原子
        for atom_data in geometry:
            element = atom_data['element']
            # 将元素符号转换为原子类型
            element_to_atomic_num = {
                'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
                'P': 15, 'S': 16, 'Cl': 17, 'Br': 35, 'I': 53
            }
            atomic_num = element_to_atomic_num.get(element, 6)  # 默认使用碳
            
            new_atom = Chem.Atom(atomic_num)
            mol.AddAtom(new_atom)
        
        # 转换为完整的分子对象
        mol = mol.GetMol()
        
        # 添加氢原子
        #mol = Chem.AddHs(mol)
        
        # 创建构象并设置坐标
        conf = Chem.Conformer(len(geometry))
        for i, atom_data in enumerate(geometry):
            conf.SetAtomPosition(i, (atom_data['x'], atom_data['y'], atom_data['z']))
        
        mol.AddConformer(conf)
        
        return mol
        
    except Exception as e:
        print(f"从几何结构创建分子时出错: {e}")
        return None

def list_scan_points(scan_points):
    """
    列出所有扫描点的摘要信息
    """
    if not scan_points:
        print("没有可用的扫描点数据")
        return
    
    print(f"扫描点总数: {len(scan_points)}")
    print("步数\t能量 (Hartree)")
    for point in scan_points:
        print(f"{point['step']}\t{point['energy']:.8f}")

def view_scan_point_interactive(scan_points, original_mol=None):
    """
    交互式查看扫描点
    """
    if not scan_points:
        print("没有可用的扫描点数据")
        return
    
    # 列出所有扫描点
    list_scan_points(scan_points)
    
    # 让用户选择扫描点
    while True:
        try:
            step_input = input("\n请输入要查看的扫描点步数 (输入 'q' 退出): ")
            if step_input.lower() == 'q':
                break
            
            step_number = int(step_input)
            if step_number < 1 or step_number > len(scan_points):
                print(f"无效的步数，请输入 1 到 {len(scan_points)} 之间的数字")
                continue
            
            # 提取并显示扫描点
            result = extract_scan_point_to_mol(scan_points, step_number, original_mol)
            if result:
                mol, energy = result
                print(f"扫描点 {step_number} 已显示")
                
        except ValueError:
            print("请输入有效的数字")
        except Exception as e:
            print(f"处理扫描点时出错: {e}")

##函数IRC
def read_irc_output(filename):
    """
    从Gaussian IRC输出文件中读取每个IRC点的结构和能量
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        irc_points = []  # 存储所有IRC点
        cur_geometry = None  # 当前几何结构
        cur_energy = None  # 当前能量
        cur_path = None    # 当前路径编号
        cur_point = None   # 当前点编号
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # 检查是否是IRC点结束标志 Point Number:  10          Path Number:   2
            if 'Point Number:' in line and 'Path Number:' in line:
                # 提取路径和点编号
                print(line)
                path_match = re.search(r'Path Number:\s+(\d+)', line)
                point_match = re.search(r'Point Number:\s+(\d+)', line)
                
                if path_match and point_match:
                    new_path = int(path_match.group(1))
                    new_point = int(point_match.group(1))
                    
                    # 如果当前有几何结构和能量，则保存为一个IRC点
                    if cur_geometry is not None and cur_energy is not None:
                        irc_points.append({
                            'path_number': cur_path,
                            'point_number': cur_point,
                            'energy': cur_energy,
                            'geometry': cur_geometry
                        })
                        print(f"找到IRC点: 路径 {cur_path}, 点 {cur_point}, 能量 = {cur_energy:.8f}")
                    
                    # 更新当前路径和点编号
                    cur_path = new_path
                    cur_point = new_point
                    # 重置当前几何结构和能量，准备接收新的点数据
                    cur_geometry = None
                    cur_energy = None
            
            # 检查是否是输入坐标部分的开头 - IRC通常使用Input orientation
            elif 'Input orientation:' in line:
                # 提取几何结构
                cur_geometry = extract_geometry_from_position(lines, i)
                # 跳过已经处理的行
                i = find_end_of_geometry(lines, i)
                continue
            
            # 检查是否是SCF能量行
            elif 'SCF Done:' in line:
                match = re.search(r'SCF Done:\s*E\([^)]+\)\s*=\s*([-\d\.]+(?:[eE][-+]?\d+)?)', line)
                if match:
                    cur_energy = float(match.group(1))
            
            # 检查IRC计算完成标志，保存最后一个点
            elif 'IRC--maximum number of cycles reached' in line or 'Reaction path following complete' in line:
                if cur_geometry is not None and cur_energy is not None and cur_path is not None and cur_point is not None:
                    irc_points.append({
                        'path_number': cur_path,
                        'point_number': cur_point,
                        'energy': cur_energy,
                        'geometry': cur_geometry
                    })
                    print(f"找到最后一个IRC点: 路径 {cur_path}, 点 {cur_point}, 能量 = {cur_energy:.8f}")
            
            i += 1
        
        # 确保保存最后一个点（如果存在）
        if cur_geometry is not None and cur_energy is not None and cur_path is not None and cur_point is not None:
            irc_points.append({
                'path_number': cur_path,
                'point_number': cur_point,
                'energy': cur_energy,
                'geometry': cur_geometry
            })
        
        print(f"成功提取 {len(irc_points)} 个IRC点的数据")
        return irc_points
        
    except FileNotFoundError:
        print(f"错误: 文件 {filename} 未找到")
        return []
    except Exception as e:
        print(f"解析IRC输出文件时出错: {e}")
        import traceback
        traceback.print_exc()
        return []

def print_irc_results(irc_points):
    """
    打印IRC结果
    """
    if not irc_points:
        print("未找到IRC点数据")
        return
    
    # 按路径分组，处理可能的None值
    paths = {}
    for point in irc_points:
        path_num = point['path_number']
        # 如果路径编号为None，分配一个默认值
        if path_num is None:
            path_num = 0
        
        if path_num not in paths:
            paths[path_num] = []
        paths[path_num].append(point)
    
    # 对每个路径的点按点编号排序
    for path_num in paths:
        # 确保点编号不为None
        paths[path_num].sort(key=lambda x: x['point_number'] if x['point_number'] is not None else 0)
    
    print(f"找到 {len(paths)} 条IRC路径:")
    
    # 对路径编号进行安全排序
    sorted_paths = sorted((p for p in paths.keys() if p is not None), key=lambda x: x if x is not None else 0)
    
    # 强行使用第一个点的能量作为参考能量（TS点）
    if irc_points and irc_points[0]['energy'] is not None:
        ref_energy = irc_points[0]['energy']
        hartree_to_kcal = 627.509  # 转换因子
        print(f"参考能量 (TS点): {ref_energy:.8f} Hartree")
    else:
        print("错误：第一个点没有能量数据")
        return
    
    for path_num in sorted_paths:
        points = paths[path_num]
        print(f"\n路径 {path_num} (共 {len(points)} 个点):")
        print("点编号\t能量 (Hartree)\t相对能量 (kcal/mol)")
        
        # 确保能量不为None
        valid_points = [p for p in points if p['energy'] is not None]
        if not valid_points:
            continue
            
        for point in valid_points:
            rel_energy = (point['energy'] - ref_energy) * hartree_to_kcal
            point_num = point['point_number'] if point['point_number'] is not None else "N/A"
            print(f"{point_num}\t{point['energy']:.8f}\t{rel_energy:+.2f}")
    
    # 打印总结信息 - 只处理有效的点
    valid_points = [point for point in irc_points if point['energy'] is not None]
    if len(valid_points) > 1:
        print("\nIRC总结:")
        all_energies = [point['energy'] for point in valid_points]
        min_energy = min(all_energies)
        max_energy = max(all_energies)
        energy_span = (max_energy - min_energy) * 627.509
        
        print(f"最低能量点: {min_energy:.8f} Hartree")
        print(f"最高能量点: {max_energy:.8f} Hartree") 
        print(f"能量跨度: {energy_span:.2f} kcal/mol")

    
def extract_irc_point_to_mol(irc_points, path_number, point_number, original_mol=None,show_3d=True):
    """
    从IRC点数据中提取特定路径和点的分子结构
    
    参数:
    irc_points: read_irc_output函数返回的IRC点列表
    path_number: 路径编号
    point_number: 点编号
    original_mol: 可选的原始分子对象（用于保持键连接信息）
    
    返回:
    (mol_data, energy) 元组，其中mol_data包含分子信息，energy是该点的能量
    如果未找到指定点，返回None
    """
    # 查找指定的IRC点
    target_point = None
    for point in irc_points:
        if point['path_number'] == path_number and point['point_number'] == point_number:
            target_point = point
            break
    
    if target_point is None:
        print(f"未找到路径 {path_number} 点 {point_number} 的数据")
        # 打印可用的路径和点以供参考
        if irc_points:
            available_points = set((p['path_number'], p['point_number']) for p in irc_points)
            print(f"可用的IRC点: {sorted(available_points)}")
        else:
            print("没有找到任何IRC点数据")
        return None
    
    energy = target_point['energy']
    geometry = target_point['geometry']

    # 如果有原始分子对象，使用它来更新坐标
    if original_mol is not None:
        # 创建一个字典格式的results，与update_mol_coordinates函数兼容
        results = {'geometry': target_point['geometry']}
        mol = update_mol_coordinates(original_mol, results)
        if mol is None:
            print("更新分子坐标失败")
            return None
    else:
        # 如果没有提供原始分子，创建一个新的分子对象
        mol = create_mol_from_geometry(target_point['geometry'])
        if mol is None:
            print("从几何结构创建分子失败")
            return None
    
    # 显示3D结构
    if show_3d:
        view = view_3d_structure(mol, f"扫描点 {path_number,point_number}")
        #if view:
        #    view.show()    
            
    # 创建通用的分子数据结构
    mol_data = {
        'path_number': path_number,
        'point_number': point_number,
        'energy': energy,
        'geometry': geometry,
        'elements': [atom['element'] for atom in geometry],
        'coordinates': [(atom['x'], atom['y'], atom['z']) for atom in geometry]
    }
    
    print(f"IRC路径 {path_number} 点 {point_number} 已提取，能量: {energy:.8f} Hartree")
    return mol, target_point['energy']
    #return mol_data, energy

#函数#connect two mol and constraint opt
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def connect_molecules_with_constrained_optimization(mol1, i1, i2, mol2, j1, j2, 
                                                  fixed_atoms=None, 
                                                  distance_constraints=None,
                                                  optimization_steps=500, 
                                                  force_field='MMFF94'):
    """
    将mol1和mol2连接起来，删除指定的H原子并在相应的C原子间形成新键，然后进行带约束的力场优化
    
    参数:
    mol1: 第一个分子
    i1: mol1中要连接的C原子索引
    i2: mol1中要删除的H原子索引  
    mol2: 第二个分子（COOH基团）
    j1: mol2中要连接的C原子索引
    j2: mol2中要删除的H原子索引
    fixed_atoms: 要固定的原子索引列表（这些原子在优化过程中位置不变）
    distance_constraints: 距离约束字典，格式为 {(atom1, atom2): target_distance}
    optimization_steps: 优化步数，默认500
    force_field: 使用的力场，'MMFF94'或'UFF'，默认MMFF94
    
    返回:
    连接并优化后的新分子
    """
    
    # 创建可编辑的分子对象
    combined = Chem.RWMol(Chem.CombineMols(mol1, mol2))
    
    # 计算mol1中的原子数，用于调整mol2的原子索引
    mol1_num_atoms = mol1.GetNumAtoms()
    
    # 调整mol2的原子索引
    j1_new = j1 + mol1_num_atoms
    j2_new = j2 + mol1_num_atoms
    
    # 在要连接的原子间添加键
    combined.AddBond(i1, j1_new, Chem.BondType.SINGLE)
    
    # 删除不需要的H原子
    combined.RemoveAtom(j2_new)  # 先删除mol2的H
    combined.RemoveAtom(i2)      # 再删除mol1的H
    
    # 转换为普通分子对象
    new_mol = combined.GetMol()
    
    # 修复化学键类型和价态问题
    new_mol = sanitize_molecule(new_mol)
    
    # 初始几何优化前的预处理
    try:
        # 确保分子有3D坐标
        if new_mol.GetNumConformers() == 0:
            new_mol = Chem.AddHs(new_mol)
            AllChem.EmbedMolecule(new_mol, randomSeed=42)
        else:
            # 如果有构象但需要添加H，先移除构象再添加H
            new_mol = Chem.AddHs(new_mol)
            AllChem.EmbedMolecule(new_mol, randomSeed=42)
        
        # 进行带约束的力场优化
        new_mol = optimize_with_constraints(
            new_mol, 
            fixed_atoms=fixed_atoms,
            distance_constraints=distance_constraints,
            max_iters=optimization_steps,
            force_field=force_field
        )
        
    except Exception as e:
        print(f"力场优化过程中出现错误: {e}")
        print("返回未优化的结构")
    
    return new_mol

def sanitize_molecule(mol):
    """
    修复分子的化学键类型和价态问题
    
    参数:
    mol: 要修复的分子
    
    返回:
    修复后的分子
    """
    try:
        # 尝试进行标准的化学检查
        Chem.SanitizeMol(mol)
        return mol
    except:
        # 如果标准检查失败，尝试手动修复
        print("标准化学检查失败，尝试手动修复...")
        
        # 创建可编辑分子
        rw_mol = Chem.RWMol(mol)
        
        # 尝试修复键序
        try:
            # 使用Kekulize处理芳香结构
            Chem.Kekulize(rw_mol, clearAromaticFlags=True)
        except:
            print("Kekulization失败，跳过此步骤")
        
        # 转换为普通分子
        mol = rw_mol.GetMol()
        
        # 尝试设置合理的键序
        try:
            # 使用感知器检测键序
            AllChem.Compute2DCoords(mol)
            Chem.AssignAtomChiralTagsFromStructure(mol)
            Chem.AssignStereochemistry(mol)
        except:
            print("键序检测失败，跳过此步骤")
        
        return mol

def optimize_with_constraints(mol, fixed_atoms=None, distance_constraints=None, 
                             max_iters=500, force_field='MMFF94'):
    """
    对分子进行带约束的几何优化
    """
    try:
        # 确保分子有3D坐标
        if mol.GetNumConformers() == 0:
            # 如果分子没有H原子，先添加H
            if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
                mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # 获取当前构象
        conf_id = mol.GetConformer().GetId()
        
        # 尝试MMFF94力场优化
        if force_field.upper() == 'MMFF94':
            try:
                # 首先尝试计算MMFF94力场参数
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
                if mmff_props is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)
                    if ff is not None:
                        # 添加固定原子约束
                        if fixed_atoms is not None:
                            print('fixed_atoms in MMFF94:',fixed_atoms)
                            for atom_idx in fixed_atoms:
                                ff.AddFixedPoint(atom_idx)
                        
                        # 添加距离约束
                        if distance_constraints is not None:
                            for (atom1, atom2), target_dist in distance_constraints.items():
                                ff.AddDistanceConstraint(atom1, atom2, target_dist, target_dist, 100.0)
                        
                        # 进行优化
                        result = ff.Minimize(maxIts=max_iters)
                        e0=ff.CalcEnergy()
                        print(f"MMFF94优化完成，能量: {e0}")
                        return mol
                    else:
                        print("无法创建MMFF94力场")
                else:
                    print("无法获取MMFF94属性")
            except Exception as e:
                print(f"MMFF94优化失败: {e}")
        
        # 如果MMFF94失败，尝试UFF力场
        print("尝试UFF力场优化...")
        try:
            # 对于UFF，先确保分子有环信息
            try:
                Chem.GetSymmSSSR(mol)
            except:
                print("环信息初始化失败，继续UFF优化")
            
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
            if ff is not None:
                # 添加固定原子约束
                if fixed_atoms is not None:
                    for atom_idx in fixed_atoms:
                        ff.AddFixedPoint(atom_idx)
                
                # 添加距离约束
                if distance_constraints is not None:
                    for (atom1, atom2), target_dist in distance_constraints.items():
                        ff.AddDistanceConstraint(atom1, atom2, target_dist, target_dist, 100.0)
                
                # 进行优化
                result = ff.Minimize(maxIts=max_iters)
                print(f"UFF优化完成，能量: {result}")
            else:
                print("无法创建UFF力场")
        except Exception as e:
            print(f"UFF优化也失败: {e}")
            
    except Exception as e:
        print(f"带约束的几何优化错误: {e}")
    
    return mol

def align_molecules_for_connection(mol1, i1, i2, mol2, j1, j2):
    """
    调整两个分子的位置和方向，为连接做准备
    
    参数:
    mol1: 第一个分子
    i1: mol1中要连接的原子索引(RDKit索引)
    i2: mol1中要删除的原子索引(RDKit索引)
    mol2: 第二个分子
    j1: mol2中要连接的原子索引(RDKit索引)
    j2: mol2中要删除的原子索引(RDKit索引)
    
    返回:
    调整后的mol1和mol2
    """
    # 复制分子以避免修改原始结构
    mol1_aligned = Chem.Mol(mol1)
    mol2_aligned = Chem.Mol(mol2)
    
    # 确保分子有构象
    if mol1_aligned.GetNumConformers() == 0:
        raise ValueError("mol1没有3D构象")
    if mol2_aligned.GetNumConformers() == 0:
        raise ValueError("mol2没有3D构象")
    
    # 获取构象
    conf1 = mol1_aligned.GetConformer()
    conf2 = mol2_aligned.GetConformer()
    
    # 1. 获取mol1中相关原子的坐标
    pos_i1 = np.array(conf1.GetAtomPosition(i1))
    pos_i2 = np.array(conf1.GetAtomPosition(i2))
    
    # 2. 将mol1平移，使i2(H原子)位于原点
    translation1 = -pos_i2
    for atom_idx in range(mol1_aligned.GetNumAtoms()):
        pos = np.array(conf1.GetAtomPosition(atom_idx))
        new_pos = pos + translation1
        conf1.SetAtomPosition(atom_idx, Chem.rdGeometry.Point3D(new_pos[0], new_pos[1], new_pos[2]))
    
    # 更新i1的位置
    pos_i1 = np.array(conf1.GetAtomPosition(i1))
    
    # 3. 旋转mol1，使i1原子位于X轴负半轴
    vec_i1_to_i2 = pos_i1  # 因为i2现在在原点
    
    # 计算旋转矩阵，将vec_i1_to_i2旋转到(-1, 0, 0)方向
    target_vec = np.array([-1.0, 0.0, 0.0])
    
    # 归一化向量
    vec_i1_to_i2_norm = vec_i1_to_i2 / np.linalg.norm(vec_i1_to_i2)
    target_vec_norm = target_vec / np.linalg.norm(target_vec)
    
    # 计算旋转轴和角度
    rotation_axis = np.cross(vec_i1_to_i2_norm, target_vec_norm)
    rotation_axis_norm = rotation_axis / np.linalg.norm(rotation_axis) if np.linalg.norm(rotation_axis) > 0 else np.array([0, 0, 1])
    
    # 计算旋转角度
    dot_product = np.dot(vec_i1_to_i2_norm, target_vec_norm)
    rotation_angle = np.arccos(np.clip(dot_product, -1.0, 1.0))
    
    # 应用旋转到mol1的所有原子
    if rotation_angle > 1e-10:  # 只有当需要旋转时才执行
        rotation_matrix = rotation_matrix_from_axis_angle(rotation_axis_norm, rotation_angle)
        for atom_idx in range(mol1_aligned.GetNumAtoms()):
            pos = np.array(conf1.GetAtomPosition(atom_idx))
            new_pos = np.dot(rotation_matrix, pos)
            conf1.SetAtomPosition(atom_idx, Chem.rdGeometry.Point3D(new_pos[0], new_pos[1], new_pos[2]))
    
    # 4. 获取mol2中相关原子的坐标
    pos_j1 = np.array(conf2.GetAtomPosition(j1))
    pos_j2 = np.array(conf2.GetAtomPosition(j2))
    
    # 5. 将mol2平移，使j2(H原子)位于原点
    translation2 = -pos_j2
    for atom_idx in range(mol2_aligned.GetNumAtoms()):
        pos = np.array(conf2.GetAtomPosition(atom_idx))
        new_pos = pos + translation2
        conf2.SetAtomPosition(atom_idx, Chem.rdGeometry.Point3D(new_pos[0], new_pos[1], new_pos[2]))
    
    # 更新j1的位置
    pos_j1 = np.array(conf2.GetAtomPosition(j1))
    
    # 6. 旋转mol2，使j1原子位于X轴正半轴
    vec_j2_to_j1 = pos_j1  # 因为j2现在在原点
    
    # 计算旋转矩阵，将vec_j2_to_j1旋转到(1, 0, 0)方向
    target_vec2 = np.array([1.0, 0.0, 0.0])
    
    # 归一化向量
    vec_j2_to_j1_norm = vec_j2_to_j1 / np.linalg.norm(vec_j2_to_j1)
    target_vec2_norm = target_vec2 / np.linalg.norm(target_vec2)
    
    # 计算旋转轴和角度
    rotation_axis2 = np.cross(vec_j2_to_j1_norm, target_vec2_norm)
    rotation_axis2_norm = rotation_axis2 / np.linalg.norm(rotation_axis2) if np.linalg.norm(rotation_axis2) > 0 else np.array([0, 0, 1])
    
    # 计算旋转角度
    dot_product2 = np.dot(vec_j2_to_j1_norm, target_vec2_norm)
    rotation_angle2 = np.arccos(np.clip(dot_product2, -1.0, 1.0))
    
    # 应用旋转到mol2的所有原子
    if rotation_angle2 > 1e-10:  # 只有当需要旋转时才执行
        rotation_matrix2 = rotation_matrix_from_axis_angle(rotation_axis2_norm, rotation_angle2)
        for atom_idx in range(mol2_aligned.GetNumAtoms()):
            pos = np.array(conf2.GetAtomPosition(atom_idx))
            new_pos = np.dot(rotation_matrix2, pos)
            conf2.SetAtomPosition(atom_idx, Chem.rdGeometry.Point3D(new_pos[0], new_pos[1], new_pos[2]))
    
    # 7. 将mol2沿X轴移动，使两个连接原子之间有合适的距离
    # 典型的C-C键长度约为1.54 Å
    bond_length = 1.54
    for atom_idx in range(mol2_aligned.GetNumAtoms()):
        pos = np.array(conf2.GetAtomPosition(atom_idx))
        new_pos = pos + np.array([bond_length, 0, 0])
        conf2.SetAtomPosition(atom_idx, Chem.rdGeometry.Point3D(new_pos[0], new_pos[1], new_pos[2]))
    
    return mol1_aligned, mol2_aligned

def rotation_matrix_from_axis_angle(axis, angle):
    """
    根据旋转轴和角度计算旋转矩阵
    """
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    one_minus_cos_a = 1 - cos_a
    
    rotation_matrix = np.array([
        [cos_a + x*x*one_minus_cos_a, x*y*one_minus_cos_a - z*sin_a, x*z*one_minus_cos_a + y*sin_a],
        [y*x*one_minus_cos_a + z*sin_a, cos_a + y*y*one_minus_cos_a, y*z*one_minus_cos_a - x*sin_a],
        [z*x*one_minus_cos_a - y*sin_a, z*y*one_minus_cos_a + x*sin_a, cos_a + z*z*one_minus_cos_a]
    ])
    
    return rotation_matrix

def convert_user_indices_to_rdkit(user_indices):
    """
    将用户索引(从1开始)转换为RDKit索引(从0开始)
    """
    if isinstance(user_indices, list):
        return [idx - 1 for idx in user_indices if idx > 0]
    elif isinstance(user_indices, dict):
        converted = {}
        for (atom1, atom2), distance in user_indices.items():
            converted[(atom1-1, atom2-1)] = distance
        return converted
    elif isinstance(user_indices, int):
        return user_indices - 1 if user_indices > 0 else user_indices
    else:
        return user_indices

def map_indices_after_connection(original_indices, deleted_indices, mol1_num_atoms, mol2_num_atoms):
    """
    在连接后重新映射原子索引
    
    参数:
    original_indices: 原始索引列表或字典 (RDKit索引)
    deleted_indices: 被删除的原子索引列表
    mol1_num_atoms: mol1的原子数量
    mol2_num_atoms: mol2的原子数量
    
    返回:
    重新映射后的索引
    """
    # 创建完整的索引映射表
    # 首先，处理mol1的原子
    index_map = {}
    new_idx = 0
    
    # 映射mol1的原子
    for i in range(mol1_num_atoms):
        if i not in deleted_indices:
            index_map[i] = new_idx
            new_idx += 1
    
    # 映射mol2的原子 (注意：mol2的原始索引需要加上mol1_num_atoms)
    for i in range(mol2_num_atoms):
        mol2_idx = i + mol1_num_atoms
        if mol2_idx not in deleted_indices:
            index_map[mol2_idx] = new_idx
            new_idx += 1
    
    # 现在映射原始索引
    if isinstance(original_indices, list):
        mapped = []
        for idx in original_indices:
            if idx in index_map:
                mapped.append(index_map[idx])
            else:
                print(f"警告: 原子 {idx} 已被删除，无法固定")
        return mapped
    elif isinstance(original_indices, dict):
        mapped = {}
        for (atom1, atom2), distance in original_indices.items():
            if atom1 in index_map and atom2 in index_map:
                mapped[(index_map[atom1], index_map[atom2])] = distance
            else:
                print(f"警告: 原子对 ({atom1}, {atom2}) 中至少有一个原子已被删除，无法应用约束")
        return mapped
    else:
        return index_map.get(original_indices, -1)



def create_precise_index_mapping(mol1, mol2, i1, i2, j1, j2, fixed_atoms_rdkit):
    """
    创建精确的索引映射表
    
    参数:
    mol1: 第一个分子
    mol2: 第二个分子
    i1, i2: mol1中要连接和删除的原子索引(RDKit索引)
    j1, j2: mol2中要连接和删除的原子索引(RDKit索引)
    fixed_atoms_rdkit: 要固定的原子索引列表(RDKit索引)
    
    返回:
    映射后的固定原子索引列表
    """
    # 创建映射表
    index_map = {}
    
    # 映射mol1的原子
    new_idx = 0
    for atom_idx in range(mol1.GetNumAtoms()):
        if atom_idx != i2:  # 跳过被删除的原子
            index_map[atom_idx] = new_idx
            new_idx += 1
    
    # 映射mol2的原子
    mol1_num_atoms = mol1.GetNumAtoms()
    for atom_idx in range(mol2.GetNumAtoms()):
        if atom_idx != j2:  # 跳过被删除的原子
            # mol2的原子在合并后的索引是 atom_idx + mol1_num_atoms
            original_combined_idx = atom_idx + mol1_num_atoms
            index_map[original_combined_idx] = new_idx
            new_idx += 1
    
    # 现在映射固定原子
    fixed_atoms_mapped = []
    for atom_idx in fixed_atoms_rdkit:
        if atom_idx in index_map:
            fixed_atoms_mapped.append(index_map[atom_idx])
        else:
            print(f"警告: 原子 {atom_idx+1} (用户索引) 已被删除，无法固定")
    
    return fixed_atoms_mapped

def simple_connect_molecules(mol1, ii1, ii2, mol2, jj1, jj2, fixed_atoms=None, distance_constraints=None):
    """
    连接两个分子，预先调整位置和方向
    
    参数:
    mol1: 第一个分子
    ii1: mol1中要连接的原子索引(用户计数，从1开始)
    ii2: mol1中要删除的原子索引(用户计数，从1开始)
    mol2: 第二个分子
    jj1: mol2中要连接的原子索引(用户计数，从1开始)
    jj2: mol2中要删除的原子索引(用户计数，从1开始)
    fixed_atoms: 要固定的原子索引列表(用户计数，从1开始)
    distance_constraints: 距离约束字典，键为原子对(用户计数，从1开始)
    
    返回:
    连接后的新分子
    """
    # 将用户索引转换为RDKit索引
    i1 = ii1 - 1
    i2 = ii2 - 1
    j1 = jj1 - 1
    j2 = jj2 - 1
    
    # 转换固定原子和距离约束的索引
    fixed_atoms_rdkit = convert_user_indices_to_rdkit(fixed_atoms) if fixed_atoms else None
    distance_constraints_rdkit = convert_user_indices_to_rdkit(distance_constraints) if distance_constraints else None
    
    print(f"用户索引转换:")
    print(f"  mol1: ii1={ii1}->i1={i1}, ii2={ii2}->i2={i2}")
    print(f"  mol2: jj1={jj1}->j1={j1}, jj2={jj2}->j2={j2}")
    
    # 验证原子索引
    if i1 >= mol1.GetNumAtoms() or i2 >= mol1.GetNumAtoms():
        raise ValueError(f"mol1原子索引超出范围: 分子有{mol1.GetNumAtoms()}个原子")
    if j1 >= mol2.GetNumAtoms() or j2 >= mol2.GetNumAtoms():
        raise ValueError(f"mol2原子索引超出范围: 分子有{mol2.GetNumAtoms()}个原子")
    
    # 调整两个分子的位置和方向
    print("调整分子位置和方向...")
    mol1_aligned, mol2_aligned = align_molecules_for_connection(mol1, i1, i2, mol2, j1, j2)
    
    # 创建可编辑的分子对象
    combined = Chem.RWMol(Chem.CombineMols(mol1_aligned, mol2_aligned))
    
    # 计算mol1中的原子数，用于调整mol2的原子索引
    mol1_num_atoms = mol1_aligned.GetNumAtoms()
    
    # 调整mol2的原子索引
    j1_new = j1 + mol1_num_atoms
    j2_new = j2 + mol1_num_atoms
    
    # 验证原子类型
    atom_i1 = combined.GetAtomWithIdx(i1)
    atom_i2 = combined.GetAtomWithIdx(i2)
    atom_j1 = combined.GetAtomWithIdx(j1_new)
    atom_j2 = combined.GetAtomWithIdx(j2_new)
    
    print(f"连接原子验证:")
    print(f"  mol1 原子{i1}: {atom_i1.GetSymbol()}")
    print(f"  mol1 原子{i2}: {atom_i2.GetSymbol()} (将被删除)")
    print(f"  mol2 原子{j1_new}: {atom_j1.GetSymbol()}")
    print(f"  mol2 原子{j2_new}: {atom_j2.GetSymbol()} (将被删除)")
    
    # 获取连接前的键
    bond_i = mol1.GetBondBetweenAtoms(i1, i2)
    bond_j = mol2.GetBondBetweenAtoms(j1, j2)
    
    # 在要连接的原子间添加键
    combined.AddBond(i1, j1_new, Chem.BondType.SINGLE)
    
    # 删除不需要的H原子
    # 注意：先删除索引较大的原子，以避免索引变化
    indices_to_remove = sorted([i2, j2_new], reverse=True)
    for idx in indices_to_remove:
        combined.RemoveAtom(idx)
    
    # 转换为普通分子对象并初始化环信息
    new_mol = combined.GetMol()
    
    # 确保分子有正确的环信息
    try:
        Chem.GetSymmSSSR(new_mol)
    except:
        print("警告：环信息初始化失败，但这可能不影响后续操作")
    
    # 为优化准备约束信息
    fixed_atoms_mapped = None
    if fixed_atoms_rdkit:
        # 创建一个简化的映射：仅保留在删除原子之前的原子
        fixed_atoms_mapped = []
        mol1_num = mol1.GetNumAtoms()
        mol2_num = mol2.GetNumAtoms()
        
        for atom_idx in fixed_atoms_rdkit:
            if atom_idx < mol1_num:  # 来自mol1的原子
                # 如果这个原子不是被删除的i2，则映射它
                if atom_idx != i2:
                    # 在删除i2后，需要调整索引
                    if atom_idx > i2:
                        fixed_atoms_mapped.append(atom_idx - 1)
                    else:
                        fixed_atoms_mapped.append(atom_idx)
            else:  # 来自mol2的原子，但需要转换为原始mol2的索引
                mol2_idx = atom_idx - mol1_num
                if mol2_idx < mol2_num and mol2_idx != j2:  # 不是被删除的j2
                    # 在删除i2后，mol2原子的新索引 = 原子在合并后的原始索引 - 1 (因为删除了i2)
                    # 但需要先计算在删除所有原子后的位置
                    # 简化处理：计算在合并后的原始位置，然后减去删除的原子数
                    original_combined_idx = mol2_idx + mol1_num
                    # 删除i2后，如果atom_idx > i2，需要减1
                    if original_combined_idx > i2:
                        adjusted_idx = original_combined_idx - 1
                    else:
                        adjusted_idx = original_combined_idx
                    # 再考虑删除j2_new的影响
                    if adjusted_idx > j2_new:
                        adjusted_idx -= 1
                    fixed_atoms_mapped.append(adjusted_idx)
        
        print(f"固定原子映射: {fixed_atoms_rdkit} -> {fixed_atoms_mapped}")
    
    # 进行优化
    if fixed_atoms_mapped:
        print("进行约束优化...")
        try:
            new_mol = optimize_with_constraints(
                new_mol, 
                fixed_atoms=fixed_atoms_mapped,
                max_iters=200,
                force_field='MMFF94'
            )
        except Exception as e:
            print(f"优化失败: {e}，返回未优化的结构")
    else:
        print("没有约束，跳过优化")
    
    return new_mol

def optimize_with_constraints(mol, fixed_atoms=None, distance_constraints=None, 
                             max_iters=500, force_field='MMFF94'):
    """
    对分子进行带约束的几何优化
    
    参数:
    mol: 要优化的分子
    fixed_atoms: 要固定的原子索引列表
    distance_constraints: 距离约束字典，格式为 {(atom1, atom2): target_distance}
    max_iters: 最大迭代次数
    force_field: 使用的力场，'MMFF94'或'UFF'
    
    返回:
    优化后的分子
    """
    try:
        # 确保分子有3D坐标
        if mol.GetNumConformers() == 0:
            # 如果分子没有H原子，先添加H
            if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
                mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # 获取当前构象
        conf_id = mol.GetConformer().GetId()
        
        # 尝试MMFF94力场优化
        if force_field.upper() == 'MMFF94':
            try:
                # 首先尝试计算MMFF94力场参数
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
                if mmff_props is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)
                    if ff is not None:
                        # 计算优化前的能量
                        initial_energy = ff.CalcEnergy()
                        print(f"MMFF94优化前能量: {initial_energy:.2f} kcal/mol")
                        
                        # 添加固定原子约束
                        if fixed_atoms is not None:
                            print('固定原子索引 (MMFF94):', fixed_atoms)
                            for atom_idx in fixed_atoms:
                                # 验证原子索引有效
                                if atom_idx < mol.GetNumAtoms():
                                    ff.AddFixedPoint(atom_idx)
                                    print(f"  固定原子 {atom_idx}: {mol.GetAtomWithIdx(atom_idx).GetSymbol()}")
                                else:
                                    print(f"  警告: 原子索引 {atom_idx} 超出范围")
                        
                        # 添加距离约束
                        if distance_constraints is not None:
                            for (atom1, atom2), target_dist in distance_constraints.items():
                                ff.AddDistanceConstraint(atom1, atom2, target_dist, target_dist, 100.0)
                        
                        # 进行优化，返回的是状态码，不是能量
                        result = ff.Minimize(maxIts=max_iters)
                        
                        # 计算优化后的能量
                        final_energy = ff.CalcEnergy()
                        print(f"MMFF94优化后能量: {final_energy:.2f} kcal/mol")
                        print(f"能量变化: {final_energy - initial_energy:.2f} kcal/mol")
                        
                        # 解释状态码
                        if result == 0:
                            print("MMFF94优化: 成功完成")
                        elif result == 1:
                            print("MMFF94优化: 达到最大迭代次数，可能未完全收敛")
                        else:
                            print(f"MMFF94优化: 返回状态码 {result}")
                        
                        return mol
                    else:
                        print("无法创建MMFF94力场")
                else:
                    print("无法获取MMFF94属性")
            except Exception as e:
                print(f"MMFF94优化失败: {e}")
        
        # 如果MMFF94失败，尝试UFF力场
        print("尝试UFF力场优化...")
        try:
            # 对于UFF，先确保分子有环信息
            try:
                Chem.GetSymmSSSR(mol)
            except:
                print("环信息初始化失败，继续UFF优化")
            
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
            if ff is not None:
                # 计算优化前的能量
                initial_energy = ff.CalcEnergy()
                print(f"UFF优化前能量: {initial_energy:.2f} kcal/mol")
                
                # 添加固定原子约束
                if fixed_atoms is not None:
                    print('固定原子索引 (UFF):', fixed_atoms)
                    for atom_idx in fixed_atoms:
                        # 验证原子索引有效
                        if atom_idx < mol.GetNumAtoms():
                            ff.AddFixedPoint(atom_idx)
                            print(f"  固定原子 {atom_idx}: {mol.GetAtomWithIdx(atom_idx).GetSymbol()}")
                        else:
                            print(f"  警告: 原子索引 {atom_idx} 超出范围")
                
                # 添加距离约束
                if distance_constraints is not None:
                    for (atom1, atom2), target_dist in distance_constraints.items():
                        ff.AddDistanceConstraint(atom1, atom2, target_dist, target_dist, 100.0)
                
                # 进行优化，返回的是状态码，不是能量
                result = ff.Minimize(maxIts=max_iters)
                
                # 计算优化后的能量
                final_energy = ff.CalcEnergy()
                print(f"UFF优化后能量: {final_energy:.2f} kcal/mol")
                print(f"能量变化: {final_energy - initial_energy:.2f} kcal/mol")
                
                # 解释状态码
                if result == 0:
                    print("UFF优化: 成功完成")
                elif result == 1:
                    print("UFF优化: 达到最大迭代次数，可能未完全收敛")
                else:
                    print(f"UFF优化: 返回状态码 {result}")
            else:
                print("无法创建UFF力场")
        except Exception as e:
            print(f"UFF优化也失败: {e}")
            
    except Exception as e:
        print(f"带约束的几何优化错误: {e}")
    
    return mol
    
def optimize_with_constraints_old(mol, fixed_atoms=None, distance_constraints=None, 
                             max_iters=500, force_field='MMFF94'):
    """
    对分子进行带约束的几何优化
    
    参数:
    mol: 要优化的分子
    fixed_atoms: 要固定的原子索引列表(RDKit索引)
    distance_constraints: 距离约束字典，键为原子对(RDKit索引)
    max_iters: 最大迭代次数
    force_field: 使用的力场
    
    返回:
    优化后的分子
    """
    try:
        # 确保分子有3D坐标
        if mol.GetNumConformers() == 0:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # 获取当前构象
        conf_id = mol.GetConformer().GetId()
        
        # UFF力场优化（更宽容）
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        
        print('fixed_atoms in UFF:',fixed_atoms)
        if ff is not None:
            # 添加固定原子约束
            valid_fixed_atoms = []
            for atom_idx in fixed_atoms:
                if atom_idx < mol.GetNumAtoms():
                    valid_fixed_atoms.append(atom_idx)
                    atom = mol.GetAtomWithIdx(atom_idx)
                    print(f"  固定原子 {atom_idx}: {atom.GetSymbol()}")
                else:
                    print(f"  警告: 原子索引 {atom_idx} 超出范围")
            #fixed_atoms = valid_fixed_atoms
            print('fixed_atoms in UFF:',valid_fixed_atoms)
       

            
            # 添加距离约束
            if distance_constraints is not None:
                for (atom1, atom2), target_dist in distance_constraints.items():
                    ff.AddDistanceConstraint(atom1, atom2, target_dist, target_dist, 100.0)
            
            # 进行优化
            ff.Minimize(maxIts=max_iters)
            print("带约束的UFF优化完成")
        else:
            print("警告: 无法创建UFF力场")
            
    except Exception as e:
        print(f"几何优化错误: {e}")
    
    return mol

from rdkit import Chem
from rdkit.Chem import AllChem

def print_molecule_info(mol, mol_name="Molecule"):
    """
    打印分子的详细信息，包括原子信息和连接表
    
    参数:
    mol: RDKit分子对象
    mol_name: 分子名称，用于标识输出
    """
    print(f"\n{'='*60}")
    print(f"{mol_name} 详细信息")
    print(f"{'='*60}")
    
    # 基本信息
    print(f"分子式: {Chem.MolToSmiles(mol)}")
    print(f"原子总数: {mol.GetNumAtoms()}")
    print(f"键总数: {mol.GetNumBonds()}")
    
    # 检查是否有3D坐标
    if mol.GetNumConformers() > 0:
        conformer = mol.GetConformer()
        print(f"有3D坐标: 是 (共 {mol.GetNumConformers()} 个构象)")
    else:
        print("有3D坐标: 否")
    
    # 打印每个原子的详细信息
    print(f"\n{'='*60}")
    print("原子详细信息")
    print(f"{'='*60}")
    
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        atomic_num = atom.GetAtomicNum()
        symbol = atom.GetSymbol()
        formal_charge = atom.GetFormalCharge()
        hybridization = atom.GetHybridization()
        is_aromatic = atom.GetIsAromatic()
        degree = atom.GetDegree()
        total_valence = atom.GetTotalValence()
        explicit_valence = atom.GetExplicitValence()
        implicit_valence = atom.GetImplicitValence()
        
        # 获取坐标（如果有）
        if mol.GetNumConformers() > 0:
            pos = conformer.GetAtomPosition(atom_idx)
            coord_str = f"({pos.x:.3f}, {pos.y:.3f}, {pos.z:.3f})"
        else:
            coord_str = "无3D坐标"
        
        print(f"原子 {atom_idx}:")
        print(f"  元素: {symbol} (原子序数: {atomic_num})")
        print(f"  坐标: {coord_str}")
        print(f"  形式电荷: {formal_charge}")
        print(f"  杂化: {hybridization}")
        print(f"  芳香性: {is_aromatic}")
        print(f"  连接度: {degree}")
        print(f"  总价电子: {total_valence}")
        print(f"  显式价电子: {explicit_valence}")
        print(f"  隐式价电子: {implicit_valence}")
        
        # 打印连接的原子和键信息
        if degree > 0:
            neighbors_info = []
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                bond_type = bond.GetBondType()
                bond_type_str = str(bond_type).split('.')[-1]
                neighbors_info.append(f"{neighbor_idx}({neighbor.GetSymbol()}, {bond_type_str})")
            
            print(f"  连接原子: {', '.join(neighbors_info)}")
        print()

def print_connection_table(mol, mol_name="Molecule"):
    """
    打印分子的连接表（邻接表）
    
    参数:
    mol: RDKit分子对象
    mol_name: 分子名称，用于标识输出
    """
    print(f"\n{'='*60}")
    print(f"{mol_name} 连接表")
    print(f"{'='*60}")
    
    # 创建连接表
    connection_table = {}
    
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        
        connections = []
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            bond_type = bond.GetBondType()
            bond_type_str = str(bond_type).split('.')[-1]
            
            connections.append((neighbor_idx, neighbor.GetSymbol(), bond_type_str))
        
        connection_table[atom_idx] = {
            'symbol': symbol,
            'connections': connections
        }
    
    # 打印连接表
    for atom_idx, info in connection_table.items():
        symbol = info['symbol']
        connections = info['connections']
        
        if connections:
            conn_str = ", ".join([f"{idx}({sym}, {bond})" for idx, sym, bond in connections])
        else:
            conn_str = "无连接"
        
        print(f"原子 {atom_idx}({symbol}): {conn_str}")

def print_bond_information(mol, mol_name="Molecule"):
    """
    打印分子中所有键的详细信息
    
    参数:
    mol: RDKit分子对象
    mol_name: 分子名称，用于标识输出
    """
    print(f"\n{'='*60}")
    print(f"{mol_name} 键信息")
    print(f"{'='*60}")
    
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        
        begin_idx = begin_atom.GetIdx()
        end_idx = end_atom.GetIdx()
        begin_symbol = begin_atom.GetSymbol()
        end_symbol = end_atom.GetSymbol()
        
        bond_type = bond.GetBondType()
        bond_type_str = str(bond_type).split('.')[-1]
        
        is_aromatic = bond.GetIsAromatic()
        is_conjugated = bond.GetIsConjugated()
        is_in_ring = bond.IsInRing()
        
        print(f"键 {begin_idx}({begin_symbol}) - {end_idx}({end_symbol}):")
        print(f"  类型: {bond_type_str}")
        print(f"  芳香性: {is_aromatic}")
        print(f"  共轭性: {is_conjugated}")
        print(f"  在环中: {is_in_ring}")
        print()

def get_atom_connections(mol, atom_idx):
    """
    获取特定原子的连接信息
    
    参数:
    mol: RDKit分子对象
    atom_idx: 原子索引
    
    返回:
    连接信息字典
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    symbol = atom.GetSymbol()
    
    connections = []
    for neighbor in atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
        bond_type = bond.GetBondType()
        bond_type_str = str(bond_type).split('.')[-1]
        
        connections.append({
            'atom_idx': neighbor_idx,
            'symbol': neighbor.GetSymbol(),
            'bond_type': bond_type_str
        })
    
    return {
        'atom_idx': atom_idx,
        'symbol': symbol,
        'connections': connections
    }

def print_detailed_molecule_analysis(mol, mol_name="Molecule"):
    """
    打印分子的完整分析信息
    
    参数:
    mol: RDKit分子对象
    mol_name: 分子名称，用于标识输出
    """
    print_molecule_info(mol, mol_name)
    print_connection_table(mol, mol_name)
    print_bond_information(mol, mol_name)
    
    # 打印总结信息
    print(f"\n{'='*60}")
    print(f"{mol_name} 总结")
    print(f"{'='*60}")
    
    # 统计不同元素的数量
    element_count = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        element_count[symbol] = element_count.get(symbol, 0) + 1
    
    print("元素组成:")
    for element, count in element_count.items():
        print(f"  {element}: {count}")
    
    # 统计不同键类型的数量
    bond_type_count = {}
    for bond in mol.GetBonds():
        bond_type = str(bond.GetBondType()).split('.')[-1]
        bond_type_count[bond_type] = bond_type_count.get(bond_type, 0) + 1
    
    print("\n键类型统计:")
    for bond_type, count in bond_type_count.items():
        print(f"  {bond_type}: {count}")