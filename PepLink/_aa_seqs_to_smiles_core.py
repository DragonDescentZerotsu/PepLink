import copy

from rdkit import Chem
from rdkit.Chem import rdmolops
from typing import Union, List, Dict
import csv
import json

class AAs(object):
    """
    氨基酸类，获得N、C端原子对象和可编辑分子对象
    可以添加dummy atom
    """
    def __init__(self, smiles, idx=None, full_name=None, one_letter_code=None):
        self.smiles = smiles
        self.idx = idx
        self.full_name = full_name
        self.one_letter_code = one_letter_code
        self.mol = Chem.MolFromSmiles(self.smiles)  # 只读
        self.N_terminal_atom = None  # 只读格式的 atom
        self.C_terminal_atom = None
        # rw 表示可编辑的，这里存储可编辑对象的分子和N、C端原子
        self.rwmol = None
        self.rw_N_terminal_atom = None
        self.rw_C_terminal_atom = None
        self.N_dummy = None
        self.C_dummy = None
        # 初始化找到N、C端原子
        self._init_aa()

    def _init_aa(self):
        """
        找到 N 端氮原子和 C 端碳原子并存储到属性里
        :return:
        """
        self.find_N_C_terminal_atoms()

    def find_N_C_terminal_atoms(self):
        """
        找到分子中符合条件的 N 端氮原子和 C 端碳原子。
        :return: N_terminal_idx, C_terminal_idx
        """
        N_terminal_atom, C_terminal_atom = self.strict_find_N_C_terminal_atoms()
        if N_terminal_atom is None or C_terminal_atom is None:
            N_terminal_atom, C_terminal_atom = self.loose_find_N_C_terminal_atoms()
        if N_terminal_atom is None and C_terminal_atom is None:
            print(f'\n{self.idx}: {self.smiles}: Failed to find N terminal atom and C terminal atom')
        elif N_terminal_atom is None and C_terminal_atom is not None:
            print(f'\n{self.idx}: {self.smiles}: Failed to find N terminal atom')
        elif N_terminal_atom is not None and C_terminal_atom is None:
            print(f'\n{self.idx}: {self.smiles}: Failed to find C terminal atom')
        self.N_terminal_atom = N_terminal_atom
        self.C_terminal_atom = C_terminal_atom

    def strict_find_N_C_terminal_atoms(self):
        """
        找到分子中符合条件的 N 端氮原子和 C 端碳原子。
        条件：
          - C 端碳的 degree 必须是 3，并且三个邻居分别是两个氧和一个碳。
          - C 端碳的两跳邻居中包含满足条件的氮原子，氮和 C 端碳之间必须通过碳相连。
          - 满足条件的 N 端氮原子必须至少连接两个氢。
        :param mol: RDKit 的分子对象
        :return: (N_terminal_idx, C_terminal_idx)
        """
        N_terminal_atom = None
        C_terminal_atom = None

        for atom in self.mol.GetAtoms():
            # 识别符合 C 端碳条件的原子：碳，degree=3，并且邻居包括两个氧和一个碳
            if atom.GetSymbol() == "C" and atom.GetDegree() == 3:
                oxygen_count = 0
                carbon_neighbor = None
                double_bond_to_oxygen = False
                _OH_exists = False

                # 遍历 C 端碳的邻居
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "O":
                        oxygen_count += 1
                        # 检查与氧的键类型是否为双键
                        bond = self.mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            double_bond_to_oxygen = True
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            if neighbor.GetDegree() == 1:
                                _OH_exists = True
                    elif neighbor.GetSymbol() == "C":
                        carbon_neighbor = neighbor

                # 确保有两个氧和一个碳邻居，且一个氧通过双键连接, 通过单键连接的 O 必须是只有 1 个neighbour的 -OH
                if oxygen_count == 2 and carbon_neighbor is not None and double_bond_to_oxygen and _OH_exists:
                    # 检查 C 端碳的两跳邻居中的氮原子
                    found_N_terminal = False
                    for two_hop_neighbor in carbon_neighbor.GetNeighbors():
                        if two_hop_neighbor.GetSymbol() == "N" and two_hop_neighbor.GetDegree() == 1:
                            # 确保 N 端氮只有一个邻居，且是 carbon_neighbor
                            if two_hop_neighbor.GetNeighbors()[0].GetIdx() == carbon_neighbor.GetIdx():
                                # 找到符合条件的 N 端氮和 C 端碳
                                C_terminal_atom = atom
                                N_terminal_atom = two_hop_neighbor
                                found_N_terminal = True
                                break

                    # 如果找到合格的 N 端氮，退出循环
                    if found_N_terminal:
                        break
        if N_terminal_atom is not None and C_terminal_atom is not None:
            N_terminal_atom.SetProp('N_terminal', 'True')
            C_terminal_atom.SetProp('C_terminal', 'True')
        return N_terminal_atom, C_terminal_atom

    def loose_find_N_C_terminal_atoms(self):
        """
        找到分子中符合条件的最接近的 N 端氮原子和 C 端碳原子。
        条件：
          - C 端碳的 degree 必须是 3，并且三个邻居分别是两个氧和一个碳。
          - 满足条件的 N 端氮原子必须 degree 为 1，且唯一的邻居是一个碳原子。

          - 没有标准的 C 端就找一个任意有 -OH 连接的 C 或者是其他原子
          - 没有标准的 N 端就找一个离 羧基 C 最近的 N
        :param mol: RDKit 的分子对象
        :return: (N_terminal_idx, C_terminal_idx)
        """
        # 存储符合条件的 C 端碳和 N 端氮的索引
        candidate_C_atoms_idx = []
        candidate_N_atoms_idx = []

        # 没有标准的 N 端就找一个离 羧基 C 最近的 N
        more_loose_candidate_N_atoms_idx = []

        # 没有标准的 C 端就找一个任意有 -OH 连接的 C 或者是其他原子
        more_loose_candidate_C_atoms_idx = []

        # 遍历分子中的所有原子，寻找符合条件的 C 端碳和 N 端氮
        for atom in self.mol.GetAtoms():
            # 符合 C 端碳的条件：degree 为 3，且有两个氧和一个碳邻居
            if atom.GetSymbol() == "C" and atom.GetDegree() == 3:
                oxygen_count = 0
                carbon_neighbor = None
                double_bond_to_oxygen = False
                _OH_exists = False

                # 遍历 C 端碳的邻居
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "O":
                        oxygen_count += 1
                        # 检查与氧的键类型是否为双键
                        bond = self.mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            double_bond_to_oxygen = True
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            if neighbor.GetDegree() == 1:
                                _OH_exists = True
                    elif neighbor.GetSymbol() == "C":
                        carbon_neighbor = neighbor

                # 确保有两个氧邻居、一个碳邻居，且一个氧通过双键连接, 通过单键连接的 O 必须是只有 1 个neighbour的 -OH
                if oxygen_count == 2 and carbon_neighbor is not None and double_bond_to_oxygen and _OH_exists:
                    candidate_C_atoms_idx.append(atom.GetIdx())
                    for two_hop_neighbor in carbon_neighbor.GetNeighbors():
                        if two_hop_neighbor.GetSymbol() == "N" and two_hop_neighbor.GetExplicitValence() < 3:
                            candidate_N_atoms_idx.append(two_hop_neighbor.GetIdx())
                            break


            # 符合 N 端氮的条件：degree 为 1，且唯一的邻居是碳原子
            elif atom.GetSymbol() == "N" and atom.GetDegree() == 1:
                neighbors = atom.GetNeighbors()
                if len(neighbors) == 1 and neighbors[0].GetSymbol() == "C":
                    candidate_N_atoms_idx.append(atom.GetIdx())

        # 初始化最短路径和对应的原子对
        shortest_distance = float('inf')
        N_terminal_atom = None
        C_terminal_atom = None

        # 遍历所有候选对，计算跳数距离
        if len(candidate_C_atoms_idx) > 0 and len(candidate_N_atoms_idx) > 0:
            for c_idx in candidate_C_atoms_idx:
                for n_idx in candidate_N_atoms_idx:
                    # 计算跳数距离
                    distance = len(rdmolops.GetShortestPath(self.mol, c_idx, n_idx))
                    if distance < shortest_distance:
                        shortest_distance = distance
                        N_terminal_atom = self.mol.GetAtomWithIdx(n_idx)
                        C_terminal_atom = self.mol.GetAtomWithIdx(c_idx)

        # 如果没有标准意义上的N端N原子, 找到所有的N原子
        elif len(candidate_C_atoms_idx) > 0 and len(candidate_N_atoms_idx) == 0:
            for atom in self.mol.GetAtoms():
                if atom.GetSymbol() == "N":
                    more_loose_candidate_N_atoms_idx.append(atom.GetIdx())
            for c_idx in candidate_C_atoms_idx:
                for n_idx in more_loose_candidate_N_atoms_idx:
                    # 计算跳数距离
                    distance = len(rdmolops.GetShortestPath(self.mol, c_idx, n_idx))
                    if distance < shortest_distance:
                        shortest_distance = distance
                        N_terminal_atom = self.mol.GetAtomWithIdx(n_idx)
                        C_terminal_atom = self.mol.GetAtomWithIdx(c_idx)

        # 如果没有标准意义上的 C 端 C 原子, 找到所有的有 -OH 连接的 C 原子或其他原子作为替代的 C 端 C 原子
        elif len(candidate_C_atoms_idx) == 0 and len(candidate_N_atoms_idx) > 0:
            for atom in self.mol.GetAtoms():
                if atom.GetSymbol() == "O" and atom.GetDegree() == 1:
                    more_loose_candidate_C_atoms_idx.append(atom.GetNeighbors()[0].GetIdx())

            # 由于一个原子上可能连接两个 -OH 所以用 set 去重
            more_loose_candidate_C_atoms_idx = list(set(more_loose_candidate_C_atoms_idx))
            for c_idx in more_loose_candidate_C_atoms_idx:
                for n_idx in candidate_N_atoms_idx:
                    # 计算跳数距离
                    distance = len(rdmolops.GetShortestPath(self.mol, c_idx, n_idx))
                    if distance < shortest_distance:
                        shortest_distance = distance
                        N_terminal_atom = self.mol.GetAtomWithIdx(n_idx)
                        C_terminal_atom = self.mol.GetAtomWithIdx(c_idx)

        # 如果没有标准意义上的 C 端 C 原子和 N 端 N 原子
        elif len(candidate_C_atoms_idx) == 0 and len(candidate_N_atoms_idx) == 0:
            for atom in self.mol.GetAtoms():
                if atom.GetSymbol() == "O" and atom.GetDegree() == 1:
                    more_loose_candidate_C_atoms_idx.append(atom.GetNeighbors()[0].GetIdx())
                if atom.GetSymbol() == "N":
                    more_loose_candidate_N_atoms_idx.append(atom.GetIdx())

            # 由于一个原子上可能连接两个 -OH 所以用 set 去重
            more_loose_candidate_C_atoms_idx = list(set(more_loose_candidate_C_atoms_idx))
            for c_idx in more_loose_candidate_C_atoms_idx:
                for n_idx in more_loose_candidate_N_atoms_idx:
                    # 计算跳数距离
                    distance = len(rdmolops.GetShortestPath(self.mol, c_idx, n_idx))
                    if distance < shortest_distance:
                        shortest_distance = distance
                        N_terminal_atom = self.mol.GetAtomWithIdx(n_idx)
                        C_terminal_atom = self.mol.GetAtomWithIdx(c_idx)

        # 到这里都还没有找到 C 和 N 只能说明是没有 C 或者没有 N
        # 这里只去找一个，即C端或者N端
        if N_terminal_atom is None and C_terminal_atom is None:
            for atom in self.mol.GetAtoms():

                # 这里 C 端的标准依旧是只要有 -OH 连接就行
                if atom.GetSymbol() == "O" and atom.GetDegree() == 1:
                    more_loose_candidate_C_atoms_idx.append(atom.GetNeighbors()[0].GetIdx())
                if atom.GetSymbol() == "N":
                    more_loose_candidate_N_atoms_idx.append(atom.GetIdx())

            # 由于一个原子上可能连接两个 -OH 所以用 set 去重
            more_loose_candidate_C_atoms_idx = list(set(more_loose_candidate_C_atoms_idx))

            def find_atom_with_lowest_volancy(idx_list, mol, find_max = False):
                """
                找到 index list 中 volancy 最小或最大的原子的 index 返回
                :param idx_list: C_termianl_idx_list 或者 N_termianl_idx_list
                :param mol: 分子对象
                :param find_max: 是否返回 valency 最大的原子
                :return: idx of min valency
                """
                min_volancy = float('inf')
                max_valency = 0
                idx_of_min_valency = 0

                for idx in idx_list:
                    atom = mol.GetAtomWithIdx(idx)
                    explicit_valence = atom.GetExplicitValence()  # 显式电子配位数
                    if find_max:
                        if explicit_valence > max_valency:
                            max_valency = explicit_valence
                            idx_of_min_valency = idx
                    else:
                        if explicit_valence < min_volancy:
                            min_volancy = explicit_valence
                            idx_of_min_valency = idx

                return idx_of_min_valency

            # 如果有 C 端原子
            if len(more_loose_candidate_C_atoms_idx) > 0:

                # 选取valency最大的作为 C 端原子
                min_valency_idx = find_atom_with_lowest_volancy(more_loose_candidate_C_atoms_idx, self.mol, find_max = True)
                C_terminal_atom = self.mol.GetAtomWithIdx(min_valency_idx)

            elif len(more_loose_candidate_N_atoms_idx) > 0:
                # 选取valency最小的作为 N 端原子
                min_valency_idx = find_atom_with_lowest_volancy(more_loose_candidate_N_atoms_idx, self.mol)
                N_terminal_atom = self.mol.GetAtomWithIdx(min_valency_idx)


        if N_terminal_atom is not None:
            N_terminal_atom.SetProp('N_terminal', 'True')
        if C_terminal_atom is not None:
            C_terminal_atom.SetProp('C_terminal', 'True')
        return N_terminal_atom, C_terminal_atom


    def add_dummy_atoms(self, is_n_terminal: bool=False, is_c_terminal: bool=False):
        """
        在指定的 N 端氮和 C 端碳上添加哑原子以便于连接，并根据氨基酸的位置删除不必要的原子。
        :param is_n_terminal: 是否为多肽链的 N 端
        :param is_c_terminal: 是否为多肽链的 C 端
        :return: 修改后的分子对象
        """
        N_dummy, C_dummy = None, None
        rw_mol = Chem.RWMol(self.mol)

        # 更新N、C端原子对象, 只有在 N，C 端原子存在的情况下才在其上添加 dummy atom
        if self.N_terminal_atom is not None:
            is_n_terminal = False
            rw_N_terminal_atom = rw_mol.GetAtomWithIdx(self.N_terminal_atom.GetIdx())
        else:
            is_n_terminal = True

        if self.C_terminal_atom is not None:
            is_c_terminal = False
            rw_C_terminal_atom = rw_mol.GetAtomWithIdx(self.C_terminal_atom.GetIdx())
        else:
            is_c_terminal = True

        # 如果不是 C 端氨基酸，则删除 C 端羧基中的 -OH，并在 C 端碳上添加哑原子
        if not is_c_terminal:
            for neighbor in rw_C_terminal_atom.GetNeighbors():
                # 这个 C 和 -OH 的 O 通过单键连接, 同时还得确保这个邻居 O 确实是 -OH， -OH 中 O 的度为1
                if neighbor.GetSymbol() == 'O' and neighbor.GetDegree() == 1 and rw_mol.GetBondBetweenAtoms(neighbor.GetIdx(), rw_C_terminal_atom.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                    rw_mol.RemoveAtom(neighbor.GetIdx())  # 删除 -OH 中的氧原子
                    break
            # 在 C 端碳上添加哑原子
            C_dummy = rw_mol.AddAtom(Chem.Atom(0))  # 哑原子
            rw_mol.AddBond(rw_C_terminal_atom.GetIdx(), C_dummy, Chem.rdchem.BondType.SINGLE)

        # 如果不是 N 端氨基酸，在 N 端氮上添加哑原子
        if not is_n_terminal:
            # 在 N 端氮上添加哑原子
            N_dummy = rw_mol.AddAtom(Chem.Atom(0))  # 哑原子
            rw_mol.AddBond(rw_N_terminal_atom.GetIdx(), N_dummy, Chem.rdchem.BondType.SINGLE)

        # 返回修改后的分子对象
        self.rwmol = rw_mol
        if self.N_terminal_atom is not None:
            self.rw_N_terminal_atom = rw_N_terminal_atom
            self.rw_N_terminal_atom.SetProp('N_terminal', 'True')
        if self.C_terminal_atom is not None:
            self.rw_C_terminal_atom = rw_C_terminal_atom
            self.rw_C_terminal_atom.SetProp('C_terminal', 'True')
        if N_dummy is not None:
            self.N_dummy = rw_mol.GetAtomWithIdx(N_dummy)
            self.N_dummy.SetProp('N_dummy', 'True')
        if C_dummy is not None:
            self.C_dummy = rw_mol.GetAtomWithIdx(C_dummy)
            self.C_dummy.SetProp('C_dummy', 'True')

    def add_position_property(self, residue: int) -> None:
        """
        在可编辑分子对象中添加 residue 属性
        方便后续添加 intrachain bond 和 interchain bond
        :param residue: 氨基酸在 peptide 中的位置，第一个位置的编号是 1 不是 0
        :return: None
        """
        if self.rwmol is None:
            self.add_dummy_atoms(is_n_terminal=False, is_c_terminal=False)
        for atom in self.rwmol.GetAtoms():
            atom.SetProp('residue', str(residue))


class Peptide(object):
    """
    拼接每一个peptide
    """
    def __init__(self,
                 aa_seqs: Union[str, List[str]],
                 aa_smiles_dict :Dict,
                 idx:int=None,
                 intrachain_bonds: List=[],
                 interchain_bonds: List=[],
                 unusual_aas: List=[],
                 cTerminus: str = None,
                 nTerminus: str = None,
                 c_terminus_modify_name_smiles: Dict = None,
                 n_terminus_modify_name_smiles: Dict = None):
        self.aa_seqs = aa_seqs
        self.AAs_seqs = []  # List
        self.idx = idx
        self.noise_data_flag = False
        self.intrachain_bonds = intrachain_bonds  # [{}, {}, {}, ...]
        self.interchain_bonds = interchain_bonds
        self.unusual_aas = unusual_aas
        self.cTerminus = cTerminus
        self.nTerminus = nTerminus
        self.mols = []  # for Multimer
        self.single_chain_smiles = []  # for Multimer
        self.mol = None
        self.smiles = None
        self.main_chain_linked_mols = []
        self.intrachainBonds_linked_mols = []
        self.cTerminus_modified_mols = []
        self.ncTerminus_modified_mols = []
        self._init_aas(aa_smiles_dict)
        self.link_main_chain()
        self.link_intrachain_bonds()
        # TODO: 这里有点问题，如果是 Multimer 的话那么 cTerminus 应该是个 list 而不是一个值？
        if self.cTerminus is not None and not self.noise_data_flag:
        # if self.cTerminus == 'OMe':
            for rw_intrachainBonds_linked_mol in self.intrachainBonds_linked_mols:
                self.c_terminus_modification(rw_intrachainBonds_linked_mol, c_terminus_modify_name_smiles)
        else:
            self.cTerminus_modified_mols = self.intrachainBonds_linked_mols

        if self.nTerminus is not None and not self.noise_data_flag:
        # if self.cTerminus == 'OMe':
            for cTerminus_modified_mol in self.cTerminus_modified_mols:
                self.n_terminus_modification(cTerminus_modified_mol, n_terminus_modify_name_smiles)
        else:
            self.ncTerminus_modified_mols = self.cTerminus_modified_mols

        if not self.noise_data_flag:
            self.remove_redundant_dummy_atoms()

    def _init_aas(self, aa_smiles_dict:Dict):
        """
        初始化 peptide，里面所有的氨基酸都变成 AA 类
        :param aa_smiles_dict: 氨基酸名称和 SMILES 对
        :return: None
        """
        def _init_single_aa_seq(seq:str, unusual_aas, AMP_id) -> List[AAs]:
            """
            初始化一条 aa_seq 为AAs类的AAs_seq
            :param seq: 氨基酸序列
            :return: list of AAs
            """
            AAs_seq = []

            # 记录处理到第几个 x 或者 X 了
            unusual_aas_count = 0
            for residue, aa in enumerate(seq.strip()):
                if aa == 'X' or aa == 'x':
                    # 提取对应的 non-canonical 氨基酸
                    try:
                        unusual_aa = unusual_aas[unusual_aas_count]
                    except IndexError:
                        # 如果有序列中 Xx 个数和 unusual_aa 个数不一样的问题，把最后一个补成 Leucine
                        print(f'\nID {AMP_id} Missing one non-canonical amino acid, replace it with Leucine(L)')
                        unusual_aa = {"position": len(seq.strip()),
                                      "modificationType":{
                                          "name": "L"
                                      }}

                    try:
                        assert residue +1 == unusual_aa['position'], 'unexpected non-canonical residue position'
                    except AssertionError:
                        # 如果 peptide 中实际的 Xx 的位置和提供的 position 匹配不上则提示
                        print('\nunexpected non-canonical residue position')
                        print(f'The id of this peptide is {AMP_id} and residue {residue+1}, with wrong position being {unusual_aa["position"]}')
                        print('Keep processing other peptides')
                    unusual_aa_name = unusual_aa['modificationType']['name']
                    AA = AAs(aa_smiles_dict[unusual_aa_name], idx=unusual_aa_name)
                    unusual_aas_count += 1
                else:
                    AA = AAs(aa_smiles_dict[aa], idx=aa, one_letter_code=aa)

                # 为每个原子添加 'residue' 属性
                AA.add_position_property(residue+1)
                AAs_seq.append(AA)
            return AAs_seq

        if isinstance(self.aa_seqs, str):

            # 把 unusual amino acids 的 position 按顺序排列好方便之后处理
            self.unusual_aas.sort(key=lambda x: x['position'])
            self.AAs_seqs.append(_init_single_aa_seq(self.aa_seqs, self.unusual_aas, self.idx))
        elif isinstance(self.aa_seqs, list):
            for seq, unusual_aas_in_this_seq in zip(self.aa_seqs, self.unusual_aas):
                unusual_aas_in_this_seq.sort(key=lambda x: x['position'])
                self.AAs_seqs.append(_init_single_aa_seq(seq, unusual_aas_in_this_seq, self.idx))
        else:
            print('aa_seqs should be str or list')
            exit(1)

    def link_main_chain(self):
        """
        把所有的AAs_seqs中的氨基酸都连接起来成为一个大的mol对象
        """
        for AAs_seq in self.AAs_seqs:
            peptide_mol = AAs_seq[0].rwmol.GetMol()
            for AA in AAs_seq[1:]:
                combined = Chem.CombineMols(peptide_mol, AA.rwmol.GetMol())
                rw_mol = Chem.RWMol(combined)

                # 记录前一个COOH的C和后一个-NH2的N
                prior_C_dummy_atom = later_N_dummy_atom = None
                C_min_residue = float('inf')
                N_max_residue = 0

                # 从 N 端到 C 端连接 peptide
                for atom in rw_mol.GetAtoms():

                    # 找到C端C原子的 dummy atom
                    if atom.GetSymbol() == '*' and atom.HasProp('C_dummy'):
                        if int(atom.GetProp('residue')) <= C_min_residue:
                            C_min_residue = int(atom.GetProp('residue'))
                            prior_C_dummy_atom = atom

                    # 找到N端N原子的 dummy atom
                    if atom.GetSymbol() == '*' and atom.HasProp('N_dummy'):
                        if int(atom.GetProp('residue')) >= N_max_residue:
                            N_max_residue = int(atom.GetProp('residue'))
                            later_N_dummy_atom = atom

                # 添加肽键
                if prior_C_dummy_atom is not None and later_N_dummy_atom is not None:
                    rw_mol.AddBond(prior_C_dummy_atom.GetNeighbors()[0].GetIdx(), later_N_dummy_atom.GetNeighbors()[0].GetIdx(), Chem.BondType.SINGLE)

                    # 移除哑原子, 按照降序删除，否则原子序号会出错
                    index_to_remove = [prior_C_dummy_atom.GetIdx(), later_N_dummy_atom.GetIdx()]
                    for idx_to_remove in sorted(index_to_remove, reverse=True):
                        rw_mol.RemoveAtom(idx_to_remove)
                    peptide_mol = rw_mol.GetMol()

            self.main_chain_linked_mols.append(peptide_mol)

    @staticmethod
    def judge_side_COOH(rw_main_chain_linked_mol:Chem.RWMol, atom_idx:int) -> bool:
        """
        判断一个原子是不是 sidechain 上的 -COOH 上的 C 原子
        :param rw_main_chain_linked_mol: 需要被检查的可编辑分子
        :param atom_idx: 要检查的原子 idx
        :return: bool
        """
        atom = rw_main_chain_linked_mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == "C" and atom.GetDegree() == 3 and atom.GetExplicitValence() == 4 and not atom.HasProp('C_terminal'):
            O_count = 0  # 氧原子计数
            double_bond_to_oxygen = False  # 是否找到双键 O
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    O_count += 1
                    neighbor_bond = rw_main_chain_linked_mol.GetBondBetweenAtoms(atom.GetIdx(),neighbor.GetIdx())
                    if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        double_bond_to_oxygen = True
                    # 记录应该被去掉的那一个 -OH
                    if neighbor_bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        oxygen_idx_to_remove = neighbor.GetIdx()

            # 确保有两个氧邻居，且一个氧通过双键连接
            if O_count == 2 and double_bond_to_oxygen:
                return True
            else:
                return False

    def link_intrachain_bonds(self):

        def link_one_chain_intrachain_bonds(rw_main_chain_linked_mol, bonds, seq_len):
            """
            连接一条单独的多肽链的 intrachainBonds
            :param rw_main_chain_linked_mol: main chain 已经连接好的可编辑分子对象
            :param bonds: DBAASP 中原始数据格式的bonds
            :param seq_len: 正在处理的这条序列的长度
            :return: 所有 intrachainBonds 连接好的分子对象
            """
            def correct_bonds_positions(position1, position2, seq_len):
                """
                有可能会有 position 是 0 或者超出 seq_len 范围的情况，修正这一点
                :param position1: bond position1
                :param position2: bond position2
                :param seq_len: length of this sequence
                :return: corrected positions
                """
                # 记录哪些位置被修改过
                corrected_positions = {}
                # 排序获得最小值和最大值
                position1, position2 = sorted([position1, position2])
                if position1 == 0:
                    corrected_positions['position1'] = position1
                    position1 = 1
                if position2 > seq_len:
                    corrected_positions['position2'] = position2
                    position2 = seq_len
                return position1, position2, corrected_positions

            def visualize_to_be_linked_atoms(mol, idx1, idx2):
                from matplotlib import pyplot as plt
                from rdkit.Chem import Draw

                img = Draw.MolToImage(mol, size=(5000, 4000),highlightAtoms=[idx1, idx2])

                # 假设 img 是 PIL.Image.Image 对象
                # 修改图片的大小为一个大尺寸
                fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）

                # 设置黑色背景
                fig.patch.set_facecolor('black')

                # 显示图片
                plt.imshow(img)
                plt.axis('off')  # 关闭坐标轴
                plt.show()

            def DSB(rw_main_chain_linked_mol, position1, position2, chainParticipating=None, seq_len:int = None):
                """
                link Disulfide Bond, -S-S-
                :param rw_main_chain_linked_mol: 可编辑主链已连接分子
                :param position1: 二硫键第一个位置
                :param position2: 二硫键第二个位置
                :param chainParticipating: str, 'SSB' or 'SMB' or 'MMB'
                :param seq_len: 正在处理的这条序列的长度
                :return: RWMol, -S-S- 已连接好
                """

                position1, position2, _ = correct_bonds_positions(position1, position2, seq_len)

                # 初始化两个要被连接起来的原子的 Idx
                link_atom_1_index = None
                link_atom_2_index = None
                for atom in rw_main_chain_linked_mol.GetAtoms():
                    if atom.GetSymbol() == 'S' and int(atom.GetProp('residue')) == position1:
                        link_atom_1_index = atom.GetIdx()
                    elif atom.GetSymbol() == 'S' and int(atom.GetProp('residue')) == position2:
                        link_atom_2_index = atom.GetIdx()
                        # 找到第二个 S 原子之后就退出寻找 S 的循环
                        break

                # 不能在同一个原子上进行连接
                if link_atom_1_index is not None and link_atom_2_index is not None:
                    # 找到的两个 S 原子之间添加单键
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, link_atom_2_index, Chem.BondType.SINGLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(f'Peptide {self.idx}: Disulfide Bond no S atom found on |position_1: {position1}| AND |position_2: {position2}|, please check')
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(f'Peptide {self.idx}: Disulfide Bond no S atom found on |position_1: {position1}|, please check')
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(f'Peptide {self.idx}: Disulfide Bond no S atom found on |position_2: {position2}|, please check')

                return rw_main_chain_linked_mol

            def AMD(rw_main_chain_linked_mol, position1:int, position2:int, chainParticipating:str=None, seq_len:int=None):
                """
                任何羧基和N形成的单键，但其实有双键的性质，R-(O=)C-N(-R)-R
                :param chainParticipating: str, 'SSB' or 'SMB' or 'MMB'
                :param seq_len: 正在处理的这条序列的长度
                :return: RWMol, Amide 已连接好
                """

                position1, position2, _ = correct_bonds_positions(position1, position2, seq_len)

                link_atom_1_index = None
                link_atom_2_index = None
                dummy_atoms_idxs_to_remove = []
                oxygen_idx_to_remove = None
                positions = [position1, position2]  # 用于防止在同一个氨基酸中既匹配 C 又匹配 N

                # Mainchain-Mainchain Bond
                if chainParticipating == 'MMB':
                    # 第一遍循环先找 C 端
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        # Mainchain-Mainchain Bond 的 Amide 的 C 端一定是有 * 相连的
                        if atom.HasProp('C_terminal'):
                            if int(atom.GetProp('residue')) in positions:
                                # 预先记录其上的 dummy atom idx 便于之后删除
                                for neighbor in atom.GetNeighbors():
                                    if neighbor.GetSymbol() == '*':
                                        dummy_atoms_idxs_to_remove.append(neighbor.GetIdx())
                                        # rw_main_chain_linked_mol.RemoveAtom(neighbor.GetIdx())
                                        # 现在可以获得要连接 C 端原子的 Idx
                                        if int(atom.GetProp('residue')) == position1:
                                            # 匹配到一个之后就把其从 positions 中去掉防止重复匹配
                                            positions.remove(position1)
                                            link_atom_1_index = atom.GetIdx()
                                        else:
                                            positions.remove(position2)
                                            link_atom_2_index = atom.GetIdx()
                                        break

                                if neighbor.GetSymbol() == '*':
                                    break
                    # 第二遍循环找 N 端
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        if atom.HasProp('N_terminal'):
                            if int(atom.GetProp('residue')) in positions:
                                # 预先记录其上的 dummy atom idx 便于之后删除
                                for neighbor in atom.GetNeighbors():
                                    if neighbor.GetSymbol() == '*':
                                        dummy_atoms_idxs_to_remove.append(neighbor.GetIdx())
                                        break

                                # N 端原子不一定要有 * 相连接
                                if int(atom.GetProp('residue')) == position1:
                                    # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                                    positions.remove(position1)
                                    link_atom_1_index = atom.GetIdx()
                                else:
                                    positions.remove(position2)
                                    link_atom_2_index = atom.GetIdx()

                # Sidechain-Sidechain Bond
                if chainParticipating == 'SSB':
                    C_found = False
                    N_found_min_valency: int = 10  # 随便初始化一个比较大的值, 一个氨基酸侧链上可能有很多 N 原子，选配位数最少的那个
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        
                        if C_found:
                            break
                        
                        if int(atom.GetProp('residue')) in positions:
                            # 侧脸羧基 -COOH 的条件，度为3，配位数为4，并且没有 C_terminal 属性，同时还要求两个 O 邻居，其中一个是双键
                            if atom.GetSymbol() == "C" and atom.GetDegree() == 3 and atom.GetExplicitValence() == 4 and not atom.HasProp('C_terminal') and not C_found:
                                O_count = 0  # 氧原子计数
                                double_bond_to_oxygen = False  # 是否找到双键 O
                                for neighbor in atom.GetNeighbors():
                                    if neighbor.GetSymbol() == 'O':
                                        O_count += 1
                                        neighbor_bond = rw_main_chain_linked_mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                                        if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                            double_bond_to_oxygen = True
                                        # 记录应该被去掉的那一个 -OH
                                        if neighbor_bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                            oxygen_idx_to_remove = neighbor.GetIdx()

                                # 确保有两个氧邻居，且一个氧通过双键连接
                                if O_count == 2 and double_bond_to_oxygen:
                                    C_found = True
                                    dummy_atoms_idxs_to_remove.append(oxygen_idx_to_remove)
                                    if int(atom.GetProp('residue')) == position1:
                                        # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                                        positions.remove(position1)
                                        link_atom_1_index = atom.GetIdx()
                                    else:
                                        positions.remove(position2)
                                        link_atom_2_index = atom.GetIdx()

                                # 如果不是符合标准的 -COOH C 原子，则清零要移除的 C
                                else:
                                    oxygen_idx_to_remove = None

                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        if int(atom.GetProp('residue')) in positions:
                            # 寻找 N 原子，标准是不能在主链上且配位数要小于3
                            if atom.GetSymbol() == "N" and atom.GetExplicitValence() < 3 and not atom.HasProp('N_terminal'):
                                # 一个氨基酸侧链上可能有很多 N 原子，选配位数最少的那个
                                if atom.GetExplicitValence() < N_found_min_valency:
                                    N_found_min_valency = atom.GetExplicitValence()
                                    if int(atom.GetProp('residue')) == position1:
                                        link_atom_1_index = atom.GetIdx()
                                    else:
                                        link_atom_2_index = atom.GetIdx()

                # Sidechain-Mainchain Bond
                if chainParticipating == 'SMB':
                    # 分几种情况：
                    # 先找 C_dummy, 找到了的话就去找另一个position上的侧链 N
                    # 如果上面失败了，没有找到 C_dummy 或者没有找到另一个 position 上的侧链 N，那么重新开始寻找侧链上的 -COOH，去连接主链上的 N

                    C_dummy_found = False
                    side_N_found = False
                    N_found_min_valency: int = 10
                    # 找 C_dummy
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        if int(atom.GetProp('residue')) in positions:
                            if atom.HasProp('C_dummy'):
                                dummy_atoms_idxs_to_remove.append(atom.GetIdx())
                                if int(atom.GetProp('residue')) == position1:
                                    positions.remove(position1)
                                    # 这里 atom 是 dummy atom，要找唯一的邻居才是 C_terminal
                                    link_atom_1_index = atom.GetNeighbors()[0].GetIdx()
                                else:
                                    positions.remove(position2)
                                    link_atom_2_index = atom.GetNeighbors()[0].GetIdx()
                                C_dummy_found = True
                                break
                    # 找侧链 N，只有在找到 C_dummy时才这么做
                    if C_dummy_found:
                        for atom in rw_main_chain_linked_mol.GetAtoms():
                            if int(atom.GetProp('residue')) in positions:
                                # 寻找 N 原子，标准是不能在主链上且配位数要小于3
                                if atom.GetSymbol() == "N" and atom.GetExplicitValence() < 3 and not atom.HasProp('N_terminal'):
                                    side_N_found = True
                                    # 一个氨基酸侧链上可能有很多 N 原子，选配位数最少的那个
                                    if atom.GetExplicitValence() < N_found_min_valency:
                                        N_found_min_valency = atom.GetExplicitValence()
                                        if int(atom.GetProp('residue')) == position1:
                                            link_atom_1_index = atom.GetIdx()
                                        else:
                                            link_atom_2_index = atom.GetIdx()
                    # 如果没有找到侧链 N，那么要从头重新寻找 侧链 -COOH
                    if not side_N_found:
                        # 复原 posisions 和 link_atom_1_index，link_atom_2_index, dummy_atoms_idxs_to_remove
                        positions = [position1, position2]
                        link_atom_1_index = None
                        link_atom_2_index = None
                        dummy_atoms_idxs_to_remove = []

                        # 初始化 flag
                        side_C_found = False
                        for atom in rw_main_chain_linked_mol.GetAtoms():
                            # 寻找侧链 -COOH
                            if int(atom.GetProp('residue')) in positions:
                                if atom.GetSymbol() == "C" and atom.GetDegree() == 3 and atom.GetExplicitValence() == 4 and not atom.HasProp('C_terminal') and not side_C_found:
                                    O_count = 0  # 氧原子计数
                                    double_bond_to_oxygen = False  # 是否找到双键 O
                                    for neighbor in atom.GetNeighbors():
                                        if neighbor.GetSymbol() == 'O':
                                            O_count += 1
                                            neighbor_bond = rw_main_chain_linked_mol.GetBondBetweenAtoms(atom.GetIdx(),neighbor.GetIdx())
                                            if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                                double_bond_to_oxygen = True
                                            # 记录应该被去掉的那一个 -OH
                                            if neighbor_bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                                oxygen_idx_to_remove = neighbor.GetIdx()

                                    # 确保有两个氧邻居，且一个氧通过双键连接
                                    if O_count == 2 and double_bond_to_oxygen:
                                        side_C_found = True
                                        dummy_atoms_idxs_to_remove.append(oxygen_idx_to_remove)
                                        if int(atom.GetProp('residue')) == position1:
                                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                                            positions.remove(position1)
                                            link_atom_1_index = atom.GetIdx()
                                        else:
                                            positions.remove(position2)
                                            link_atom_2_index = atom.GetIdx()
                                        break

                                    # 如果不是符合标准的 -COOH C 原子，则清零要移除的 C
                                    else:
                                        oxygen_idx_to_remove = None

                        # 找到侧链 C 之后找主链的 N
                        if side_C_found:
                            for atom in rw_main_chain_linked_mol.GetAtoms():
                                if atom.HasProp('N_terminal'):
                                    if int(atom.GetProp('residue')) in positions:
                                        # 预先记录其上的 dummy atom idx 便于之后删除
                                        for neighbor in atom.GetNeighbors():
                                            if neighbor.GetSymbol() == '*':
                                                dummy_atoms_idxs_to_remove.append(neighbor.GetIdx())
                                                break

                                        # N 端原子不一定要有 * 相连接
                                        if int(atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                                            positions.remove(position1)
                                            link_atom_1_index = atom.GetIdx()
                                        else:
                                            positions.remove(position2)
                                            link_atom_2_index = atom.GetIdx()
                                        break

                # 所有情况都处理完之后连接两个需要连接的原子
                if link_atom_1_index is not None and link_atom_2_index is not None:
                    # 找到的两个 S 原子之间添加单键
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, link_atom_2_index, Chem.BondType.SINGLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(f'Peptide {self.idx}: Amide Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}| AND |position_2: {position2}|, please check')
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(f'Peptide {self.idx}: Amide Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_2_index).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}|, please check')
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(f'Peptide {self.idx}: Amide Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_1_index).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}|, please check')

                # 降序删除多余的 dummy atoms，否则会出错 Range Error
                if len(dummy_atoms_idxs_to_remove) > 0:
                    for dummy_atom_idx in sorted(dummy_atoms_idxs_to_remove, reverse=True):
                        rw_main_chain_linked_mol.RemoveAtom(dummy_atom_idx)

                # 去掉连接侧链 -COOH 应该删除的 -OH
                # if oxygen_idx_to_remove is not None:
                #     rw_main_chain_linked_mol.RemoveAtom(oxygen_idx_to_remove)

                return rw_main_chain_linked_mol

            def TIE(rw_main_chain_linked_mol, position1:int, position2:int, chainParticipating:str=None, seq_len:int=None):
                """
                一般结构形式为 R–S–R’，由 cysteine 和 dehydrated serine 生成，serine 上有一个 -OH
                :param rw_main_chain_linked_mol: 可编辑主链已连接分子
                :param position1: TIE第一个位置
                :param position2: TIE第二个位置
                :param chainParticipating: str, 'SSB' or 'SMB' or 'MMB'
                :param seq_len: 正在处理的这条序列的长度
                :return: RWMol, -R-S-R'- 已连接好的可编辑分子
                """
                # 纠正可能的 position 错误
                position1, position2, correct_positions = correct_bonds_positions(position1, position2, seq_len)

                # 初始化两个要被连接起来的原子的 Idx
                link_atom_1_index = None
                link_atom_2_index = None
                positions = [position1, position2]
                # 可能这里并不会使用
                dummy_atoms_idxs_to_remove = []
                # 先找 S 原子
                for atom in rw_main_chain_linked_mol.GetAtoms():
                    if atom.GetSymbol() == 'S' and int(atom.GetProp('residue')) in positions:
                        if int(atom.GetProp('residue')) == position1:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            link_atom_1_index = atom.GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = atom.GetIdx()
                        break

                # 找到 S 再找 -OH，要尽量排除 -COOH 上的 -OH 和 =O
                hydroxyl_idx_list =[]
                side_COOH_OH_idx_list = []
                oxygen_idx_to_remove = None
                for atom in rw_main_chain_linked_mol.GetAtoms():
                    # 注意这里一定要是 atom.GetExplicitValence() == 1 而不是 degree=1，因为双键O也是只有一个邻居
                    if atom.GetSymbol() == 'O' and atom.GetExplicitValence() == 1 and int(atom.GetProp('residue')) in positions:
                        hydroxyl_idx_list.append(atom.GetIdx())

                # 如果不止一个 -OH ，那可能有 -COOH 或者不止一个 -OH, 尝试去掉其中的 -COOH 中的 -OH
                if len(hydroxyl_idx_list) > 1:
                    # 利用deepcopy进行循环这样就不会改变原始的 list 导致循环出错
                    for hydroxyl_idx in copy.deepcopy(hydroxyl_idx_list):
                        oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(hydroxyl_idx)
                        oxygen_neighbor_atom = oxygen_atom.GetNeighbors()[0]
                        # 检查这个邻居是否是 -COOH 中的 C 原子
                        if oxygen_neighbor_atom.GetSymbol() == "C" and oxygen_neighbor_atom.GetDegree() == 3 and oxygen_neighbor_atom.GetExplicitValence() == 4 and not oxygen_neighbor_atom.HasProp('C_terminal'):
                            O_count = 0  # 氧原子计数
                            double_bond_to_oxygen = False  # 是否找到双键 O
                            for neighbor in oxygen_neighbor_atom.GetNeighbors():
                                if neighbor.GetSymbol() == 'O':
                                    O_count += 1
                                    neighbor_bond = rw_main_chain_linked_mol.GetBondBetweenAtoms(oxygen_neighbor_atom.GetIdx(), neighbor.GetIdx())
                                    if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                        double_bond_to_oxygen = True

                            # 确保有两个氧邻居，且一个氧通过双键连接，那就是找到了
                            if O_count == 2 and double_bond_to_oxygen:
                                # 去掉那些和 —COOH 相连的 -OH
                                hydroxyl_idx_list.remove(hydroxyl_idx)
                                side_COOH_OH_idx_list.append(hydroxyl_idx)


                # 尝试去掉 -COOH 中的 -OH 之后如果还有不止一个 -OH，那就随便选择第一个作为侧链 -OH
                if len(hydroxyl_idx_list) >= 1:
                    oxygen_idx_to_remove = hydroxyl_idx_list[0]
                    # 获得这个 O 原子对象
                    oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                    if int(oxygen_atom.GetProp('residue')) == position1:
                        # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                        positions.remove(position1)
                        # 这个 -OH 是脱水过程中需要被去掉的部分，他的邻居才是要链接的原子
                        link_atom_1_index = oxygen_atom.GetNeighbors()[0].GetIdx()
                    else:
                        positions.remove(position2)
                        link_atom_2_index = oxygen_atom.GetNeighbors()[0].GetIdx()

                # 如果没有单独的 -OH 被找到，那可能是在 -COOH 上
                if len(hydroxyl_idx_list) == 0 and len(side_COOH_OH_idx_list) > 0:
                    oxygen_idx_to_remove = side_COOH_OH_idx_list[0]
                    # 获得这个 O 原子对象
                    oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                    if int(oxygen_atom.GetProp('residue')) == position1:
                        # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                        positions.remove(position1)
                        # 这个 -OH 是脱水过程中需要被去掉的部分，他的邻居才是要链接的原子
                        link_atom_1_index = oxygen_atom.GetNeighbors()[0].GetIdx()
                    else:
                        positions.remove(position2)
                        link_atom_2_index = oxygen_atom.GetNeighbors()[0].GetIdx()

                # 如果连 -COOH 上的 -OH 也没有找到，那可能是 SMB，找 C dummy atom
                if len(hydroxyl_idx_list) == 0 and len(side_COOH_OH_idx_list) == 0:
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        # 注意这里一定要是 atom.GetExplicitValence() == 1 而不是 degree=1，因为双键O也是只有一个邻居
                        if atom.GetSymbol() == '*' and atom.HasProp('C_dummy') and int(atom.GetProp('residue')) in positions:
                            dummy_atoms_idxs_to_remove.append(atom.GetIdx())
                            if int(atom.GetProp('residue')) == position1:
                                # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                                positions.remove(position1)
                                # 这个 -OH 是脱水过程中需要被去掉的部分，他的邻居才是要链接的原子
                                link_atom_1_index = atom.GetNeighbors()[0].GetIdx()
                            else:
                                positions.remove(position2)
                                link_atom_2_index = atom.GetNeighbors()[0].GetIdx()
                            break


                # 可视化 debug 代码，运行时记得注释掉
                # if link_atom_1_index is not None and link_atom_2_index is not None:
                #     visualize_to_be_linked_atoms(rw_main_chain_linked_mol, link_atom_1_index, link_atom_2_index)

                # 连接新键并检查是否有错
                # 所有情况都处理完之后连接两个需要连接的原子
                if link_atom_1_index is not None and link_atom_2_index is not None:
                    # 找到的两个 S 原子之间添加单键
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, link_atom_2_index, Chem.BondType.SINGLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(f"Peptide {self.idx}: Thioether Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}，original {correct_positions.get('position1')}| AND |position_2: {position2}，original {correct_positions.get('position2')}|, please check")
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(f"Peptide {self.idx}: Thioether Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_2_index).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}, original {correct_positions.get('position1')}|, please check")
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(f"Peptide {self.idx}: Thioether Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_1_index).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}, original {correct_positions.get('position2')}|, please check")

                if oxygen_idx_to_remove is not None:
                    dummy_atoms_idxs_to_remove.append(oxygen_idx_to_remove)
                # 降序删除多余的 dummy atoms，否则会出错 Range Error
                if len(dummy_atoms_idxs_to_remove) > 0:
                    for dummy_atom_idx in sorted(dummy_atoms_idxs_to_remove, reverse=True):
                        rw_main_chain_linked_mol.RemoveAtom(dummy_atom_idx)

                # 去掉连接侧链 -COOH 应该删除的 -OH
                # if oxygen_idx_to_remove is not None:
                #     rw_main_chain_linked_mol.RemoveAtom(oxygen_idx_to_remove)

                return rw_main_chain_linked_mol

            def DCB(rw_main_chain_linked_mol, position1:int, position2:int, chainParticipating:str=None, seq_len:int=None):
                """
                Dicarbon bond
                目前看来大部分的双键都是靠 -C=C 和 C=C- 形成的，变成 -C=C-, 失去一个 -C=C-
                :param rw_main_chain_linked_mol:
                :param position1:
                :param position2:
                :param chainParticipating:
                :param seq_len:
                :return:
                """

                position1, position2, correct_positions = correct_bonds_positions(position1, position2, seq_len)

                # 初始化两个要被连接起来的原子的 Idx
                link_atom_1_index = None
                link_atom_2_index = None
                positions = [position1, position2]

                carbon_idx_to_remove = []

                for atom in rw_main_chain_linked_mol.GetAtoms():
                    # 判断最外面双键 =C 的规则很简单：只有一个邻居，并且 valence 是2
                    if atom.GetSymbol() == 'C' and atom.GetDegree() == 1 and atom.GetExplicitValence() == 2 and int(atom.GetProp('residue')) in positions:
                        if int(atom.GetProp('residue')) == position1:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            # 这个 =C 是脱水过程中需要被去掉的部分，他的邻居才是要链接的原子
                            link_atom_1_index = atom.GetNeighbors()[0].GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = atom.GetNeighbors()[0].GetIdx()

                        # 这个 =C 需要被移除
                        carbon_idx_to_remove.append(atom.GetIdx())

                        # 两个要连接的原子都找到之后就直接退出循环
                        if link_atom_1_index is not None and link_atom_2_index is not None:
                            break

                # 连接新键并检查是否有错
                # 所有情况都处理完之后连接两个需要连接的原子
                if link_atom_1_index is not None and link_atom_2_index is not None:
                    # 找到的两个 =C 原子之间添加双键 =
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, link_atom_2_index,Chem.BondType.DOUBLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Dicarbon Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}，original {correct_positions.get('position1')}| AND |position_2: {position2}，original {correct_positions.get('position2')}|, please check")
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Dicarbon Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_2_index).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}, original {correct_positions.get('position1')}|, please check")
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Dicarbon Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_1_index).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}, original {correct_positions.get('position2')}|, please check")

                # 这都是画图用的部分
                # img = Draw.MolToImage(rw_main_chain_linked_mol, size=(5000, 4000))

                # 假设 img 是 PIL.Image.Image 对象
                # 修改图片的大小为一个大尺寸
                # fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）

                # 设置黑色背景
                # fig.patch.set_facecolor('black')

                # 显示图片
                # plt.imshow(img)
                # plt.axis('off')  # 关闭坐标轴
                # plt.show()

                # 降序删除多余的 =C atoms，否则会出错 Range Error
                if len(carbon_idx_to_remove) > 0:
                    for carbon_atom_idx in sorted(carbon_idx_to_remove, reverse=True):
                        rw_main_chain_linked_mol.RemoveAtom(carbon_atom_idx)

                return rw_main_chain_linked_mol


            def EST(rw_main_chain_linked_mol, position1:int, position2:int, chainParticipating:str=None, seq_len:int=None):
                """
                Ester bond, 由酸 -COOH 和醇 -OH 生成
                DBAASP 中全部都是 Side-Main 或者 Main-Main Bond
                :param rw_main_chain_linked_mol:
                :param position1:
                :param position2:
                :param chainParticipating:
                :param seq_len:
                :return:
                """
                position1, position2, correct_positions = correct_bonds_positions(position1, position2, seq_len)

                # 初始化两个要被连接起来的原子的 Idx
                link_atom_1_index = None
                link_atom_2_index = None
                positions = [position1, position2]
                dummy_atoms_idxs_to_remove = []

                if chainParticipating == 'SMB' or chainParticipating == 'MMB':
                    # 先找 C_dummy, 找到了的话就去找另一个 position 上的侧链 -OH

                    C_dummy_found = False
                    side_OH_found = False
                    # 找 C_dummy
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        if int(atom.GetProp('residue')) in positions:
                            if atom.HasProp('C_dummy'):
                                dummy_atoms_idxs_to_remove.append(atom.GetIdx())
                                if int(atom.GetProp('residue')) == position1:
                                    positions.remove(position1)
                                    # 这里 atom 是 dummy atom，要找唯一的邻居才是 C_terminal
                                    link_atom_1_index = atom.GetNeighbors()[0].GetIdx()
                                else:
                                    positions.remove(position2)
                                    link_atom_2_index = atom.GetNeighbors()[0].GetIdx()
                                C_dummy_found = True
                                break

                    # 找侧链的 -OH
                    # 找到 C_dummy 再找 -OH，要尽量排除 -COOH 上的 -OH 和 =O
                    hydroxyl_idx_list = []
                    side_COOH_OH_idx_list = []
                    oxygen_idx_to_remove = None
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        # 注意这里一定要是 atom.GetExplicitValence() == 1 而不是 degree=1，因为双键O也是只有一个邻居
                        if atom.GetSymbol() == 'O' and atom.GetExplicitValence() == 1 and int(atom.GetProp('residue')) in positions:
                            hydroxyl_idx_list.append(atom.GetIdx())

                    # 如果不止一个 -OH ，那可能有 -COOH 或者不止一个 -OH, 尝试去掉其中的 -COOH 中的 -OH
                    if len(hydroxyl_idx_list) > 1:
                        # 利用deepcopy进行循环这样就不会改变原始的 list 导致循环出错
                        for hydroxyl_idx in copy.deepcopy(hydroxyl_idx_list):
                            oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(hydroxyl_idx)
                            oxygen_neighbor_atom = oxygen_atom.GetNeighbors()[0]
                            # 检查这个邻居是否是 -COOH 中的 C 原子
                            if oxygen_neighbor_atom.GetSymbol() == "C" and oxygen_neighbor_atom.GetDegree() == 3 and oxygen_neighbor_atom.GetExplicitValence() == 4 and not oxygen_neighbor_atom.HasProp('C_terminal'):
                                O_count = 0  # 氧原子计数
                                double_bond_to_oxygen = False  # 是否找到双键 O
                                for neighbor in oxygen_neighbor_atom.GetNeighbors():
                                    if neighbor.GetSymbol() == 'O':
                                        O_count += 1
                                        neighbor_bond = rw_main_chain_linked_mol.GetBondBetweenAtoms(oxygen_neighbor_atom.GetIdx(), neighbor.GetIdx())
                                        if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                            double_bond_to_oxygen = True

                                # 确保有两个氧邻居，且一个氧通过双键连接，那就是找到了
                                if O_count == 2 and double_bond_to_oxygen:
                                    # 去掉那些和 —COOH 相连的 -OH
                                    hydroxyl_idx_list.remove(hydroxyl_idx)
                                    side_COOH_OH_idx_list.append(hydroxyl_idx)

                    # 尝试去掉 -COOH 中的 -OH 之后如果还有不止一个 -OH，那就随便选择第一个作为侧链 -OH
                    if len(hydroxyl_idx_list) >= 1:
                        # 这个侧链上的 -OH 形成 ester 时不用被去掉
                        oxygen_idx_to_remove = hydroxyl_idx_list[0]
                        # 获得这个 O 原子对象
                        oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                        if int(oxygen_atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            # 这个 -OH 脱水过程中不用被去掉
                            link_atom_1_index = oxygen_atom.GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = oxygen_atom.GetIdx()

                    # 如果没有单独的 -OH 被找到，那可能是在 -COOH 上
                    if len(hydroxyl_idx_list) == 0 and len(side_COOH_OH_idx_list) > 0:
                        oxygen_idx_to_remove = side_COOH_OH_idx_list[0]
                        # 获得这个 O 原子对象
                        oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                        if int(oxygen_atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            # 这个 -OH 脱水过程中不用被去掉
                            link_atom_1_index = oxygen_atom.GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = oxygen_atom.GetIdx()

                # if chainParticipating == 'MMB':
                #     print('stop')
                #     # 这都是画图用的部分
                #     img = Draw.MolToImage(rw_main_chain_linked_mol, size=(5000, 4000))
                #     fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）
                #     fig.patch.set_facecolor('black')  # 设置黑色背景
                #     plt.imshow(img)  # 显示图片
                #     plt.axis('off')  # 关闭坐标轴
                #     plt.show()

                # 连接新键并检查是否有错
                # 所有情况都处理完之后连接两个需要连接的原子
                if link_atom_1_index is not None and link_atom_2_index is not None:
                    # 找到的两个 =C 原子之间添加双键 =
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, link_atom_2_index, Chem.BondType.SINGLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Ester Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}，original {correct_positions.get('position1')}| AND |position_2: {position2}，original {correct_positions.get('position2')}|, please check")
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Ester Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_2_index).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}, original {correct_positions.get('position1')}|, please check")
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Ester Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_1_index).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}, original {correct_positions.get('position2')}|, please check")

                # 降序删除多余的 dummy atoms，否则会出错 Range Error
                if link_atom_1_index is not None and link_atom_2_index is not None:
                    if len(dummy_atoms_idxs_to_remove) > 0:
                        for dummy_atom_idx in sorted(dummy_atoms_idxs_to_remove, reverse=True):
                            rw_main_chain_linked_mol.RemoveAtom(dummy_atom_idx)

                # 这都是画图用的部分
                # img = Draw.MolToImage(rw_main_chain_linked_mol, size=(5000, 4000))
                # fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）
                # fig.patch.set_facecolor('black')  # 设置黑色背景
                # plt.imshow(img)  # 显示图片
                # plt.axis('off')  # 关闭坐标轴
                # plt.show()

                return rw_main_chain_linked_mol

            def p_XylB(rw_main_chain_linked_mol, position1: int, position2: int, chainParticipating: str = None,seq_len: int = None):
                """
                -S-C-C₆-C-S-, para-Xylene thioether bridge connecting a pair of CYS residues
                :param rw_main_chain_linked_mol:
                :param position1:
                :param position2:
                :param chainParticipating:
                :param seq_len:
                :return:
                """
                position1, position2, correct_positions = correct_bonds_positions(position1, position2, seq_len)

                p_Xylene_smiles = "CC1=CC=C(C)C=C1"
                p_Xylene_mol = Chem.MolFromSmiles(p_Xylene_smiles)
                rw_p_Xylene_mol = Chem.RWMol(p_Xylene_mol)

                # look for -CH3
                CH3_idxs = []
                for atom in rw_p_Xylene_mol.GetAtoms():
                    # 单独添加的分子都是 residue = 0
                    atom.SetProp('residue', '0')
                    if atom.GetExplicitValence() == 1:
                        CH3_idxs.append(atom.GetIdx())
                p_Xulene_link_atom_idx_1, p_Xulene_link_atom_idx_2 = CH3_idxs

                # 初始化两个要被连接起来的原子的 Idx
                link_atom_1_index = None
                link_atom_2_index = None
                positions = [position1, position2]

                # 先找 S 原子
                for atom in rw_main_chain_linked_mol.GetAtoms():
                    if atom.GetSymbol() == 'S' and int(atom.GetProp('residue')) in positions:
                        if int(atom.GetProp('residue')) == position1:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            link_atom_1_index = atom.GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = atom.GetIdx()

                # 合并分子，这会带来 atom idx 的平移
                combined_mol = rdmolops.CombineMols(rw_p_Xylene_mol, rw_main_chain_linked_mol)
                rw_main_chain_linked_mol = Chem.RWMol(combined_mol)

                # 计算 offset 抵消平移
                offset = rw_p_Xylene_mol.GetNumAtoms()
                if link_atom_1_index is not None:
                    link_atom_1_index += offset
                if link_atom_2_index is not None:
                    link_atom_2_index += offset

                # 连接新键并检查是否有错
                # 所有情况都处理完之后连接两个需要连接的原子
                if link_atom_1_index is not None and link_atom_2_index is not None:
                    # 添加 -C-S- 单键
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, p_Xulene_link_atom_idx_1, Chem.BondType.SINGLE)
                    rw_main_chain_linked_mol.AddBond(link_atom_2_index, p_Xulene_link_atom_idx_2, Chem.BondType.SINGLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: -S-C-C₆-C-S- Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}，original {correct_positions.get('position1')}| AND |position_2: {position2}，original {correct_positions.get('position2')}|, please check")
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: -S-C-C₆-C-S- Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_2_index).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}, original {correct_positions.get('position1')}|, please check")
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: -S-C-C₆-C-S- Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_1_index).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}, original {correct_positions.get('position2')}|, please check")

                # 这都是画图用的部分
                # img = Draw.MolToImage(rw_main_chain_linked_mol, size=(5000, 4000))
                # fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）
                # fig.patch.set_facecolor('black')  # 设置黑色背景
                # plt.imshow(img)  # 显示图片
                # plt.axis('off')  # 关闭坐标轴
                # plt.show()

                return rw_main_chain_linked_mol

            def TRZB(rw_main_chain_linked_mol, position1: int, position2: int, chainParticipating: str = None,seq_len: int = None):
                """
                Sidechain-Sidechain Bond linked with Triazolic bridge
                :param rw_main_chain_linked_mol:
                :param position1:
                :param position2:
                :param chainParticipating:
                :param seq_len:
                :return:
                """

                position1, position2, correct_positions = correct_bonds_positions(position1, position2, seq_len)
                positions = [position1, position2]
                rw_main_chain_linked_mol.UpdatePropertyCache(strict=False)

                # 找 -C#C （alkyne）
                alkyne_idx_near_C = None  # 离 C 更近的炔的那个 C 原子的 idx
                alkyne_idx_far_C = None
                for atom in rw_main_chain_linked_mol.GetAtoms():
                    if int(atom.GetProp('residue')) in positions and atom.GetSymbol() == 'C' and atom.GetDegree() == 2 and atom.GetExplicitValence() == 4:
                        alkyne_idx_near_C = atom.GetIdx()
                        for neighbor in atom.GetNeighbors():
                            bond = rw_main_chain_linked_mol.GetBondBetweenAtoms(neighbor.GetIdx(), alkyne_idx_near_C)
                            if bond.GetBondType() == Chem.BondType.TRIPLE:
                                alkyne_idx_far_C = neighbor.GetIdx()
                                rw_main_chain_linked_mol.RemoveBond(alkyne_idx_near_C, alkyne_idx_far_C)  # 移除三键
                                rw_main_chain_linked_mol.AddBond(alkyne_idx_near_C, alkyne_idx_far_C, Chem.rdchem.BondType.DOUBLE)  # 添加双键
                                if int(atom.GetProp('residue')) == position1:
                                    # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                                    positions.remove(position1)
                                else:
                                    positions.remove(position2)
                                break

                rw_main_chain_linked_mol.UpdatePropertyCache(strict=False)

                # 找 -N=N#N 中间的那个 N
                N3_idx_mid = None
                N3_idx_near_C = None  # 离 C 更近的炔的那个 N 原子的 idx
                N3_idx_far_C = None
                for atom in rw_main_chain_linked_mol.GetAtoms():
                    if int(atom.GetProp('residue')) in positions and atom.GetSymbol() == 'N' and atom.GetDegree() == 2 and atom.GetExplicitValence() == 4:
                        N3_idx_mid = atom.GetIdx()
                        # 去掉电荷
                        atom.SetFormalCharge(0)
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetDegree() == 1:
                                N3_idx_far_C = neighbor.GetIdx()
                                # 去掉电荷
                                neighbor.SetFormalCharge(0)
                            if neighbor.GetDegree() == 2:
                                N3_idx_near_C = neighbor.GetIdx()
                        rw_main_chain_linked_mol.RemoveBond(N3_idx_mid, N3_idx_near_C)  # 移除双键
                        rw_main_chain_linked_mol.AddBond(N3_idx_mid, N3_idx_near_C, Chem.rdchem.BondType.SINGLE)  # 添加单键

                        rw_main_chain_linked_mol.RemoveBond(N3_idx_mid, N3_idx_far_C)  # 移除三键
                        rw_main_chain_linked_mol.AddBond(N3_idx_mid, N3_idx_far_C, Chem.rdchem.BondType.DOUBLE)  # 添加双键
                        break

                        # 连接新键并检查是否有错
                        # 所有情况都处理完之后连接两个需要连接的原子
                if None not in [alkyne_idx_near_C, alkyne_idx_far_C, N3_idx_near_C, N3_idx_far_C, N3_idx_mid]:
                    # 添加 -C-S- 单键
                    rw_main_chain_linked_mol.AddBond(N3_idx_near_C, alkyne_idx_far_C, Chem.BondType.SINGLE)
                    rw_main_chain_linked_mol.AddBond(N3_idx_far_C, alkyne_idx_near_C, Chem.BondType.SINGLE)
                elif alkyne_idx_near_C is None and None in [N3_idx_near_C, N3_idx_far_C, N3_idx_mid]:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Triazolic bridge Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}，original {correct_positions.get('position1')}| AND |position_2: {position2}，original {correct_positions.get('position2')}|, please check")
                elif alkyne_idx_near_C is None and None not in [N3_idx_near_C, N3_idx_far_C, N3_idx_mid]:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Triazolic bridge Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(N3_idx_mid).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}, original {correct_positions.get('position1')}|, please check")
                elif alkyne_idx_near_C is not None and None in [N3_idx_near_C, N3_idx_far_C, N3_idx_mid]:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Triazolic bridge Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(alkyne_idx_near_C).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}, original {correct_positions.get('position2')}|, please check")

                # 这都是画图用的部分
                # img = Draw.MolToImage(rw_main_chain_linked_mol, size=(5000, 4000))
                # fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）
                # fig.patch.set_facecolor('black')  # 设置黑色背景
                # plt.imshow(img)  # 显示图片
                # plt.axis('off')  # 关闭坐标轴
                # plt.show()

                return rw_main_chain_linked_mol

            def E_but_2_enyl_B(rw_main_chain_linked_mol, position1: int, position2: int, chainParticipating: str = None, seq_len: int = None):
                """
                -C-C=C-C-
                Sidechain-Sidechain C=C bond cross-linked with (E)-but-2-enyl
                :param rw_main_chain_linked_mol:
                :param position1:
                :param position2:
                :param chainParticipating:
                :param seq_len:
                :return:
                """

                position1, position2, correct_positions = correct_bonds_positions(position1, position2, seq_len)

                E_but_2_enyl_smiles = "CC=CC"
                E_but_2_enyl_mol = Chem.MolFromSmiles(E_but_2_enyl_smiles)
                rw_E_but_2_enyl_mol = Chem.RWMol(E_but_2_enyl_mol)

                # look for -CH3
                CH3_idxs = []
                for atom in rw_E_but_2_enyl_mol.GetAtoms():
                    # 单独添加的分子都是 residue = 0
                    atom.SetProp('residue', '0')
                    if atom.GetExplicitValence() == 1:
                        CH3_idxs.append(atom.GetIdx())
                E_but_2_enyl_link_atom_idx_1, E_but_2_enyl_link_atom_idx_2 = CH3_idxs

                # 初始化两个要被连接起来的原子的 Idx
                link_atom_1_index = None
                link_atom_2_index = None
                positions = [position1, position2]

                # 找两个 N 原子
                for atom in rw_main_chain_linked_mol.GetAtoms():
                    if atom.GetSymbol() == 'N' and int(atom.GetProp('residue')) in positions and not atom.HasProp('N_terminal'):
                        if int(atom.GetProp('residue')) == position1:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            link_atom_1_index = atom.GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = atom.GetIdx()

                    # 两个 N 都找到了就退出循环
                    if None not in [link_atom_1_index, link_atom_2_index]:
                        break

                # 合并分子，这会带来 atom idx 的平移
                combined_mol = rdmolops.CombineMols(rw_E_but_2_enyl_mol, rw_main_chain_linked_mol)
                rw_main_chain_linked_mol = Chem.RWMol(combined_mol)

                # 计算 offset 抵消平移
                offset = rw_E_but_2_enyl_mol.GetNumAtoms()
                if link_atom_1_index is not None:
                    link_atom_1_index += offset
                if link_atom_2_index is not None:
                    link_atom_2_index += offset

                # 连接新键并检查是否有错
                # 所有情况都处理完之后连接两个需要连接的原子
                if None not in [link_atom_1_index, link_atom_2_index]:
                    # 添加 -C-S- 单键
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, E_but_2_enyl_link_atom_idx_1, Chem.BondType.SINGLE)
                    rw_main_chain_linked_mol.AddBond(link_atom_2_index, E_but_2_enyl_link_atom_idx_2,Chem.BondType.SINGLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: (E)-but-2-enyl Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}，original {correct_positions.get('position1')}| AND |position_2: {position2}，original {correct_positions.get('position2')}|, please check")
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: (E)-but-2-enyl Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_2_index).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}, original {correct_positions.get('position1')}|, please check")
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: (E)-but-2-enyl Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_1_index).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}, original {correct_positions.get('position2')}|, please check")

                # 这都是画图用的部分
                # img = Draw.MolToImage(rw_main_chain_linked_mol, size=(5000, 4000))
                # fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）
                # fig.patch.set_facecolor('black')  # 设置黑色背景
                # plt.imshow(img)  # 显示图片
                # plt.axis('off')  # 关闭坐标轴
                # plt.show()

                return rw_main_chain_linked_mol

            def BisMeBn_B(rw_main_chain_linked_mol, position1: int, position2: int, chainParticipating: str = None, seq_len: int = None):
                """
                -C-C6-C-
                侧链-侧链连接是通过二亚甲基苯（Bismethylenebenzene）实现的。也就是说，这两个侧链是通过一个带有两个甲基的苯环（类似于苯桥结构）交联在一起。
                和上面Sidechain-Sidechain C=C bond cross-linked with (E)-but-2-enyl形成机制一样
                :param rw_main_chain_linked_mol:
                :param position1:
                :param position2:
                :param chainParticipating:
                :param seq_len:
                :return:
                """

                position1, position2, correct_positions = correct_bonds_positions(position1, position2, seq_len)

                E_but_2_enyl_smiles = "c1cc(C)ccc1C"
                E_but_2_enyl_mol = Chem.MolFromSmiles(E_but_2_enyl_smiles)
                rw_E_but_2_enyl_mol = Chem.RWMol(E_but_2_enyl_mol)

                # look for -CH3
                CH3_idxs = []
                for atom in rw_E_but_2_enyl_mol.GetAtoms():
                    # 单独添加的分子都是 residue = 0
                    atom.SetProp('residue', '0')
                    if atom.GetExplicitValence() == 1:
                        CH3_idxs.append(atom.GetIdx())
                E_but_2_enyl_link_atom_idx_1, E_but_2_enyl_link_atom_idx_2 = CH3_idxs

                # 初始化两个要被连接起来的原子的 Idx
                link_atom_1_index = None
                link_atom_2_index = None
                positions = [position1, position2]

                # 找两个 N 原子
                for atom in rw_main_chain_linked_mol.GetAtoms():
                    if atom.GetSymbol() == 'N' and int(atom.GetProp('residue')) in positions and not atom.HasProp('N_terminal'):
                        if int(atom.GetProp('residue')) == position1:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            link_atom_1_index = atom.GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = atom.GetIdx()

                    # 两个 N 都找到了就退出循环
                    if None not in [link_atom_1_index, link_atom_2_index]:
                        break

                # 合并分子，这会带来 atom idx 的平移
                combined_mol = rdmolops.CombineMols(rw_E_but_2_enyl_mol, rw_main_chain_linked_mol)
                rw_main_chain_linked_mol = Chem.RWMol(combined_mol)

                # 计算 offset 抵消平移
                offset = rw_E_but_2_enyl_mol.GetNumAtoms()
                if link_atom_1_index is not None:
                    link_atom_1_index += offset
                if link_atom_2_index is not None:
                    link_atom_2_index += offset

                # 连接新键并检查是否有错
                # 所有情况都处理完之后连接两个需要连接的原子
                if None not in [link_atom_1_index, link_atom_2_index]:
                    # 添加 -C-S- 单键
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, E_but_2_enyl_link_atom_idx_1, Chem.BondType.SINGLE)
                    rw_main_chain_linked_mol.AddBond(link_atom_2_index, E_but_2_enyl_link_atom_idx_2,Chem.BondType.SINGLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: BisMeBn-B Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}，original {correct_positions.get('position1')}| AND |position_2: {position2}，original {correct_positions.get('position2')}|, please check")
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: BisMeBn-B Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_2_index).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}, original {correct_positions.get('position1')}|, please check")
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: BisMeBn-B Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_1_index).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}, original {correct_positions.get('position2')}|, please check")

                # 这都是画图用的部分
                # img = Draw.MolToImage(rw_main_chain_linked_mol, size=(5000, 4000))
                # fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）
                # fig.patch.set_facecolor('black')  # 设置黑色背景
                # plt.imshow(img)  # 显示图片
                # plt.axis('off')  # 关闭坐标轴
                # plt.show()

                return rw_main_chain_linked_mol

            def but_2_ynyl_B(rw_main_chain_linked_mol, position1: int, position2: int, chainParticipating: str = None, seq_len: int = None):
                """
                -C-C6-C-
                侧链-侧链连接是通过二亚甲基苯（Bismethylenebenzene）实现的。也就是说，这两个侧链是通过一个带有两个甲基的苯环（类似于苯桥结构）交联在一起。
                和上面Sidechain-Sidechain C=C bond cross-linked with (E)-but-2-enyl形成机制一样
                :param rw_main_chain_linked_mol:
                :param position1:
                :param position2:
                :param chainParticipating:
                :param seq_len:
                :return:
                """

                position1, position2, correct_positions = correct_bonds_positions(position1, position2, seq_len)

                E_but_2_enyl_smiles = "CC#CC"
                E_but_2_enyl_mol = Chem.MolFromSmiles(E_but_2_enyl_smiles)
                rw_E_but_2_enyl_mol = Chem.RWMol(E_but_2_enyl_mol)

                # look for -CH3
                CH3_idxs = []
                for atom in rw_E_but_2_enyl_mol.GetAtoms():
                    # 单独添加的分子都是 residue = 0
                    atom.SetProp('residue', '0')
                    if atom.GetExplicitValence() == 1:
                        CH3_idxs.append(atom.GetIdx())
                E_but_2_enyl_link_atom_idx_1, E_but_2_enyl_link_atom_idx_2 = CH3_idxs

                # 初始化两个要被连接起来的原子的 Idx
                link_atom_1_index = None
                link_atom_2_index = None
                positions = [position1, position2]

                # 找两个 N 原子
                for atom in rw_main_chain_linked_mol.GetAtoms():
                    if atom.GetSymbol() == 'N' and int(atom.GetProp('residue')) in positions and not atom.HasProp('N_terminal'):
                        if int(atom.GetProp('residue')) == position1:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            link_atom_1_index = atom.GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = atom.GetIdx()

                    # 两个 N 都找到了就退出循环
                    if None not in [link_atom_1_index, link_atom_2_index]:
                        break

                # 合并分子，这会带来 atom idx 的平移
                combined_mol = rdmolops.CombineMols(rw_E_but_2_enyl_mol, rw_main_chain_linked_mol)
                rw_main_chain_linked_mol = Chem.RWMol(combined_mol)

                # 计算 offset 抵消平移
                offset = rw_E_but_2_enyl_mol.GetNumAtoms()
                if link_atom_1_index is not None:
                    link_atom_1_index += offset
                if link_atom_2_index is not None:
                    link_atom_2_index += offset

                # 连接新键并检查是否有错
                # 所有情况都处理完之后连接两个需要连接的原子
                if None not in [link_atom_1_index, link_atom_2_index]:
                    # 添加 -C-S- 单键
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, E_but_2_enyl_link_atom_idx_1, Chem.BondType.SINGLE)
                    rw_main_chain_linked_mol.AddBond(link_atom_2_index, E_but_2_enyl_link_atom_idx_2,Chem.BondType.SINGLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: but_2_ynyl Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}，original {correct_positions.get('position1')}| AND |position_2: {position2}，original {correct_positions.get('position2')}|, please check")
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: but_2_ynyl Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_2_index).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}, original {correct_positions.get('position1')}|, please check")
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: but_2_ynyl Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_1_index).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}, original {correct_positions.get('position2')}|, please check")

                # 这都是画图用的部分
                # img = Draw.MolToImage(rw_main_chain_linked_mol, size=(5000, 4000))
                # fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）
                # fig.patch.set_facecolor('black')  # 设置黑色背景
                # plt.imshow(img)  # 显示图片
                # plt.axis('off')  # 关闭坐标轴
                # plt.show()

                return rw_main_chain_linked_mol

            def AMN(rw_main_chain_linked_mol, position1: int, position2: int, chainParticipating: str = None, seq_len: int = None):
                """
                通常通过 -OH 和 -NHn 形成
                :param rw_main_chain_linked_mol:
                :param position1:
                :param position2:
                :param chainParticipating:
                :param seq_len:
                :return:
                """

                position1, position2, correct_positions = correct_bonds_positions(position1, position2, seq_len)

                # 初始化两个要被连接起来的原子的 Idx
                link_atom_1_index = None
                link_atom_2_index = None
                positions = [position1, position2]
                dummy_atoms_idxs_to_remove = []
                oxygen_idx_to_remove = None
                N_idx_to_remove = None

                if chainParticipating == 'SSB':

                    # 先找侧链上的 N
                    N_min_valence_idx = None
                    N_min_valence = 10
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        if int(atom.GetProp('residue')) in positions and not atom.HasProp(
                                'N_terminal') and atom.GetSymbol() == 'N' and atom.GetExplicitValence() < 3:
                            if atom.GetExplicitValence() < N_min_valence:
                                N_min_valence_idx = atom.GetIdx()
                                N_min_valence = atom.GetExplicitValence()

                    # 如果找到了侧链上的 N
                    if N_min_valence_idx is not None:
                        side_N_found = True
                        atom = rw_main_chain_linked_mol.GetAtomWithIdx(N_min_valence_idx)
                        if int(atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                            positions.remove(position1)
                            link_atom_1_index = N_min_valence_idx
                        else:
                            positions.remove(position2)
                            link_atom_2_index = N_min_valence_idx

                    # 再找侧链上的 -OH
                    # 找到 N_dummy 再找 -OH，要尽量排除 -COOH 上的 -OH 和 =O
                    hydroxyl_idx_list = []
                    side_COOH_OH_idx_list = []
                    oxygen_idx_to_remove = None
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        # 注意这里一定要是 atom.GetExplicitValence() == 1 而不是 degree=1，因为双键O也是只有一个邻居
                        if atom.GetSymbol() == 'O' and atom.GetExplicitValence() == 1 and int(
                                atom.GetProp('residue')) in positions:
                            hydroxyl_idx_list.append(atom.GetIdx())

                    # 如果不止一个 -OH ，那可能有 -COOH 或者不止一个 -OH, 尝试去掉其中的 -COOH 中的 -OH
                    if len(hydroxyl_idx_list) > 1:
                        # 利用deepcopy进行循环这样就不会改变原始的 list 导致循环出错
                        for hydroxyl_idx in copy.deepcopy(hydroxyl_idx_list):
                            oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(hydroxyl_idx)
                            oxygen_neighbor_atom = oxygen_atom.GetNeighbors()[0]
                            # 检查这个邻居是否是 -COOH 中的 C 原子
                            if oxygen_neighbor_atom.GetSymbol() == "C" and oxygen_neighbor_atom.GetDegree() == 3 and oxygen_neighbor_atom.GetExplicitValence() == 4 and not oxygen_neighbor_atom.HasProp(
                                    'C_terminal'):
                                O_count = 0  # 氧原子计数
                                double_bond_to_oxygen = False  # 是否找到双键 O
                                for neighbor in oxygen_neighbor_atom.GetNeighbors():
                                    if neighbor.GetSymbol() == 'O':
                                        O_count += 1
                                        neighbor_bond = rw_main_chain_linked_mol.GetBondBetweenAtoms(
                                            oxygen_neighbor_atom.GetIdx(), neighbor.GetIdx())
                                        if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                            double_bond_to_oxygen = True

                                # 确保有两个氧邻居，且一个氧通过双键连接，那就是找到了
                                if O_count == 2 and double_bond_to_oxygen:
                                    # 去掉那些和 —COOH 相连的 -OH
                                    hydroxyl_idx_list.remove(hydroxyl_idx)
                                    side_COOH_OH_idx_list.append(hydroxyl_idx)

                    # 尝试去掉 -COOH 中的 -OH 之后如果还有不止一个 -OH，那就随便选择第一个作为侧链 -OH
                    if len(hydroxyl_idx_list) >= 1:
                        # 这个侧链上的 -OH 形成 ester 时不用被去掉
                        oxygen_idx_to_remove = hydroxyl_idx_list[0]
                        # 获得这个 O 原子对象
                        oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                        if int(oxygen_atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            # 这个 -OH 形成 Amine 的脱水过程中需要被去掉
                            link_atom_1_index = oxygen_atom.GetNeighbors()[0].GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = oxygen_atom.GetNeighbors()[0].GetIdx()

                    # 如果没有单独的 -OH 被找到，那可能是在 -COOH 上
                    if len(hydroxyl_idx_list) == 0 and len(side_COOH_OH_idx_list) > 0:
                        oxygen_idx_to_remove = side_COOH_OH_idx_list[0]
                        # 获得这个 O 原子对象
                        oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                        if int(oxygen_atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                            # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                            positions.remove(position1)
                            # 这个 -OH 形成 Amine 的脱水过程中需要被去掉
                            link_atom_1_index = oxygen_atom.GetNeighbors()[0].GetIdx()
                        else:
                            positions.remove(position2)
                            link_atom_2_index = oxygen_atom.GetNeighbors()[0].GetIdx()

                if chainParticipating == 'SMB':

                    C_dummy_found = False
                    side_N_found = False
                    # 找 C_dummy
                    for atom in rw_main_chain_linked_mol.GetAtoms():
                        if int(atom.GetProp('residue')) in positions:
                            if atom.HasProp('C_dummy'):
                                dummy_atoms_idxs_to_remove.append(atom.GetIdx())
                                if int(atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                                    positions.remove(position1)
                                    # 这里 atom 是 dummy atom，要找唯一的邻居才是 C_terminal
                                    link_atom_1_index = atom.GetNeighbors()[0].GetIdx()
                                else:
                                    positions.remove(position2)
                                    link_atom_2_index = atom.GetNeighbors()[0].GetIdx()
                                C_dummy_found = True
                                break

                    # 如果找到 C_dummy 了就去侧链上找 N
                    if C_dummy_found:
                        N_min_valence_idx = None
                        N_min_valence = 10
                        for atom in rw_main_chain_linked_mol.GetAtoms():
                            if int(atom.GetProp('residue')) in positions and not atom.HasProp('N_terminal') and atom.GetSymbol() == 'N' and atom.GetExplicitValence()<3:
                                if atom.GetExplicitValence()<N_min_valence:
                                    N_min_valence_idx = atom.GetIdx()
                                    N_min_valence = atom.GetExplicitValence()

                        # 如果找到了侧链上的 N
                        if N_min_valence_idx is not None:
                            side_N_found = True
                            atom = rw_main_chain_linked_mol.GetAtomWithIdx(N_min_valence_idx)
                            if int(atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                                positions.remove(position1)
                                link_atom_1_index = N_min_valence_idx
                            else:
                                positions.remove(position2)
                                link_atom_2_index = N_min_valence_idx

                        # 如果找到了 C_dummy 但是没找到 侧链上的 N
                        else:
                            link_atom_1_index = None
                            link_atom_2_index = None
                            positions = [position1, position2]
                            dummy_atoms_idxs_to_remove = []

                    # 如果其实没找到 C_dummy 或者没找到侧链上的 N，那就复原参数
                    if not C_dummy_found or not side_N_found:
                        link_atom_1_index = None
                        link_atom_2_index = None
                        positions = [position1, position2]
                        dummy_atoms_idxs_to_remove = []

                        N_dummy_found = False
                        side_OH_found = False
                        # 找 N_dummy
                        for atom in rw_main_chain_linked_mol.GetAtoms():
                            if int(atom.GetProp('residue')) in positions:
                                if atom.HasProp('N_dummy'):
                                    dummy_atoms_idxs_to_remove.append(atom.GetIdx())
                                    if int(atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                                        positions.remove(position1)
                                        # 这里 atom 是 dummy atom，要找唯一的邻居才是 C_terminal
                                        link_atom_1_index = atom.GetNeighbors()[0].GetIdx()
                                    else:
                                        positions.remove(position2)
                                        link_atom_2_index = atom.GetNeighbors()[0].GetIdx()
                                    N_dummy_found = True
                                    break

                        # 找到 N dummy 就去找侧链上的 -OH
                        if N_dummy_found:
                            # 找侧链的 -OH
                            # 找到 N_dummy 再找 -OH，要尽量排除 -COOH 上的 -OH 和 =O
                            hydroxyl_idx_list = []
                            side_COOH_OH_idx_list = []
                            oxygen_idx_to_remove = None
                            for atom in rw_main_chain_linked_mol.GetAtoms():
                                # 注意这里一定要是 atom.GetExplicitValence() == 1 而不是 degree=1，因为双键O也是只有一个邻居
                                if atom.GetSymbol() == 'O' and atom.GetExplicitValence() == 1 and int(atom.GetProp('residue')) in positions:
                                    hydroxyl_idx_list.append(atom.GetIdx())

                            # 如果不止一个 -OH ，那可能有 -COOH 或者不止一个 -OH, 尝试去掉其中的 -COOH 中的 -OH
                            if len(hydroxyl_idx_list) > 1:
                                # 利用deepcopy进行循环这样就不会改变原始的 list 导致循环出错
                                for hydroxyl_idx in copy.deepcopy(hydroxyl_idx_list):
                                    oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(hydroxyl_idx)
                                    oxygen_neighbor_atom = oxygen_atom.GetNeighbors()[0]
                                    # 检查这个邻居是否是 -COOH 中的 C 原子
                                    if oxygen_neighbor_atom.GetSymbol() == "C" and oxygen_neighbor_atom.GetDegree() == 3 and oxygen_neighbor_atom.GetExplicitValence() == 4 and not oxygen_neighbor_atom.HasProp(
                                            'C_terminal'):
                                        O_count = 0  # 氧原子计数
                                        double_bond_to_oxygen = False  # 是否找到双键 O
                                        for neighbor in oxygen_neighbor_atom.GetNeighbors():
                                            if neighbor.GetSymbol() == 'O':
                                                O_count += 1
                                                neighbor_bond = rw_main_chain_linked_mol.GetBondBetweenAtoms(
                                                    oxygen_neighbor_atom.GetIdx(), neighbor.GetIdx())
                                                if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                                    double_bond_to_oxygen = True

                                        # 确保有两个氧邻居，且一个氧通过双键连接，那就是找到了
                                        if O_count == 2 and double_bond_to_oxygen:
                                            # 去掉那些和 —COOH 相连的 -OH
                                            hydroxyl_idx_list.remove(hydroxyl_idx)
                                            side_COOH_OH_idx_list.append(hydroxyl_idx)

                            # 尝试去掉 -COOH 中的 -OH 之后如果还有不止一个 -OH，那就随便选择第一个作为侧链 -OH
                            if len(hydroxyl_idx_list) >= 1:
                                # 这个侧链上的 -OH 形成 ester 时不用被去掉
                                oxygen_idx_to_remove = hydroxyl_idx_list[0]
                                # 获得这个 O 原子对象
                                oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                                if int(oxygen_atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                                    # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                                    positions.remove(position1)
                                    # 这个 -OH 形成 Amine 的脱水过程中需要被去掉
                                    link_atom_1_index = oxygen_atom.GetNeighbors()[0].GetIdx()
                                else:
                                    positions.remove(position2)
                                    link_atom_2_index = oxygen_atom.GetNeighbors()[0].GetIdx()

                            # 如果没有单独的 -OH 被找到，那可能是在 -COOH 上
                            if len(hydroxyl_idx_list) == 0 and len(side_COOH_OH_idx_list) > 0:
                                oxygen_idx_to_remove = side_COOH_OH_idx_list[0]
                                # 获得这个 O 原子对象
                                oxygen_atom = rw_main_chain_linked_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                                if int(oxygen_atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                                    # 匹配到一个之后就把 position1 或者 position2 删掉防止重复匹配
                                    positions.remove(position1)
                                    # 这个 -OH 形成 Amine 的脱水过程中需要被去掉
                                    link_atom_1_index = oxygen_atom.GetNeighbors()[0].GetIdx()
                                else:
                                    positions.remove(position2)
                                    link_atom_2_index = oxygen_atom.GetNeighbors()[0].GetIdx()

                    # 如果到这里都还没找到，拿其实是一个侧链 N 和 主链 N 融合了 https://www.mdpi.com/2079-6382/11/8/1080
                    if None in [link_atom_1_index, link_atom_2_index]:
                        # 重新开始一边找 主链 N 一边找侧链N，先找 侧链 N
                        N_min_valence_idx = None
                        N_min_valence = 10
                        side_N_found = False
                        N_idx_to_remove = None
                        for atom in rw_main_chain_linked_mol.GetAtoms():
                            if int(atom.GetProp('residue')) in positions and not atom.HasProp(
                                    'N_terminal') and atom.GetSymbol() == 'N' and atom.GetExplicitValence() < 3:
                                if atom.GetExplicitValence() < N_min_valence:
                                    N_min_valence_idx = atom.GetIdx()
                                    N_min_valence = atom.GetExplicitValence()

                        # 如果找到了侧链上的 N
                        if N_min_valence_idx is not None:
                            side_N_found = True
                            atom = rw_main_chain_linked_mol.GetAtomWithIdx(N_min_valence_idx)
                            if int(atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                                positions.remove(position1)
                                link_atom_1_index = atom.GetNeighbors()[0].GetIdx()
                            else:
                                positions.remove(position2)
                                link_atom_2_index = atom.GetNeighbors()[0].GetIdx()
                            N_idx_to_remove = N_min_valence_idx

                        if side_N_found:
                            # 开始寻找主链的 N
                            for atom in rw_main_chain_linked_mol.GetAtoms():
                                if int(atom.GetProp('residue')) in positions and atom.HasProp('N_terminal'):
                                    if int(atom.GetProp('residue')) == position1 and link_atom_1_index is None:
                                        positions.remove(position1)
                                        link_atom_1_index = atom.GetIdx()
                                    else:
                                        positions.remove(position2)
                                        link_atom_2_index = atom.GetIdx()
                                    break

                # 连接新键并检查是否有错
                # 所有情况都处理完之后连接两个需要连接的原子
                if link_atom_1_index is not None and link_atom_2_index is not None:
                    # 找到的两个 =C 原子之间添加单键
                    rw_main_chain_linked_mol.AddBond(link_atom_1_index, link_atom_2_index,Chem.BondType.SINGLE)
                elif link_atom_1_index is None and link_atom_2_index is None:
                    # 这是都是数据有误的情况
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Amine Bond, chainParticipating: {chainParticipating} no match on |position_1: {position1}，original {correct_positions.get('position1')}| AND |position_2: {position2}，original {correct_positions.get('position2')}|, please check")
                elif link_atom_1_index is None and link_atom_2_index is not None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Amine Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_2_index).GetSymbol()} atom on |position_2: {position2}|, no match on |position_1: {position1}, original {correct_positions.get('position1')}|, please check")
                elif link_atom_1_index is not None and link_atom_2_index is None:
                    self.noise_data_flag = True
                    print(
                        f"Peptide {self.idx}: Amine Bond, chainParticipating: {chainParticipating} found {rw_main_chain_linked_mol.GetAtomWithIdx(link_atom_1_index).GetSymbol()} atom on |position_1: {position1}|, no match on |position_2: {position2}, original {correct_positions.get('position2')}|, please check")

                idx_to_remove = []

                if len(dummy_atoms_idxs_to_remove) > 0:
                    idx_to_remove.extend(dummy_atoms_idxs_to_remove)
                if oxygen_idx_to_remove is not None:
                    idx_to_remove.append(oxygen_idx_to_remove)
                if N_idx_to_remove is not None:
                    idx_to_remove.append(N_idx_to_remove)

                # 降序删除多余的 dummy atoms，否则会出错 Range Error
                if link_atom_1_index is not None and link_atom_2_index is not None:
                    if len(idx_to_remove) > 0:
                        for atom_idx in sorted(idx_to_remove, reverse=True):
                            rw_main_chain_linked_mol.RemoveAtom(atom_idx)

                # 这都是画图用的部分
                # img = Draw.MolToImage(rw_main_chain_linked_mol, size=(5000, 4000))
                # fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）
                # fig.patch.set_facecolor('black')  # 设置黑色背景
                # plt.imshow(img)  # 显示图片
                # plt.axis('off')  # 关闭坐标轴
                # plt.show()

                return rw_main_chain_linked_mol

            # 建立 bond 名到对应函数的映射字典
            bond_link_func_dict = {'DSB': DSB, 'AMD': AMD, 'TIE':TIE, 'DCB':DCB, 'EST':EST, 'p-XylB':p_XylB, 'TRZB':TRZB, '(E)-but-2-enyl-B':E_but_2_enyl_B, 'BisMeBn-B':BisMeBn_B, 'but-2-ynyl-B':but_2_ynyl_B, 'AMN':AMN}


            # 处理 bonds 中的每一个 bond，bonds 的格式和 DBAASP 里面原始的 intrachainBonds 是一样的
            for bond in bonds:

                # 测试一种 intrachainBond 连接的代码
                # if bond['type']['name'] not in ['AMN']: # or bond['chainParticipating']['name'] != 'MMB':
                if bond['type']['name'] not in list(bond_link_func_dict.keys()):
                    # 如果有没实现的键连接，抛弃这个序列
                    self.noise_data_flag = True
                    continue
                rw_main_chain_linked_mol = bond_link_func_dict[bond['type']['name']](rw_main_chain_linked_mol, bond['position1'], bond['position2'], bond.get('chainParticipating', {}).get('name') if bond['chainParticipating'] is not None else 'SSB', seq_len)

            return rw_main_chain_linked_mol

        # 把这些连接好主链的都先转换成可编辑的分子
        rw_main_chain_linked_mols = []
        for main_chain_linked_mol in self.main_chain_linked_mols:
            rw_main_chain_linked_mols.append(Chem.RWMol(main_chain_linked_mol))

        # 如果只有一条序列，那么 self.intrachain_bonds 的格式是：[{}, {}, {},...]
        if isinstance(self.aa_seqs, str):
            self.intrachainBonds_linked_mols.append(link_one_chain_intrachain_bonds(rw_main_chain_linked_mols[0], self.intrachain_bonds, len(self.aa_seqs)))
        if isinstance(self.aa_seqs, list):
            # TODO：multimer
            pass

    def link_interchain_bonds(self):
        # TODO：multimer
        pass

    def c_terminus_modification(self, rw_intrachainBonds_linked_mol, c_terminus_modify_name_smiles_dict):

        # 有效的modification才继续
        if type(c_terminus_modify_name_smiles_dict[self.cTerminus]) is float:
            self.noise_data_flag = True
            return None
        # 加载分子
        modification_mol = Chem.MolFromSmiles(c_terminus_modify_name_smiles_dict[self.cTerminus])
        rw_modification_mol = Chem.RWMol(modification_mol)

        # 初始化一些变量
        atom_idx_to_remove = []

        # 先寻找 modification 上的 dummy atom ，指定好的方便连接
        modi_dummy_atom_idx = None
        modi_atom_to_be_linked_idx = None
        for atom in rw_modification_mol.GetAtoms():
            if atom.GetSymbol() == '*':
                modi_dummy_atom_idx = atom.GetIdx()
                atom_idx_to_remove.append(modi_dummy_atom_idx)
                modi_atom_to_be_linked_idx = atom.GetNeighbors()[0].GetIdx()
                break

        # 没找到修饰分子的 dummy atom 就去找 N 原子
        if modi_atom_to_be_linked_idx is None:
            # 找 modification_mol valence 最小的 N
            N_min_valence_idx = None
            N_min_valence = 10
            for atom in rw_modification_mol.GetAtoms():
                if atom.GetSymbol() == 'N' and atom.GetExplicitValence() < 3:
                    if atom.GetExplicitValence() < N_min_valence:
                        N_min_valence_idx = atom.GetIdx()
                        N_min_valence = atom.GetExplicitValence()

                        modi_atom_to_be_linked_idx = N_min_valence_idx

        # 如果 dummy atom 和 N 原子都没找到，那就是和 -OH 发生反应
        if modi_atom_to_be_linked_idx is None:
            # 再找侧链上的 -OH
            # 找到 N_dummy 再找 -OH，要尽量排除 -COOH 上的 -OH 和 =O
            hydroxyl_idx_list = []
            side_COOH_OH_idx_list = []
            oxygen_idx_to_remove = None
            for atom in rw_modification_mol.GetAtoms():
                # 注意这里一定要是 atom.GetExplicitValence() == 1 而不是 degree=1，因为双键O也是只有一个邻居
                if atom.GetSymbol() == 'O' and atom.GetExplicitValence() == 1:
                    hydroxyl_idx_list.append(atom.GetIdx())

            # 如果不止一个 -OH ，那可能有 -COOH 或者不止一个 -OH, 尝试去掉其中的 -COOH 中的 -OH
            if len(hydroxyl_idx_list) > 1:
                # 利用deepcopy进行循环这样就不会改变原始的 list 导致循环出错
                for hydroxyl_idx in copy.deepcopy(hydroxyl_idx_list):
                    oxygen_atom = rw_modification_mol.GetAtomWithIdx(hydroxyl_idx)
                    oxygen_neighbor_atom = oxygen_atom.GetNeighbors()[0]
                    # 检查这个邻居是否是 -COOH 中的 C 原子
                    if oxygen_neighbor_atom.GetSymbol() == "C" and oxygen_neighbor_atom.GetDegree() == 3 and oxygen_neighbor_atom.GetExplicitValence() == 4:
                        O_count = 0  # 氧原子计数
                        double_bond_to_oxygen = False  # 是否找到双键 O
                        for neighbor in oxygen_neighbor_atom.GetNeighbors():
                            if neighbor.GetSymbol() == 'O':
                                O_count += 1
                                neighbor_bond = rw_modification_mol.GetBondBetweenAtoms(oxygen_neighbor_atom.GetIdx(), neighbor.GetIdx())
                                if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    double_bond_to_oxygen = True

                        # 确保有两个氧邻居，且一个氧通过双键连接，那就是找到了
                        if O_count == 2 and double_bond_to_oxygen:
                            # 去掉那些和 —COOH 相连的 -OH
                            hydroxyl_idx_list.remove(hydroxyl_idx)
                            side_COOH_OH_idx_list.append(hydroxyl_idx)

            # 尝试去掉 -COOH 中的 -OH 之后如果还有不止一个 -OH，那就随便选择第一个作为侧链 -OH
            if len(hydroxyl_idx_list) >= 1:
                # 这个侧链上的 -OH 形成 ester 时不用被去掉
                oxygen_idx_to_remove = hydroxyl_idx_list[0]
                # atom_idx_to_remove.append(oxygen_idx_to_remove)
                # 获得这个 O 原子对象
                # oxygen_atom = rw_modification_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                modi_atom_to_be_linked_idx = oxygen_idx_to_remove

            # 如果没有单独的 -OH 被找到，那可能是在 -COOH 上
            if len(hydroxyl_idx_list) == 0 and len(side_COOH_OH_idx_list) > 0:
                oxygen_idx_to_remove = side_COOH_OH_idx_list[0]
                # atom_idx_to_remove.append(oxygen_idx_to_remove)
                # 获得这个 O 原子对象
                # oxygen_atom = rw_modification_mol.GetAtomWithIdx(oxygen_idx_to_remove)
                modi_atom_to_be_linked_idx = oxygen_idx_to_remove

        # 如果所有可以和 -COOH 反应的基团都没有找到
        if modi_atom_to_be_linked_idx is None:
            self.noise_data_flag = True
            print(f'C terminal modification error: no * / -NH2 / -OH found in modification mol: {self.cTerminus}, AMP id: {self.idx}')

        # 如果找到了一个就去找主链的 C_dummy
        else:
            # 找主链的C_dummy
            main_chain_C_dummy_idx = None
            main_chain_C_idx = None
            for atom in rw_intrachainBonds_linked_mol.GetAtoms():
                if atom.GetSymbol() == '*' and atom.HasProp('C_dummy'):
                    main_chain_C_dummy_idx = atom.GetIdx()
                    main_chain_C_idx = atom.GetNeighbors()[0].GetIdx()
                    break
            combined_mol = rdmolops.CombineMols(rw_modification_mol, rw_intrachainBonds_linked_mol)
            rw_c_terminus_modified_mol = Chem.RWMol(combined_mol)

            offset = rw_modification_mol.GetNumAtoms()
            if main_chain_C_idx is not None:
                main_chain_C_idx += offset
                main_chain_C_dummy_idx += offset
                atom_idx_to_remove.append(main_chain_C_dummy_idx)
                # 添加键连接
                rw_c_terminus_modified_mol.AddBond(main_chain_C_idx, modi_atom_to_be_linked_idx, Chem.rdchem.BondType.SINGLE)
                # 删除多余的 atoms
                for atom_idx in sorted(atom_idx_to_remove, reverse=True):
                    rw_c_terminus_modified_mol.RemoveAtom(atom_idx)
                self.cTerminus_modified_mols.append(rw_c_terminus_modified_mol)
            else:
                print(f'C terminal modification error: no C_dummy atom found in the main chain of AMP {self.idx} {self.cTerminus}')
                self.noise_data_flag = True

    def n_terminus_modification(self, cTerminus_modified_mol, n_terminus_modify_name_smiles_dict):
        # 加载分子
        modification_mol = Chem.MolFromSmiles(n_terminus_modify_name_smiles_dict[self.nTerminus])
        rw_modification_mol = Chem.RWMol(modification_mol)

        # 初始化一些变量
        atom_idx_to_remove = []

        # 先寻找 modification 上的 dummy atom ，指定好的方便连接
        modi_dummy_atom_idx = None
        # dummy atom 所有的邻居都要和主链连接然后将 dummy atom 删去
        modi_atom_to_be_linked_idxs = []
        for atom in rw_modification_mol.GetAtoms():
            if atom.GetSymbol() == '*':
                modi_dummy_atom_idx = atom.GetIdx()
                atom_idx_to_remove.append(modi_dummy_atom_idx)
                for neighbor in atom.GetNeighbors():
                    modi_atom_to_be_linked_idxs.append(neighbor.GetIdx())
                break

        # 没找到修饰分子的 dummy atom 就去找 -COOH 上的 -OH
        if len(modi_atom_to_be_linked_idxs) == 0:
            for atom in rw_modification_mol.GetAtoms():
                # 寻找 -COOH
                if atom.GetSymbol() == "C" and atom.GetDegree() == 3 and atom.GetExplicitValence() == 4:
                    O_count = 0  # 氧原子计数
                    double_bond_to_oxygen = False  # 是否找到双键 O
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'O':
                            O_count += 1
                            neighbor_bond = rw_modification_mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                            if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                double_bond_to_oxygen = True
                            # 记录应该被去掉的那一个 -OH
                            if neighbor_bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                oxygen_idx_to_remove = neighbor.GetIdx()

                    # 确保有两个氧邻居，且一个氧通过双键连接
                    if O_count == 2 and double_bond_to_oxygen:
                        atom_idx_to_remove.append(oxygen_idx_to_remove)
                        modi_atom_to_be_linked_idxs.append(atom.GetIdx())
                        break

        # 如果 dummy atom 和 -COOH 上的 -OH 都没找到，那就找 -OH
        if len(modi_atom_to_be_linked_idxs) == 0:
            for atom in rw_modification_mol.GetAtoms():
                if atom.GetSymbol() == "O" and atom.GetExplicitValence() == 1:
                    atom_idx_to_remove.append(atom.GetIdx())
                    modi_atom_to_be_linked_idxs.append(atom.GetNeighbors()[0].GetIdx())
                    break

        # 如果所有可以和 -COOH 反应的基团都没有找到
        if len(modi_atom_to_be_linked_idxs) == 0:
            self.noise_data_flag = True
            print(f'N terminal modification error: no * / -COOH / -OH found in modification mol: {self.nTerminus}, AMP id: {self.idx}')

        else:
            # 找主链的 N_dummy
            main_chain_N_dummy_idx = None
            main_chain_N_idx = None
            for atom in cTerminus_modified_mol.GetAtoms():
                if atom.GetSymbol() == '*' and atom.HasProp('N_dummy'):
                    main_chain_N_dummy_idx = atom.GetIdx()
                    main_chain_N_idx = atom.GetNeighbors()[0].GetIdx()
                    break
            combined_mol = rdmolops.CombineMols(rw_modification_mol, cTerminus_modified_mol)
            rw_n_terminus_modified_mol = Chem.RWMol(combined_mol)

            offset = rw_modification_mol.GetNumAtoms()
            if main_chain_N_idx is not None:
                main_chain_N_idx += offset
                main_chain_N_dummy_idx += offset
                atom_idx_to_remove.append(main_chain_N_dummy_idx)
                # 添加键连接
                for modi_atom_to_be_linked_idx in modi_atom_to_be_linked_idxs:
                    rw_n_terminus_modified_mol.AddBond(main_chain_N_idx, modi_atom_to_be_linked_idx, Chem.rdchem.BondType.SINGLE)

                # 如果 N 连接的原子数超过允许的范围了就加电荷平衡 valence
                main_N_atom = rw_n_terminus_modified_mol.GetAtomWithIdx(main_chain_N_idx)


                # 删除多余的 atoms
                for atom_idx in sorted(atom_idx_to_remove, reverse=True):
                    rw_n_terminus_modified_mol.RemoveAtom(atom_idx)

                rw_n_terminus_modified_mol.UpdatePropertyCache(strict=False)
                if main_N_atom.GetDegree() > 3:
                    main_N_atom.SetFormalCharge(main_N_atom.GetDegree() - 3)
                self.ncTerminus_modified_mols.append(rw_n_terminus_modified_mol)
            else:
                print(f'N terminal modification error: no N_dummy atom found in the main chain of AMP {self.idx} {self.nTerminus}')
                self.noise_data_flag = True



    def remove_redundant_dummy_atoms(self):
        """
        去除连接完之后多余不必要的 dummy atoms
        应该在所有键，包括 interchain Bonds 都连接完之后再做，但是这里先测试一下连完intrachain bonds的部分
        :return: None
        """
        C_dummy_atom_idx = None
        N_dummy_atom_idx = None
        for atom in self.ncTerminus_modified_mols[0].GetAtoms():
            if atom.HasProp('C_dummy'):
                COOH_C_atom = atom.GetNeighbors()[0]
                C_dummy_atom_idx = atom.GetIdx()
            if atom.HasProp('N_dummy'):
                N_dummy_atom_idx = atom.GetIdx()
        if N_dummy_atom_idx is not None or C_dummy_atom_idx is not None:
            idx_to_remove = [idx for idx in [N_dummy_atom_idx, C_dummy_atom_idx] if idx is not None]
            for atom_idx in sorted(idx_to_remove, reverse=True):
                self.ncTerminus_modified_mols[0].RemoveAtom(atom_idx)
        if C_dummy_atom_idx is not None:
            oxygen = Chem.MolFromSmiles("[OH]")
            rw_oxygen = Chem.RWMol(oxygen)
            oxygen_idx = self.ncTerminus_modified_mols[0].AddAtom(rw_oxygen.GetAtomWithIdx(0))
            self.ncTerminus_modified_mols[0].AddBond(COOH_C_atom.GetIdx(), oxygen_idx, Chem.BondType.SINGLE)

        # 这都是画图用的部分
        # img = Draw.MolToImage(self.intrachainBonds_linked_mols[0], size=(5000, 4000))
        # fig = plt.figure(figsize=(20, 20))  # figsize 参数控制图的尺寸（单位：英寸）
        # fig.patch.set_facecolor('black')  # 设置黑色背景
        # plt.imshow(img)  # 显示图片
        # plt.axis('off')  # 关闭坐标轴
        # plt.show()


def get_aa_smiles_dict(smiles_path:str) -> Dict:
    """
    读取氨基酸对应的smiles并且加上g的
    :param smiles_path: 到存储 smiles 的路径, 文件名应该是 all_aa_smiles.csv
    :return: aa_smiles_dict
    """
    with open(smiles_path, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        aa_smiles_dict = {row['aa']: row['SMILES'] for row in reader}

    # G和g的SMILES应该一样
    aa_smiles_dict['g'] = aa_smiles_dict.get('G')
    return aa_smiles_dict

def load_DBAASP_data(wo_PubChem_path:str):
    """
    返回一个list，里面都是处理好的 [[aa_seqs], [intrachain bonds dict], [interchain bonds dicy]]
    :param wo_PubChem_path: path to peptides_wo_PubChem_SMILES_data.json
    :return: [[aa_seqs], [intrachain bonds dict], [interchain bonds dicy]]
    """
    with open(wo_PubChem_path, 'r', encoding='utf-8') as file:
        data = json.load(file)  # 将JSON内容加载为Python字典

    aa_seqs = []
    intrachain_bonds = []
    interchain_bonds = []
    for AMP in data:
        if AMP['complexity']['name'] == 'Monomer':
            aa_seqs.append(AMP['sequence'])
            intrachain_bonds.append(AMP['intrachainBonds'])
            interchain_bonds.append(AMP['interchainBonds'])


if __name__ == '__main__':
    from PepLink.data import load_default_aa_smiles

    aa_smiles_dict = load_default_aa_smiles()
    print(f"Loaded {len(aa_smiles_dict)} amino-acid mappings from PepLink/data.")
