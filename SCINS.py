import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold



class SCINS_generator:
    '''

    Class to calculate SCINS of a given SMILES string

    '''

    def __init__(self,smi):
        self._smi = smi
        self._mol = Chem.MolFromSmiles(smi)
        self._scaf = MurckoScaffold.GetScaffoldForMol(MurckoScaffold.MakeScaffoldGeneric(self._mol))
        self._scaf_atoms = self._scaf.GetAtoms()
        self._scaf_bonds = self._scaf.GetBonds()
        self._scaf_smi = Chem.MolToSmiles(self._scaf)
        self._ring_system = self.GetRingSystemsofscaf()
        self._ring_system_count = self.count_ring_systems()
        self._bin_values = [1,2,3,4,7]

        # Linkers: [direct bond between rings, linear chain between rings, branched chain between rings]
        self._linkers = [0,0,0]
        self._chain_binning = [0,0,0,0,0]

    def count_ring_systems(self):
        '''
        Count Ring systems and separate them into exactly one ring, two rings or three or more
        :return:
        '''
        # rings: 1, 2, >3
        at_in_ring = [a.IsInRing() for a in self._scaf_atoms]
        ring_comp = [0,0,0]
        for system in self._ring_system:
            more_than_one = 0
            for i in system:
                n_atoms = self._scaf.GetAtomWithIdx(i).GetNeighbors()
                # if molecule has more than 2 neighbors in ringsystem count additional ring
                num_ri = 0
                for at in n_atoms:
                    idx = at.GetIdx()
                    if at_in_ring[idx] and idx in system:
                        num_ri = num_ri + 1

                if num_ri > 2:
                    more_than_one = more_than_one + 1
            if more_than_one == 0:
                ring_comp[0] = ring_comp[0] + 1
            elif more_than_one == 2:
                ring_comp[1] = ring_comp[1] + 1
            else:
                ring_comp[2] = ring_comp[2] + 1

        return ring_comp

    def GetRingSystems(self,mol,includeSpiro=False):
        '''
        function to get Ring Systems
        :param mol:     RDKit Molecule, molecule for which information is seeked
        :param includeSpiro: Should be False -- from RDKit standard documentation
        :return: Ring systems
        '''
        ri = mol.GetRingInfo()
        systems = []
        for ring in ri.AtomRings():
            ringAts = set(ring)
            nSystems = []
            for system in systems:
                nInCommon = len(ringAts.intersection(system))
                if nInCommon and (includeSpiro or nInCommon > 1):
                    ringAts = ringAts.union(system)
                else:
                    nSystems.append(system)
            nSystems.append(ringAts)
            systems = nSystems
        return systems

    def GetRingSystemsofscaf(self):
        return self.GetRingSystems(self._scaf)

    def explore_bonds(self,bond,bond_indices,branched,at_in_ring):
        '''
        Explore bonds in recursive fashion and return indices of linker and if branched or not
        :param bond:           RDKit Bond, Bond to be analyze
        :param bond_indices:   vector Bool, if bond idx is set to True bond has been considered for linker calculation
        :param branched:       Bool, if True Linker is branched
        :param at_in_ring:     vector of Bool, if True atom is in ring
        :return:
        '''

        # Get begining and end atom index
        b_at = bond.GetBeginAtomIdx()
        e_at = bond.GetEndAtomIdx()

        bond_indices.append(bond.GetIdx())

        if at_in_ring[b_at] and at_in_ring[e_at]:
            return bond_indices,branched

        # if not in ring explore
        if not at_in_ring[b_at]:
            # set to in ring for marking considered
            at_in_ring[b_at] = True
            atom = self._scaf.GetAtomWithIdx(b_at)
            bonds = atom.GetBonds()
            bond_idx = [b.GetIdx() for b in bonds]

            # Check if branched
            if len(bond_idx) > 2:
                branched = True

            for b_idx in bond_idx:
                # Check if allready explored
                if b_idx not in bond_indices:
                    bond = self._scaf.GetBondWithIdx(b_idx)
                    # Explore new bond
                    new_indices, new_branched = self.explore_bonds(bond,bond_indices,branched,at_in_ring)
                    bond_indices = bond_indices + new_indices
                    bond_indices = list(set(bond_indices))

                    if new_branched == True:
                        branched = True

        if not at_in_ring[e_at]:
            # set to in ring for marking considered
            at_in_ring[e_at] = True
            atom = self._scaf.GetAtomWithIdx(e_at)
            bonds = atom.GetBonds()
            bond_idx = [b.GetIdx() for b in bonds]

            # Check if branched
            if len(bond_idx) > 2:
                branched = True

            for b_idx in bond_idx:
                # Check if allready explored
                if b_idx not in bond_indices:
                    bond = self._scaf.GetBondWithIdx(b_idx)
                    # Explore new bond
                    new_indices, new_branched = self.explore_bonds(bond,bond_indices,branched,at_in_ring)
                    bond_indices = bond_indices + new_indices
                    bond_indices = list(set(bond_indices))

                    if new_branched == True:
                        branched = True

        return bond_indices,branched

    def calc_fragments(self):
        '''
        Function to calculate linkers, and bin them by length
        :return:
        '''
        # set linkers to zero
        self._linkers[0] = 0
        self._linkers[1] = 0
        self._linkers[2] = 0
        self._chain_binning = [0,0,0,0,0]

        # Calculate Linkers
        bond_in_ring = [b.IsInRing() for b in self._scaf_bonds]
        at_in_ring = [a.IsInRing() for a in self._scaf_atoms]

        # Go through bonds
        b = 0
        while 0 in bond_in_ring:
            bond = self._scaf.GetBondWithIdx(b)

            # only consider bonds not in rings
            if not bond_in_ring[b]:
                b_at = bond.GetBeginAtomIdx()
                e_at = bond.GetEndAtomIdx()

                # if both connected atoms are in rings linker between atoms
                if at_in_ring[b_at] and at_in_ring[e_at]:
                    self._linkers[0] = self._linkers[0] + 1
                    bond_in_ring[b] = True
                    self._chain_binning[0] = self._chain_binning[0] + 1

                else:
                    # Explore linker
                    indices, branched = self.explore_bonds(bond=bond,bond_indices=[b],branched=False,at_in_ring=at_in_ring)

                    for idx in indices:
                        bond_in_ring[idx] = True
                    if branched:
                        self._linkers[2] = self._linkers[2] + 1
                    else:
                        self._linkers[1] = self._linkers[1] + 1

                        # Bin lengths
                        chain_len = len(indices)

                        if chain_len == 2:
                            self._chain_binning[1] = self._chain_binning[1] + 1
                        if chain_len == 3 or chain_len == 4:
                            self._chain_binning[2] = self._chain_binning[2] + 1
                        if chain_len == 5 or chain_len == 6:
                            self._chain_binning[3] = self._chain_binning[3] + 1
                        if chain_len >= 7:
                            self._chain_binning[4] = self._chain_binning[4] + 1

            b = b + 1


    def calc_A(self):
        '''
        Number of Chain Assemblies. Chain assemblies are contiguous linkers between ring assemblies. They are uncovered by removing all ring bonds in the molecule
        :return:    int, number of chain assemblies in molecule
        '''

        return self._linkers[0] + self._linkers[1] + self._linkers[2]

    def calc_B(self):
        '''
        Chains are all unbranched linkers needed to cover all nonring bonds in the molecule.
        :return:
        '''

        return self._linkers[0] + self._linkers[1]

    def calc_C(self):
        '''
        Number of rings
        :return:    int, number of rings in molecule
        '''

        return Chem.rdMolDescriptors.CalcNumRings(self._scaf)

    def calc_D(self):
        '''
        Number of Ring Assemblies. Ring assemblies are fragments remaining when all acyclic bonds have been removed.
        :return:    int, number of Ring assemblies
        '''

        return len(self._ring_system)

    def calc_E(self):
        '''
        Number of Bridge Bonds. A contiguous path of more than one bond shared between more than one rings counts as bridge bond.
        Here interpreted by amount of Bridges counted by RDKit
        :return:
        '''

        return int(Chem.rdMolDescriptors.CalcNumBridgeheadAtoms(self._scaf)/2)

    def calc_F(self):
        '''
        Number of Ring Assemblies Consisting of Exactly One Ring.
        :return:
        '''

        return self._ring_system_count[0]

    def calc_G(self):
        '''
        Number of Ring Assemblies Consisting of Exactly Two Rings.
        :return:
        '''

        return self._ring_system_count[1]

    def calc_H(self):
        '''
        Number of Ring Assemblies Consisting of more than 3 Ring.
        :return:
        '''

        return self._ring_system_count[2]

    def calc_I(self):
        '''
        Number of Macrocylces defined as at least 12 atoms
        :return:
        '''
        macro = Chem.MolFromSmarts("[r{12-}]")
        macro_t = self._scaf.GetSubstructMatches(macro)
        macro_id = [macro_t[i][0] for i in range(len(macro_t))]
        cycles = 0

        if str(macro_t) == "()":
            return 0
        else:
            cylce_id = 0
            for system in self._ring_system:
                cylce_id = cylce_id + 1
                system_list = list(system)
                # Check if part of cycle
                for i, m in enumerate(macro_id):
                    if m in system_list:
                        macro_id[i] = -cylce_id

            # Return unique number of macrocycles
            return len(np.unique(macro_id))

    def calc_lowest_bin(self):
        '''
        Calculate the bin value of the shortest chain still binned
        :return:
        '''
        # Delete entry in bin
        for i, c_b in enumerate(self._chain_binning):
            if c_b != 0:
                self._chain_binning[i] = self._chain_binning[i] -1
                return self._bin_values[i]
        return 0

    def calc_J(self):
        '''
        Chains are all unbranched linkers needed to cover all nonring bonds in the molecule
        Binned Length of Shortest Shortest Chain. If the binned length of the shortest chain exists, it is used; otherwise, it is zero
        :return:
        '''

        return self.calc_lowest_bin()

    def calc_K(self):
        '''
        WARNING ALLWAYS CALL calc_J before!
        Chains are all unbranched linkers needed to cover all nonring bonds in the molecule
        Binned Length of Second Shortest Chain. If the binned length of the second shortest chain exists, it is used; otherwise, it is zero.
        :return:
        '''

        return self.calc_lowest_bin()

    def calc_L(self):
        '''
        WARNING ALLWAYS CALL calc_J before!
        Chains are all unbranched linkers needed to cover all nonring bonds in the molecule
        Binned Length of Third Shortest Chain. If the binned length of the second shortest chain exists, it is used; otherwise, it is zero.
        :return:
        '''

        return self.calc_lowest_bin()

    def calc_M(self):
        '''
        WARNING ALLWAYS CALL calc_J before!
        Chains are all unbranched linkers needed to cover all nonring bonds in the molecule
        Binned Length of Fourth Shortest Chain. If the binned length of the second shortest chain exists, it is used; otherwise, it is zero.
        :return:
        '''

        return self.calc_lowest_bin()

    def Calculate_SCINS(self,type="vec"):
        '''
        Calculate SCINS value
        :param type:    str, Output format
        :return:        vector of int, SCINS value
        '''
        self.calc_fragments()
        SCINS = [self.calc_A(),self.calc_B(),self.calc_C(),self.calc_D(),self.calc_E(),self.calc_F(),self.calc_G(),self.calc_H(),self.calc_I(),self.calc_J(),self.calc_K(),self.calc_L(),self.calc_M()]
        if type == "vec":
            return SCINS
        elif type == "str":
            strings = [str(integer) for integer in SCINS]
            return "".join(strings)
        elif type == "code":
            strings = [str(integer) for integer in SCINS]
            string = "".join(strings)
            return string[:5] + "-" + string[5:9] + "-" +string[9:]
        else:
            exit("TYPE is not supported! please use vec, str or code")