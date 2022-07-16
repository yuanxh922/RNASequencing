"""module to create sequences
"""
import sys
import random
import numpy as np
from math import log
from loguru import logger
from .bases import Bases
from .bases import Modifications
from .compounds import Compound
from .adducts import Adducts

class CompoundFeature:
    """represents the feature of Compound dataset
    """
    def __init__(self):
        self.min_mass = sys.maxsize
        self.max_mass = 0
        self.min_rt = sys.maxsize
        self.max_rt = 0
        self.min_vol = sys.maxsize
        self.max_vol = 0

class SeqManager:
    """sequences managemant
    """

    def __init__(self):
        """initial the class
        """
        self._bases = Bases()
        self._modifications = Modifications()
        self._sequences = list()

    @property
    def sequences(self):
        return self._sequences

    def make_cpds(self, orientation, length, width, N, noises, seqs):
        """create a dataset, which contains a sequnce with <length> nodes
        """
        tail = False
        anchor = 0
        rev_anchor = 0
        if orientation == 3:
            anchor = 694.2397
            #anchor = 826.3184
            rev_anchor = 18.0106
        elif orientation == 30:
            rev_anchor = 18.0106
            anchor = -61.9557
        elif orientation == 5:
            tail = True
            anchor = 593.1743
            rev_anchor = -61.9557
        elif orientation == 50:
            anchor = 18.0106
            rev_anchor = -61.9557
        elif orientation == 9:
            anchor = 18.0106 + 79.9663
            rev_anchor = -61.9557

        symbol_seqs = list()
        core_cpds = None
        if seqs:
            symbol_seqs, core_cpds = self.generate_multi_cpds(seqs, anchor, rev_anchor, tail)
        else:
            symbol_seqs, core_cpds = self.make_multi_core_seq(N, anchor, rev_anchor, length, width, tail)
            # print('len of symbol_seqs {}'.format(len(symbol_seqs)))
        #core_cpds.append(self.make_anchor_point(anchor))
        adducts = Adducts()
        salt_cpds = self.make_salts(adducts.combs(), core_cpds)
        effective_noises = self.make_effective_noises(salt_cpds)
        noises_inside = self.make_ineffective_noises_inside_zone(effective_noises)
        noises_lower_zone = self.make_ineffective_noises_in_lower_zone(effective_noises)
        noises_upper_zone = self.make_ineffective_noises_in_upper_zone(effective_noises)

        final_points = list()
        final_points.extend(core_cpds)
        if noises:
            final_points.extend(salt_cpds)
            final_points.extend(effective_noises)
            final_points.extend(noises_inside)
            final_points.extend(noises_lower_zone)
            final_points.extend(noises_upper_zone)

            noises_left_zone = self.make_ineffective_noises_in_left_zone(final_points)
            final_points.extend(noises_left_zone)

        final_points = self.make_extra_columns(final_points)
        return symbol_seqs, final_points

    def make_extra_columns(self, points):
        points = self._make_col_width(points)
        points = self._make_col_qs(points)
        return points

    def _is_symbol_occupied(self, seqs, seq):
        """we'd keep
        """
        if not seqs:
            return False
        exist_symbols = [seq[:3] for seq in seqs]
        if seq[:3] in exist_symbols:
            return True
        return False

    def generate_multi_cpds(self, seqs, anchor, rev_anchor, tail=False):
        """create compounds for seqs
        """
        symbol_size = len(self._bases.symbol_names)
        multi_cpds = list()
        symbol_seqs = list()
        for symbol_seq in seqs:
            _, cpds = self.make_core_seq(anchor, symbol_seq, tail)
            #multi_cpds.extend(cpds)
            #symbol_seqs.append(symbol_seq)

            symbol_seq_list = list(symbol_seq)
            rev_symbol_seq, rev_core_cpds = self.make_reversed_core_seq(rev_anchor, symbol_seq_list, not tail)
            multi_cpds.extend(rev_core_cpds)
            symbol_seqs.append(rev_symbol_seq)

        multi_cpds.sort(key = lambda x: x.mass)
        return symbol_seqs, multi_cpds

    def generate_cpds_5p(self, seq, anchor):
        seq = list(seq)
        return self.make_core_seq(anchor, seq, False)

    def generate_cpds_3p(self, seq, anchor):
        symbol_seq_list = list(seq)
        return self.make_reversed_core_seq(anchor, symbol_seq_list, True)

    def make_multi_core_seq(self, N, anchor, rev_anchor, length, width, tail=False):
        """create N seqs, each one starts at anchor, has <length+bias> items
        """
        symbol_size = len(self._bases.symbol_names)
        multi_cpds = list()
        symbol_seqs = list()
        while True:
            bias = 0
            if width > 0:
                bias = random.randint(1, width)
            symbol_seq = self.make_symbol_seq(length+bias)
            symbol_seq, cpds = self.make_core_seq(anchor, symbol_seq, tail)
            if self._is_symbol_occupied(symbol_seqs, symbol_seq):
                continue
            multi_cpds.extend(cpds)
            symbol_seqs.append(symbol_seq)

            rev_symbol_seq, rev_core_cpds = self.make_reversed_core_seq(rev_anchor, symbol_seq, not tail)
            multi_cpds.extend(rev_core_cpds)
            symbol_seqs.append(rev_symbol_seq)

            if len(symbol_seqs) == 2*N:
                break
        multi_cpds.sort(key = lambda x: x.mass)
        return symbol_seqs, multi_cpds

    def make_anchor_point(self, anchor):
        _rt = log(anchor, 1.5)
        _vol = random.randint(15e5, 30e5)
        cpd = Compound()
        cpd.mass = anchor
        cpd.rt = _rt
        cpd.vol = _vol
        return cpd

    def make_reversed_core_seq(self, anchor, seq, tail=False):
        seq.reverse()
        return self.make_seq_points(anchor, seq, tail)

    def make_symbol_seq(self, length):
        symbol_names = self._bases.symbol_names
        modification_names = self._modifications.symbol_names
        symbol_seq = self._make_symbol_seq_with_mod(length, symbol_names, modification_names)
        return symbol_seq

    #def make_core_seq(self, anchor, length, tail=False):
    def make_core_seq(self, anchor, symbol_seq, tail=False):
        """create a sequence with <length> nodes
        """
        #symbol_names = self._bases.symbol_names
        #modification_names = self._modifications.symbol_names
        #symbol_seq = self._make_symbol_seq_with_mod(length, symbol_names, modification_names)
        return self.make_seq_points(anchor, symbol_seq, tail)

    def make_seq_points(self, anchor, symbol_seq, tail=False):
        mass_seq = self._make_mass_seq(anchor, symbol_seq, tail)
        rt_seq = self._make_rt_seq(mass_seq)
        length = len(symbol_seq)
        vol_seq = self._make_vol_seq(length)

        data = dict()
        data["symbol_seq"] = "".join(symbol_seq)
        mass_seq_with_anchor = list()
        mass_seq_with_anchor.append(anchor)
        mass_seq_with_anchor.extend(mass_seq)
        data["mass_seq"] = mass_seq_with_anchor
        self._sequences.append(data)

        # print("{} {}".format("".join(symbol_seq), mass_seq))
        cpds = list()
        print("Nucle\tMass\t\tRT\tVol")
        mass_of_name = self._bases.mass_of_name
        for i in range(length):
            print("{}\t{:.4f}\t{:.2f}\t{:.4f}".format(symbol_seq[i], mass_seq[i], rt_seq[i], vol_seq[i]))
            cpd = Compound()
            cpd.mass = mass_seq[i]
            cpd.rt = rt_seq[i]
            cpd.vol = vol_seq[i]
            cpds.append(cpd)
        #seq = zip(mass_seq, rt_seq, vol_seq)
        #return list(seq)
        return symbol_seq, cpds

    def _make_symbol_seq(self, length, nodes):
        """create <length> length sequence with rand nodes from <nodes>
        :param nodes: a list of names which represents the final sequences
        """
        return [random.choice(nodes) for _ in range(length)]

    def _make_symbol_seq_with_mod(self, length, nodes, modifications):
        seq = self._make_symbol_seq(length, nodes)
        num_modif = int(length / 5.0)
        chose_idxs = random.choices(range(length)[1:-1], k=num_modif)
        for idx in chose_idxs:
            seq[idx] = random.choice(modifications)
        return seq


    def _make_mass_seq(self, anchor, symbol_seq, tail=False):
        """make mass sequences, start fro anchor, add mass for each symbol
        """
        def mass_of_name(symbol):
            mass = self._bases.mass_of_name(symbol)
            if not mass:
                mass = self._modifications.mass_of_name(symbol)
            return mass

        def mass_of_tail(symbol):
            mass = self._bases.mass_of_tail(symbol)
            if not mass:
                mass = self._modifications.mass_of_tail(symbol)
            return mass
        #mass_of_name = self._bases.mass_of_name
        #mass_of_tail = self._bases.mass_of_tail
        base_seq = [mass_of_name(symbol) for symbol in symbol_seq[:-1]]
        # print("base_seq {}".format(base_seq))
        if tail:
            base_seq.append(mass_of_tail(symbol_seq[-1]))
        else:
            base_seq.append(mass_of_name(symbol_seq[-1]))
        mass_seq = list()
        _mass = anchor
        for item in base_seq:
            _mass += item
            stringency = self._bases.ppm2dm(_mass, 1)
            noise = random.uniform(-1 * stringency, stringency)
            mass_seq.append(round(_mass + noise, 4))

        return mass_seq

    def _make_rt_seq(self, mass_seq):
        """make rt sequences, by math.log(), 
        taking mass as input x, rt as output y
        """
        n = random.uniform(1.2, 1.6)
        return [log(x, n) for x in mass_seq]

    def _make_vol_seq(self, length):
        avg_vol = random.randint(15e5, 30e5)
        sum_vol = avg_vol * length
        rng = np.random.default_rng()
        percent = random.uniform(0.5, 0.8)
        vol_seq = rng.multinomial(percent*sum_vol, [1/(1.0*(length-1))]*(length-1))
        vol_seq = list(vol_seq)
        vol_seq = [random.uniform(0, 1) + item for item in vol_seq]
        vol_seq.append((1-percent)*sum_vol)
        return vol_seq

    def make_salts(self, base_salts, cpds):
        salt_cpds = list()
        for cpd in cpds:
            salt_num = random.randint(0, 10)
            for _ in range(salt_num):
                salt_cpd = Compound()
                salt_cpd.mass = cpd.mass + random.choice(base_salts) + random.uniform(10e-4, 10e-3)
                salt_cpd.rt = cpd.rt + random.uniform(-1, 1)
                salt_cpd.vol = cpd.vol * random.uniform(0.8, 1.2)
                salt_cpds.append(salt_cpd)

        return salt_cpds

    def make_effective_noises(self, cpds):
        """make effective noises.
        generate N(between [0, 50]) nodes
        mass between ppm [-50, 50]
        RT between [-0.2, 0.2]
        Vol between [80%, 120%]
        """
        effective_cpds = list()
        for cpd in cpds:
            effective_num = random.randint(1, 5)
            for _ in range(effective_num):
                effective_cpd = Compound()
                stringency = self._bases.ppm2dm(cpd.mass, 50)
                noise = random.uniform(-1 * stringency, stringency)
                effective_cpd.mass = cpd.mass + noise
                effective_cpd.rt = cpd.rt + random.uniform(-0.2, 0.2)
                effective_cpd.vol = cpd.vol * random.uniform(0.8, 1.2)
                effective_cpds.append(effective_cpd)

        return effective_cpds

    def _extrem_values(self, cpds):
        """ find out the extrem values from compounds datasets.
        max and min values of mass, RT, Vol, etc.
        """
        feature = CompoundFeature()
        feature.min_mass = sys.maxsize
        feature.max_mass = 0
        feature.min_rt = sys.maxsize
        feature.max_rt = 0
        feature.min_vol = sys.maxsize
        feature.max_vol = 0

        for cpd in cpds:
            if cpd.mass < feature.min_mass:
                feature.min_mass = cpd.mass
            if cpd.mass > feature.max_mass:
                feature.max_mass = cpd.mass

            if cpd.rt < feature.min_rt:
                feature.min_rt = cpd.rt
            if cpd.rt > feature.max_rt:
                feature.max_rt = cpd.rt

            if cpd.vol < feature.min_vol:
                feature.min_vol = cpd.vol
            if cpd.vol > feature.max_vol:
                feature.max_vol = cpd.vol

        return feature

    def _make_random_compounds(self, N, feature):
        """ generate N nodes
        mass [0.9 * min_mass, 1.1 * max_mass]
        RT [0.9 * min_rt, 1.1 * max_rt]
        Vol between [10w, 400w]
        """
        ineffective_cpds = list()
        for _ in range(N):
            ineffective_cpd = Compound()
            ineffective_cpd.mass = random.uniform(feature.min_mass * 0.9, feature.max_mass * 1.1)
            ineffective_cpd.rt = random.uniform(feature.min_rt * 0.9, feature.max_rt * 1.1)
            ineffective_cpd.vol = random.uniform(10e4, 400e4)
            ineffective_cpds.append(ineffective_cpd)

        return ineffective_cpds

    def make_ineffective_noises_inside_zone(self, cpds):
        """make ineffective noises, actually its random scatters.
        generate N(between [50, 200]) nodes
        mass [0.9 * min_mass, 1.1 * max_mass]
        RT [0.9 * min_rt, 1.1 * max_rt]
        Vol between [10w, 400w]
        """
        feature = self._extrem_values(cpds)
        ms_list = [cpd.mass for cpd in cpds]
        def func_high(x):
            return log(x, 1.15)
        def func_low(x):
            return log(x, 1.7)
        ineffective_num = random.randint(50, 200)
        ineffective_cpds = list()
        for _ in range(ineffective_num):
            ineffective_cpd = Compound()
            ineffective_cpd.mass = random.uniform(feature.min_mass * 0.9, feature.max_mass * 1.1)
            ineffective_cpd.rt = random.uniform(func_low(ineffective_cpd.mass) * 0.9, func_high(ineffective_cpd.mass) * 1.1)
            ineffective_cpd.vol = random.uniform(10e4, 400e4)
            ineffective_cpds.append(ineffective_cpd)

        return ineffective_cpds

    def make_ineffective_noises_in_lower_zone(self, cpds):
        """make ineffective noises, actually its random scatters.
        generate N(between [100, 300]) nodes
        """
        feature = self._extrem_values(cpds)
        feature.max_rt = 0.35 * (feature.max_rt - feature.min_rt) + feature.min_rt
        ineffective_num = random.randint(100, 300)
        ineffective_cpds = list()
        return self._make_random_compounds(ineffective_num, feature)

    def make_ineffective_noises_in_upper_zone(self, cpds):
        """make ineffective noises, actually its random scatters.
        generate N(between [50, 300]) nodes
        """
        feature = self._extrem_values(cpds)
        feature.min_rt = 0.5* (feature.max_rt - feature.min_rt) + feature.min_rt
        ineffective_num = random.randint(50, 300)
        return self._make_random_compounds(ineffective_num, feature)

    def make_ineffective_noises_in_left_zone(self, cpds):
        """make ineffective noises, actually its random scatters.
        generate N(between [100, 300]) nodes
        """
        feature = self._extrem_values(cpds)
        feature.max_mass = feature.min_mass
        feature.min_mass = 0
        ineffective_num = random.randint(100, 300)
        return self._make_random_compounds(ineffective_num, feature)

    def _make_col_width(self, cpds):
        def func(x):
            y = 0.08 * x + 0.09
            return y

        for cpd in cpds:
            cpd.width = func(cpd.rt)

        return cpds

    def _make_col_qs(self, cpds):
        for cpd in cpds:
            cpd.qs = random.randint(50, 100)
            if cpd.vol > 5E6:
                cpd.qs = 100

        return cpds

    def to_excel(self, cpds, fpath):
        df = Compound.cpds2df(cpds)
        df.to_excel(fpath)
