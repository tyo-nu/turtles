# Tests for functions in turtles_utils.py
#
# Sequences for test files taken from Cobalt experiments (first 1000 seqs
# from each file) for 0 Control, 20 min timepoint, and 1 Control conditions


import os
from math import isclose

import numpy as np
import pandas as pd
import pytest

from turtles.turtles_utils import (
    bin_seq, calc_aitchison_distance, calc_cum_length_dists,
    calc_norm_len_base_pcts, calc_signal, calc_switch_bin, calc_switch_bins,
    calc_switch_time_w_rates, calc_switch_times, closure, clr, cutoff_float,
    enumerate_indices_from_counts, generate_aitch_df, get_length_dists,
    get_norm_len_base_counts, get_seq_intervals, get_total_base_pcts,
    parse_fastq_lengths, read_seqs)


tests_path = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(tests_path, 'data/')

zero_control_conds = ['test_C0_rep1_reads', 'test_C0_rep2_reads']
one_control_conds = ['test_C1_rep1_reads', 'test_C1_rep2_reads']

condition_dict = {'test_C01_rep1_reads': '01 at 20 min',
                  'test_C01_rep2_reads': '01 at 20 min',
                  'test_C0_rep1_reads': '0 Control',
                  'test_C0_rep2_reads': '0 Control',
                  'test_C1_rep1_reads': '1 Control',
                  'test_C1_rep2_reads': '1 Control'}

rep_dict = {'test_C01_rep1_reads': 1,
            'test_C01_rep2_reads': 2,
            'test_C0_rep1_reads': 1,
            'test_C0_rep2_reads': 2,
            'test_C1_rep1_reads': 1,
            'test_C1_rep2_reads': 2}


@pytest.fixture(scope='module')
def default_seqs():
    """Default sequence input used in all tests requiring sequence data.
    Exceptions include tests in the TestSeqIO class."""

    seqs = read_seqs(data_dir, filename_end='trimmed.fq', cutoff=5.8)

    return seqs


@pytest.fixture(scope='module')
def default_counts_dict(default_seqs):
    """Default binned base counts. Generated from defaults_seqs."""

    num_bins = 1000

    counts_dict = get_norm_len_base_counts(default_seqs, num_bins=num_bins)

    return counts_dict


@pytest.fixture(scope='module')
def default_pcts_dict(default_counts_dict):
    """Default base count frequencies. Generated from default_counts_dict."""

    pcts_dict = calc_norm_len_base_pcts(default_counts_dict)

    return pcts_dict


@pytest.fixture(scope='module')
def default_aitch_df(default_pcts_dict):
    """Default DataFrame with calculated Aitchison distances."""

    clr_data = calc_aitchison_distance(default_pcts_dict)

    df = generate_aitch_df(default_pcts_dict, clr_data, condition_dict,
                           rep_dict, zero_control_conds, one_control_conds)

    return df


@pytest.fixture(scope='module')
def default_signal_df(default_aitch_df):
    """Default DataFrame with calculated signal values."""

    df = calc_signal(default_aitch_df, zero_control_conds, one_control_conds)

    return df


@pytest.fixture(scope='module')
def default_len_dists(default_seqs):
    """Default length distributions."""

    len_dists = get_length_dists(default_seqs)

    return len_dists


@pytest.fixture(scope='module')
def default_averages_df(default_len_dists):
    """Default DataFrame with average information on length distributions (such
    as mean length, length standard deviation, rate, etc.)."""

    averages = parse_fastq_lengths(default_len_dists)

    return averages


class TestSeqIO(object):
    """Tests for reading NGS sequences using turtles_utils.read_seqs() with
    various parameters. Also tests turtles_utils.get_seq_intervals()."""

    @pytest.fixture(scope="class")
    def seqs(self):
        return read_seqs(data_dir, filename_end='trimmed.fq')

    def test_all_dirs_found(self, seqs):
        assert len(seqs) == 6

    def test_seqs_read(self, seqs):
        for condition in seqs:
            assert seqs[condition]

    def test_all_seqs_read(self, seqs):
        for condition in seqs:
            assert len(seqs[condition]) == 1000

    def test_first_seqs_correct(self, seqs):
        first_seqs = ['CTAGGCAC',
                      'GTAGGAGAATGAG',
                      'GAGGTAAAATAGCTGGAGGATC',
                      'AAGTGACGGGGGATGCAGAGGGGTTGTCC',
                      'GATCGAGTGAAAGGATAGGAAGGTCA',
                      'CGGGTGAGTCGAGACGAGAAGGATA']

        for i, condition in enumerate(seqs):
            first_seq = list(first_seqs[i])
            assert seqs[condition][0] == first_seq

    def test_last_seqs_correct(self, seqs):
        last_seqs = ['CGAGCCTAAAAC',
                     'CGGGCGGGTAGCGAA',
                     'CTGAGACGGGAGGGAGG',
                     'CCGCGCTAGAGCG',
                     'AAGGGATAGCGCTAAA',
                     'AAAGCAGATGACGAGATAAG']

        for i, condition in enumerate(seqs):
            last_seq = list(last_seqs[i])
            assert seqs[condition][-1] == last_seq

    def test_degen(self):
        seqs_degen = read_seqs(data_dir, filename_end='trimmed.fq', degen=5)

        first_seqs = ['CAC',
                      'AGAATGAG',
                      'AAAATAGCTGGAGGATC',
                      'ACGGGGGATGCAGAGGGGTTGTCC',
                      'AGTGAAAGGATAGGAAGGTCA',
                      'GAGTCGAGACGAGAAGGATA']

        for i, condition in enumerate(seqs_degen):
            first_seq = list(first_seqs[i])
            assert seqs_degen[condition][0] == first_seq

    def test_cutoff_int(self):
        seqs_cutoff_int = read_seqs(data_dir, filename_end='trimmed.fq',
                                    cutoff=5.0)

        first_seqs = ['CTA',
                      'GTAGGAGA',
                      'GAGGTAAAATAGCTGGA',
                      'AAGTGACGGGGGATGCAGAGGGGT',
                      'GATCGAGTGAAAGGATAGGAA',
                      'CGGGTGAGTCGAGACGAGAA']

        for i, condition in enumerate(seqs_cutoff_int):
            first_seq = list(first_seqs[i])
            assert seqs_cutoff_int[condition][0] == first_seq

    def test_cutoff_float(self):
        seq = 'ACGTACGT'
        seq = cutoff_float(seq, 1.5)

        possibilities = ['ACGTACG', 'ACGTAC']

        assert seq in possibilities

    def test_cutoff_float_all(self, default_seqs):
        first_seqs = [('CTA', 'CT'),
                      ('GTAGGAGA', 'GTAGGAG'),
                      ('GAGGTAAAATAGCTGGA', 'GAGGTAAAATAGCTGG'),
                      ('AAGTGACGGGGGATGCAGAGGGGT', 'AAGTGACGGGGGATGCAGAGGGG'),
                      ('GATCGAGTGAAAGGATAGGAA', 'GATCGAGTGAAAGGATAGGA'),
                      ('CGGGTGAGTCGAGACGAGAA', 'CGGGTGAGTCGAGACGAGA')]

        for i, condition in enumerate(default_seqs):
            first_seq_possibilities = (list(seq) for seq in first_seqs[i])
            assert default_seqs[condition][0] in first_seq_possibilities

    def test_degen_cutoff(self):
        seqs_degen_cutoff = read_seqs(data_dir, filename_end='trimmed.fq',
                                      degen=5, cutoff=5.8)

        first_seqs = [('G', 'GGAAACGGGAAAGCTAAATCAAGAGA',
                       'GGAAACGGGAAAGCTAAATCAAGAG'),  # first seq too short
                      ('AGA', 'AG'),
                      ('AAAATAGCTGGA', 'AAAATAGCTGG'),
                      ('ACGGGGGATGCAGAGGGGT', 'ACGGGGGATGCAGAGGGG'),
                      ('AGTGAAAGGATAGGAA', 'AGTGAAAGGATAGGA'),
                      ('GAGTCGAGACGAGAA', 'GAGTCGAGACGAGA')]

        for i, condition in enumerate(seqs_degen_cutoff):
            first_seq_possibilities = (list(seq) for seq in first_seqs[i])
            assert seqs_degen_cutoff[condition][0] in first_seq_possibilities

    def test_exact_length(self):
        seqs_len10 = read_seqs(data_dir, filename_end='trimmed.fq',
                               seq_len_req=10)

        for condition in seqs_len10:
            for seq in seqs_len10[condition]:
                assert len(seq) == 10

    def test_data_random_sampling(self):
        seqs = read_seqs(data_dir, filename_end='trimmed.fq', p_discard=0.5)

        for condition in seqs:
            # With 1000 seqs, chance of having less than 400 or more than 600
            # seqs by chance is low (at least six standard devs from mean).
            assert len(seqs[condition]) > 400
            assert len(seqs[condition]) < 600

    def test_get_seqs_for_one_specific_condition(self):
        seqs = read_seqs(data_dir, filename_end='trimmed.fq',
                         cond_text='C0_rep1')

        cond_name = 'test_C0_rep1_reads'

        assert len(seqs) == 1
        assert cond_name in seqs

    def test_get_seqs_for_multiple_specific_conditions(self):
        seqs = read_seqs(data_dir, filename_end='trimmed.fq',
                         cond_text=['C0_rep1', 'C1_rep1'])

        cond1_name = 'test_C0_rep1_reads'
        cond2_name = 'test_C1_rep1_reads'

        assert len(seqs) == 2
        assert cond1_name in seqs
        assert cond2_name in seqs

    def test_neg_degen_raises_ValueError(self):
        with pytest.raises(ValueError):
            read_seqs(data_dir, filename_end='trimmed.fq', degen=-1)

    def test_neg_cutoff_raises_ValueError(self):
        with pytest.raises(ValueError):
            read_seqs(data_dir, filename_end='trimmed.fq', cutoff=-1)

    def test_neg_p_discard_raises_ValueError(self):
        with pytest.raises(ValueError):
            read_seqs(data_dir, filename_end='trimmed.fq', p_discard=-0.5)

    def test_p_discard_greater_than_1_raises_ValueError(self):
        with pytest.raises(ValueError):
            read_seqs(data_dir, filename_end='trimmed.fq', p_discard=1.5)

    def test_get_seq_intervals_single(self):
        seq_dict = {'test': ['ACGTACGTACGTACG']}

        seq_parsed_begin = get_seq_intervals(seq_dict, [(1, 5)])[0]['test'][0]
        seq_parsed_mid = get_seq_intervals(seq_dict, [(6, 10)])[0]['test'][0]
        seq_parsed_end = get_seq_intervals(seq_dict, [(11, 15)])[0]['test'][0]

        assert seq_parsed_begin == 'ACGTA'
        assert seq_parsed_mid == 'CGTAC'
        assert seq_parsed_end == 'GTACG'

    def test_get_seq_intervals_multiple(self):
        seq_dict = {'test': ['ACGTACGTACGTACG']}

        seq_parsed_begin_end = get_seq_intervals(seq_dict, [(1, 5), (11, 15)])

        assert seq_parsed_begin_end == [{'test': ['ACGTA']},
                                        {'test': ['GTACG']}]

    def test_get_seq_intervals_pos_int(self):
        seq_dict = {'test': ['ACGTACGTACGTACG']}

        seq_parsed = get_seq_intervals(seq_dict, [11])[0]['test'][0]

        assert seq_parsed == 'GTACG'

    def test_get_seq_intervals_neg_int(self):
        seq_dict = {'test': ['ACGTACGTACGTACG']}

        seq_parsed = get_seq_intervals(seq_dict, [-5])[0]['test'][0]

        assert seq_parsed == 'GTACG'

    def test_get_seq_intervals_combined(self):
        seq_dict = {'test': ['ACGTACGTACGTACG']}

        seq_parsed_begin_end = get_seq_intervals(seq_dict,
                                                 [(1, 5), (11, 15), 11, -5])

        assert seq_parsed_begin_end == [{'test': ['ACGTA']},
                                        {'test': ['GTACG']},
                                        {'test': ['GTACG']},
                                        {'test': ['GTACG']}]


class TestSeqBinning(object):
    """Tests for functions in turtles_utils.py that bin sequence data."""

    def test_binning_divisible(self):
        seq = list('ACAAGATAAA')
        n_bins = 20

        bins = np.arange(0, 1, 1 / n_bins)

        a_inds = [(0, 2), (4, 8), (10, 12), (14, 20)]
        c_inds = [(2, 4)]
        g_inds = [(8, 10)]
        t_inds = [(12, 14)]

        binned_seq_inds = bin_seq(seq, bins)

        assert binned_seq_inds['A'] == a_inds
        assert binned_seq_inds['C'] == c_inds
        assert binned_seq_inds['G'] == g_inds
        assert binned_seq_inds['T'] == t_inds

    def test_binning_indivisible(self):
        seq = list('AACGT')
        n_bins = 8

        bins = np.arange(0, 1, 1 / n_bins)

        # Each base gets binned into each bin with a left edge in its
        # normalized interval.
        # For example, A would have a normalized interval of (0, 0.4) since the
        # first two bases out of five are A. Because this includes the values
        # (0/8, 1/8, 2/8, 3/8) = (0, 0.125, 0.25, 0.375), A gets the first four
        # bins.
        a_inds = [(0, 4)]
        c_inds = [(4, 5)]
        g_inds = [(5, 7)]
        t_inds = [(7, 8)]

        binned_seq_inds = bin_seq(seq, bins)

        assert binned_seq_inds['A'] == a_inds
        assert binned_seq_inds['C'] == c_inds
        assert binned_seq_inds['G'] == g_inds
        assert binned_seq_inds['T'] == t_inds

    def test_binning_seq_len_1(self):
        seq = list('A')
        n_bins = 10

        bins = np.arange(0, 1, 1 / n_bins)

        a_inds = [(0, 10)]
        c_inds = []
        g_inds = []
        t_inds = []

        binned_seq_inds = bin_seq(seq, bins)

        assert binned_seq_inds['A'] == a_inds
        assert binned_seq_inds['C'] == c_inds
        assert binned_seq_inds['G'] == g_inds
        assert binned_seq_inds['T'] == t_inds

    def test_bin_counts_all_seqs(self, default_seqs, default_counts_dict):
        num_bins = 1000

        assert default_counts_dict

        for condition in default_counts_dict:
            a_counts = default_counts_dict[condition]['A']
            c_counts = default_counts_dict[condition]['C']
            g_counts = default_counts_dict[condition]['G']
            t_counts = default_counts_dict[condition]['T']
            n_counts = default_counts_dict[condition]['N']

            num_seqs = len(default_seqs[condition])

            assert len(a_counts) == num_bins
            assert len(c_counts) == num_bins
            assert len(g_counts) == num_bins
            assert len(t_counts) == num_bins
            assert len(n_counts) == num_bins

            total_counts = sum([sum(a_counts), sum(c_counts), sum(g_counts),
                                sum(t_counts), sum(n_counts)])

            assert total_counts == num_bins * num_seqs

    def test_bin_counts_sample_seqs(self):
        sample_seqs = ['CAGG', 'TGGCAA', 'AAACCCGGGTTT', 'TGGT']
        seqs_dict = {'test': sample_seqs}

        counts_dict = get_norm_len_base_counts(seqs_dict, num_bins=12)

        a_counts = [1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1]
        c_counts = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
        g_counts = [0, 0, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1]
        t_counts = [2, 2, 1, 0, 0, 0, 0, 0, 0, 2, 2, 2]

        assert list(counts_dict['test']['A']) == a_counts
        assert list(counts_dict['test']['C']) == c_counts
        assert list(counts_dict['test']['G']) == g_counts
        assert list(counts_dict['test']['T']) == t_counts


class TestBaseFrequencyCalcs(object):
    """Tests for functions in turtles.turtles_utils.py that calculate base
    frequency."""

    def test_all_percents_add_to_100(self, default_pcts_dict):

        for condition in default_pcts_dict:
            a_pcts = default_pcts_dict[condition]['A']
            c_pcts = default_pcts_dict[condition]['C']
            g_pcts = default_pcts_dict[condition]['G']
            t_pcts = default_pcts_dict[condition]['T']

            n_bins = len(a_pcts)

            for i in range(n_bins):
                total = sum([a_pcts[i], c_pcts[i], g_pcts[i], t_pcts[i]])

                # May not be exactly 100% (e.g. 99.99999999%) due to precision
                # of floating point numbers
                assert isclose(total, 1.00, rel_tol=0.000001)

    def test_calculated_percents_correct(self):
        counts_dict = {'test': {'A': [12, 15, 13], 'C': [20, 11, 2],
                                'G': [15, 14, 5], 'T': [3, 0, 5]}}

        correct_pcts = {'A': [0.24, 0.375, 0.52], 'C': [0.40, 0.275, 0.08],
                        'G': [0.30, 0.35, 0.20], 'T': [0.06, 0.0, 0.20]}

        pcts_dict = calc_norm_len_base_pcts(counts_dict)

        assert list(pcts_dict['test']['A']) == correct_pcts['A']
        assert list(pcts_dict['test']['C']) == correct_pcts['C']
        assert list(pcts_dict['test']['G']) == correct_pcts['G']
        assert list(pcts_dict['test']['T']) == correct_pcts['T']

    def test_calculated_percent_correct_with_n(self):
        counts_dict = {'test': {'A': [12, 15, 6], 'C': [20, 11, 4],
                                'G': [15, 14, 7], 'T': [3, 0, 3],
                                'N': [0, 0, 5]}}

        correct_n_pcts = [0.0, 0.0, 0.20]

        pcts_dict = calc_norm_len_base_pcts(counts_dict, exclude_n=False)

        assert list(pcts_dict['test']['N'] == correct_n_pcts)

    def test_all_total_percents_add_to_100(self, default_seqs):
        pcts_dict = get_total_base_pcts(default_seqs)

        for condition in pcts_dict:
            a_pct = pcts_dict[condition]['A']
            c_pct = pcts_dict[condition]['C']
            g_pct = pcts_dict[condition]['G']
            t_pct = pcts_dict[condition]['T']

            total = sum([a_pct, c_pct, g_pct, t_pct])

            assert isclose(total, 1.00, rel_tol=0.000001)

    def test_total_percents_correct(self):
        seqs = {'test': ['ACGT', 'AAAA', 'GAGA', 'TATA', 'CCCC']}

        correct_pcts = {'A': 0.45, 'C': 0.25, 'G': 0.15, 'T': 0.15}

        pcts_dict = get_total_base_pcts(seqs)
        total_pcts = pcts_dict['test']

        assert total_pcts['A'] == correct_pcts['A']
        assert total_pcts['C'] == correct_pcts['C']
        assert total_pcts['G'] == correct_pcts['G']
        assert total_pcts['T'] == correct_pcts['T']

    def test_all_total_percents_add_to_100_with_n(self, default_seqs):
        pcts_dict = get_total_base_pcts(default_seqs, exclude_n=False)

        for condition in pcts_dict:
            a_pct = pcts_dict[condition]['A']
            c_pct = pcts_dict[condition]['C']
            g_pct = pcts_dict[condition]['G']
            t_pct = pcts_dict[condition]['T']
            n_pct = pcts_dict[condition]['N']

            total = sum([a_pct, c_pct, g_pct, t_pct, n_pct])

            assert isclose(total, 1.00, rel_tol=0.000001)

    def test_total_percents_correct_with_n(self):
        seqs = {'test': ['ACGT', 'AAAA', 'GAGA', 'TNTA', 'CNNC']}

        correct_pcts = {'A': 0.40, 'C': 0.15, 'G': 0.15, 'T': 0.15, 'N': 0.15}

        pcts_dict = get_total_base_pcts(seqs, exclude_n=False)
        total_pcts = pcts_dict['test']

        assert total_pcts['A'] == correct_pcts['A']
        assert total_pcts['C'] == correct_pcts['C']
        assert total_pcts['G'] == correct_pcts['G']
        assert total_pcts['T'] == correct_pcts['T']
        assert total_pcts['N'] == correct_pcts['N']


class TestAitchisonDistanceCalcs(object):
    """Tests for calculation of Aitchison distance."""

    def test_closure(self):
        x = [3, 2, 5]

        ans = closure(x)

        correct_ans = [0.3, 0.2, 0.5]

        assert list(ans) == correct_ans

    def test_clr(self):
        x = np.array([[3, 2, 5],
                      [2, 1, 2]])

        clr_data = clr(x)

        correct_vals = [[-0.0351, -0.4405, 0.4757],
                        [0.2310, -0.4621, 0.2310]]

        for i, composition_vals in enumerate(clr_data):
            for j, val in enumerate(composition_vals):
                assert isclose(val, correct_vals[i][j], abs_tol=0.001)

    def test_aitch_dist(self):
        pcts_dict = {'test': {'A': [0.24, 0.375, 0.52],
                              'C': [0.40, 0.175, 0.08],
                              'G': [0.30, 0.35, 0.20],
                              'T': [0.06, 0.10, 0.20]}}

        clr_data = calc_aitchison_distance(pcts_dict)
        clr_data = clr_data['test']

        correct_clr = [[0.1630813, 0.53822219, 0.94570627],
                       [0.67390692, -0.22391786, -0.92609591],
                       [0.38622485, 0.46922932, -0.00980518],
                       [-1.22321306, -0.78353365, -0.00980518]]

        for i, base in enumerate(clr_data):
            base_clr = list(clr_data[base])

            for j in range(len(base_clr)):
                assert isclose(base_clr[j], correct_clr[i][j], abs_tol=0.001)


class TestSignalDataFrameGeneration(object):
    """Tests on the dataframe generated from signal calculations."""

    @pytest.fixture(scope='class')
    def default_aitch_df(self, default_pcts_dict):
        clr_data = calc_aitchison_distance(default_pcts_dict)

        df = generate_aitch_df(default_pcts_dict, clr_data, condition_dict,
                               rep_dict, zero_control_conds, one_control_conds)

        return df

    @pytest.fixture(scope='class')
    def default_signal_df(self, default_aitch_df):
        df = calc_signal(default_aitch_df, zero_control_conds,
                         one_control_conds)

        return df

    def test_aitch_df_columns(self, default_aitch_df):
        columns = ['Directory', 'Condition', 'Replicate', 'Bin Number',
                   'Aitch Dist (from 0)', 'Aitch Dist (from 1)', 'A % Aitch',
                   'C % Aitch', 'G % Aitch', 'T % Aitch',
                   'A % Aitch Diff from 0', 'C % Aitch Diff from 0',
                   'G % Aitch Diff from 0', 'T % Aitch Diff from 0',
                   'A % Aitch Diff from 1', 'C % Aitch Diff from 1',
                   'G % Aitch Diff from 1', 'T % Aitch Diff from 1', 'A %',
                   'C %', 'G %', 'T %']

        assert list(default_aitch_df.columns) == columns

    def test_aitch_df_length(self, default_pcts_dict, default_aitch_df):
        correct_df_length = len(default_pcts_dict) * 1000

        assert len(default_aitch_df) == correct_df_length

    def test_aitch_df_01_control_first_row_vals(self, default_aitch_df):
        first_row = default_aitch_df.loc[0, :]
        last_row = default_aitch_df.loc[1999, :]

        assert first_row['Condition'] == '01 at 20 min'
        assert first_row['Directory'] == 'test_C01_rep1_reads'
        assert first_row['Replicate'] == 1
        assert first_row['Bin Number'] == 1

        assert last_row['Replicate'] == 2
        assert last_row['Bin Number'] == 1000

        assert first_row['Aitch Dist (from 0)'] \
            < first_row['Aitch Dist (from 1)']

        assert first_row['Aitch Dist (from 1)'] \
            > last_row['Aitch Dist (from 1)']

    def test_aitch_df_0_control_first_row_vals(self, default_aitch_df):
        first_row = default_aitch_df.loc[2000, :]

        assert first_row['Condition'] == '0 Control'
        assert first_row['Directory'] == 'test_C0_rep1_reads'
        assert first_row['Replicate'] == 1
        assert first_row['Bin Number'] == 1
        assert first_row['Aitch Dist (from 0)'] \
            < first_row['Aitch Dist (from 1)']

    def test_aitch_df_1_control_first_row_vals(self, default_aitch_df):
        first_row = default_aitch_df.loc[4000, :]

        assert first_row['Condition'] == '1 Control'
        assert first_row['Directory'] == 'test_C1_rep1_reads'
        assert first_row['Replicate'] == 1
        assert first_row['Bin Number'] == 1
        assert first_row['Aitch Dist (from 0)'] \
            > first_row['Aitch Dist (from 1)']

    def test_signal_df_columns(self, default_signal_df):
        columns = ['Directory', 'Condition', 'Replicate', 'Bin Number',
                   'Aitch Dist (from 0)', 'Aitch Dist (from 1)', 'A % Aitch',
                   'C % Aitch', 'G % Aitch', 'T % Aitch',
                   'A % Aitch Diff from 0', 'C % Aitch Diff from 0',
                   'G % Aitch Diff from 0', 'T % Aitch Diff from 0',
                   'A % Aitch Diff from 1', 'C % Aitch Diff from 1',
                   'G % Aitch Diff from 1', 'T % Aitch Diff from 1', 'A %',
                   'C %', 'G %', 'T %', 'Aitch Fraction', 'Signal']

        assert list(default_signal_df.columns) == columns

    def test_signal_df_signal_vals(self, default_signal_df):
        for directory in default_signal_df.Directory.unique():
            df = default_signal_df[default_signal_df.Directory == directory]
            signal_col = np.array(df.loc[:, 'Signal'])

            if condition_dict[directory] == '0 Control':
                for signal in signal_col:
                    assert isclose(signal, 0, abs_tol=0.20)

            elif condition_dict[directory] == '1 Control':
                for signal in signal_col:
                    assert isclose(signal, 1, abs_tol=0.20)

    def test_signal_calc(self):
        df = pd.DataFrame()

        dir_col = ['test0'] * 5 + ['test01'] * 5 + ['test1'] * 5
        cond_col = ['0 Control'] * 5 + ['01 at 30 min'] * 5 + ['1 Control'] * 5
        rep_col = [1] * 15
        aitch0_col = [0, 0, 0, 0, 0, 0.0, 0.24, 0.50, 0.75, 1.0, 1, 1, 1, 1, 1]
        aitch1_col = [1, 1, 1, 1, 1, 1.0, 0.72, 0.50, 0.25, 0.0, 0, 0, 0, 0, 0]

        df['Directory'] = dir_col
        df['Condition'] = cond_col
        df['Replicate'] = rep_col
        df['Aitch Dist (from 0)'] = aitch0_col
        df['Aitch Dist (from 1)'] = aitch1_col

        df = calc_signal(df, zero_control_conds, one_control_conds)

        correct_signal = aitch0_col  # aitch0 + aitch1 = 1
        correct_signal[6] = 0.25     # except for at position 6

        assert list(df['Signal']) == correct_signal


class TestLengthDistCalcs(object):
    """Tests for functions calculating and processing sequence length
    distributions."""

    def test_get_length_dists_all_seqs(self, default_seqs, default_len_dists):
        assert len(default_len_dists) == 6

        for condition in default_len_dists:
            len_dists = default_len_dists[condition]

            assert len(len_dists) == 201  # default is 0 to 200 max length

            total_seqs_correct = len(default_seqs[condition])
            total_seqs_from_lens = sum(len_dists)

            assert total_seqs_from_lens == total_seqs_correct

    def test_get_length_dists_calc(self):
        seqs = {'test': ['A', 'C', 'T', 'AA', 'AC', 'G', 'TTT', 'A', '', 'AA']}

        correct_len_dist = [1, 5, 3, 1]

        len_dist = list(get_length_dists(seqs, max_len=3)['test'])

        assert len(len_dist) == 4
        assert len_dist == correct_len_dist

    def test_cum_length_dists_calc(self):
        seqs = {'test': ['A', 'C', 'T', 'AA', 'AC', 'G', 'TTT', 'A', '', 'AA']}

        len_dists = get_length_dists(seqs, max_len=3)

        correct_cum_len_dist = [0.1, 0.6, 0.9, 1.0]

        cum_len_dist = list(calc_cum_length_dists(len_dists)['test'])

        for val, correct_val in zip(cum_len_dist, correct_cum_len_dist):
            assert isclose(val, correct_val, rel_tol=0.000001)

    def test_parse_fastq_lengths_all_seqs(self, default_len_dists,
                                          default_averages_df):
        correct_columns = ['Directory', 'Mean', 'Std Devs', 'Rate (nt/min)']

        columns = list(default_averages_df.columns)

        assert columns == correct_columns

        num_dirs = len(default_len_dists)
        num_rows = num_dirs

        assert len(default_averages_df) == num_rows

    def test_parse_fastq_lengths_calcs(self):
        len_dists = {'test': [1, 5, 3, 1]}

        correct_mean = 1.4
        correct_std = 0.8
        correct_rate = 0.07

        averages = parse_fastq_lengths(len_dists, expt_time=20)

        mean = averages['Mean'].values[0]
        std = averages['Std Devs'].values[0]
        rate = averages['Rate (nt/min)'].values[0]

        assert mean == correct_mean
        assert std == correct_std
        assert isclose(rate, correct_rate)

    def test_enumerate_indices_from_counts(self):
        counts = [1, 5, 3, 1]

        correct_inds = [0, 1, 1, 1, 1, 1, 2, 2, 2, 3]

        inds = enumerate_indices_from_counts(counts)

        assert inds == correct_inds


class TestTimepointCalcs(object):
    """Makes sure timepoint calculations are correct."""

    def test_calc_switch_bins_all_seqs(self, default_signal_df,
                                       default_averages_df):
        averages = calc_switch_bins(default_averages_df, default_signal_df)

        assert 'Switch Bin' in list(averages.columns)

        switch_bins = averages['Switch Bin'].values

        sb_rep1 = switch_bins[0]
        sb_rep2 = switch_bins[1]

        # For switching at 20 min, switch bin should be around 300ish
        assert 100 < sb_rep1
        assert 100 < sb_rep2
        assert sb_rep1 < 500
        assert sb_rep2 < 500

    def test_calc_switch_bin_01(self):
        signals = [0, 0.12, 0.24, 0.31, 0.42, 0.52, 0.66, 0.77, 0.89, 0.94, 1]

        correct_switch_bin = 5.8

        switch_bin = calc_switch_bin(signals)

        assert switch_bin == correct_switch_bin

    def test_calc_switch_bin_10(self):
        signals = [1, 0.91, 0.8, 0.6, 0.44, 0.4, 0.31, 0.29, 0.2, 0.12, 0]

        correct_switch_bin = 4.625

        switch_bin = calc_switch_bin(signals, mode='10')

        assert switch_bin == correct_switch_bin

    def test_calc_switch_times_all_seqs(self, default_averages_df,
                                        default_signal_df):
        averages = calc_switch_bins(default_averages_df, default_signal_df)
        averages = calc_switch_times(default_averages_df,
                                     start_control_conds=zero_control_conds,
                                     end_control_conds=one_control_conds)

        assert 'Switch Time' in list(averages.columns)

        times = averages['Switch Time'].values

        time_rep1 = times[0]
        time_rep2 = times[1]

        # Predicted times should be within 5 minutes of 20 min
        assert abs(time_rep1 - 20) < 5
        assert abs(time_rep2 - 20) < 5

    def test_calc_switch_time_w_rates_calcs(self):
        switch_bins = pd.Series([100, 250, 400])
        num_bins = 1000
        t_expt = 60
        r_start = 1
        r_end = 2

        switch_times = calc_switch_time_w_rates(switch_bins, num_bins, t_expt,
                                                r_start, r_end)
        switch_times = list(switch_times)

        correct_switch_times = [3.16, 8.57, 15]

        for i in range(len(switch_times)):
            assert isclose(switch_times[i], correct_switch_times[i],
                           abs_tol=0.01)
