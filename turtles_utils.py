# Utility functions for analysis of TdT NGS data for TURTLES system
#
# Author: Jonathan Strutz (jonathanstrutz2021@u.northwestern.edu)
#
# Note: I copied two functions, closure and clr, from skbio.stats.composition
# to this file since I had issues installing skbio on Windows (I tried both
# conda and pip).

import os
import re

import numpy as np
import pandas as pd
from Bio import SeqIO


def get_pct_len_base_pair_counts(data_dir, num_bins=1000, seq_len_req=None,
                                 filename_end='trimmed.fq'):
    """Read R1 fastq files in subdirectories in data_dir. In each file, count
    the number of A, C, G, and T at each position. Then normalize counts by
    the sequence's length by binning.

    Parameters
    ----------
    data_dir : str
        Filepath to folder which contains subfolders which each have an R1
        fastq file.
    num_bins : int (default: 1000)
        The number of bins to divide each sequence up into.
    seq_len_req : int (default: 1000)
        If provided, only sequences of length seq_len_req will be processed.
    filename_end : str (default: 'trimmed.fq')
        Suffix of the fastq filenames to read. Should be the same for every
        condition.

    Returns
    -------
    counts_dict : dict
        Contains an array of counts at each bin (indexed) for each condition.
        Structure is {'directory_name': {'Base': [#, #, ..., #]}}.
    """

    counts_dict = {}

    bins = np.arange(0, 1 + 1 / num_bins, 1 / num_bins)

    for directory in os.listdir(data_dir):
        for filename in os.listdir(data_dir + directory):
            if filename.endswith(filename_end) and 'R1' in filename:
                data_file = filename

        filepath = data_dir + directory + '/' + data_file
        if 'R1_Trimmed_Reads' not in filepath and 'No_TdT' not in filepath:

            counts_dict[directory] = {base: np.zeros(num_bins)
                                      for base in ['A', 'C', 'G', 'N', 'T']}

            print('Loading', directory, '\n')
            with open(filepath) as infile:
                fasta_sequences = SeqIO.parse(infile, 'fastq')

                for fasta in fasta_sequences:
                    seq = list(str(fasta.seq))
                    seq_len = len(seq)
                    seq_fracs = [i / seq_len for i in range(seq_len)]

                    if not seq_len_req or seq_len == seq_len_req:

                        bin_ind = np.digitize(seq_fracs, bins) - 1
                        bin_ind = np.append(bin_ind, num_bins - 1)

                        inds = {base: [] for base in ['A', 'C', 'G', 'T', 'N']}
                        prev_base = seq[0]
                        i_init = 0
                        if seq_len > 1:
                            for i, base in enumerate(seq[1:]):
                                i = i + 1
                                if base != prev_base:
                                    inds[prev_base].append((bin_ind[i_init],
                                                            bin_ind[i]))
                                    i_init = i

                                if i == len(seq) - 1:
                                    inds[base].append((bin_ind[i_init],
                                                       num_bins))

                                prev_base = seq[i]

                        else:
                            inds[seq[0]].append((0, num_bins))

                        for base in ['A', 'C', 'G', 'T', 'N']:
                            for (i1, i2) in inds[base]:
                                counts_dict[directory][base][i1:i2] += 1

    return counts_dict


def read_seqs_position(data_dir, positions, filename_end='trimmed.fq',
                       seq_min=1):
    """Read sequences at given base position ranges.

    Parameters
    ----------
    data_dir : str
        Filepath to folder which contains subfolders which each have an R1
        fastq file.
    positions : list of int tuples and/or negative ints
        Base positions to parse each seq. For example, [(1, 5), (11, 15)]
        would return a list of two dicts, one for seqs from bp 1-5, the other
        for seqs from bp 11-15. If just an int, goes from that position to the
        end of the sequence, e.g. -5 would get the last 5 bases.
    filename_end : str (default: 'trimmed.fq')
        Suffix of the fastq filenames to read. Should be the same for every
        condition.
    seq_min : int (default: 1)
        Ignore sequences with length less than seq_min.

    Returns
    -------
    seqs : list of dicts
        Key is condition, value is list of new seqs."""

    seqs = [{} for _ in range(len(positions))]

    for directory in os.listdir(data_dir):
        for filename in os.listdir(os.path.join(data_dir, directory)):
            if filename.endswith(filename_end) and 'R1' in filename:
                data_file = filename

        fastq_filepath = os.path.join(data_dir, directory, data_file)
        if 'R1_Trimmed_Reads' not in fastq_filepath:

            with open(fastq_filepath) as infile:
                print('Loading', data_file, '\n')
                fasta_records = SeqIO.parse(infile, 'fastq')

                fasta_seqs = [str(rec.seq) for rec in fasta_records]
                for i, pos in enumerate(positions):
                    if isinstance(pos, tuple):
                        seqs[i][directory] = [seq[(pos[0] - 1): pos[1]]
                                              for seq in fasta_seqs if
                                              len(seq) >= seq_min]
                    elif isinstance(pos, int) and pos < 0:
                        seqs[i][directory] = [seq[(pos):]
                                              for seq in fasta_seqs if
                                              len(seq) >= seq_min]
                    else:
                        raise Exception('Positions should be negative ints '
                                        'and/or tuples of positive integers.')

    return seqs


def read_seqs_condition(data_dir, cond_text, filename_end='trimmed.fq'):
    """Read sequences for given condition(s).

    Parameters
    ----------
    data_dir : str
        Filepath to folder which contains subfolders which each have an R1
        fastq file.
    cond_text : str
        Text within condition name to search for. Conditions with this text
        will be used. Searches within directory name for condition.
    filename_end : str (default: 'trimmed.fq')
        Suffix of the fastq filenames to read. Should be the same for every
        condition.

    Returns
    -------
    fasta_seqs : list
        List of seqs in given condition."""

    for directory in os.listdir(data_dir):
        if cond_text in directory:
            for filename in os.listdir(os.path.join(data_dir, directory)):
                if filename.endswith(filename_end) and 'R1' in filename:
                    data_file = filename

            fastq_filepath = os.path.join(data_dir, directory, data_file)
            if 'R1_Trimmed_Reads' not in fastq_filepath:

                with open(fastq_filepath) as infile:
                    print('Loading', data_file, '\n')
                    fasta_records = SeqIO.parse(infile, 'fastq')

                    fasta_seqs = [str(rec.seq) for rec in fasta_records]

    return fasta_seqs


def parse_fastq_lengths(data_dir, expt_time=60, min_len=1):
    """For each condition in the overall data_dir, look at each R1 fastq file
    (all sequences, not deduplicated or trimmed based on length), and tally
    the lengths of all the sequences. Then calculate mean, std length.

    Parameters
    ----------
    data_dir : str
        directory of all condition subfolders (which each contain fastq files
        and such).
    expt_time : numeric
        Length of the experiment (in minutes).
    min_len : int
        Ignore sequences with length lower than min_len.

    Returns
    -------
    averages : pd.DataFrame
        Pandas DataFrame with fastq length information for each condition.
        Columns are 'Condition', 'Mean' (length), 'Std Devs' (length), and
        'Rate (nt/min)'.
    """

    conds = []
    means = []
    stds = []

    for cond_folder in os.listdir(data_dir):
        for datafile in os.listdir(os.path.join(data_dir, cond_folder)):
            if datafile.endswith('_trimmed.fq') and 'R1' in datafile:

                datafile_path = os.path.join(data_dir, cond_folder, datafile)
                with open(datafile_path) as infile:
                    fastq = SeqIO.parse(infile, "fastq")

                    lengths = []
                    for rec in fastq:
                        seq = str(rec.seq)
                        if len(seq) >= min_len:
                            lengths.append(len(seq))
                    lengths = np.array(lengths)

                    conds.append(cond_folder)
                    means.append(np.mean(lengths))
                    stds.append(np.std(lengths))

    averages = pd.DataFrame()
    averages['Condition'] = conds
    averages['Mean'] = means
    averages['Std Devs'] = stds

    averages['Rate (nt/min)'] = [mean / expt_time for mean in means]

    return averages


def calc_switch_bins(data_dir, data_df, condition_dict, expt_time=60,
                     mode='01', rep_ind=None):
    """Assumes hat we want to use the normalized 0 --> 1 signal in data_df
    rather than the norm % difference. Assumes 0 --> 1 signal is precalculated
    and tabulated in data_df in 'Signal' column.

    Parameters
    ----------
    data_dir : str
        Path to overall directory containing condition folders.
    data_df : pd.DataFrame
        Should contain columns 'Condition' and 'Signal'.
    condition_dict : dict
        Structure is {directory: condition}
    expt_time : numeric
        Length of the experiment (in minutes).
    mode : str
        Determines which way we are switching (from 0 to 1 or from 1 to 0).
        Also, will look for this string to mark switch conditions. For example,
        if mode == '01', then all rows in data_df with '01' in Condition will
        be assumed to be switching experiments.
    rep_ind : int
        If given, get replicate number by doing condition[rep_ind]. Otherwise,
        uses a regular expression.

    Returns
    -------
    averages : pd.DataFrame
        A DataFrame with columns: Condition, Mean, Std Devs, Switch_Bin.
    """

    averages = parse_fastq_lengths(data_dir, expt_time)

    test_conditions = list(averages.Condition)

    signal_dict = {}
    for condition in test_conditions:
        cond = condition_dict[condition]
        if rep_ind:
            rep = condition[rep_ind]
        else:
            rep = condition[3:][re.search("[0-9]+", condition[3:]).start()]
        sub_data_df = data_df[data_df.Condition == cond]
        sub_data_df = sub_data_df[data_df.Replicate == rep]
        signals = sub_data_df['Signal']
        signal_dict[condition] = np.array(signals)

    averages = calc_switch_bin(averages, test_conditions, signal_dict,
                               mode=mode)

    return averages


def calc_switch_bin(averages, test_conditions, signal_dict, mode='01'):
    """Find two bin interval which contains 0.50 signal. Then interpolate
    between those two bins to find the exact switch bin.

    Parameters
    ----------
    averages : pandas.DataFrame
        Just needs to have all conditions listed in 'Condition' column.
    test_conditions : list
        List of non-control conditions.
    signal_dict : dict
        Structure is {condition: [array of norm diffs]}.
    mode : str
        Determines which way we are switching (from 0 to 1 or from 1 to 0).

    Returns
    -------
    averages : pandas.DataFrame
        Updated with switch base in 'Switch_Bin' column.
    """

    for cond in test_conditions:
        for index, diff in enumerate(signal_dict[cond]):
            curr_diffs = signal_dict[cond]
            if mode == '01' and diff > 0.5:
                switch_bin_upper = index
                switch_bin_lower = index - 1
                m = curr_diffs[switch_bin_upper] - curr_diffs[switch_bin_lower]
                x = switch_bin_lower + (0.5 - curr_diffs[switch_bin_lower]) / m
                switch_bin = x + 1
                averages.loc[averages.Condition == cond, 'Switch_Bin'] \
                    = switch_bin
                break
            elif mode == '10' and diff < 0.5:
                switch_bin_upper = index - 1
                switch_bin_lower = index
                m = curr_diffs[switch_bin_lower] - curr_diffs[switch_bin_upper]
                x = switch_bin_upper + (0.5 - curr_diffs[switch_bin_upper]) / m
                switch_bin = x + 1
                averages.loc[averages.Condition == cond, 'Switch_Bin'] \
                    = switch_bin
                break

    return averages


def calc_switch_times(averages, start_control_conds=None,
                      end_control_conds=None, t_expt=60, num_bins=1000):
    """Calculate switch times for 01 or 10 (does not support more than one
    switch). If no controls for start and end are specified, assumes constant
    rate throughout. However, if controls are specified, assumes constant 0
    rate and constant (but potentially different) 1 rate.

    Parameters
    ----------
    averages : pd.DataFrame
        A DataFrame with columns: Condition, Mean, Std Devs, Switch_Bin.
    start_control_conds : list (default: None)
        Name of the conditions in averages.Condition that are the 1 control if
        switching from 1 --> 0 or the 0 control if switching from 0 --> 1.
    end_control_conds : list (default: None)
        Name of the conditions in averages.Condition that are the 0 control if
        switching from 1 --> 0 or the 1 control if switching from 0 --> 1.
    t_expt : numeric (default: 60)
        Length of the experiment (in minutes).
    num_bins : int (default: 1000)
        The number of bins each sequence is divided up into.

    Returns
    -------
    averages : pd.DataFrame
        A DataFrame with columns: Condition, Mean, Std Devs, Switch_Bin,
        Switch_Time.
    """

    if start_control_conds and end_control_conds:
        r_start_df = averages[averages['Condition'].isin(start_control_conds)]
        r_start = r_start_df['Rate (nt/min)'].mean()

        r_end_df = averages[averages['Condition'].isin(end_control_conds)]
        r_end = r_end_df['Rate (nt/min)'].mean()

        alpha = r_start / r_end

        denom = (num_bins / averages['Switch_Bin'] + alpha - 1)
        averages['Switch_Time'] = (alpha * t_expt) / denom

    else:
        averages['Switch_Time'] = averages['Switch_Bin'] / num_bins * t_expt

    return averages


def calc_0_control_means(pct_dict, zero_control_conds, max_base_n=40):
    """Calculate the average percent dNTP incorporation for the 0 control
    conditions listed in control_conds. Does this for each base A, C, G, and T.
    Averages over time.

    Parameters
    ----------
    pct_dict : dict
        Structure is {condition: {base: pct_list}} where pct_list is a list of
        percent dNTP incorporation by position.
    zero_control_conds : list
        A list of keys (0 conditions) in pct_dict which go into the average
        calculation.
    max_base_n : int
        Base location to calculate up to. If 40, will calculate means based on
        the first 40 bases, for example.

    Returns
    -------
    zero_control_means : dict
        Structure is {base: average_pct}.
    """

    zero_control_means = \
        {base: np.mean([pct_dict[cond][base] for cond in zero_control_conds],
                       axis=0)[0:max_base_n] for base in ['A', 'C', 'G', 'T']}

    return zero_control_means


def calc_cum_length_dist(len_data):
    """Takes length data dataframe and calculates the cumulative length
    distribution from the regular length distribution. Takes into account the
    different distribution for each condition, replicate pair.

    Parameters
    ----------
    len_data : pd.DataFrame
        Should have columns 'Condition', 'Replicate', and '% Count'.

    Returns
    -------
    cum_count_col : list
        A list of cumulative counts. This can then be assigned to a column in
        length data (will have same length as length data).
    """

    cum_count_col = []

    for cond in len_data.Condition.unique():
        for rep in len_data[len_data.Condition == cond].Replicate.unique():

            sub_len_data = len_data[len_data.Condition == cond]
            sub_len_data = sub_len_data[len_data.Replicate == rep]
            pct_counts = sub_len_data['% Count']

            cum_counts = []
            cum_val = 0
            for pct in pct_counts:
                cum_val += pct
                cum_counts.append(cum_val)

            cum_count_col += cum_counts

    return cum_count_col


def closure(mat):
    """
    NOTE:
    Copied from skbio.stats.composition.closure since I had issues installing
    the skbio package on Windows.

    Performs closure to ensure that all elements add up to 1.
    Parameters
    ----------
    mat : array_like
       a matrix of proportions where
       rows = compositions
       columns = components
    Returns
    -------
    array_like, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1
    Raises
    ------
    ValueError
       Raises an error if any values are negative.
    ValueError
       Raises an error if the matrix has more than 2 dimension.
    ValueError
       Raises an error if there is a row that has all zeros.
    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import closure
    >>> X = np.array([[2, 2, 6], [4, 4, 2]])
    >>> closure(X)
    array([[ 0.2,  0.2,  0.6],
           [ 0.4,  0.4,  0.2]])
    """
    mat = np.atleast_2d(mat)
    if np.any(mat < 0):
        raise ValueError("Cannot have negative proportions")
    if mat.ndim > 2:
        raise ValueError("Input matrix can only have two dimensions or less")
    if np.all(mat == 0, axis=1).sum() > 0:
        raise ValueError("Input matrix cannot have rows with all zeros")
    mat = mat / mat.sum(axis=1, keepdims=True)
    return mat.squeeze()


def clr(mat):
    r"""
    NOTE:
    Copied from skbio.stats.composition.clr since I had issues installing the
    skbio package on Windows.

    Performs centre log ratio transformation.
    This function transforms compositions from Aitchison geometry to
    the real space. The :math:`clr` transform is both an isometry and an
    isomorphism defined on the following spaces
    :math:`clr: S^D \rightarrow U`
    where :math:`U=
    \{x :\sum\limits_{i=1}^D x = 0 \; \forall x \in \mathbb{R}^D\}`
    It is defined for a composition :math:`x` as follows:
    .. math::
        clr(x) = \ln\left[\frac{x_1}{g_m(x)}, \ldots, \frac{x_D}{g_m(x)}\right]
    where :math:`g_m(x) = (\prod\limits_{i=1}^{D} x_i)^{1/D}` is the geometric
    mean of :math:`x`.
    Parameters
    ----------
    mat : array_like, float
       a matrix of proportions where
       rows = compositions and
       columns = components
    Returns
    -------
    numpy.ndarray
         clr transformed matrix
    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import clr
    >>> x = np.array([.1, .3, .4, .2])
    >>> clr(x)
    array([-0.79451346,  0.30409883,  0.5917809 , -0.10136628])
    """
    mat = closure(mat)
    lmat = np.log(mat)
    gm = lmat.mean(axis=-1, keepdims=True)
    return (lmat - gm).squeeze()
