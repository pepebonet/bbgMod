import os
import sys
import queue
import click
import functools
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from time import sleep
from scipy import stats
from multiprocessing import Process, Queue, Pipe

from tombo import tombo_helper as th
from tombo.tombo_helper import TomboReads


def window(a, w=5, o=1):
    sh = (a.size - w + 1, w); st = a.strides * 2
    return np.lib.stride_tricks.as_strided(
        a, strides = st, shape = sh)[0::o].copy()


def calc_window_fishers_method(pvals, lag, method='fisher'):
    """Compute Fisher's Method over a moving window across a set of p-values
    """
    f_pvals = np.empty(pvals.shape); f_pvals[:] = np.NAN
    result = window(pvals)
    f_pvals[...,lag:-lag] = [stats.combine_pvalues(i, 'stouffer')[1] for i in result]
    return f_pvals


def compute_ks_tests(samp_base_levels, ctrl_base_levels):
    samp_valid_indices = np.logical_not(np.isnan(samp_base_levels))
    ctrl_valid_indices = np.logical_not(np.isnan(ctrl_base_levels))

    stats_ks =  np.array([stats.ks_2samp(
        np.sort(pos_samp_levels[samp_valid_indices[i]]),
        np.sort(pos_ctrl_levels[ctrl_valid_indices[i]]))
                     for i, (pos_samp_levels, pos_ctrl_levels) in enumerate(zip(
                             samp_base_levels, ctrl_base_levels))])
    return stats_ks[:, 1]


def get_stats(samp_base_levels, ctrl_base_levels, cov_start, cov_end, fm_offset):
    cov_reg_stats = compute_ks_tests(
        samp_base_levels[cov_start:cov_end],
        ctrl_base_levels[cov_start:cov_end])

    cov_reg_stats = calc_window_fishers_method(
        cov_reg_stats, fm_offset)

    indeces = np.argwhere(
        np.logical_not(np.isnan(cov_reg_stats)) == True).flatten()
    return cov_reg_stats, indeces


def get_reg_output(reg_data, cov_reg_stats, samp_cov, ctrl_cov, 
    reg_stats, reg_poss, reg_cov, reg_ctrl_cov, sequences, cov_start, 
    cov_end, fm_offset, indeces):

    reg_stats.append(cov_reg_stats[indeces])
    reg_poss.append(np.arange(reg_data.start - fm_offset + cov_start,
        reg_data.start - fm_offset + cov_end)[indeces])
    reg_cov.append(samp_cov[cov_start:cov_end][indeces])
    reg_ctrl_cov.append(ctrl_cov[cov_start:cov_end][indeces])
    sequences.append([reg_data.seq[i] for i in indeces - fm_offset + cov_start])

    return reg_stats, reg_poss, reg_cov, reg_ctrl_cov, sequences


def get_raw_currents(reg_data, ctrl_reg_data, fm_offset):
    samp_base_levels = reg_data.copy().update(
        start=reg_data.start - fm_offset,
        end=reg_data.end + fm_offset).get_base_levels()

    ctrl_base_levels = ctrl_reg_data.copy().update(
        start=ctrl_reg_data.start - fm_offset,
        end=ctrl_reg_data.end + fm_offset).get_base_levels()

    return samp_base_levels, ctrl_base_levels


def get_reg_cov(samp_base_levels, ctrl_base_levels, min_test_reads):

    samp_cov = np.logical_not(np.isnan(samp_base_levels)).sum(axis=1)
    ctrl_cov = np.logical_not(np.isnan(ctrl_base_levels)).sum(axis=1)
    cov_regs = np.where(np.diff(np.concatenate([[False,], np.logical_and(
        np.greater_equal(samp_cov, min_test_reads),
        np.greater_equal(ctrl_cov, min_test_reads)), [False,]])))[0]

    return samp_cov, ctrl_cov, cov_regs


def output_as_df(reg_data, reg_stats, reg_poss, reg_cov, reg_ctrl_cov, sequences):
    chrom = np.repeat(reg_data.chrm, len(np.concatenate(reg_poss)))
    strand = np.repeat(reg_data.strand, len(np.concatenate(reg_poss)))
    df = pd.DataFrame({
        'chrm': chrom, 
        'pos': np.concatenate(reg_poss),
        'strand': strand,
        'seq': np.concatenate(sequences),
        'cov_treat': np.concatenate(reg_cov), 
        'cov_untreat': np.concatenate(reg_ctrl_cov),
        'pval': np.concatenate(reg_stats)
    })
    return df


def compute_group_reg_stats(
        reg_data_naive, reads_index, ctrl_reads_index, fm_offset, min_test_reads):
    
    ctrl_reg_data = reg_data_naive.copy().add_reads(ctrl_reads_index).add_seq()
    reg_data = reg_data_naive.add_reads(reads_index).add_seq()

    samp_base_levels, ctrl_base_levels = get_raw_currents(
        reg_data, ctrl_reg_data, fm_offset
    )
    samp_cov, ctrl_cov, cov_regs = get_reg_cov(
        samp_base_levels, ctrl_base_levels, min_test_reads
    )
    if len(cov_regs) == 0:
        return pd.DataFrame()

    reg_stats, reg_poss, reg_cov, reg_ctrl_cov, sequences = [], [], [], [], []
    for cov_start, cov_end in zip(cov_regs[:-1:2], cov_regs[1::2]):

        if cov_end - cov_start < (fm_offset * 2) + 1: continue
        cov_reg_stats, indeces = get_stats(
            samp_base_levels, ctrl_base_levels, cov_start, cov_end, fm_offset
        )
        reg_stats, reg_poss, reg_cov, reg_ctrl_cov, sequences = get_reg_output(
            reg_data, cov_reg_stats, samp_cov, ctrl_cov, reg_stats, reg_poss, 
            reg_cov, reg_ctrl_cov, sequences, cov_start, cov_end, fm_offset, indeces
        )

    if len(reg_stats) == 0:
        return pd.DataFrame()

    return output_as_df(reg_data, reg_stats, reg_poss, 
            reg_cov, reg_ctrl_cov, sequences)


def _test_signif_worker(
        region_q, stats_q, progress_q, reads_index, fm_offset,
        min_test_reads, ctrl_reads_index):
    while True:
        try:
            reg_data = region_q.get(block=False)
        except queue.Empty:
            # sometimes throws false empty error with get(block=False)
            sleep(0.01)
            if not region_q.empty(): continue
            break

        try:
            stat_type_reg_stats = compute_group_reg_stats(
                reg_data, reads_index, ctrl_reads_index, fm_offset, min_test_reads
            )
            
            stats_q.put(stat_type_reg_stats)
            progress_q.put(1)

        except:
            progress_q.put(1)
            continue

    return


def test_significance(reads_index, region_size, cpus, output,
        min_test_reads, fm_offset, ctrl_reads_index):
    """Test for significant shifted signal in mutliprocessed batches
    """
    #TODO <JB> add progress by filling progress_q and adding tqdm
    region_q = Queue()
    stats_q = Queue()
    progress_q = Queue()

    # split chromosomes into separate regions to process independently
    num_regions = 0

    for chrm, strand, reg_start in reads_index.iter_cov_regs(
            1, region_size, ctrl_reads_index):
        region_q.put(th.intervalData(
            chrm=chrm, start=reg_start, end=reg_start + region_size,
            strand=strand))
        num_regions += 1
    # wait for queue items to register in queue and avoid queue appearing empty
    sleep(0.1)

    test_ps = []
    for p_id in range(cpus):
        p = Process(target=_test_signif_worker, args=(region_q, stats_q, 
        progress_q, reads_index, fm_offset, min_test_reads, ctrl_reads_index))
        p.start()
        test_ps.append(p)

    main_prog_conn, prog_conn = Pipe()
    prog_p = Process(target=_get_progress_queue,
                        args=(progress_q, prog_conn, num_regions))
    prog_p.daemon = True
    prog_p.start()

    # main region stats queue getter
    write_folder = os.path.join(output, 'tmp')
    main_stats_conn, stats_conn = Pipe()
    stats_p = Process(target=_get_stats_queue, args=(
        stats_q, stats_conn, write_folder))
    stats_p.daemon = True
    stats_p.start()

    # wait for test processes to finish
    for test_p in test_ps:
        test_p.join()

    # in a very unlikely case the progress queue could die while the
    # main process remains active and thus we would have a deadlock here
    if prog_p.is_alive():
        # send signal to getter queue to finish and return results
        main_prog_conn.send(True)
        # returns total number of processed reads if that is needed
        main_prog_conn.recv()

    main_stats_conn.send(True)
    main_stats_conn.recv()

    return write_folder


def _get_progress_queue(progress_q, prog_conn, num_regions):
    th.status_message(
        'Performing modified base detection across genomic regions.')
    bar = tqdm(total=num_regions, smoothing=0)

    tot_num_rec_proc = 0
    while True:
        try:
            iter_val = progress_q.get(block=False)
            tot_num_rec_proc += iter_val
            bar.update(iter_val)
        except queue.Empty:
            if prog_conn.poll():
                break
            sleep(0.1)
            continue

    bar.close()
    prog_conn.send(tot_num_rec_proc)

    return

def _get_stats_queue(stats_q, stats_conn, output):
    count = 0; os.mkdir(output)
    while True:
        try:
            stats = stats_q.get(block=False)
            out_file = os.path.join(output, 'stats_modifications_{}.tsv'.format(count))
            stats.to_csv(out_file, sep='\t', index=None)
            count += 1
        except queue.Empty:
            # wait for main process to send indicator that all regions
            # have been processed
            if stats_conn.poll():
                sleep(0.1)
                break
            sleep(0.1)
            continue

    # Clear leftover values from queues
    while not stats_q.empty():
        stats = stats_q.get(block=False)
        out_file = os.path.join(output, 'stats_modifications_{}.tsv'.format(count))
        stats.to_csv(out_file, sep='\t', index=None)
        count += 1

    stats_conn.send(True)

    return


def arrange_final_output(write_folder):
    header_file = os.path.join(write_folder, 'stats_modifications_0.tsv')
    all_files = os.path.join(write_folder, '*.tsv')
    out_file = os.path.join(write_folder.rsplit('/', 1)[0], 'modifications.tsv')
    cmd = 'head -1 {} > {}; tail -n +2 -q {} >> {}'. format(
        header_file, out_file, all_files, out_file)
    subprocess.call(cmd, shell=True)
    subprocess.call('rm -r {}'.format(write_folder), shell=True)


def main(fast5_basedirs, alternate_fast5_basedirs, cpus, 
    corrected_group, basecall_subgroups, output):
    print('Performing two-sample group comparison significance testing.')
    reads_index = th.TomboReads(
            [fast5_basedirs], corrected_group, basecall_subgroups)
    
    ctrl_reads_index = th.TomboReads(
        [alternate_fast5_basedirs], corrected_group,
        basecall_subgroups)

    write_folder = test_significance(reads_index,
        region_size=1000, min_test_reads=30, fm_offset=2,
        ctrl_reads_index=ctrl_reads_index, cpus=cpus, output=output)

    arrange_final_output(write_folder)
    