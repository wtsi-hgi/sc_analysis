#-- import modules --#
import io
import os
import sys
import argparse
import re
import gc
import subprocess
import numpy as np
import polars as pl
from itertools import product
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice
from collections import Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def pigz_open(path: str):
    """
    Open a gzip file using pigz for faster decompression
    Parameters:
        -- path: file path to the gzip file
    Returns:
        -- io_wrapper: a TextIOWrapper for reading the decompressed file
    """
    return subprocess.Popen(["pigz", "-dc", path], stdout = subprocess.PIPE)

def fastq_iter_pe(handle1, handle2):
    """
    FASTQ paired-end parser yielding ((header1, seq1, qual1), (header2, seq2, qual2))
    Parameters:
        -- handle1: file handle for read 1 FASTQ
        -- handle2: file handle for read 2 FASTQ
    Yields:
        -- ((header1, seq1, qual1), (header2, seq2, qual2))
    """
    while True:
        header1 = handle1.readline()
        header2 = handle2.readline()
        if not header1 or not header2:
            break

        seq1 = handle1.readline()
        seq2 = handle2.readline()
        plus1 = handle1.readline()
        plus2 = handle2.readline()
        qual1 = handle1.readline()
        qual2 = handle2.readline()

        # check for truncated records
        if not (seq1 and plus1 and qual1 and seq2 and plus2 and qual2):
            raise ValueError("Truncated FASTQ record in one of the pairs")

        yield ((header1.rstrip("\n"), seq1.rstrip("\n"), qual1.rstrip("\n")), 
               (header2.rstrip("\n"), seq2.rstrip("\n"), qual2.rstrip("\n")))
        
def init_worker():
    """
    Initilize worker
    """
    global dict_tf_barcodes
    global dict_cell_barcodes

def match_hamming_numpy(seq: str, pattern: str, max_mismatches: int) -> int:
    """
    Find approximate match of pattern in seq by hamming distance allowing max_mismatches
    Parameters:
        -- seq: the target sequence
        -- pattern: the pattern to match
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- int: the Hamming distance, or the maximum length if they differ in length
    """
    k = len(pattern)
    n = len(seq)
    if k > n:
        return -1

    # convert the sequence to np.uint8 arrays of ASCII codes
    seq_arr = np.frombuffer(seq.encode("ascii"), dtype = np.uint8)
    pat_arr = np.frombuffer(pattern.encode("ascii"), dtype = np.uint8)

    # sliding window: create a 2D view of seq_arr of shape (n-k+1, k)
    windows = np.lib.stride_tricks.sliding_window_view(seq_arr, window_shape = k)

    # calculate hamming distances vectorized
    mismatches = np.sum(windows != pat_arr, axis = 1)
    matches = np.where(mismatches <= max_mismatches)[0]
    return int(matches[0]) if matches.size > 0 else -1

def build_barcode_lookup(barcode_dict, max_mismatch, alphabet = "ACGT"):
    """
    Build a fast lookup table for cell barcodes allowing up to `max_mismatch` mismatches.
    Parameters:
        -- barcode_dict: Dictionary where keys are valid barcode sequences (strings)
        -- max_mismatch: Number of mismatches to tolerate
        -- alphabet: Alphabet used for barcodes (default = "ACGT").
    Returns:
    -------
        -- dict: Lookup dictionary mapping any possible mutated barcode sequence
    """
    lookup = {}
    for cb in barcode_dict:
        lookup[cb] = cb

        if max_mismatch >= 1:
            for i in range(len(cb)):
                original_base = cb[i]
                for base in alphabet:
                    if base == original_base:
                        continue
                    mutated = cb[:i] + base + cb[i+1:]
                    lookup.setdefault(mutated, cb)

    return lookup

def process_pe_pair(read_pair: tuple) -> list:
    """
    Detect canonical splicing event in a pair of reads
    Parameters:
        -- read: tuple (read1, read2) and read1, read2 are tuple (header, sequence, quality)
    Returns:
        -- tuple: tuple: (read1, read2, barcode) 
    """
    read1, read2 = read_pair
    read1_seq = read1[1]
    read2_seq = read2[1]

    stats = { 'n_processed_reads':        1,
              'n_failed_reads':           0,
              'n_cell_barcode_not_found': 0,
              'n_marker_not_found':       0,
              'n_tf_barcode_not_matched': 0,
              'n_valid':                  0 }

    if "N" in read1_seq:
        stats["n_failed_reads"] += 1
        return {}, stats
    
    read1_cb = read1_seq[:16]
    if args.has_umi:
        read1_umi = read1_seq[16:26]
    else:
        read1_umi = "no_umi"

    if read1_cb not in dict_cell_barcodes:
        stats["n_cell_barcode_not_found"] += 1
        return {}, stats
    found_cb = dict_cell_barcodes[read1_cb]

    marker_region = read2_seq[args.marker_start : args.marker_end + 1]
    tf_found = match_hamming_numpy(marker_region, args.marker_seq, args.max_mismatch)
    if tf_found == -1:
        stats["n_marker_not_found"] += 1
        return {}, stats

    tf_barcode_start = args.marker_start + tf_found + len(args.marker_seq)
    tf_barcode = read2_seq[tf_barcode_start : tf_barcode_start + args.tf_barcode_len]
    if tf_barcode not in dict_tf_barcodes:
        stats["n_tf_barcode_not_matched"] += 1
        return {}, stats

    dict_out = { "cell_barcode" : found_cb,
                 "read_umi"     : read1_umi,
                 "tf_barcode"   : tf_barcode,
                 "tf_name"      : dict_tf_names[dict_tf_barcodes[tf_barcode]] }
    stats["n_valid"] += 1
    return dict_out, stats

def batch_process_pe_pairs(batch_reads: list) -> list:
    """
    Process a batch of read pairs to extract variants and barcodes.
    Parameters:
        -- batch_reads
    Returns:
        -- list of tuples
    """
    results = []
    for read_pair in batch_reads:
        result = process_pe_pair(read_pair)
        results.append(result)
    return results

def function_processpool_pe(args):
    """
    Wrapper function for process pool as ProcessPoolExecutor expects a function rather than returned results.
    """
    return batch_process_pe_pairs(args)

def process_pe_pairs_in_chunk(path_read1, path_read2):
    """
    Read paired-end FASTQ files in chunks and process reads in parallel
    Parameters:
        -- path_read1: path to paired-end read1 FASTQ file
        -- path_read2: path to paired-end read2 FASTQ file
    Yields:
        -- DataFrame: barcode
    """
    fh_read1 = io.TextIOWrapper(pigz_open(path_read1).stdout) if path_read1.endswith(".gz") else open(path_read1)
    fh_read2 = io.TextIOWrapper(pigz_open(path_read2).stdout) if path_read2.endswith(".gz") else open(path_read2)
    read_iter = fastq_iter_pe(fh_read1, fh_read2)

    with ProcessPoolExecutor(max_workers = args.threads, initializer = init_worker) as executor:
        while True:
            read_chunk = list(islice(read_iter, args.chunk_size))
            if not read_chunk:
                break

            # Divide chunk into batches
            # if process_long_read is very fast, we can use a larger batch size to make better use of CPU resources
            # if process_long_read is very slow, we can use a smaller batch size to make better use of CPU resources
            batch_size = min(args.chunk_size, 40000)
            read_batches = [
                read_chunk[i:i+batch_size]
                for i in range(0, len(read_chunk), batch_size)
            ]

            barcode_counts = Counter()
            chunk_stats = { 'n_processed_reads':        0,
                            'n_failed_reads':           0,
                            'n_cell_barcode_not_found': 0,
                            'n_marker_not_found':       0,
                            'n_tf_barcode_not_matched': 0,
                            'n_valid':                  0 }
            # executor.map() --> memeory may grow accumulatively
            # 1. creates all tasks upfront for the entire input list (read_batches)
            # 2. stores them all in an internal queue inside the Executor.
            # 3. waits for results in order — not as they finish.
            # 
            # as_completed()
            # 1. Yields each Future as soon as it finishes, regardless of submission order
            # 2. You can process and discard the result immediately.
            futures = [ executor.submit(function_processpool_pe, batch) for batch in read_batches ]
            for future in as_completed(futures):
                batch_result = future.result()
                for dict_out, stats in batch_result:
                    if dict_out:
                        key = (dict_out["cell_barcode"], 
                               dict_out["read_umi"],
                               dict_out["tf_barcode"],
                               dict_out["tf_name"])
                        barcode_counts[key] += 1
                    for k in chunk_stats:
                        chunk_stats[k] += stats[k]
                # -- free memory -- #
                del batch_result
                gc.collect()

            if barcode_counts:
                df_yield = pl.DataFrame({"cell_barcode" : [k[0] for k in barcode_counts.keys()],
                                         "read_umi"     : [k[1] for k in barcode_counts.keys()],
                                         "tf_barcode"   : [k[2] for k in barcode_counts.keys()],
                                         "tf_name"      : [k[3] for k in barcode_counts.keys()],
                                         "count"        : [v for v in barcode_counts.values()]})
            else:
                df_yield = pl.DataFrame(schema={"cell_barcode" : pl.Utf8, 
                                                "read_umi"     : pl.Utf8, 
                                                "tf_barcode"   : pl.Utf8, 
                                                "tf_name"      : pl.Utf8, 
                                                "count"        : pl.Int64})

            # -- free memory -- #
            del read_chunk, read_batches, futures, barcode_counts
            gc.collect()

            yield df_yield, chunk_stats
    fh_read1.close()
    fh_read2.close()


#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process a bwa bam file for canonical splicing events.", allow_abbrev = False)
    parser.add_argument("--reads",          type = str, required = True,        help = "fastq file(s), eg: se_reads.fq.gz or pe_r1.fq.gz,pe_r2.fq.gz")
    parser.add_argument("--tf_barcode",     type = str, required = True,        help = "csv file of tf barcode, per barcode per line")
    parser.add_argument("--cell_barcode",   type = str, required = True,        help = "csv file of cell barcode, per barcode per line")
    parser.add_argument("--tf_barcode_len", type = int, default = 24,           help = "the length of the tf barcode")
    parser.add_argument("--marker_seq",     type = str, default = "GAAAGGACGA", help = "sequence of the marker, identifies sequence before barcode to determine barcode position")
    parser.add_argument("--marker_start",   type = int, default = 25,           help = "start index of key region where the marker sequence expects to be found")
    parser.add_argument("--marker_end",     type = int, default = 50,           help = "end index of key region where the marker sequence expects to be found")
    parser.add_argument("--max_mismatch",   type = int, default = 1,            help = "max mismatches allowed in barcode matches")
    parser.add_argument("--has_umi",        action = "store_true",              help = "whether R1 reads have umi sequences")
    parser.add_argument("--output_dir",     type = str, default = os.getcwd(),  help = "output directory")
    parser.add_argument("--output_prefix",  type = str, default = '',           help = "output prefix")
    parser.add_argument("--chunk_size",     type = int, default = 100000,       help = "chunk size for processing reads")
    parser.add_argument("--threads",        type = int, default = 40,           help = "number of threads")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    input_path = args.reads.strip()
    input_path_parts = [r.strip() for r in input_path.split(",") if r.strip()]
    if len(input_path_parts) != 2:
        raise ValueError(f"Paired-end reads expects 2 comma-separated files, but got {len(input_path_parts)}: {input_path}")
    path_read1, path_read2 = input_path_parts

    if args.output_prefix == '':
        output_prefix = os.path.splitext(os.path.basename(input_path_parts[0]))[0]
    else:
        output_prefix = args.output_prefix

    # -- read input files -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Reading files, please wait...", flush = True)
    df_tf_barcodes = pl.read_csv(args.tf_barcode, separator = ",", has_header = True)
    dict_tf_names = { row["barcode"]: row["tf_name"] for row in df_tf_barcodes.iter_rows(named = True) }
    dict_tf_barcodes = build_barcode_lookup(dict_tf_names, max_mismatch = args.max_mismatch)

    df_cell_barcodes = pl.read_csv(args.cell_barcode, separator = ",", has_header = False)
    dict_cell_barcodes = { row[0]: True for row in df_cell_barcodes.iter_rows(named = False) }
    dict_cell_barcodes = build_barcode_lookup(dict_cell_barcodes, max_mismatch = args.max_mismatch)

    # -- free memory -- #
    del df_tf_barcodes, df_cell_barcodes
    gc.collect()

    # -- prepare output files -- #
    os.makedirs(args.output_dir, exist_ok = True)
    os.chdir(args.output_dir)

    tf_stats = f"{output_prefix}.tf_stats.tsv"
    if os.path.exists(tf_stats):
        os.remove(tf_stats)

    tf_barcodes = f"{output_prefix}.tf_barcodes.tsv"
    if os.path.exists(tf_barcodes):
        os.remove(tf_barcodes)

    # -- parallel processing -- #
    total_stats = { 'n_processed_reads':        0,
                    'n_failed_reads':           0,
                    'n_cell_barcode_not_found': 0,
                    'n_marker_not_found':       0,
                    'n_tf_barcode_not_matched': 0,
                    'n_valid':                  0 }

    list_results = []
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing fastq file, please wait...", flush = True)
    for i, (chunk_result, chunk_stats) in enumerate(process_pe_pairs_in_chunk(path_read1, path_read2)):
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed chunk {i+1} with {args.chunk_size} reads", flush = True)
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}      No. of valid reads: {chunk_stats['n_valid']}", flush = True)          
        if not chunk_result.is_empty():
            list_results.append(chunk_result)
        for k in chunk_stats:
            total_stats[k] += chunk_stats[k]
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Finished.", flush = True)

     # -- clean and format the extracted barcodes from reads -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Generating results, please wait...", flush = True)
    list_results_filtered = [df for df in list_results if df.height > 0]
    if list_results_filtered:
        df_barcode = pl.concat(list_results_filtered, how = "vertical")
        df_barcode_counts = ( df_barcode.group_by(["cell_barcode", "read_umi", "tf_barcode", "tf_name"])
                                        .agg(pl.sum("count").alias("count"))
                                        .sort(["cell_barcode", "tf_name"]) )
        df_barcode_counts.write_csv(tf_barcodes, separator = "\t", null_value = "NA")
    else:
        with open(tf_barcodes, "w") as f:
            f.write("no barcode found in the reads, please check your barcode marker or positions!\n")

    with open(tf_stats, "w") as fh:
        fh.write(f"Number of processed reads\t{total_stats['n_processed_reads']}\n")
        fh.write(f"Number of failed reads\t{total_stats['n_failed_reads']}\n")
        fh.write(f"Number of reads with cell barcode not found\t{total_stats['n_cell_barcode_not_found']}\n")
        fh.write(f"Number of reads which the marker was not found\t{total_stats['n_marker_not_found']}\n")
        fh.write(f"Number of reads with unmatched barcode\t{total_stats['n_tf_barcode_not_matched']}\n")
        fh.write(f"Number of valid reads\t{total_stats['n_valid']}\n")

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Done.", flush = True)