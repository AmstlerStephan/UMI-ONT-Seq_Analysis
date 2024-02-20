import argparse
import logging
import os

import pysam
import pandas as pd
import sys

def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Command line interface to telemap"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-l",
        "--log",
        dest="log",
        choices=[
            "DEBUG",
            "INFO",
            "WARNING",
            "ERROR",
            "CRITICAL",
            "debug",
            "info",
            "warning",
            "error",
            "critical",
        ],
        default="INFO",
        help="Print debug information",
    )
    parser.add_argument(
        "--bam_file", 
        dest="BAM_FILE", 
        type=str, 
        required=True, 
        help="input bam file to extract haplotypes from"
    )
    parser.add_argument(
        "--cluster_stats", 
        dest="CLUSTER_STATS", 
        type=str, 
        help="cluster stats file"
    )
    parser.add_argument(
        "-o", 
        "--output", 
        dest="OUTPUT", 
        required=True, 
        help="Output folder"
    )
    parser.add_argument(
        "--min_reads_per_cluster", 
        dest="MIN_CLUSTER_SIZE", 
        type=int,
        help="Minimal cluster size"
    )
    parser.add_argument(
        "--max_reads_per_cluster", 
        dest="MAX_CLUSTER_SIZE", 
        type=int,
        help="Maximal cluster size"
    )
    
    args = parser.parse_args(argv)

    return args

def filter_bam(args):
    bam_file = args.BAM_FILE
    cluster_stats = args.CLUSTER_STATS
    min_cluster_size = args.MIN_CLUSTER_SIZE
    max_cluster_size = args.MAX_CLUSTER_SIZE
    output = args.OUTPUT
    outfile = os.path.join(output, "extracted_STR.tsv")  
    STR_start = 2472
    STR_end = 2506       
    STR_stats = list()
    clusters_above_threshold = get_clusters(cluster_stats, min_cluster_size, max_cluster_size)    

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            name = read.query_name[:-2]
            if name in clusters_above_threshold:
                STR_length = read.get_overlap(STR_start, STR_end)
                STR_start_pos, STR_end_pos = get_STR_pos(read, STR_start, STR_end)
                STR_sequence = read.get_forward_sequence()[STR_start_pos:STR_end_pos]
                STR_stats.append((name, STR_length, STR_sequence))
        
        write_STR_stats(outfile, STR_stats)

def write_STR_stats(outfile, STR_stats):
    with open(outfile, "w") as out:
        print("read_name\tSTR_length\tSTR_sequence", file=out)
        for read_name, STR_length, STR_sequence in STR_stats:
            print("{}\t{}\t{}".format(read_name, STR_length, STR_sequence), file=out)
        
def get_STR_pos(read, STR_start, STR_end):
    pairs = read.get_aligned_pairs(with_seq=True)
    STR_start_pos = None
    STR_end_pos = None
    i = 0
    # pairs_subset = pairs[STR_start - 200 : STR_end + 300]
    pairs_subset = pairs
    
    while STR_start_pos == None:
        for query_pos, ref_pos, base in pairs_subset:
            if ref_pos == STR_start - i: 
                STR_start_pos = query_pos
            if i > 200:
                break
        i = i + 1
    
    i = 0
    while STR_end_pos == None:
        for query_pos, ref_pos, base in pairs_subset:
            if ref_pos == STR_end + i:
                STR_end_pos = query_pos
            if i > 200:
                break
        i = i + 1

    return STR_start_pos, STR_end_pos                      

def get_clusters(cluster_stats_file, min_cluster_size, max_cluster_size):
    cluster_stats = pd.read_csv(cluster_stats_file, sep = "\t")
    
    cluster_stats_filtered = cluster_stats[(cluster_stats["cluster_written"] == 1) & 
                                           (cluster_stats["reads_written_fwd"] + cluster_stats["reads_written_rev"] >= min_cluster_size) &
                                           (cluster_stats["reads_found"] <= max_cluster_size)]
    
    return cluster_stats_filtered["cluster_id"].to_list()    

def main(argv=sys.argv[1:]):
    """
    Basic command line interface to telemap.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log.upper())
    logging.basicConfig(level=numeric_level, format="%(message)s")

    filter_bam(args)


if __name__ == "__main__":
    main()