import argparse
import logging
import os
import sys

import pysam
import pandas as pd

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
        "-t", "--threads", dest="THREADS", type=int, default=1, help="Number of threads."
    )
    parser.add_argument(
        "--min_cluster_size", dest="MIN_CLUSTER_SIZE", type=int, default = 2, help="Minimal number of reads per cluster"
    )
    parser.add_argument(
        "--max_cluster_size", dest="MAX_CLUSTER_SIZE", type=int, default = 100, help="Maximal number of reads per cluster"
    )
    parser.add_argument(
        "--mode", dest="MODE", type=str, default="range", help="Mode of cluster filtering [single_size | range]"
    )
    parser.add_argument(
        "-o", "--output", dest="OUT", type=str, required=False, help="Output folder"
    )
    parser.add_argument(
        "--bam_file", dest="BAM", type=str, help="BAM file"
    )
    parser.add_argument(
        "--cluster_stats", dest="CLUSTER_STATS", type=str, help="cluster stats file"
    )
    args = parser.parse_args(argv)

    return args

def filter_reads(args):
    cluster_stats = args.CLUSTER_STATS
    bam_file = args.BAM
    output = args.OUT
    min_cluster_size = args.MIN_CLUSTER_SIZE 
    max_cluster_size = args.MAX_CLUSTER_SIZE 
    mode = args.MODE
    
    for cluster_threshold in range(min_cluster_size, max_cluster_size, 1):
        clusters_above_threshold = get_clusters(cluster_stats, cluster_threshold, mode)    
        outfile_name = get_outfile_name(output, cluster_threshold, max_cluster_size, mode)
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            outfile = pysam.AlignmentFile(outfile_name, "w", template = bam)
            for read in bam.fetch(until_eof=True):
                name = read.query_name[:-2]
                if name in clusters_above_threshold:
                    outfile.write(read)              

def get_outfile_name(output, cluster_threshold, max_cluster_size, mode):
    if mode == "single_size":
        outfile_name = os.path.join(output, "clustersize_{}.bam".format(cluster_threshold) )          
    elif mode == "range":
        outfile_name = os.path.join(output, "clustersize_{}_to_{}.bam".format(cluster_threshold, max_cluster_size) )         
    return outfile_name
    
def get_clusters(cluster_stats_file, cluster_threshold, mode):
    cluster_stats = pd.read_csv(cluster_stats_file, sep = "\t")
    if mode == "single_size":
        cluster_stats_filtered = cluster_stats[(cluster_stats["cluster_written"] == 1) & ( cluster_stats["reads_found"] == cluster_threshold)]
    elif mode == "range":
        cluster_stats_filtered = cluster_stats[(cluster_stats["cluster_written"] == 1) & ( cluster_stats["reads_found"] >= cluster_threshold)]
    else:
        raise TypeError("Mode not available. Choose 'single_size' or 'range'")
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

    filter_reads(args)


if __name__ == "__main__":
    main()
