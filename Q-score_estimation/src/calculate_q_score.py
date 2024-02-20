import argparse
import logging
import os
import sys
import json
import math
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
        "--mutation_classification", dest="MUTATION_CLASSIFICATION", type=str, help="Mutation classification of the plasmids"
    )
    parser.add_argument(
        "--barcode", dest="BARCODE", type=str, help="Barcode of the sample"
    )
    parser.add_argument(
        "--sample_sheet", dest="SAMPLE_SHEET", type=str, help="Sample sheet containing barcode sample information"
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

def calculate_q_score(args):
    mutation_classification = args.MUTATION_CLASSIFICATION
    bam_file = args.BAM
    output = args.OUT
    sample_sheet = args.SAMPLE_SHEET
    barcode = args.BARCODE
    
    sample = get_sample(sample_sheet, barcode)
    sample_type = "TypeA" if sample[4:].startswith("1000") else "TypeB"
    fragment = int(sample[-4:])
    
    mutations = get_mutations(mutation_classification, fragment)
    errors = get_errors(bam_file, sample_type, fragment, mutations)
    write_error_stats(errors, fragment, sample_type, output)
    write_qscore_stats(errors, fragment, sample_type, output)

def write_qscore_stats(errors, fragment, sample_type, output):
    outfile = os.path.join(output, "cluster_qscore_stats.tsv")
    with open(outfile, "w") as out:
        print("cluster\tqscore\tn_errors\terror_probability\tfragment\tsample_type", file=out)
        for read, info in errors.items():
            qscore = info["qscore"]
            n_errors = info["n_errors"]
            error_probability = info["error_probability"]
            print("{}\t{}\t{}\t{}\t{}\t{}".format(read, qscore, n_errors, error_probability, fragment, sample_type), file = out)    
    
def write_error_stats(errors, fragment, sample_type, output):
    outfile = os.path.join(output, "cluster_error_stats.tsv")
    with open(outfile, "w") as out:
        print("cluster\tqscore\tn_errors\terror_probability\tposition\tread_position\treference_base\tread_base\terror_type\tfwd_context\trev_context\tfragment\tsample_type", file=out)
        for read, info in errors.items():
            qscore = info["qscore"]
            n_errors = info["n_errors"]
            error_probability = info["error_probability"]
            if len(info["errors"]) > 0:
                for error in info["errors"]:
                    position = error["position"]
                    read_position = error["read_position"]
                    reference_base = error["reference_base"]
                    read_base = error["read_base"]
                    error_type = error["error_type"]
                    fwd_context = error["fwd_context"]
                    rev_context = error["rev_context"]
                    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                        read, 
                        qscore, 
                        n_errors, 
                        error_probability, 
                        position, 
                        read_position, 
                        reference_base, 
                        read_base, 
                        error_type,
                        fwd_context,
                        rev_context,
                        fragment, 
                        sample_type
                        ), file = out)    
            else:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    read, 
                    qscore, 
                    n_errors, 
                    error_probability, 
                    "-", 
                    "-", 
                    "-", 
                    "-", 
                    "-", 
                    "-", 
                    "-", 
                    fragment, 
                    sample_type
                    ), file = out)    
                
def get_q_measures(errors, fragment):
    n_errors = len(errors)
    error_probability = n_errors / fragment
    if n_errors == 0:
        q_score = 60
    else:
        q_score = -10 * math.log10(error_probability)
    return q_score, n_errors, error_probability
    

def get_errors(bam_file, sample_type, fragment, mutations):
    errors = dict(dict(list(dict())))

    with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam.fetch(until_eof = True):
                errors[read.query_name] = dict()
                errors[read.query_name]["errors"] = list(dict())
                read_sequence = read.query_sequence
                pairs = get_read_pairs(read)
                last_ref_pos = 0
                last_read_pos = 0
                for read_pos, ref_pos, ref_base in pairs:
                    is_insertion = ref_pos is None
                    is_del = read_pos is None
                    
                    if not (is_del or is_insertion):
                        last_ref_pos = ref_pos
                        last_read_pos = read_pos
                        read_base = read_sequence[read_pos]
                        if read_base.upper() == ref_base.upper():
                            continue
                        else:
                            error_type = "mismatch"
                    else:
                        if is_del:
                            last_ref_pos = ref_pos
                            read_base = "D" 
                            error_type = "deletion"
                            read_pos = round(last_read_pos + 0.1, 2)
                        elif is_insertion:
                            read_base = read_sequence[read_pos]
                            error_type = "insertion"
                            ref_base = ":"
                            ref_pos = round(last_ref_pos + 0.1, 2)
                                                
                    if is_error(ref_pos, read_base, mutations, sample_type):
                        fwd_context_length = 5    
                        rev_context_length = 5
                        last_read_pos_rounded = round(last_read_pos)
                        residual_bases = len(read_sequence) - last_read_pos_rounded
                        if residual_bases <= 5:
                            fwd_context_length = residual_bases
                        if last_read_pos <= 5:
                            rev_context = last_read_pos_rounded
                        fwd_context = read_sequence[last_read_pos_rounded:last_read_pos_rounded + fwd_context_length]
                        rev_context = read_sequence[last_read_pos_rounded - rev_context_length:last_read_pos_rounded]
                        errors[read.query_name]["errors"].append(
                            {"position" : ref_pos + 1,
                            "read_position" : read_pos + 1,
                            "reference_base" : ref_base, 
                            "read_base" : read_base,
                            "error_type" : error_type, 
                            "fwd_context" : fwd_context, 
                            "rev_context" : rev_context})
                    last_ref_pos = ref_pos
                    last_read_pos = read_pos
                
                q_score, n_errors, error_probability = get_q_measures(errors[read.query_name]["errors"], fragment)
                errors[read.query_name]["qscore"] = q_score
                errors[read.query_name]["n_errors"] = n_errors
                errors[read.query_name]["error_probability"] = error_probability    
    return errors 

def get_read_pairs(read):
    query_start = read.query_alignment_start
    query_end = read.query_alignment_end
    pairs = read.get_aligned_pairs(with_seq=True)
    return pairs[query_start:query_end]

def is_error(ref_pos, read_base, mutations, sample_type):
    # python is 0-based
    ref_pos = ref_pos + 1
    is_error = True
    read_base = read_base.upper()
    
    if ref_pos in mutations.index:
        mutation_type = mutations["type_annot"].get(ref_pos)
        mutation = mutations[sample_type].get(ref_pos)
        if (mutation_type == "STR") or (mutation == read_base):
            is_error = False
    
    return is_error

def get_mutations(mutation_classification, fragment):
    mutations = pd.read_csv(mutation_classification)
    mutations_filtered = mutations[mutations["fragment"] == fragment]
    mutations_filtered.set_index(["Position"], inplace = True)
    return mutations_filtered
    
def get_sample(sample_sheet_path, barcode):
    nanopore_barcode = "NB{}".format(barcode[-2:])
    sample_sheet = json.load(open(sample_sheet_path, "r"))
    for sample in sample_sheet:
        if sample["Barcode"] == nanopore_barcode:
            return sample["Sample"]


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

    calculate_q_score(args)


if __name__ == "__main__":
    main()
