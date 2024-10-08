#! usr/bin/python3
import pysam
import sys
import argparse

parser = argparse.ArgumentParser()

#Script arguments
parser.add_argument("--input_bam", help="Input BAM file to filter, it should be indexed first.")
parser.add_argument("--output_bam", help="Output filtered BAM file.")
parser.add_argument("--min_length", help="Read minimum length filter.")
parser.add_argument("--min_qual", help="Read minimum base quality score filter.")

args=parser.parse_args()

input_bam = args.input_bam
output_bam = args.output_bam
min_length = args.min_length
min_qual = args.min_qual

#Print given arguments
print("List of given arguments:")
print("input_bam: " + input_bam)
print("output_bam: " + output_bam)
print("min_length: " + min_length)
print("min_qual: " + min_qual)

#Read input and output BAM files
samfile = pysam.AlignmentFile(input_bam, "rb", check_sq = False)
filtered_samfile = pysam.AlignmentFile(output_bam, "wb", template = samfile)

#Metrics variables
nb_reads_before = 0
nb_reads_after = 0

print("\nReading BAM file...")
for read in samfile.fetch(until_eof=True):
    nb_reads_before += 1
    if (read.query_length >= int(min_length) and read.get_tag("qs") >= int(min_qual)):
        #print("Keeping read " + read.query_name + " to filtered BAM...")
        nb_reads_after += 1
        filtered_samfile.write(read)

print("Nb of reads before filtering: " + str(nb_reads_before))
print("Nb of reads after filtering: " + str(nb_reads_after))

samfile.close()
filtered_samfile.close()

print("\nFinished!")