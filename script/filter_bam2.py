#! usr/bin/python3
import pysam
import sys

# bam files needs to be indexed first

input_file=sys.argv[1]
output_file=sys.argv[2]
samfile=pysam.AlignmentFile(input_file,"rb",check_sq=False)

filtered = pysam.AlignmentFile(output_file, "wb", template=samfile)
for read in samfile.fetch(until_eof=True):
	if (read.query_length >= 200 and read.get_tag("qs") >= 10):
		filtered.write(read)

samfile.close()
filtered.close()
