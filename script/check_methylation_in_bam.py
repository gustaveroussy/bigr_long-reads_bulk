#! usr/bin/python3
import pysam
import sys

#BAM files must be indexed first !

input_file=sys.argv[1]
samfile=pysam.AlignmentFile(input_file,"rb",check_sq=False)

for read in samfile.fetch(until_eof=True):
	break

try :
	read.get_tag("MM")
except KeyError:
	mm_exist = False
else:
	mm_exist = True

print(mm_exist)
