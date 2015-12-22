# Add depth entries to a VCF

import gzip
import sys

if len(sys.argv) != 3:
	sys.exit('Usage: python add-depth.py <vcf.gz> <depth_prefix>')


infile_vcf = gzip.GzipFile(sys.argv[1], 'rt')

# Two stages: 
#  1) Modify VCF header to add new fields
#  2) Add depth to data rows
#
# Achieve with an FSM.
# States:
#   0   In preamble, before FORMAT section
#	1	In preamble, during FORMAT section
#	2	In preamble, after FORMAT section
#	3	In header (line starts with #)
#	4	In data
#
# Transitions:
#   0 -> 1	When ##FORMAT is encountered
#   1 -> 2	When ##FORMAT is not encountered, but the line still starts with ##
#   1 -> 3	When the line starts with #, but not ##
#	2 -> 3	When the line starts with #, but not ##
#   3 -> 4	Automatically
#
# Actions on lines:
#   0: Emit line
#   1: Emit line
#   2: Emit line
#   3: Check only one sample.  Emit line.
#   4: Add depth data and emit line.
#
# Actions on transitions:
#   1 -> 2: Emit depth specific format fields
#   1 -> 3: Emit depth specific format fields

state = 0
this_chrom = ''

for line in infile_vcf:
	if state == 0:
		if line.startswith('##FORMAT=<'):
			state = 1
	
	if state == 1:
		if not line.startswith('##FORMAT=<'):
			# Emit the new format entries
			sys.stdout.write('##FORMAT=<ID=KCCG_PERF_DP_MIN,Number=1,Type=Integer,Description="Minimum read depth">\n')
			sys.stdout.write('##FORMAT=<ID=KCCG_PERF_DP_MOD,Number=1,Type=Integer,Description="Modal read depth">\n')
			sys.stdout.write('##FORMAT=<ID=KCCG_PERF_DP_MAX,Number=1,Type=Integer,Description="Maximum read depth">\n')
			if line.startswith('##'):
				state = 2
			else:
				state = 3

	if state == 2:
		if not line.startswith('##'):
			state = 3

	if state != 4:
		sys.stdout.write(line)
	else:
		lineparts = line.rstrip().split('\t')

		# Is this line's chromosome different from the last?  If so,
		# load the new chromosome's depths
		if lineparts[0] != this_chrom:
			this_chrom = lineparts[0]
			this_chrom_depths = {}
			depth_file = open(sys.argv[2] + this_chrom + '.depth', 'rt')
			for depth_line in depth_file:
				depth_line_entries = [int(x) for x in depth_line.rstrip().split('\t')[1:]]
				this_chrom_depths[(depth_line_entries[0], depth_line_entries[1])] = tuple(depth_line_entries[2:])

		# Fetch the depth for this line.
		line_depths = this_chrom_depths[(int(lineparts[1]) - 1, int(lineparts[1]) + len(lineparts[3]) - 1)]

		# Modify line in place and emit
		lineparts[8] += ':KCCG_PERF_DP_MIN:KCCG_PERF_DP_MOD:KCCG_PERF_DP_MAX'
		lineparts[9] += ':{}:{}:{}'.format(*line_depths)
		sys.stdout.write('\t'.join(lineparts) + '\n')

	if state == 3:
		# Check that the header has only one sample
		assert len(line.rstrip().split('\t')) == 10, 'Single-sample VCFs only supported.'
		state = 4
