import sys
import subprocess
import collections

if len(sys.argv) != 7:
	sys.exit('Usage: depth.py <bamfile> <bedfile> <chrom> <minBaseQ> <minMapQ> <outfile>')

bamfile = sys.argv[1]
bedfile = open(sys.argv[2], 'rt')
chrom = sys.argv[3]
min_base_q = sys.argv[4]
min_map_q = sys.argv[5]
outfile = open(sys.argv[6], 'wt')

# samtools depth has the annoying property of omitting
# output for positions lacking any reads.  We would like
# these positions to be included in the depths counter, 
# with depth 0.  The following code reads the output of 
# samtools depth as it is emitted, and adds in any missing 
# zero depth entries.

for line in bedfile:
	lineparts = line.rstrip().split('\t')
	this_chrom = lineparts[0]
	if this_chrom != chrom:
		continue
	start_bed = int(lineparts[1])
	end_bed = int(lineparts[2])

	proc = subprocess.Popen(['samtools', 'depth', '-r', '{}:{}-{}'.format(chrom, start_bed+1, end_bed), '-q', min_base_q, '-Q', min_map_q, bamfile], stdout = subprocess.PIPE)

	depths = collections.Counter()

	# next_pos is the position (in 1-based / samtools 
	# coordinates) of the next base that we expect to see in
	# the samtools depth output
	next_pos = start_bed + 1
	for line in proc.stdout:
		# Get the next line from samtools depth
		line_chrom, line_pos, line_depth = line.strip().split('\t')
		line_pos = int(line_pos)

		# If line_pos > next_pos, then there have been positions
		# skipped by samtools depth.  Add these as zero depth.
		depths[0] += line_pos - next_pos

		# Add the samtools depth output
		depths[line_depth] += 1

		# Update the next expected position
		next_pos = line_pos + 1

	# samtools depth may have finished early, before depth for 
	# the whole requested interval has been emitted.  In this case,
	# add in the missing zero-depth entries at the end.
	depths[0] += end_bed + 1 - next_pos

	depths += collections.Counter()

	# Extract the minimum, maximum, and modal depths in the region
	depth_mode = depths.most_common(1)[0][0]
	depth_min = min((item[0] for item in depths.items()))
	depth_max = max((item[0] for item in depths.items()))

	outfile.write('\t'.join((chrom, str(start_bed), str(end_bed), str(depth_min), str(depth_mode), str(depth_max))) + '\n')
