BEGIN {
	FS=""
	CURRENT_CHROM=""
	CURRENT_LENGTH=-1
}

{
	if (substr($0, 1, 1) == ">")
	{
		if (CURRENT_CHROM != "")
			print CURRENT_CHROM, 0, CURRENT_LENGTH
		gsub(/ .*/, "", $0)
		CURRENT_CHROM = substr($0, 2)
		CURRENT_LENGTH = 0		
	}
	else
	{
		gsub(/[ \t]+/, "" $0)
		CURRENT_LENGTH += length($0)
	}
}

END {
	if (CURRENT_CHROM != "")
		print CURRENT_CHROM, 0, CURRENT_LENGTH
}

