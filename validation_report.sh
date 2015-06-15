#!/bin/bash
set -e

#####################################################################
# VERSION
#####################################################################
VERSION="20150615-1"

#####################################################################
# SOFTWARE AND DATA LOCATIONS
#####################################################################
# Head location for resources bundle
RESOURCES_HEAD="/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources"

# Software
RSCRIPT="/home/marpin/bin/Rscript"
JAVA="/usr/java/latest/bin/java"
RTG_CORE="${RESOURCES_HEAD}/rtg-core/rtg-core.jar"
RTG_THREADS=4
RTG_VCFEVAL="${JAVA} -Xmx8G -jar ${RTG_CORE} vcfeval -T ${RTG_THREADS}"
GIT="/home/marpin/bin/git"

# Data
GOLD_CALLS_VCFGZ="${RESOURCES_HEAD}/gold_standard/calls.vcf.gz"
GOLD_CALLS_VCFGZTBI="${RESOURCES_HEAD}/gold_standard/calls.vcf.gz.tbi"
GOLD_HARDMASK_VALID_REGIONS_BEDGZ="${RESOURCES_HEAD}/gold_standard/valid_regions.bed.gz"
REFERENCE_SDF="${RESOURCES_HEAD}/reference/ref.sdf/"
FUNCTIONAL_REGIONS_BEDGZ_PREFIX="${RESOURCES_HEAD}/functional_regions/"
MASK_REGIONS_BEDGZ_PREFIX="${RESOURCES_HEAD}/mask_regions/"
REFERENCE_BSGENOME="BSgenome.HSapiens.1000g.37d5"		# This is a custom package, available at /share/ClusterShare/biodata/contrib/marpin/reference/hs37d5/build/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz

# Scratch space location (default; this is also set by a command
# line option, to allow rudimentary restarting of runs)
SCRATCH_DEFAULT="/directflow/ClinicalGenomicsPipeline/tmp"

EXEC_DIR=$(pwd)

#####################################################################
# COMMAND LINE PARSING
#####################################################################
print_usage() {
cat << EOF
Usage: ${0##*/} [-d CHROM] [-t TMPDIR] [-r] [-f] [-x] [-o OUTFILE] <INFILE>

Create a WGS validation report.

    INFILE       Input NA12878 genotype calls, in vcf.gz format.
    -o OUTFILE   Write the report to OUTFILE (default: report.pdf)
    -d CHROM     Debug mode.  Currently, adds additional debug 
                 information to the report, and examines chromosome
                 CHROM only.
    -t TMPDIR    Path to a temporary scratch space directory.  If 
                 missing, mktemp will be used to create one.
    -f           Force noninteractive mode; all Y/N prompts will be
                 automatically answered Y.  Default is to run in 
                 interactive mode.
    -r           Resume flag.  If set, and TMPDIR points to a 
                 partially-completed run, will attempt to resume
                 this run.  Default is to not resume; if a partial
                 run is found in TMPDIR and -f is not set, query
                 the user as to what to do.  With -f set, any partial
                 run data in TMPDIR is automatically cleared, and
                 the validation report generation started from the 
                 beginning.  Regardless of the resume setting, the
                 final report generation is always repeated with
                 every invocation.
    -x 			 Generate an extended report, with threshold and
    			 score diagnostics.
    -h           Display this help and exit.

Version ${VERSION}

Mark Pinese
EOF
}

OPTIND=1
input_vcfgz_path=""
debug=0
debug_chrom="-"
noninteractive=0
output_pdf_path="${EXEC_DIR}/report.pdf"
resume=0
temp_supplied=0
extended=0

while getopts "d:o:t:hfrx" opt; do
	case "$opt" in
		h)
			print_usage
			exit 0
			;;
		d)
			debug=1
			debug_chrom=$OPTARG
			;;
		o)
			output_pdf_path=$(readlink -f $OPTARG)
			;;
		f)
			noninteractive=1
			;;
		t)
			temp_supplied=1
			SCRATCH=$OPTARG
			;;
		r)
			resume=1
			;;
		x)
			extended=1
			;;
		'?')
			print_usage >&2
			exit 1
			;;
	esac
done

shift $((OPTIND-1))

if [ $# -ne 1 ]; then
	print_usage >&2
	exit 1
fi

input_vcfgz_path=$1


#####################################################################
# PARAMETER CHECKING
#####################################################################
if [ ! -e ${input_vcfgz_path} ]; then
	echo -e "\033[0;31mError: Input file ${input_vcfgz_path} not found.\033[0m"
	exit 2
fi


if [ -e ${output_pdf_path} ]; then
	if [ ${resume} -eq 0 ]; then
		echo -e "\033[1;33mOutput file ${output_pdf_path} already exists, and the resume flag is not set; if you continue this will be overwritten.\033[0m"
		echo -en "\033[1;33mDo you wish to continue? [YES to continue; any other string to cancel] \033[0m"
		prompt="YES"
		if [ ${noninteractive} -eq 0 ]; then
			read prompt
		fi
		if [ "${prompt}" == "YES" ]; then
			echo -e "\033[0;32mRemoving output file and continuing...\033[0m"
			rm -f ${output_pdf_path}
		else
			echo -e "\033[0;31mCancelled.\033[0m"
			exit 3
		fi
	else
		echo -e "\033[1;33mOutput file ${output_pdf_path} already exists, and resume flag is set.  Attempting to resume run...\033[0m"
	fi
fi

if [ ${temp_supplied} -eq 0 ]; then
	SCRATCH=$(mktemp -d --tmpdir=${SCRATCH_DEFAULT} valrept.XXXXXXXXXX)
fi

RTG_OVERLAP_SCRATCH="${SCRATCH}/overlap"
KNITR_SCRATCH="${SCRATCH}/knitr"


#####################################################################
# SOFTWARE VERSIONING
#####################################################################
VERSION_GIT_BRANCH=$(${GIT} rev-parse --abbrev-ref HEAD)
VERSION_GIT_COMMIT=$(${GIT} rev-parse --verify HEAD)
VERSION_EXEC_HOST=$(uname -a)


#####################################################################
# DEBUG REPORTING
#####################################################################
if [ $debug -eq 1 ]; then
	echo -e "\033[1;36mDebug enabled.\033[0m"
	echo -e "\033[1;36mVariables:\033[0m"
	echo -e "\033[1;36m  RSCRIPT=                          \"${RSCRIPT}\"\033[0m"
	echo -e "\033[1;36m  JAVA=                             \"${JAVA}\"\033[0m"
	echo 
	echo -e "\033[1;36m  RTG_CORE=                         \"${RTG_CORE}\"\033[0m"
	echo -e "\033[1;36m  RTG_THREADS=                      \"${RTG_THREADS}\"\033[0m"
	echo -e "\033[1;36m  RTG_VCFEVAL=                      \"${RTG_VCFEVAL}\"\033[0m"
	echo 
	echo -e "\033[1;36m  RESOURCES_HEAD=                   \"${RESOURCES_HEAD}\"\033[0m"
	echo -e "\033[1;36m  GOLD_CALLS_VCFGZ=                 \"${GOLD_CALLS_VCFGZ}\"\033[0m"
	echo -e "\033[1;36m  GOLD_CALLS_VCFGZTBI=              \"${GOLD_CALLS_VCFGZTBI}\"\033[0m"
	echo -e "\033[1;36m  GOLD_HARDMASK_VALID_REGIONS_BEDGZ=\"${GOLD_HARDMASK_VALID_REGIONS_BEDGZ}\"\033[0m"
	echo -e "\033[1;36m  KCCG_HARDMASK_CALLABLE_REGIONS=   \"${KCCG_HARDMASK_CALLABLE_REGIONS}\"\033[0m"
	echo -e "\033[1;36m  REFERENCE_SDF=                    \"${REFERENCE_SDF}\"\033[0m"
	echo 
	echo -e "\033[1;36m  REFERENCE_BSGENOME=               \"${REFERENCE_BSGENOME}\"\033[0m"
	echo 
	echo -e "\033[1;36m  SCRATCH=                          \"${SCRATCH}\"\033[0m"
	echo -e "\033[1;36m  RTG_OVERLAP_SCRATCH=              \"${RTG_OVERLAP_SCRATCH}\"\033[0m"
	echo -e "\033[1;36m  KNITR_SCRATCH=                    \"${KNITR_SCRATCH}\"\033[0m"
	echo
	echo -e "\033[1;36m  VERSION=                          \"${VERSION}\"\033[0m"
fi


#####################################################################
# CREATE TEMPORARY DIRECTORIES
#####################################################################
mkdir -p ${SCRATCH}
mkdir -p ${KNITR_SCRATCH}



#####################################################################
# VCF OVERLAP EVALUATION
#####################################################################
echo -e "\033[0;32mComputing VCF overlaps...\033[0m"
if [ -e ${RTG_OVERLAP_SCRATCH} ]; then
	if [ ${resume} -eq 0 ]; then
		echo -e "\033[1;33mOverlap scratch directory ${RTG_OVERLAP_SCRATCH} already exists, and the resume flag is not set; if you continue this will be overwritten.\033[0m"
		echo -en "\033[1;33mDo you wish to continue? [YES to continue; any other string to cancel] \033[0m"
		prompt="YES"
		if [ ${noninteractive} -eq 0 ]; then
			read prompt
		fi
		if [ "${prompt}" == "YES" ]; then
			echo -e "\033[0;32mClearing scratch directory and continuing...\033[0m"
			rm -rf ${RTG_OVERLAP_SCRATCH}
		else
			echo -e "\033[0;31mCancelled.\033[0m"
			exit 3
		fi
	fi
fi
# Note on the use of --all-records below: if this is absent, then
# the performance of a given cutoff (for example, GQ > thresh) is
# actually measuring the more complex condition of GQ > thresh AND
# FILTER = PASS.  With --all-records, the performance measures just
# that of the cutoff (eg GQ > thresh) alone.  Considerg including 
# a script option to control this behaviour, because I can think 
# of valid use cases for both configurations.  Alternatively, leave
# --all-records set, and perform some magic in R to move things
# about appropriately.  Specifically, the absence of --all-records
# can be simulated by moving sites in R, as follows:
#   TP for which FILT != PASS  -->  FN
#   FP for which FILT != PASS  -->  TN
#   FN for which FILT != PASS  -->  FN (no change)
if [ ${resume} -eq 0 ] || [ ! -e ${RTG_OVERLAP_SCRATCH}/tp.vcf.gz ] || [ ! -e ${RTG_OVERLAP_SCRATCH}/fp.vcf.gz ] || [ ! -e ${RTG_OVERLAP_SCRATCH}/fn.vcf.gz ]; then
	${RTG_VCFEVAL} --all-records -b ${GOLD_CALLS_VCFGZ} -c ${input_vcfgz_path} -t ${REFERENCE_SDF} -o ${RTG_OVERLAP_SCRATCH} > /dev/null 2>&1
elif [ ${resume} -eq 1 ]; then
	echo -e "\033[0;32mOverlap output files found and resume flag set; not recomputing.\033[0m"
fi

# TODO: Parse ${RTG_OVERLAP_SCRATCH}/vcfeval.log to identify regions to exclude
# eg Evaluation too complex (5001 unresolved paths, 18033 iterations) at reference region 2:105849275-105849281. Variants in this region will not be included in results.
# This will be required to remove the TNs in this region, but it's polish.


#####################################################################
# CALCULATIONS FOR REPORT
#####################################################################
echo -e "\033[0;32mPerforming calculations for report...\033[0m"
mkdir -p ${KNITR_SCRATCH}
# knitr doesn't play well with building knits outside of its working
# directory.  Currently we get around this with a bit of a kludge, 
# by copying the report files to ${KNITR_SCRATCH}, then executing in 
# that directory.
cp -f report.Rnw ${KNITR_SCRATCH}
cp -f report_functions.R ${KNITR_SCRATCH}
cp -f report_debug.Rnw ${KNITR_SCRATCH}
cp -f report_calculations.R ${KNITR_SCRATCH}
cd ${KNITR_SCRATCH}

# Run the script
# SECURITY WARNING: Code injection possible in VERSION_ variables.
# Ensure that git branch and execution host names can not be under
# malicious control.
if [ ${resume} -eq 0 ] || [ ! -e report_data.rda ]; then
	${RSCRIPT} --vanilla report_calculations.R ${debug} ${debug_chrom} ${extended} ${input_vcfgz_path} ${RTG_OVERLAP_SCRATCH}/tp.vcf.gz ${RTG_OVERLAP_SCRATCH}/fp.vcf.gz ${RTG_OVERLAP_SCRATCH}/fn.vcf.gz ${GOLD_CALLS_VCFGZ} ${REFERENCE_BSGENOME} ${GOLD_HARDMASK_VALID_REGIONS_BEDGZ} ${FUNCTIONAL_REGIONS_BEDGZ_PREFIX} ${MASK_REGIONS_BEDGZ_PREFIX} "'""${VERSION}""'" "${VERSION_GIT_BRANCH}" "${VERSION_GIT_COMMIT}" "'""${VERSION_EXEC_HOST}""'"
elif [ ${resume} -eq 1 ]; then
	echo -e "\033[0;32mReport calculation output files found and resume flag set; not recomputing.\033[0m"
fi


#####################################################################
# REPORT GENERATION
#####################################################################
echo -e "\033[0;32mGenerating report...\033[0m"

${RSCRIPT} --vanilla -e "library(knitr); knit('report.Rnw', output = 'report.tex')"

# Latex often 'fails' (returns a nonzero exit status), but still 
# generates a report.  Keep going when this happens, and test for
# failure explicitly later.
set +e

# Remove the report.pdf that may be present in the scratch directory,
# so we can later check whether pdflatex successfully built a report
# or not.
rm -f ${KNITR_SCRATCH}/report.pdf

# Run pdflatex
pdflatex -interaction nonstopmode report.tex

# Check  whether the report.pdf was generated
if [ ! -e ${KNITR_SCRATCH}/report.pdf ]; then
	echo -e "\033[0;31mError: pdflatex did not successfully generate a report.pdf.\033[0m"
	echo -e "\033[0;31mCheck ${KNITR_SCRATCH}/report.tex and the latex log for errors.\033[0m"
	exit 4
fi

# Copy the completed report to the final destination
cp "${KNITR_SCRATCH}/report.pdf" "${output_pdf_path}"

echo -e "\033[0;32mReport generated successfully.\033[0m"

# Clean up scratch if everything worked OK, and we're not a debug
# run.
# if [ $debug -eq 0 ]; then
# 	echo -e "\033[0;32mClearing scratch space...\033[0m"
# 	rm -rf ${SCRATCH}
# else
# 	echo -e "\033[1;36mDebug run; not clearing scratch space.\033[0m"
# 	echo -e "\033[1;36mScratch can be found at: \"${SCRATCH}\"\033[0m"
# fi

echo -e "\033[0;32mDone.\033[0m"
