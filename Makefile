#####################################################################
# 
#  Makefile overseeing the full validation reporter pipeline
#  
#  Execute as:
#    make INPUT=<path to input vcf.gz>
#  which will create a final version report in report.pdf.
#  
#  To create a report with additional debugging output, use
#    make INPUT=<path to input.vcf.gz> DEBUG=TRUE
#  which will create the debug report in report_debug.pdf.
#  
#  Additional targets of interest:
#    clean: Remove latex intermediate files
#    ultraclean: clean, and additionally remove the output 
#      and overlap working files.
#  
#  
#  Mark Pinese
#  25 May 2015	MP 	Started writing
#  
#####################################################################

.PHONY: all clean ultraclean report
.SECONDARY: 


RSCRIPT=~/bin/Rscript
JAVA=/usr/java/latest/bin/java
RTG_CORE=~/software/rtg-core/build/rtg-core.jar
RTG_THREADS=4
RTG_VCFEVAL=$(JAVA) -Xmx4G -jar $(RTG_CORE) vcfeval -T $(RTG_THREADS)


all: report


clean: 
	rm -f report.tex report.aux report_debug.tex report_debug.aux figure/* cache/*


ultraclean: clean
	rm -f report.pdf report_debug.pdf scratch/overlap/*


ifeq ($(DEBUG), TRUE)
report: report_debug.pdf
else
report: report.pdf
endif


report.pdf: report.Rnw 
	$(RSCRIPT) -e "knitr::knit2pdf('report.Rnw', 'report.pdf')" <calls>


report_debug.pdf: report.Rnw 
	$(RSCRIPT) -e "knitr::knit2pdf('report.Rnw', report_debug.pdf)" -d <calls>


scratch/overlap/fn.vcf.gz scratch/overlap/fp.vcf.gz scratch/overlap/tp.vcf.gz: data/giab.vcf.gz data/1000g_v37_phase2.sdf/* $(INPUT)
	mkdir -p scratch/overlap && \
    $(RTG_VCFEVAL) -b data/giab.vcf.gz -c $(INPUT) -t data/1000g_v37_phase2.sdf/ -o scratch/overlap


