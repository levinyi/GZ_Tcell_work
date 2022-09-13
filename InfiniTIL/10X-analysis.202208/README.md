# 10X-analysis
Contains functions and pipelines for analyzing 10X single-cell data.
This includes sets of wrappers to make analyzing data with Seurat simpler to
follow and less verbose (at least when using similar parameters and settings)
to our "standard".\
In addition, this repository includes a set of tools to help automate this
"standard" analysis.\
As an additional note, the goal of this analysis is to perform basic,
relatively "generic" calculations that are performed on all 10X data. Further
characterization and analysis specific to each dataset are not implemented here.
The pipeline implemented in this repository does, however, convert the final
output to a format easily loadable by BioTuring Browser, facilitating manual
and visual exploration without knowledge of any programming language.

Built and maintained by: Seon Kinrot (RootPath Genomics)

Current modules are:
1. R files defining helper functions and wrappers (under src/)
   1. SeuratWrappers
   2. PlotFunctions
   3. AutomationFunctions
2. Notebooks and example scripts (under Examples/)
   1. Example_notebook (needs update)
3. ShinyQC
4. R scripts defining CLIs and automation pipelines:
   1. ShinyQC.R - a CLI interface for calling ShinyQC
   2. BatchCR.R - a CLI interface to automate CellRanger runs on multiple files
   and/ or samples
   3. BatchSeurat.R - a CLI interface for automating Seurat analysis of data
   after CellRanger. One key feature is being able to specify individual
   samples with positional arguments, rather than having to run on all data in
   a folder

## SeuratWrappers
A collection of R functions to streamline 10X data analysis in Seurat.
Contains functions to load, integrate and merge data from single-cell GEX,
VDJ and CSP data, as well as to automate some "standard" workflows
(e.g. dimensionality reduction and clustering) with less manual steps.

## PlotFunctions
A collection of R functions to be used in conjunction with SeuratWrappers.
The functions in this file automate the process of generating some routine
plots from 10X data. They also contain some manually defined sets of genes
and signatures, avoiding the need to copy and paste them to every individual
analysis.

## AutomationFunctions
Contains a set of wrapper functions for specific portions of the automation
pipeline in BatchAnalysis and BatchCR. This helps in breaking down the pipeline,
de-cluttering the main script and making it more readable.

## Example_notebook
An R Markdown notebook demonstrating the use of SeuratWrappers and
PlotFunctions to analyze a 10X dataset. Specifically, this is based on a
CITE-seq experiment where 8 individual samples were integrated with GEX, VDJ
and CSP. Of note, this is a bit more cumbersome than using automated scripts
(such as BatchAnalysis), but has the advantage of more intermediate QC
opportunities, more flexibility (e.g. in setting up subgroups for merging of
VDJ data) and the possibility of additional "non-standard" analyses.

** This notebook has not been updated in a while and may contain some
deprecated behaviors and notations

## ShinyQC
A Shiny app that takes in a list of Seurat objects and lets the user go through
them one by one and perform QC filtering interactively.\
Example screenshot:\
!["Example screenshot from app"](./ShinyQC/www/Screenshot.png
                                 "Screenshot from app")

Can also be run as a standalone module by saving a list of Seurat object in
ShinyQC/data/input.rds and calling ShinyQC.R. Post-QC data will be saved under
ShinyQC/data/output.rds. Usage:
```
>>Rscript ShinyQC.R --help

Usage: ShinyQC [-h] [-p PORT]
 Command line interface for automated 10X QC filtering

Options:
	-h, --help
		Show this help message and exit

	-p PORT, --port=PORT
		Set the port to use for the Shiny App (default: 8886)

For support, contact Seon Kinrot (RootPath Genomics)
```

## BatchCR
An R script defining a CLI interface, meant to help automate CellRanger runs.
This module is essentially identical to the CellRanger portion of BatchAnalysis.
Usage and help text:
```
>> Rscript BatchCR.R --help

Usage: BatchCR [-h] [-u CPUS] [-x MAX_MEM]
[-g] [-c | -s] [-m] [-v]
-t TRANSCRIPTOME [-j VDJ_REF]
[-f FEATURE_REF] [-i INPUT]
[sample 1] [sample 2]...
 Command line interface for automated CellRanger count runs

Options:
	-h, --help
		Show this help message and exit

	-u CPUS, --cpus=CPUS
		The number of cores to use for CellRanger (default 24)

	-x MAX_MEMORY, --max_memory=MAX_MEMORY
		Max memory usage in CellRanger, in GB (default 256

	-g, --gex_only
		Only analyze GEX data, ignoring VDJ.
		Note this has no bearing on CSP flags

	-c, --cite_seq
		Specify this data includes CSP.
		By default, assumes CSP and GEX are fused
		(i.e. 'sep_csp' is False)

	-s, --sep_csp
		Specify that GEX and CSP matrices are separate.
		If this flag is passed, 'cite_seq' is set to True

	-m, --mouse
		Specify this is a mouse dataset [default is human]

	-v, --verbose
		Print status updates as analysis advances
		[default is not to print updates]

	-t TRANSCRIPTOME, --transcriptome=TRANSCRIPTOME
		The path to the transcriptome reference
		file. *This is a required argument*

	-j VDJ_REF, --vdj_ref=VDJ_REF
		The path to the VDJ reference file
		(required unless the '-g' flag is passed

	-f FEATURE_REF, --feature_ref=FEATURE_REF
		The path to the feature reference file
		(required if passing the '-c' or '-s' flags)

	-i INPUT, --input=INPUT
		Set the path to the input directory containing the 10X analysis
		[default is the current working directory]

For support, contact Seon Kinrot (RootPath Genomics)
```

## BatchSeurat
An R script defining a CLI interface, meant to help automate Seurat analysis on
post-Cellranger data.
This module is essentially identical to the Seurat portion of BatchAnalysis.
Usage and help text:
```
>> Rscript BatchSeurat.R --help
Usage: BatchSeurat [-h] [-g] [-c | -s] [-m] [-v]
[-i INPUT] [-o OUTPUT]
[sample 1] [sample 2]...
 Command line interface for automated CellRanger count runs

Options:
	-h, --help
		Show this help message and exit

	-g, --gex_only
		Only analyze GEX data, ignoring VDJ.
		Note this has no bearing on CSP flags

	-c, --cite_seq
		Specify this data includes CSP.
		By default, assumes CSP and GEX are fused
		(i.e. 'sep_csp' is False)

	-s, --sep_csp
		Specify that GEX and CSP matrices are separate.
		If this flag is passed, 'cite_seq' is set to True

	-m, --mouse
		Specify this is a mouse dataset [default is human]

	-v, --verbose
		Print status updates as analysis advances
		[default is not to print updates]

	-i INPUT, --input=INPUT
		Set the path to the input directory containing the 10X analysis
		[default is the current working directory]

	-o OUTPUT, --output=OUTPUT
		Set the path to which output files are
		saved [default is 'input'/analysis]

For support, contact Seon Kinrot (RootPath Genomics)
```
