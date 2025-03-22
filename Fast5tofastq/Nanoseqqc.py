##A Snakemake workflow to process Nanopore fastq files.

import glob
import os
from os.path import join, basename, dirname
import pandas as pd
from snakemake.io import expand
from snakemake.utils import R
from snakemake.io import glob_wildcards
import re
from os.path import join




#################
HOST = "jiyang.jiang"
HOMEDIR = "/jiyang.jiang"

################
SANDBOX = "lustre/scratch/ad_app_jiyang/sandbox"
GUID = "6fd9bff1-4cb1-43a7-a88f-a8c116c8fdeb"
INPUTDIR = SANDBOX + "/" + GUID + "/" + "Input"
OUTPUTDIR = SANDBOX + "/" + GUID + "/" + "Output"

GPUHOST = HOST + "@gpuserver"
GPUWORKDIR = "/home/" + HOST + "/" + "sandbox"

FAST5DIR = INPUTDIR + "/" + "fast5"
RAWFASTQDIR = INPUTDIR + "/" + "fastq"
DMP_FASTQDIR = INPUTDIR + "/" + "demultiplexed"

EXP_FASTQDIR = OUTPUTDIR + "/" + "Fastq"
FastqQCDIR  = OUTPUTDIR+ "/fastqc"
MinIONQCDIR = OUTPUTDIR + "/" + "MinIONQC" + "/fastq"

SMAPLESHEET = "SampleSheet.csv"

################Extract nanopore run information

#porerefiner_ver,1.0.0,,,,,
#library_id,TEST_TEST,,,,,
#sequencing_kit,TEST_KIT,,,,,
#sample_id,accession,barcode_id,organism,extraction_kit,comment,user
#TEST01,ACC_TEST_01,1,Salmonella coli,TEST_KIT_KIT,"a comment, with a comma in it",jiyang.jiang@gmail.com
#TEST02,ACC_TEST_02,2,Escherichia enteriditis,TEST_KIT_KIT,,jiyang.jiang@gmail.com
#TEST03,ACC_TEST_03,3,Campy stampy,TEST_KIT_KIT,,jiyang.jiang@gmail.com
#TEST04,ACC_TEST_04,1,Salmonella coli,TEST_KIT_KIT,"a comment, with a comma in it",jiyang.jiang@gmail.com
#TEST05,ACC_TEST_05,2,Escherichia enteriditis,TEST_KIT_KIT,,justin.payne@fda.hhs.gov
#TEST06,ACC_TEST_06,3,Campy stampy,TEST_KIT_KIT,,jiyang.jiang@gmail.com
#TEST07,ACC_TEST_07,1,Salmonella coli,TEST_KIT_KIT,"a comment, with a comma in it",jiyang.jiang@gmail.com
#TEST08,ACC_TEST_08,2,Escherichia enteriditis,TEST_KIT_KIT,,jiyang.jiang@gmail.com
#TEST09,ACC_TEST_09,3,Campy stampy,TEST_KIT_KIT,,jiyang.jiang@gmail.com
#TEST10,ACC_TEST_10,1,Salmonella coli,TEST_KIT_KIT,"a comment, with a comma in it",jiyang.jiang@gmail.com
#TEST11,ACC_TEST_11,2,Escherichia enteriditis,TEST_KIT_KIT,,jiyang.jiang@gmail.com
#TEST12,ACC_TEST_12,3,Campy stampy,TEST_KIT_KIT,,jiyang.jiang@gmail.com
#
#


data1 = pd.read_csv("/" + INPUTDIR + "/" + SMAPLESHEET, sep=",", nrows=2)
LIBID = str(data1.iloc[:,1][0])
print (LIBID)
SEQKIT = str(data1.iloc[:,1][1])
print (SEQKIT)
data2 = pd.read_csv("/" + INPUTDIR + "/" + SMAPLESHEET, sep=",", skiprows=3)
SAMPLESID = list(data2.sample_id)
print (SAMPLESID)
NUMSAM = len(SAMPLESID)
print (NUMSAM)
CFSANACC = list(data2.accession)
print (CFSANACC)
ORGANISM = list(data2.organism)
print (ORGANISM)
EXTRAKIT = list(set(data2.extraction_kit))
print (EXTRAKIT)
USER = list(set(data2.user))
print (USER)
BARCODES = ['barcode' + "{:02}".format(num) for num in data2.barcode_id]
print (BARCODES)
EXTRAKIT = set(data2.extraction_kit)
print (EXTRAKIT)


#zip CFSAN accession number with barcodes based on sample sheet.
CFSAN2BAR= dict(zip(CFSANACC, BARCODES))

#Extract Run id from raw fastq directory
#file = [f for f in glob.glob("/"+ RAWFASTQDIR + "/" + "*.fastq.gz")][0]
#if os.path.exists (file):
#    runid = file.split('_')[-3]
#else:
#    print ("File not exist")

#print (runid)

sam2bar = dict(zip(CFSANACC,BARCODES))

print (sam2bar.keys())
print (sam2bar.values())

ruleorder: demultiplex > MinIONQC_per_run > nanoplot_run> merge_fastq

rule all:
    input:
        DMP_FASTQDIR,
        "/" + MinIONQCDIR + "/fastq",
        expand("/" + OUTPUTDIR + "/" + "{sample}" + "_" + "nanomerge.fastq", sample=CFSANACC),
        "/" + OUTPUTDIR + "/" + "NanoPlot_PreRun.done",
        expand("/" + OUTPUTDIR + "/" + "{sample}_nanoplotqc", sample=CFSANACC),


rule demultiplex:
    input:
        glob.glob("/" + RAWFASTQDIR + "/*fastq.gz"),
    output:
        directory(DMP_FASTQDIR),
    shell:
        """
        module load guppy/3.3.3
        guppy_barcoder --input_path /{RAWFASTQDIR}  --save_path {output} -r
        module unload guppy/3.3.3
        """

rule MinIONQC_per_run:
    input:
        "/" + RAWFASTQDIR
    output:
        directory("/" + MinIONQCDIR + "/fastq"),
        touch("/" + OUTPUTDIR + "/" + "minionqc.done")
    shell:
        """
        module load MinIONQC/1.4.1
        mkdir -p /{MinIONQCDIR}/fastq
        MinIONQC.R -i {input} -o /{MinIONQCDIR}
        module unload MinIONQC/1.4.1
        """

rule nanoplot_run:
    input:
         data = "/"+ RAWFASTQDIR + "/" + "sequencing_summary.txt"
    output:
         touch("/" + OUTPUTDIR + "/" + "NanoPlot_PreRun.done"),
         data = directory("/" + OUTPUTDIR + "/" + "summary_prerun")
    shell:
         """
         module load NanoPlot/1.28.1
         NanoPlot --summary {input.data} --loglength -o {output.data}
         module unload NanoPlot/1.28.1
         """

rule mergefastq&nanoplotqc:
    input:
        fastq = lambda wc: glob.glob("/" + DMP_FASTQDIR + "/%s/*fastq" % sam2bar[wc.sample])
    output:
         fastq = "/" + OUTPUTDIR + "/" + "{sample}" + "_" + "nanomerge.fastq",
         dir = directory("/" + OUTPUTDIR + "/" + "{sample}" + "_" + "nanoplotqc"),
    shell:
        """
        cat {input.fastq} > {output.fastq}
        module load NanoPlot/1.28.1
        NanoPlot -t 2 --fastq {input.fastq} --plot hex dot -o {output.dir}
        module unload NanoPlot/1.28.1
        """

