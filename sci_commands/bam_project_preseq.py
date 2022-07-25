#!~/miniconda2/envs/python3.7/bin/python3

from pprint import isrecursive
import pandas as pd
import sys
import argparse
from multiprocessing import Pool, Process
from subprocess import Popen, PIPE
import shlex
import os
import numpy as np
import pysam
from collections import defaultdict
import argparse
import re
from time import time
from tqdm import tqdm #For progress bar use
import seaborn as sns
from collections.abc import Mapping, Container
from sys import getsizeof
import psutil
import matplotlib.pyplot as plt
import random

### Current plan is as follows: 
'''Args:
- BAM file (name sorted probably)
- Complexity.txt file
- [probably] the cutoff used for filtering out real cells, no reason to get estimates on junk
- ?

Method: read in the complexity.txt file and use the cutoff to only keep real cells,
then use pysam with the bam file to make a ldict of datafames(?) of read locations (chr start end UID strand) for each barcode
could store as {chr_start_end_strand: count} format defaultdict
using the above information, generate a histogram of nCopies \t count (i.e. how many reads are unique, doubled etc) for each barcode
From there use os to run something along the lines of 

def estimate_total_complexity(barcode):
    <write a string barcodeHist of all the values in the histogram starting a singlets, one entry per \n>
    command = shlex.split('{} pop_size -H /dev/fd/0'.format(options.preseq))
    p1 = Popen(command, stdin=PIPE, stdout=PIPE)
    stdout = p1.communicate(input=[barcodeHist])[0]
    <store some version of this stdout, there should be 3 values mean, lower c.i. higher c.i.>
    p1.wait() ### Moved in test
    p1.stdout.close()
    #return stdout? or update a global dict (lol) with the values

with Pool(options.threads) as p:
    p.map(estimate_total_complexity, barcodeDict [,chunksize]) #Probably want a different version of map, see https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.map with a chunksize >1
    #Note that starmap or starmap_async might be good options if I end up using multiple variables to estimate_total_complexity
    # Note also, that if I wanted it to maybe write to a file in realtime, map_async would probably be what I want.
    ## You can't yield a mapresult
    p.close()
    p.join()

do something with the resulting dict of values, probably spit it out to a file and plot it was a violin plot with seaborn
'''

def parse_args(args):
    parser = argparse.ArgumentParser(add_help=False, description='Estimates the total complexity of reads per cell in a BAM file, using preseq pop_size.  Requires a bam file and either a complexity.txt file with a cutoff value, or else an annotation file')
    fileGroup=parser.add_argument_group("File Options:")
    fileGroup.add_argument('-b', '--bam', type=str, nargs='+', default='', help='BAM or CRAM file, can be either name or position sorted. If CRAM, needs --reference_filename also.')
    fileGroup.add_argument('-c', '--complexity' , nargs = '?', type =str, default='', help='.complexity.txt file, used with -f')
    fileGroup.add_argument('-r', '-g', '-T', '--reference_filename', type=str, help='Reference genome path/file.fa for CRAM files')
    fileGroup.add_argument('-o','-O', '--output', default='', type=str, help='output prefix to use (defaults to bam prefix)')
    #parser.add_argument('-a', '--annot' , nargs='?', type=str, default='', help='[Not implemented] Annotation file.  If given, will use all cells in this file, skipping the complexity and cutoff files')
    settings=parser.add_argument_group("Input Settings:")
    settings.add_argument('-f', '-n', '--cutoff', type=int, default=1000, help='complexity cutoff (used with complexity.txt file)')
    settings.add_argument('-t', '-@', '--threads', type=int, default=1, help='threads to use')
    settings.add_argument('-S', default=1.0, type=float, help='Same as samtools view -S, subsamples the whole library randomly to that fraction [default=1.0 for all reads]. Consider using this for 20x10^9 read icell8 libraries. The more deeply a library is sequenced, the lower the value of S that will work.')
    settings.add_argument('-M', action='store_true', help='Flag for sciMet Library (skips certian filters during read counting)')    
    settings.add_argument('-q','--mapq', type=int, default=10, help='Mapping quality score threshold, do not consider reads under this threshold of mapping.')
    optional=parser.add_argument_group('Optional Arguments:')
    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
    optional.add_argument('-p', '--preseq' , nargs='?', type=str, default='preseq',help='path to preseq')
    optional.add_argument('-v', action="store_true", help='Verbose output [default OFF]')
    return parser.parse_args()

def debug(info):
	if options.v==True:
		sys.stderr.write('{}'.format(info))
    
def process_barcode_list():
    '''Processes the barcodes in the '''
    barcodeDict={} ### Moved here to make non-global
    total_reads=0 ###Counter for later use
    # if options.annot=='': 
    with open (options.complexity, 'r') as bFile:
        for line in bFile:
            lD=line.split('\t')
            barcode = lD[1]
            if int(lD[3])>=options.cutoff:
                barcodeDict[barcode]= defaultdict(int)
            total_reads+=int(lD[2])
    return int(total_reads/0.70), barcodeDict #Seems to be more accurate for non-mapping stuff

def get_barcode_hist(samfile, mapq, barcodeDict, total_reads):
    '''Just gets the barcode histograms matching the cells'''
    ### Populate the barcode dictionary with unique read counts
    # histDict={}
    RAM1=psutil.Process().memory_info().rss / (1024*1024*1024)
    if options.S==1.0 and options.M!=True: #Using all reads, keep it simple:
        for entry in tqdm(samfile, miniters=1000, total=total_reads):
            if int(entry.mapping_quality)>=mapq  and entry.is_proper_pair and entry.query_name.split(':')[0] in barcodeDict and not entry.is_secondary:
                #removed: and entry.is_read1
                try:
                    barcodeDict[entry.query_name.split(':')[0]]['{}\t{}\t{}1'.format(entry.reference_name.strip('chr'), entry.reference_start, entry.reference_end)]+=1#, lambda entry: '-' if entry.is_reverse==True else '+')]+=1 #taking this out for now
                    ### Note: I took out the strand specification two reads should not be identically positioned in a sciATAC library.  May want to change this, or else add a -S option
                except KeyError: ### if the barcode is not within our acceptable quality parameters
                    continue
    elif options.S==1.0 and options.M==True: #Using all reads, with sciMET (skipping proper pair checks)
        debug('Methylation Library:\n')
        for entry in tqdm(samfile, miniters=1000, total=total_reads):
            # print(entry)
            if int(entry.mapping_quality)>=mapq  and entry.query_name.split(':')[0] in barcodeDict and not entry.is_secondary:
                #removed: and entry.is_read1
                try:
                    barcodeDict[entry.query_name.split(':')[0]]['{}\t{}\t{}1'.format(entry.reference_name.strip('chr'), entry.reference_start, entry.reference_end)]+=1#, lambda entry: '-' if entry.is_reverse==True else '+')]+=1 #taking this out for now
                    ### Note: I took out the strand specification two reads should not be identically positioned in a sciATAC library.  May want to change this, or else add a -S option
                except KeyError: ### if the barcode is not within our acceptable quality parameters
                    continue
    elif options.S!=1.0 and options.M==True: ### ALSO NEEDS WORK
        for entry in tqdm(samfile, miniters=1000, total=total_reads):
            if random.uniform(0,1)<=options.S : #Uniformly subsample reads, do this before evaluating any other conditions to save time.
                if int(entry.mapping_quality)>=mapq  and entry.query_name.split(':')[0] in barcodeDict and not entry.is_secondary: ### Versus what I thought was optimized, this is 20% faster! Should probably run some more optimization tests here
                    # removed: and entry.is_read1
                    try:
                        barcodeDict[entry.query_name.split(':')[0]]['{}\t{}\t{}1'.format(entry.reference_name.strip('chr'), entry.reference_start, entry.reference_end)]+=1#, lambda entry: '-' if entry.is_reverse==True else '+')]+=1 #taking this out for now
                        ### Note: I took out the strand specification two reads should not be identically positioned in a sciATAC library.  May want to change this, or else add a -S option
                    except KeyError: ### if the barcode is not within our acceptable quality parameters
                        continue
            else:
                    continue
    else:
        for entry in tqdm(samfile, miniters=1000, total=total_reads):
            if random.uniform(0,1)<=options.S : #Uniformly subsample reads, do this before evaluating any other conditions to save time.
                if entry.is_proper_pair and int(entry.mapping_quality)>=mapq  and entry.query_name.split(':')[0] in barcodeDict and not entry.is_secondary: ### Versus what I thought was optimized, this is 20% faster! Should probably run some more optimization tests here
                    # removed: and entry.is_read1
                    try:
                        barcodeDict[entry.query_name.split(':')[0]]['{}\t{}\t{}1'.format(entry.reference_name.strip('chr'), entry.reference_start, entry.reference_end)]+=1#, lambda entry: '-' if entry.is_reverse==True else '+')]+=1 #taking this out for now
                        ### Note: I took out the strand specification two reads should not be identically positioned in a sciATAC library.  May want to change this, or else add a -S option
                    except KeyError: ### if the barcode is not within our acceptable quality parameters
                        continue
            else:
                    continue
    samfile.close() ### To hopefully help free up memory
    debug('Finished with reading BAM file, barcodeDict uses {0:.5f}GB RAM... Starting on histogram generation\n'.format(psutil.Process().memory_info().rss / (1024*1024*1024)-RAM1))
    for barcode in tqdm(barcodeDict, total=len(barcodeDict)):
        tH=defaultdict(int) #temp dict to make histogram
        for entry in barcodeDict[barcode]:
            tH[barcodeDict[barcode][entry]]+=1
        barcodeHist = ''
        for t in sorted(tH):
            barcodeHist+=''.join('{}\t{}\n'.format(t, tH[t]))
        if barcodeHist!='':
            histDict[barcode]=barcodeHist 
        else:
            continue 
    return histDict
    # del(barcodeDict)

def estimate_total_complexity(barcode):
    '''estimates the total complexity for a single barcode'''
    # print(barcode, type(barcode))
    command = shlex.split('{} pop_size -H /dev/fd/0'.format(options.preseq))
    p1 = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout = p1.communicate(input=histDict[barcode].encode())[0]
    if stdout!=b'':
        #Gets the pop size estimate: e.g. ['pop_size_estimate', 'lower_ci', 'upper_ci', '2363040.2', '1933653.9', '4193506.7', '']
        totalUnique = float(re.split('\t|\n',stdout.decode("utf-8"))[3])
        endState=p1.wait()
        p1.stdout.close()
        # output[barcode]=totalUnique ### May want to instead update a structure with total unique so as to save time.
        return (barcode,totalUnique)

def run_process(samfile):
    '''runs the histogram generation to allow the barcodeDict to be dumped from memory'''
    if options.bam[0].endswith('.bam'):
        samfile = pysam.AlignmentFile(options.bam[0], 'rb')
    elif options.bam[0].endswith('.cram'):
        samfile = pysam.AlignmentFile(options.bam[0], 'rc', reference_filename=options.reference_filename)
    total_reads, barcodeDict = process_barcode_list()
    ### TODO: Add in a function to subsample the barcodeDict for fewer cells (specifically for use with iCell8 libraries)
    debug('Finished processing barcode list using {0:.5f}GB memory... Generating complexity histograms.\n'.format(psutil.Process().memory_info().rss / (1024*1024*1024)))
    # debug("Barcode List:\n{}".format(barcodeDict.keys()))
    histDict = get_barcode_hist(samfile, options.mapq, barcodeDict, total_reads) ###Trying Adding in barcodedict for a test
    # print(histDict)
    return [list(barcodeDict.keys()), histDict]

def main(argv):
    '''options is self explanatory.
        The histDict is a dictionary of barcodes:histogram strings
        barcodeDict is where the main data is stored'''
    global options 
    options = parse_args(sys.argv)
    if options.output=='':
        options.output=options.bam[0].strip('.bam')
    global histDict
    histDict={}
    global output
    output=defaultdict(float)
    with Pool(1) as tmp:
        data = tmp.map(run_process, options.bam)
        tmp.close()
        tmp.join()
    debug('Made barcode complexity distributions. Currently using {0:.5f}GB Memory. Starting on estimations\n'.format(psutil.Process().memory_info().rss / (1024*1024*1024)))
    # print(barcodeList[0])
    barcodeList=[]
    histDict = data[0][1]
    for entry in data[0][0]: #Gets rid of barcodes with no reads in histDict
            if entry in histDict:
                    barcodeList.append(entry)
    # print(histDict, barcodeList)
    # print(len(barcodeList), barcodeList)
    with Pool(options.threads) as p:
        o = p.map(estimate_total_complexity, barcodeList)
        p.close()
        p.join()
    
    filename = options.output + ".barcode_preseq_complexity.txt"
    with open(filename, 'a') as OUTFILE:
        for entry in o:
            if entry!= None:
                outstr = entry[0] + "\t" + str(entry[1]) + "\n"
                OUTFILE.write(outstr)

    complexityData = pd.DataFrame(o, columns =['Barcode', 'EstimatedUnique'])# note: we can add the upper and lower confidence intervals if we want to. 
    complexityData2 = complexityData.apply(lambda x: np.log10(x) if x.name == 'EstimatedUnique' else x) #So we get log10(Projected complexity)
    vPlot=sns.violinplot(y=complexityData2['EstimatedUnique'])
    vPlot.set_ylabel("log10 Estimated Total Unique", fontsize=10)
    fig=vPlot.get_figure()
    fig.savefig('{}.estimated.unique.pdf'.format(options.output))

start = time()
if __name__ == "__main__":
	main(sys.argv);
	stop = time()
	m, s = divmod(stop-start, 60)
	h, m = divmod(m, 60)
	sys.stderr.write('Processing completed in {} hours, {} minutes, {} seconds\n'.format(int(h),int(m),int(s)))
	raise SystemExit
