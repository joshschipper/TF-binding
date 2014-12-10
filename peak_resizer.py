#!/usr/bin/python

print '''\nResizing peaks based on the center, or the summit, if supplied'''


import argparse, os, string, sys

#===============================================================================
#  Command Line Arguments

parser = argparse.ArgumentParser(description = 'E2F Binding Model')
parser.add_argument('-t', metavar = 'CenterType', 
                    help = 'Get sequence around peak "center" (the middle of the start and stop position), or "summit" (if available, as the last column).',
                    dest = 'centertype' ,
                    choices = ['center','summit'] ,
                    required=True)
parser.add_argument('-i', metavar = 'SequenceFile', 
                    help = 'Sequence file (without header) where the first three columns are "Chromosome Name", "Start Position", and "Stop Position" (i.e. .bed format). If centering by peak summit, an additional last column for distance from left of start position is required.' , 
                    dest = 'seqfile' , 
                    required=True)
parser.add_argument('-s', metavar = 'newsize', 
                    help = 'The distance upstream and downstream of the chosen center type; i.e. "75" will result in sequences 75bp up and downstream of the center, for 151bp total length.' , 
                    dest = 'size' , 
                    required=True)
parser.add_argument('-o', metavar = 'OutFilePrefix', 
                    help = 'Optional, the prefix that all output files will be based on (do not include file extension).', 
                    dest = 'outprefix')
args = parser.parse_args()
centertype = args.centertype
seqfile = args.seqfile
size = int(args.size)
if args.outprefix: outprefix = args.outprefix
else: outprefix = os.path.splitext(seqfile)[0] + '_' + size + 'bp_around_peak_' + centertype

    
#===============================================================================

def read_data(filename):
    ''' Creates an array for every cell in a tab separated text file'''
    data = []
    f = open(filename, 'r+') #opens the file as "f"
    for line in f: #for each line in the file
        l = string.split(line.strip(), '\t') #removes any carriage returns, then splits the tab separated line into columns
        data.append(l) #Add each line to the array
    f.close() #closing the file
    return data

def get_peak_central_seq(chipfile,n,centertype):
    '''Takes a chip file (.bed or .pk format, with first 3 columns Chr, Start, and End, and last column peak summit offset
    if using peak summits) and adjusts the numbers to be a sequence with +/- n bases around center of the peak or peak summit.
    "centertype" is either "center", or "summit".'''
    data = read_data(chipfile)
    outfile = outprefix + '.txt'
    if centertype == "summit":
        if len(data[0]) < 4 or not isinstance( data[0][-1], int ): sys.exit("Cannot find the peak summits in this file.")
        print "Getting peak summit +/-", n, "bases for", chipfile
        f=open(outfile,'w',1)
        for n1 in range(len(data)): 
            chrom = data[n1][0]
            start = int(data[n1][1])
            end = int(data[n1][2])
            summit = int(data[n1][-1])
            if summit >= start and summit <= end: 
                center = summit + start #for centering around the peak summit
            else: sys.exit("The value in the last column is not between the start and stop positions, and so is probably not the peak summit.")
            newstart = center - (n)
            newend = center + (n) +1
            newline = [chrom,newstart,newend]
            print >>f, "\t".join(map(str,newline))
        f.close()
    elif centertype == "center":
        print "Getting peak center +/-", n, "bases for", chipfile
        f=open(outfile,'w',1)
        for n1 in range(0,len(data)): #If file has a header line, change 0 to 1 here.
            chrom = data[n1][0]
            start = int(data[n1][1])
            end = int(data[n1][2])
            center = ((end - start)/2)+start #for if we want to senter around the middle of the whole sequence
            newstart = center - (n)
            newend = center + (n) +1
            newline = [chrom,newstart,newend]
            print >>f, "\t".join(map(str,newline))
        f.close()
    print "Writing new sequence to", outfile

get_peak_central_seq(seqfile,size,centertype)
