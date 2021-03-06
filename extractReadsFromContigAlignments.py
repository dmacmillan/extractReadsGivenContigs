__version__ = '0.0.2'

import argparse,os,sys,logging
try:
    import pysam
except ImportError:
    sys.exit("Pysam is a required module for this script, please install Pysam in order to continue.")

parser = argparse.ArgumentParser(description='Given a .bam or .sam alignment file containing contigs aligned to some reference, and .bam or .sam files containing reads aligned to said contigs, return all reads aligned to contigs in .bam or .sam format')
parser.add_argument('c2g',metavar='<contig-alignments>',help='The .bam or .sam file containing your contigs of interest aligned to some reference. Make sure you have the appropriate header for whatever reference you used.')
parser.add_argument('r2cs',metavar='<read-to-contig-alignments>',nargs='+',help='One or more .bam or .sam files containing reads aligned to the contigs in the previous argument.')
parser.add_argument('-o','--outfile',default="./r2c_sorted",help='The file to output the results to. Default is "./r2c_sorted".')
parser.add_argument('-f','--output-format',default='bam',choices=['bam','sam'],help='Specify what format (bam or sam) to output the results in. Default is "bam".')
parser.add_argument('-l','--log-level',choices=['debug','info','warning'],help='Set the level of log messages to see. Default is "warning".')


args = parser.parse_args()

# Determine writing mode
if (args.output_format == 'bam'):
    mode = 'wb'
elif (args.output_format == 'sam'):
    mode = 'wh'

# Set the logging level
if args.log_level:
    logging.basicConfig(filename='debug.log',filemode='w',level=getattr(logging,args.log_level.upper()))

# Open the c2g alignments
try:
    logging.debug('c2g = pysam.AlignmentFile("{}")'.format(args.c2g))
    c2g = pysam.AlignmentFile(args.c2g)
except Exception as e:
    sys.exit('{},{}'.format(type(e),e))

# Create header base
header = {'HD':{'VN':'1.0'},'SQ':[]}

# Dictionary to contain results
logging.debug('results = {}')
results = {}

fasta_file = open('merged.fa','w')

for i,alignment in enumerate(c2g):
    fcode, tcode, lib, name = alignment.query_name.split('_')
    header['SQ'].append({'LN': alignment.reference_length,'SN': alignment.query_name})
    # For each read file
    results[name] = {'name': alignment.query_name,'fcode':fcode,'tcode':tcode,'lib':lib,'reads':[],'index':i}
    fasta_file.write('>{}\n{}\n'.format(alignment.query_name,alignment.seq))
c2g.close()

logging.debug('results = {}'.format(results))
logging.debug('header = {}'.format(header))

#logging.debug('header: {}'.format(header))

#header = None

contigs_to_add = {}
headerlen = len(header['SQ'])

for readfile in args.r2cs:
    try:
        logging.debug('r2c = pysam.AlignmentFile("{}")'.format(readfile))
        r2c = pysam.AlignmentFile(readfile)
    except Exception as e:
        sys.exit('{},{}'.format(type(e),e))
    # For each read in the file matching the contig
    for contig in results:
        for read in r2c.fetch(contig):
            read.query_name = '_'.join(results[contig]['name'].split('_')[:-1]) + '_' + read.query_name
            # Alter the read reference id
            refname = r2c.getrname(read.reference_id)
            read.reference_id = results[refname]['index']
            read.next_reference_id = -1
            results[contig]['reads'].append(read)
    r2c.close()
    
# The reads-to-contigs output
logging.debug('outfile = pysam.AlignmentFile("{}","{}",header={})'.format(args.outfile,mode,header))
outfile = pysam.AlignmentFile(args.outfile + '.' + args.output_format, mode, header=header)
for contig in results:
    for read in results[contig]['reads']:
#        logging.info('outfile.write({})'.format(read.query_name))
#        logging.debug('outfile.write({})'.format(read))
        try:
            outfile.write(read)
        except Exception as e:
            sys.exit('{},{}'.format(type(e),e))

outfile.close()
