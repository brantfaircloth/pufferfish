#!/usr/bin/env python
# encoding: utf-8
"""
makeBedFromRegion.py

Created by Brant Faircloth on 2010-02-05.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

Run with:  python makeBedFromRegion.py -c db.conf

Purpose:  create a bed file from the records in the region table:

[Database]
DATABASE = db_name
USER = username
PASSWORD = password

"""

import pdb
import os
import sys
import oursql
import ConfigParser
import optparse


def interface():
    '''Get the starting parameters from a configuration file'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', '-c', dest = 'conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
metavar='FILE')
    p.add_option('--out', '-o', dest = 'out', action='store', \
type='string', default = None, help='The path to the output file.', \
metavar='FILE')
    p.add_option('--advanced', '-a', action='store_true', dest='advanced', \
default=False, help='Produce an advanced bed file')

    (options,arg) = p.parse_args()
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg


def simple(outp, data):
    for record in data:
        #pdb.set_trace()
        outp.write('%s %s %s region.id_%s\n' % (record[3], record[4], record[5], record[0]))

def advanced(cur, outp, gene, mrna):
    # get gene specific data
    query = ('''SELECT starts_exon, chromo, start, end from exons where mrna_id = %s order by start''' % (mrna))
    cur.execute(query)
    data = cur.fetchall()
    blockCount = len(data)
    chromStart = data[0][2]
    chromEnd = data[-1][3]
    blockSizes, blockStarts = [], []
    for each in data:
        exon, chrom, s, e = each
        bsz = e - s
        rst = s - chromStart
        blockSizes.append(str(bsz))
        blockStarts.append(str(rst))
    #pdb.set_trace()
    blockSizes = ','.join(blockSizes)
    blockStarts = ','.join(blockStarts)
    outp.write('%s %s %s gene_%s_mrna_%s 1000 + %s %s 0 %s %s, %s\n' % (chrom, chromStart, chromEnd, gene, mrna, chromStart, chromEnd, blockCount, blockSizes, blockStarts))

def main():
    conf = ConfigParser.ConfigParser()
    options, arg = interface()
    conf.read(options.conf)
    conn = oursql.connect(
            user=conf.get('Database','USER'), 
            passwd=conf.get('Database','PASSWORD'), 
            db=conf.get('Database','DATABASE')
            )
    cur = conn.cursor()
    outp = open(options.out,'w')
    outp.write('''track name=geneRegions description="Tetraodon Exons" useScore=0\n''')
    if not options.advanced:
        cur.execute('''SELECT * FROM exons''')
        data = cur.fetchall()
        simple(outp, data)
    if options.advanced:
        # get list of "genes"
        cur.execute('''SELECT id, mrna_id from genes''')
        for gene, mrna in cur.fetchall():
            advanced(cur, outp, gene, mrna)
    outp.close()
    cur.close()
    conn.close()


if __name__ == '__main__':
    main()