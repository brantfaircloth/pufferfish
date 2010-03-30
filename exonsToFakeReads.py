#!/usr/bin/env python
# encoding: utf-8
"""
exonsToFakeReads.py

Created by Brant Faircloth on 2010-03-09.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

For the tetraodon data, we now know what the exons are that likely make up
the genes of the organism.  But, we need to smash these exons together into
mRNA-like reads, so we can map them over to Danio, in order to get some idea
of the genes with which we will be working.  This program does that.

"""

import sys
import os
import oursql
import bx.seq.twobit


def interface():
    '''Get the starting parameters from a configuration file'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--twobit', dest = 'twobit', action='store', \
        type='string', default = None, \
        help='The path to the 2bit file.', metavar='FILE')
    
    p.add_option('--conf', dest = 'conf', action='store', \
        type='string', default = None, \
        help='The path to the configuration file.')
        
    p.add_option('--table', dest = 'table', action='store', \
        type='string', default = None, \
        help='The table containing the exons.')
        
    (options,arg) = p.parse_args()
    
    return options, arg


def getExons(cur, gene):
    '''get the exons associated with a particular gene region'''
    query = '''SELECT starts_exon, chromo, start, end from %s where \
        mrna_id = %s order by start ASC''' % (options.table, gene)
    cur.execute(query)
    return cur.fetchall()

def exonStitcher(cur, exons):
    '''for a given gene, stitch the sequence of the exons together into a
    cumulative whole that represents a quasi-mRNA'''
    exons = getExons(cur, gene)
    for exon in exons:
        se, chromo, start, end = exon
        # tb[name1][preceding:end1]

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
    tb      = bx.seq.twobit.TwoBitFile(file(options.twobit))
    # get list of genes
    query = '''SELECT distinct mrna_id from %s''' % options.table
    cur.execute(query)
    genes = cur.fetchall()
    for gene in genes:
        exons   = getExons(cur, gene)
        rna     = exonstitcher(cur, exons)
    # table containing exons


if __name__ == '__main__':
    main()

