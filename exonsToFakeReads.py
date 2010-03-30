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

import pdb
import os
import sys
import oursql
import bx.seq.twobit
import ConfigParser
import optparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def interface():
    '''Get the starting parameters from a configuration file'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--twobit', dest = 'twobit', action='store', \
        type='string', default = None, \
        help='The path to the 2bit file.', metavar='FILE')
    
    p.add_option('--configuration', dest = 'conf', action='store', \
        type='string', default = None, \
        help='The path to the configuration file.')
        
    p.add_option('--exons', dest = 'exons', action='store', \
        type='string', default = None, \
        help='The table containing the exons.')
    
    p.add_option('--genes', dest = 'genes', action='store', \
        type='string', default = None, \
        help='The table containing the genes.')

    p.add_option('--output', dest = 'output', action='store', \
        type='string', default = None, \
        help='The path to the output file.', metavar='FILE')
        
    (options,arg) = p.parse_args()
    
    if not options.conf:
        p.print_help()
        sys.exit(2)
    
    return options, arg


def getExons(cur, mrna_id, options):
    '''get the exons associated with a particular gene region'''
    query = '''SELECT starts_exon, chromo, start, end from %s where \
        mrna_id = %s order by start ASC''' % (options.exons, mrna_id)
    cur.execute(query)
    return cur.fetchall()

def exonStitcher(cur, gene, mrna_id, exons, tb):
    '''for a given gene, stitch the sequence of the exons together into a
    cumulative whole that represents a quasi-mRNA'''
    sequence = None
    exonIds = ''
    for exon in exons:
        exonIds += '%s,' % exon[0]
        se, chromo, start, end = exon
        if not sequence:
            sequence = Seq(tb[chromo][start:end], IUPAC.unambiguous_dna)
        else:
            sequence += Seq(tb[chromo][start:end], IUPAC.unambiguous_dna)
    record = SeqRecord(sequence)
    record.id = 'Tetraodon_Gene_%s' % (gene)
    record.name = record.id
    record.description = 'Tetraodon putative gene %s, mrna_id = %s, exons %s' % (gene, mrna_id, exonIds)
    return record

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
    tb  = bx.seq.twobit.TwoBitFile(file(os.path.abspath(options.twobit)))
    # get list of genes
    query = '''SELECT id, mrna_id from %s''' % options.genes
    cur.execute(query)
    genes = cur.fetchall()
    holder = []
    for gene, mrna_id in genes:
        exons   = getExons(cur, mrna_id, options)
        rna     = exonStitcher(cur, gene, mrna_id, exons, tb)
        holder.append(rna)
    SeqIO.write(holder, open(options.output, 'w'), 'fasta')


if __name__ == '__main__':
    main()

