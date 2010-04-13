#!/usr/bin/env python
# encoding: utf-8

"""
fuguGeneRegionChecker.py

Created by Brant Faircloth on 05 April 2010 15:55 PDT (-0700).
Copyright (c) 2010 Brant Faircloth. All rights reserved.

This program checks the exons associated with genes in Tetraodon to make sure
that the same exons (e.g. the gene) were located in Fugu.  In essence, this
program determined the genes present in Tetraodon are also present in Fugu and
it returns those genes present in both.
"""

import pdb
import os
import sys
import oursql
import optparse
import ConfigParser


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', dest = 'conf', action='store', 
type='string', default = None, help='The path to the configuration file.', 
metavar='FILE')

    (options,arg) = p.parse_args()
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg 


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
    cur.execute('''SELECT id, mrna_id, start, end FROM genes''')
    geneList = []
    for row in cur.fetchall():
        # tetraodon-specific info.
        tet_genes = set()
        tet_gene_id, tet_mrna_id, tet_start, tet_end = row
        tet_gene_span = tet_end - tet_start
        # fugu-specific data
        fugu_matches = set()
        fugu_chromo = set()
        fugu_start = []
        fugu_end = []
        fugu_coverage = []
        fugu_identity = []
        # get the exon ids associated with tetraodon mrna
        cur.execute('''SELECT id FROM exons WHERE mrna_id = ?''', (tet_mrna_id,))
        result = cur.fetchall()
        for r in result:
            tet_genes.add(r[0])
        # get the exons found by lastz
        cur.execute('''SELECT exons_id, name1, zstart1, end1, percent_coverage, percent_identity FROM fugumatch WHERE mrna_id = ? AND duplicate = 0''', (tet_mrna_id,))
        result = cur.fetchall()
        if result:
            for r in result:
                fugu_matches.add(r[0])
                fugu_chromo.add(r[1])
                fugu_start.append(r[2])
                fugu_end.append(r[3])
                fugu_min = min(fugu_start)
                fugu_max = max(fugu_end)
                fugu_span = fugu_max - fugu_min
                fugu_coverage.append(float(r[4]))
                fugu_identity.append(float(r[5]))
                span_diff = abs(tet_gene_span - fugu_span)
            #pdb.set_trace()
            if tet_genes == fugu_matches:
                coverage = sum(fugu_coverage)/len(fugu_coverage)
                identity = sum(fugu_identity)/len(fugu_identity)
                geneList.append({'gene_id':tet_gene_id, 'mrna_id':tet_mrna_id, 'start':fugu_min, 'end':fugu_max, 'span_diff':span_diff, 'chromo':list(fugu_chromo)[0], 'coverage':coverage, 'identity':identity})
    for record in geneList:
        cur.execute('''INSERT INTO fugugenes VALUES (?,?,?,?,?,?,?,?)''', (record['gene_id'], record['mrna_id'], record['chromo'], record['start'], record['end'], record['span_diff'], record['coverage'], record['identity']))
    conn.close()

if __name__ == '__main__':
    main()