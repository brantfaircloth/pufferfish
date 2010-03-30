#!/usr/bin/env python
# encoding: utf-8
"""
geneMixMax.py

Created by Brant Faircloth on 2010-02-05.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

Run with:  python geneMixMax.py -c db.conf

Purpose:  find the start and end bases of a putative gene in the tetnig tables
and output those starts and ends to the tetnig.genes table

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

    p.add_option('--configuration', dest='conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
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


def minMax(cur, mrna_id):
    s = ()
    e = ()
    cur.execute('''SELECT chromo, start, end from exons where mrna_id = ?''', (mrna_id,))
    #pdb.set_trace()
    for row in cur.fetchall():
        chromo, start, end = row
        s += (start,)
        e += (end, )
    return chromo, min(s), max(e)

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
    cur.execute('''SELECT id, mrna_id from genes''')
    genes = cur.fetchall()
    for gene, mrna_id in genes:
        chromo, start, end = minMax(cur, mrna_id)
        #pdb.set_trace()
        cur.execute('''update genes set chromo = ?, start= ?, end = ? where id = ?''', (chromo, start, end, gene))
    cur.close()
    conn.commit()
    conn.close()
    
if __name__ == '__main__':
    main()     