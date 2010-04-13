#!/usr/bin/env python
# encoding: utf-8

"""
lastzInserter.py

Created by Brant Faircloth on 05 April 2010 14:07 PDT (-0700).
Copyright (c) 2010 Brant Faircloth. All rights reserved.

Insert alignment records created by lastz alignments into a database using a
configuration file.
"""

import pdb
import os
import sys
import optparse
import oursql
import ConfigParser 


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', dest = 'conf', action='store', 
type='string', default = None, help='The path to the configuration file.', 
metavar='FILE')

    p.add_option('--input', dest = 'input', action='store', 
type='string', default = None, help='The path to the lastz input.', 
metavar='FILE')

    (options,arg) = p.parse_args()
    if not options.conf or not options.input:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf) or not os.path.isfile(options.input):
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
    for line in open(options.input, 'rU'):
        sline = line.strip('\n').split('\t')
        starts_exon, mrna_id = sline[6].split('_')[4], sline[6].split('_')[6]
        cur.execute('''SELECT id FROM exons WHERE mrna_id = ? AND starts_exon = ?''', (mrna_id, starts_exon))
        data = cur.fetchall()
        assert len(data) == 1, 'More than one record returned'
        exon_id = data[0][0]
        values = tuple([exon_id, mrna_id] + sline)
        #pdb.set_trace()
        cur.execute('''INSERT INTO fuguMatch VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', values)
    cur.close()
    conn.close()


if __name__ == '__main__':
    main()