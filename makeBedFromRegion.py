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
    outp = open(options.out,'w')
    outp.write('''track name=geneRegions description="Fugu Exons" useScore=0\n''')
    cur.execute('''SELECT * FROM region''')
    data = cur.fetchall()
    for record in data:
        #pdb.set_trace()
        outp.write('chrUn %s %s region.id_%s\n' % (record[1], record[2], record[0]))
    outp.close()
    cur.close()
    conn.close()


if __name__ == '__main__':
    main()