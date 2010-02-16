#!/usr/bin/env python
# encoding: utf-8
"""
splitStart.py

Created by Brant Faircloth on 2010-02-05.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

Run with:  python splitStart.py -c db.conf

Purpose:  to parse the tStart and blockSize columns of a UCSC EST or spliceEST
table out into their own table with a many to one reference back to the `id`
of the original EST table (which we add through autoincrement).  Requires a
db.conf containing (to keep these off the repos.):

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

    (options,arg) = p.parse_args()
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg

def createStartsTable(cur):
    '''create the starts table which will hold the split data from the tStarts
    column of the UCSC EST sql dump.  The UCSC est_all or intronEst table 
    should have an added id auto_increment column added as a P_KEY'''
    try:
        cur.execute('''DROP TABLE starts''')
    except:
        pass
    cur.execute('''
        CREATE TABLE `starts` (
        `id` int(10) unsigned NOT NULL auto_increment,
        `splice_id` int(10) unsigned NOT NULL,
        `exon` int(10) unsigned NOT NULL,
        `start` int(10) unsigned NOT NULL default '0',
        `length` mediumint(8) unsigned NOT NULL default '0', 
        `end` int(10) unsigned NOT NULL default '0',
        PRIMARY KEY (`id`),
        KEY `splice_id` (`splice_id`),
        KEY `exon` (`exon`),
        KEY `start` (`start`),
        KEY `end` (`end`),
        FOREIGN KEY (`splice_id`) REFERENCES `splice` (`id`)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8 
    ''')

def enterStartsData(cur, splice_id, counter, start, length, end):
    '''enter the data from the split qStart and block columns to the Starts
    table'''
    try:
        cur.execute('''INSERT INTO starts (splice_id, exon, 
        start,length, end) values (?,?,?,?,?)''', (splice_id, counter, start, \
        length, end))
    except:
        pdb.set_trace()

def getAndSplitStarts(cur):
    '''get est data from UCSC table, split the data on the commas in the 
    tStarts column and the block sizes (length) columns, determine the end
    position and enter all that into the starts table, referencing everything
    to each records (in splice table) P_KEY'''
    cur.execute('''SELECT id, blockSizes, tStarts from splice''')
    rows = cur.fetchall()
    for each in rows:
        id, blocks, tstarts = each
        # split the data for each block
        blocks = blocks.split(',')
        tstarts = tstarts.split(',')
        for l in (blocks,tstarts):
            if '' in l:
                l.remove('')
        counter = 0
        for pos, bval in enumerate(blocks):
            #pdb.set_trace()
            qend =  int(bval) + int(tstarts[pos])
            enterStartsData(cur, id, counter, int(tstarts[pos]), int(bval), qend)
            counter += 1


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
    createStartsTable(cur)
    getAndSplitStarts(cur)
    conn.commit()
    cur.close()
    conn.close()


if __name__ == '__main__':
    main()