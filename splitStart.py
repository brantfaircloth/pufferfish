#!/usr/bin/env python
# encoding: utf-8
"""
splitStart.py

Created by Brant Faircloth on 2010-02-05.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

Run with:  python splitStart.py -c db.conf

Purpose:  to parse the tStart and blockSize columns of a UCSC EST/mRNA
table out into their own table with a many to one reference back to the `id`
of the original EST/mRNA table (which we add through autoincrement).  

Requires a db.conf containing (to keep these off the repos.):

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
    p.add_option('--tablename', '-t', dest = 'table', action='store', \
type='string', default = None, help='The name of the output table', \
metavar='FILE')
    p.add_option('--reftablename', '-r', dest = 'ref', action='store', \
type='string', default = None, help='The name of the table to reference', \
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

def createStartsTable(cur, table, ref):
    '''create the starts table which will hold the split data from the tStarts
    column of the UCSC EST/mRNA sql dump.  We need to add to the UCSC table 
    added `id` auto_increment column added as a P_KEY and the ENGINE set to
    InnoDB'''
    try:
        query = '''DROP TABLE %s''' % table
        cur.execute(query)
    except:
        pass
    query ='''
        CREATE TABLE `%s` (
        `id` int(10) unsigned NOT NULL auto_increment,
        `%s_id` int(10) unsigned NOT NULL,
        `exon` int(10) unsigned NOT NULL,
        `chromo` varchar(15) NOT NULL default '',
        `start` int(10) unsigned NOT NULL default '0',
        `length` mediumint(8) unsigned NOT NULL default '0', 
        `end` int(10) unsigned NOT NULL default '0',
        PRIMARY KEY (`id`),
        KEY `%s_id` (`%s_id`),
        KEY `exon` (`exon`),
        KEY `chromo` (`chromo`),
        KEY `start` (`start`),
        KEY `end` (`end`),
        FOREIGN KEY (`%s_id`) REFERENCES `%s` (`id`)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8 
    ''' % (table, ref, ref, ref, ref, ref)
    cur.execute(query)

def enterStartsData(cur, table, ref, splice_id, chromo, counter, start, length, end):
    '''enter the data from the split qStart and block columns to the new
    table'''
    try:
        query = '''INSERT INTO %s (%s_id, exon, chromo, start, length, end) 
        values (%s,%s,'%s',%s,%s,%s)''' % (table, ref, splice_id, \
        counter, chromo, start, length, end)
        cur.execute(query)
    except:
        pdb.set_trace()

def getAndSplitStarts(cur, table, ref):
    '''get data from ref table, split the data on the commas in the 
    tStarts column and the block sizes (length) columns, determine the end
    position and enter all that into the new table, referencing everything
    to each records (in ref table) P_KEY'''
    query = '''SELECT id, tname, blockSizes, tStarts from %s''' % ref
    cur.execute(query)
    rows = cur.fetchall()
    for each in rows:
        id, chromo, blocks, tstarts = each
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
            enterStartsData(cur, table, ref, id, chromo, counter, int(tstarts[pos]), int(bval), qend)
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
    createStartsTable(cur, options.table, options.ref)
    getAndSplitStarts(cur, options.table, options.ref)
    conn.commit()
    cur.close()
    conn.close()


if __name__ == '__main__':
    main()