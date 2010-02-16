#!/usr/bin/env python
# encoding: utf-8
"""
splitCluster.py

Created by Brant Faircloth on 2010-02-05.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

Run with:  python splitCluser.py -c db.conf -d distance

Purpose:  ... 
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
    p.add_option('--distance', '-d', dest = 'dist', action='store', \
type='int', default = 0, help='The distance within which we want to search')

    (options,arg) = p.parse_args()
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg
    

def checkClusters(cur, used, each, dist, index, min_overlaps = 3):
    ''' Stuff '''
    id, splice_id, exon, start, end = each
    # this should give all regions with any sort of overlap - not including
    # abutting fragments, but that's okay
    cur.execute('''SELECT id, splice_id, exon, start, end from starts where 
        splice_id != ? and ((? > start) AND (? < end))''', \
        (splice_id, end+dist, start-dist))
    rows = cur.fetchall()
    if rows and len(rows) >= min_overlaps:
        temp = set()
        # add main record to the cluster
        rows.append(each)
        #pdb.set_trace()
        temp.update(rows)
        # check to see if there are any common elements between sets
        if len(used & temp) == 0:
            #pdb.set_trace()
            print ('%s,%s') % (splice_id, exon), rows
            # get the min start of the cluster
            mn = min([r[3] for r in rows])
            # get the max start of the cluster
            mx = max([r[4] for r in rows])
            # record the ids of fragments used so as not to dupe
            used.update(rows)
            # update the region table with pertinent starts and ends
            cur.execute('''INSERT INTO region values (?, ?, ?)''', (index, mn, mx))
            # update the region components table
            for elem in rows:
                cur.execute('''INSERT INTO region_components values (?, ?)''', (index, elem[0]))
            #pdb.set_trace()
    return used

def createRegionTable(cur):
    '''create the region table which will hold the information on DNA 
    `regions` of interest.  DNA regions are composed of exons from ESTs that 
    overlap at least 3 other exon regions from ESTs'''
    try:
        cur.execute('''DROP TABLE region_components''')
        cur.execute('''DROP TABLE region''')
    except:
        pass
    cur.execute('''
        CREATE TABLE `region` (
        `id` int(10) unsigned NOT NULL,
        `start` int(10) unsigned NOT NULL default '0',
        `end` int(10) unsigned NOT NULL default '0',
        PRIMARY KEY  (`id`)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8 
        ''')


def createRegionComponentsTable(cur):
    '''create a table holding the individual, exon components used to build
    each `region`.  This table essentially ties the region data to the data
    in the `splice` and `starts` tables.'''
    try:
        cur.execute('''DROP TABLE region_components''')
    except:
        pass
    cur.execute('''
        CREATE TABLE `region_components` (
        `region_id` int(10) unsigned NOT NULL,
        `starts_id` int(10) unsigned NOT NULL,
        KEY `region_id` (`region_id`),
        KEY `starts_id` (`starts_id`),
        FOREIGN KEY (`region_id`) REFERENCES `region` (`id`),
        FOREIGN KEY (`starts_id`) REFERENCES `starts` (`id`)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8''')

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
    #cur.execute('''SELECT splice_id, id, start, end from starts''')
    cur.execute('''SELECT id, splice_id, exon, start, end from starts''')
    rows = cur.fetchall()
    # we're going to use a set to hold our introns that are already used.
    # this saves us from having to create a table to hold these data and
    # then run a query against that table, matching on several fields.
    used = set()
    i = 0
    # create region table
    createRegionTable(cur)
    # create region_components table
    createRegionComponentsTable(cur)
    for each in rows:
        used = checkClusters(cur, used, each, options.dist, i)
        i += 1
        #if i == 20:
        #    pdb.set_trace()
    conn.commit()
    cur.close()
    conn.close()


if __name__ == '__main__':
    main() 