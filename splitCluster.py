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
    p.add_option('--tablename', '-t', dest = 'table', action='store', \
type='string', default = None, help='The name of the output table')   
    p.add_option('--reftablename', '-r', dest = 'ref', action='store', \
type='string', default = None, help='The name of the table to reference')
    p.add_option('--minoverlaps', '-m', dest = 'min', action='store', \
type='int', default = 3, help='The minimum number of overlaps we require')
    (options,arg) = p.parse_args()
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg
    


def clusterer(used, known, cur, table, chromo, mn, mx, dist):
    #pdb.set_trace()
    # find all overlapping
    query = '''SELECT id, mrna_id, exon, chromo, start, end from %s where 
        chromo = '%s' AND ((%s > start) AND (%s < end))''' % \
        (table, chromo, mx+dist, mn-dist)
    cur.execute(query)
    results = cur.fetchall()
    if results:
        overlapping = set(results)
        # are there any new overlapping
        new = overlapping.difference(known)
        if new:
            # do any remain after checking for used
            new = new.difference(used)
            known.update(new)
            if new:
                # make them used
                used.update(new)
                # get the min start of the cluster
                tmn = min([n[4] for n in new])
                if tmn < mn:
                    mn = tmn
                # get the max start of the cluster
                tmx = max([n[5] for n in new])
                if tmx > mx:
                    mx = tmx
                used, known, mn, mx = clusterer(used, known, cur, table, chromo, mn, mx, dist)
    return used, known, mn, mx 

def checkClusters(cur, ref, table, used, each, dist, index, min_overlaps, chromosome):
    ''' Stuff '''
    id, splice_id, exon, chromo, start, end = each
    if not chromosome:
        chromosome = chromo
    elif chromosome != chromo:
        # reset used
        used = set()
        chromosome = chromo
        print "Chromosome = ", chromosome
    known = set()
    # add each to the set
    known.add(each)
    # add record to used
    used.add(each)
    used, known, mn, mx = clusterer(used, known, cur, table, chromo, start, end, dist)
    cur.execute('''INSERT INTO region values (?, ?, ?, ?, ?, ?)''', (index, splice_id, exon, chromo, mn, mx))
    # update the region components table
    for elem in known:
        cur.execute('''INSERT INTO region_components values (?, ?)''', (index, elem[0]))
    return used, chromosome

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
        `mrna_id` int(10) unsigned NOT NULL,
        `starts_exon` int(10) unsigned NOT NULL,
        `chromo` varchar(15) NOT NULL default '',
        `start` int(10) unsigned NOT NULL default '0',
        `end` int(10) unsigned NOT NULL default '0',
        PRIMARY KEY  (`id`),
        KEY `mrna_id` (`mrna_id`),
        KEY `starts_exon` (`starts_exon`)
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
    query = '''SELECT id, %s_id, exon, chromo, start, end from %s''' % (options.ref, options.table)
    cur.execute(query)
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
    chromo = None
    for each in rows:
        #pdb.set_trace()
        if each in used:
            pass
        else:
            used, chromo = checkClusters(cur, options.ref, options.table, used, each, options.dist, i, options.min, chromo)
        i += 1
        #if i == 20:
        #    pdb.set_trace()
        if i%1000 == 0:
            print ('Rowcount = %s, Position = %s' % (len(rows), i))
    conn.commit()
    cur.close()
    conn.close()


if __name__ == '__main__':
    main() 