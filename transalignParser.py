#!/usr/bin/env python
# encoding: utf-8
"""
transalignParser.py

Created by Brant Faircloth on 2010-03-22.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

Parse output from a bed file intersection between our putative gene regions
and those genes aligning from transalignRefSeq.
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
    
    p.add_option('--input', dest='input', action='store', \
type='string', default = None, help='The path to the input BED file.', \
metavar='FILE')
    
    p.add_option('--configuration', dest='conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
metavar='FILE')

    (options,arg) = p.parse_args()
    
    if not options.conf or not options.input:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg


def createAccessionTable(cur):
    try:
        cur.execute('''DROP TABLE annotation''')
    except:
        pass
    cur.execute('''CREATE table annotation (id int unsigned not null 
        auto_increment, gene_id int unsigned not null, long_accession varchar(30) 
        not null, accession varchar(20) not null, refseq_id int unsigned not null default 0,
        primary key (id),
        index(refseq_id),
        index(accession),
        foreign key (gene_id) references genes (id), 
        index(gene_id)) ENGINE=InnoDB charset=utf8''')


def getGeneOverlap(cur, chromo, start, stop):
    cur.execute('''SELECT id FROM genes WHERE chromo = ? AND 
        (((start between ? and ?) OR (end between ? and ?)) OR 
        ((? between start and end) OR (? between start and end)))''', \
        (chromo, start, stop, start, stop, start, stop))
    return cur.fetchall()

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
    createAccessionTable(cur)
    conn.commit()
    #counter = 0
    for line in open(options.input, 'rU').readlines():
        if '\t' in line:
            line = line.split('\t')
            chromo, start, stop = line[0:3]
            acc = line[3].split('.')[0]
            overlaps = getGeneOverlap(cur, chromo, start, stop)
            if overlaps:
                for gene in overlaps:
                    #pdb.set_trace()
                    cur.execute('''INSERT INTO annotation (gene_id, long_accession, accession) VALUES (?,?, ?)''', (gene[0], line[3], acc))
            else:
                pdb.set_trace()
        #counter += 1
        #print line
    


if __name__ == '__main__':
    main()