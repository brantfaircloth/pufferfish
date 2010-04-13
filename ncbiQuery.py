#!/usr/bin/env python
# encoding: utf-8
"""
ncbiQuery.py

Created by Brant Faircloth on 2010-03-22.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

This program goes out to NCBI and gets the GO data (and some other data) for 
genes in our organism of choice.  It uses the accession numbers 
(annotation.accession) of reads to determine the geneId of a region 
(from ncbi nucleotide), then goes to ncbi gene to get all the other annotation 
goodies.  This program is a pain because ncbi allows you only to run it btw. 
9:00 PM and 5:00 AM EST/EDT and they shut down their servers occasionally 
at that time, which kills the script.
"""

import pdb
import os
import sys
import time
import oursql
import optparse
import ConfigParser
from Bio import Entrez

def interface():
    '''Get the starting parameters from a configuration file'''
    usage = "usage: %prog [options]"
    
    p = optparse.OptionParser(usage)
    
    p.add_option('--configuration', dest='conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
metavar='FILE')
    p.add_option('--restart', action='store_true', dest='restart', \
default=False, help='Restart the run (do not delete tables)')
    
    (options,arg) = p.parse_args()
    
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg

def searchNcbiNucleotide(acc):
    source = None
    symbol = None
    name   = None
    geneid = None
    #pdb.set_trace()
    # get the entrez nucleotide record
    try:
        nucXml = Entrez.efetch(db='nucleotide', id=acc, retmode='xml')
    except:
        pdb.set_trace()
    nuc     = Entrez.read(nucXml)
    nucXml.close()
    # this is the critter the gene was located in
    #source  = nuc[0]['GBSeq_source']
    # iterate over values of nuc[0]['GBSeq_feature-table'][1]['GBFeature_quals']
    # to get what we want
    name = nuc[0]['GBSeq_definition']
    for item in nuc[0]['GBSeq_feature-table']:
        element = item['GBFeature_quals']
        for bit in element:
            if bit['GBQualifier_name'] == 'organism':
                source = bit['GBQualifier_value']
            elif bit['GBQualifier_name'] == 'gene':
                symbol = bit['GBQualifier_value']
            #elif bit['GBQualifier_name'] == 'note':
            #    name = bit['GBQualifier_value']
            elif bit['GBQualifier_name'] == 'db_xref':
                if bit['GBQualifier_value'].split(':')[0] == 'GeneID':
                    geneid = bit['GBQualifier_value'].split(':')[1]
                else:
                    pass
    #pdb.set_trace()
    return source, symbol, name, geneid

def getNcbiGeneGoData(geneid):
    try:
        geneXml = Entrez.efetch(db='gene', id=geneid, retmode='xml')
    except:
        pdb.set_trace()
    gene    = Entrez.read(geneXml)
    geneXml.close()
    try:
        goData  = gene[0]['Entrezgene_properties'][1]['Gene-commentary_comment']
    except KeyError:
        # if there is no go data
        goData  = None
    except IndexError:
        goData  = None
    #finally:
    #    pdb.set_trace()
    return goData


def goParser(cur, goData):
    # this will iterate over `Function`,`Process`,`Component`, etc.
    goNums = []
    for category in goData:
        catName = category['Gene-commentary_label']
        for record in category['Gene-commentary_comment']:
            #pdb.set_trace()
            goDesc  = record['Gene-commentary_source'][0]['Other-source_anchor']
            goNum   = record['Gene-commentary_source'][0]['Other-source_src']\
                            ['Dbtag']['Dbtag_tag']['Object-id']['Object-id_id']
            # keep the go num so we can return
            goNums.append(goNum)
            # check to see if in GO table
            cur.execute('''SELECT id from go where id = ?''', (goNum,))
            if not cur.fetchall():
                cur.execute('''INSERT INTO go (id, category, description)
                    values (?,?,?)''', (goNum, catName, goDesc))
    return goNums

def createGoTable(cur):
    try:
        cur.execute('''DROP TABLE gorefs''')
        cur.execute('''DROP TABLE go''')
        cur.execute('''DROP TABLE refseq''')
        cur.execute('''UPDATE annotation set refseq_id = 0''')
    except:
        pass
    cur.execute('''CREATE table go (id int unsigned not null,
        category varchar(30) not null, description text not null, primary key (id),
        index(description(15))) ENGINE=InnoDB charset=utf8''')
    cur.execute('''CREATE table refseq (id int unsigned not null,
        symbol varchar(30) not null, name text, source text, primary key (id),
        index(symbol)) ENGINE=InnoDB charset=utf8''')
    cur.execute('''CREATE TABLE gorefs (refseq_id int unsigned not null,
        go_id int unsigned,
        index(refseq_id),
        index(go_id),
        foreign key (refseq_id) references refseq (id),
        foreign key (go_id) references go (id)
        ) ENGINE=InnoDB charset=utf8''')


def refSeqInserter(cur, acc, nucleotideSearch):
    #pdb.set_trace()
    source, symbol, name, geneid = nucleotideSearch
    # check to see if geneid in refseq table
    cur.execute('''SELECT id from refseq where id = ?''', (geneid,))
    if not cur.fetchall():
        cur.execute('''INSERT INTO refseq (id, symbol, name, source)
                values (?,?,?,?)''', (geneid, symbol, name, source))

def goRefInserter(cur, refseq_id, goNums=False):
    if goNums:
        for num in goNums:
            cur.execute('''INSERT INTO gorefs (refseq_id, go_id) values (?,?)''', (refseq_id, num))
    else:
        cur.execute('''INSERT INTO gorefs (refseq_id, go_id) values (?,?)''', (refseq_id, None))

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
    #pdb.set_trace()
    if not options.restart:
        createGoTable(cur)
    cur.execute('''SELECT distinct(accession) FROM annotation where refseq_id = 0''')
    accessions = cur.fetchall()
    end = None
    print 'num lookups = ', len(accessions)
    start = time.time()
    count = 0
    for acc in accessions:
        # 3 requests/sec enforced by biopython
        nucleotideSearch = searchNcbiNucleotide(acc)
        refSeqInserter(cur, acc, nucleotideSearch)
        goData = getNcbiGeneGoData(nucleotideSearch[3])
        if goData:
            goNums = goParser(cur, goData)
            goRefInserter(cur, nucleotideSearch[3], goNums)
        else:
            goRefInserter(cur, nucleotideSearch[3])
        # update the annotation table with the gene for the read
        # this will update all the records, which is what we want
        # since we're working from a distinct list
        cur.execute('''UPDATE annotation SET refseq_id = ? where
            accession = ?''', (nucleotideSearch[3], acc[0]))
        #pdb.set_trace()
        if count%100 == 0:
            print '\tIteration = ', count
            temp = time.time() - start
            print '\t\tTime: ', temp/60., ' min'
        count += 1


if __name__ == '__main__':
    Entrez.email = "brant.faircloth@gmail.com"
    main()