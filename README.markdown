Title:      README for fugu prep  
Version:    0.3 (23 Mar 2010)  
Author:     Brant C. Faircloth  
Create:     22 Mar 2010  
Edit:       23 Mar 2010  

### Get data from UCSC

    wget http://hgdownload.cse.ucsc.edu/goldenPath/tetNig2/database/all_mrna.sql
    wget http://hgdownload.cse.ucsc.edu/goldenPath/tetNig2/database/all_mrna.txt.gz

### create table as in

    http://hgdownload.cse.ucsc.edu/goldenPath/tetNig2/database/all_mrna.sql

### insert data into table

    mv all_mrna.txt all_mrna.22-Mar-2010.txt  
    
    mysql> load data infile '/Users/bcf/git/fugu/data/tetNig2/all_mrna.22-Mar-2010.txt' \
    into table x_mrna_all fields terminated by '\t';

### Use `show create table est` to show the table structure for the est data modify this create table statement to add an autoincrementing primary key we need this key later to maintain referential integrity

    CREATE TABLE `mrna` (
    `id` int unsigned NOT NULL auto_increment,
    `bin` smallint(5) unsigned NOT NULL default '0',
    `matches` int(10) unsigned NOT NULL default '0',
    `misMatches` int(10) unsigned NOT NULL default '0',
    `repMatches` int(10) unsigned NOT NULL default '0',
    `nCount` int(10) unsigned NOT NULL default '0',
    `qNumInsert` int(10) unsigned NOT NULL default '0',
    `qBaseInsert` int(10) unsigned NOT NULL default '0',
    `tNumInsert` int(10) unsigned NOT NULL default '0',
    `tBaseInsert` int(10) unsigned NOT NULL default '0',
    `strand` char(2) NOT NULL default '',
    `qName` varchar(255) NOT NULL default '',
    `qSize` int(10) unsigned NOT NULL default '0',
    `qStart` int(10) unsigned NOT NULL default '0',
    `qEnd` int(10) unsigned NOT NULL default '0',
    `tName` varchar(255) NOT NULL default '',
    `tSize` int(10) unsigned NOT NULL default '0',
    `tStart` int(10) unsigned NOT NULL default '0',
    `tEnd` int(10) unsigned NOT NULL default '0',
    `blockCount` int(10) unsigned NOT NULL default '0',
    `blockSizes` longblob NOT NULL,
    `qStarts` longblob NOT NULL,
    `tStarts` longblob NOT NULL,
    PRIMARY KEY (id),
    KEY `tName` (`tName`(8),`bin`),
    KEY `tName_2` (`tName`(8),`tStart`),
    KEY `qName` (`qName`(12)),
    KEY `tName_3` (`tName`(8),`tEnd`)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8;

### insert the data from the UCSC table

    insert into mrna (bin, matches, misMatches, repMatches, nCount,
    qNumInsert, qBaseInsert, tNumInsert, tBaseInsert, strand, qName, qSize,
    qStart, qEnd, tName, tSize, tStart, tEnd, blockCount, blockSizes, qStarts, 
    tStarts) select * from x_mrna_all;

### since these tables are duplicates, drop the original table without the primary key

    drop table x_mrna_all;

### parse the blockSize and tStart columns out into their own rows using:

    python splitStart.py --configuration=data/tetNig2/db.conf \
    --tablename=starts --reftablename=mrna

### using the data in the newly created `starts` table (from splitStart.py above), go through regions, compact that all done into a single `region` (composed of elements in `region_components` and record the area covered by the `region`

    python splitCluster.py --configuration=data/tetNig2/db.conf --distance=0 \
    --minoverlaps=0 --tablename=starts --reftablename=mrna

### rename files in mysql because it's easier than changing the program

    rename table region to exons;
    rename table region_components to exon_components;

### setup some foreign keys if not exists

    alter table exons add foreign key (mrna_id) references mrna (id);

### create the genes table

    create table genes (id int unsigned not null auto_increment, mrna_id int    
    unsigned not null, chromo varchar(15) default NULL, start int unsigned not 
    null, end int unsigned not null, primary key (id), index(mrna_id), foreign 
    key (mrna_id) references mrna(id)) ENGINE=InnoDB charset=utf8;
    insert into genes (mrna_id) select distinct(mrna_id) from exons;

### add foreign key to genes if not exists

    alter table genes add foreign key (mrna_id) references mrna (id);

### produce a bed file of the gene-region

    python makeBedFromRegion.py --configuration=data/tetNig2/db.conf \
    --out=data/tetNig2/tetNig2_putativeGenes.bed --advanced

### get the intersection between putative genes and the transalign refseq

    (do this via UCSC; overlap >= 80%)

### run the geneMinMax.py script to get ranges of putative gene regions

    python geneMinMax.py --configuration=data/tetNig2/db.conf

### parse the intersection data:

    python transalignParser.py --configuration=data/tetNig2/db.conf \
    --input=data/tetNig2/tetNig2PutativeGeneIntersecetTransmapReverse.bed

### get the annotation data from NCBI

    python ncbiQuery.py --configuration=data/tetNig2/db.conf

### run a per-species query to get go-term counts on a per-species basis

#### Humans (Function)

    select category, description, count(refseq.id) from go, refseq, gorefs 
    where category = 'Function' and refseq.source = 'Homo sapiens' and 
    gorefs.refseq_id = refseq.id and go.id = gorefs.go_id group by 
    go.description order by count(refseq.id) desc into outfile 
    '~/Database/tetNig2/tetNigHomoSapGo.Function.Terms.txt' fields terminated by '\t' ;
    
#### Humans (Process)

    select category, description, count(refseq.id) from go, refseq, gorefs 
    where category = 'Process' and refseq.source = 'Homo sapiens' and 
    gorefs.refseq_id = refseq.id and go.id = gorefs.go_id group by 
    go.description order by count(refseq.id) desc into outfile 
    '~/Database/tetNig2/tetNigHomoSapGo.Process.Terms.txt' fields terminated by '\t' ;

#### Mus musculus (Function)

    select category, description, count(refseq.id) from go, refseq, gorefs 
    where category = 'Function' and refseq.source = 'Mus musculus' and 
    gorefs.refseq_id = refseq.id and go.id = gorefs.go_id group by 
    go.description order by count(refseq.id) desc into outfile 
    '~/Database/tetNig2/tetNigMusMus.Function.Terms.txt' fields terminated by '\t' ;
    
#### Mus musculus (Process)

    select category, description, count(refseq.id) from go, refseq, gorefs 
    where category = 'Process' and refseq.source = 'Mus musculus' and 
    gorefs.refseq_id = refseq.id and go.id = gorefs.go_id group by 
    go.description order by count(refseq.id) desc into outfile 
    '~/Database/tetNig2/tetNigMusMus.Process.Terms.txt' fields terminated by '\t' ;

#### Danio rerio (Function)

    select category, description, count(refseq.id) from go, refseq, gorefs 
    where category = 'Function' and refseq.source = 'Danio rerio' and 
    gorefs.refseq_id = refseq.id and go.id = gorefs.go_id group by 
    go.description order by count(refseq.id) desc into outfile 
    '~/Database/tetNig2/tetNigDanRerGo.Function.Terms.txt' fields terminated by '\t' ;
    
#### Danio rerio (Process)

    select category, description, count(refseq.id) from go, refseq, gorefs 
    where category = 'Process' and refseq.source = 'Danio rerio' and 
    gorefs.refseq_id = refseq.id and go.id = gorefs.go_id group by 
    go.description order by count(refseq.id) desc into outfile 
    '~/Database/tetNig2/tetNigDanRerGo.Process.Terms.txt' fields terminated by '\t' ;

#### Oryzias latipes (Function)

    (No hits)
    
#### Oryzias latipes (Process)

    (No hits)