select refseq.id, refseq.symbol, gorefs.go_id, go.description from refseq, gorefs, go where gorefs.refseq_id = refseq.id and go.id = gorefs.go_id and refseq.id = 229 and go.category = 'function';

# get count of reads by go terms
select category, description, count(refseq.id) from go, refseq, gorefs where category = 'Process' and refseq.source = 'Homo sapiens' and gorefs.refseq_id = refseq.id and go.id = gorefs.go_id group by go.description order by count(refseq.id) desc;

# distinct genes represented in matches of tetraodon to fugu (may not be complete matches - may only be partial)
select count(distinct(genes.id)) from genes, fugumatch where genes.mrna_id = fugumatch.mrna_id limit 20;

#

select fugugenes.genes_id, category, description from annotation, go, refseq, gorefs, fugugenes where category = 'Function' and refseq.source = 'Homo sapiens' and gorefs.refseq_id = annotation.id and go.id = gorefs.go_id and fugugenes.genes_id = annotation.gene_id limit 10;