select refseq.id, refseq.symbol, gorefs.go_id, go.description from refseq, gorefs, go where gorefs.refseq_id = refseq.id and go.id = gorefs.go_id and refseq.id = 229 and go.category = 'function';

# get count of reads by go terms
select category, description, count(refseq.id) from go, refseq, gorefs where category = 'Function' and gorefs.refseq_id = refseq.id and go.id = gorefs.go_id group by go.description order by count(refseq.id) desc;