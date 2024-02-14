from collections import defaultdict 


gid2gene_name = {}
with open('braker_lsat_HC_gene_tagged_fixed-rc-ser.gff3') as inp:
    for line in inp:
        if line.startswith('#'):
           continue

        A = line.strip().split('\t') #Lschr1  AUGUSTUS        gene    124107  141508  .       +       .       ID=g5952;Name=Lsat_1G0000000100;conf=HC
        if A[2] == 'gene':
           id_tag, name_tag  = A[-1].split(';')[:2]
           assert 'ID=' in  id_tag
           assert 'Name=' in  name_tag

           gene_id = id_tag[3:]
           gene_name = name_tag[5:]
           if '.' in gene_id:
              gene_id = gene_id.split('.')[0]

           gid2gene_name[gene_id] = gene_name  # g5952 : Lsat_1G0000000100


#gid2gene_name['g4957']= gid2gene_name['g4956']
#gid2gene_name['g4751']= gid2gene_name['g4750']
#gid2gene_name['g4548']= gid2gene_name['g4547']


#print len(gid2gene_name) # 31719



with open('gid2gene_name.csv', 'w') as outf:
   for gid, name in gid2gene_name.items():
       outf.write(','.join([gid, name]) + '\n')  


head = ['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs']

not_found = set()
with open('queries.fasta.tsv') as inp, open('lathyrus_interproscan.csv', 'w') as outf:
    for line in inp:
        A = line.strip().split('\t') 
        
        isoform =  A[0]
        assert '.' in isoform
        gid, iso = isoform.split('.')  
        #print isoform
        if gid not in gid2gene_name: # 89 gene_ids
           not_found.add(gid)
           gid2gene_name[gid] = gid2gene_name['g' + str(int(gid[1:]) -1) ]  # must be .m1, .m2, .m3, ...,  records
        
        lathy_isoform_name = gid2gene_name[gid] + '.' + iso[1:]  
        outf.write(','.join([lathy_isoform_name] + A)  + '\n' )

print 'not_found = ', len(not_found) #   77
print 'not_found_m = ',not_found #   77

"""
not_found =  77
not_found_m =  set(['g21522', 'g1494', 'g3206', 'g9884', 'g26149', 'g15723', 'g8038', 'g23632', 'g3049', 'g23170', 'g29050', 'g836', 'g4751', 'g24043', 'g26148', 'g15847', 'g22854', 'g13073', 'g23075', 'g410', 'g4957', 'g13313', 'g6059', 'g17580', 'g25456', 'g17921', 'g25760', 'g30440', 'g31629', 'g11194', 'g12034', 'g31607', 'g4979', 'g31482', 'g6227', 'g22038', 'g2312', 'g12548', 'g1587', 'g13780', 'g10256', 'g9292', 'g27195', 'g24489', 'g28066', 'g7517', 'g12036', 'g7537', 'g23081', 'g25530', 'g25698', 'g31206', 'g30126', 'g1013', 'g19616', 'g24547', 'g11796', 'g31484', 'g9887', 'g27157', 'g16918', 'g14860', 'g30518', 'g10561', 'g29391', 'g8109', 'g1438', 'g21930', 'g29704', 'g7308', 'g29699', 'g20288', 'g29072', 'g9133', 'g1662', 'g15542', 'g4548'])

"""
          
 
