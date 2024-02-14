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
with open('out.emapper.annotations') as inp, open('lathy.emapper.annotations.tab', 'w') as outf:
    outf.write('\t'.join(['lathy_transcript_name'] + head)  + '\n' )
    for line in inp:
        if line.startswith('#'):
           continue
        A = line.strip().split('\t') 
        
        isoform =  A[0]
        assert '.' in isoform
        gid, iso = isoform.split('.')  
        #print isoform
        if gid not in gid2gene_name: # 89 gene_ids
           not_found.add(gid)
           gid2gene_name[gid] = gid2gene_name['g' + str(int(gid[1:]) -1) ]  # must be .m1, .m2, .m3, ...,  records
        
        lathy_isoform_name = gid2gene_name[gid] + '.' + iso[1:]  
        outf.write('\t'.join([lathy_isoform_name] + A)  + '\n' )

print 'not_found_m = ', len(not_found) # 89

"""
not_found_m =  set(['g11031', 'g29050', 'g4751', 'g24434', 'g8038', 'g410', 'g31629', 'g6227', 'g9292', 'g25698', 'g17313', 'g17312', 'g9884', 'g20782', 'g8109', 'g1438', 'g29699', 'g23632', 'g9133', 'g12034', 'g12036', 'g31482', 'g22038', 'g2312', 'g31484', 'g1587', 'g13780', 'g28066', 'g26148', 'g26149', 'g23081', 'g29391', 'g24547', 'g27157', 'g21930', 'g7308', 'g29076', 'g11194', 'g29072', 'g18599', 'g4548', 'g3206', 'g3049', 'g836', 'g21522', 'g4957', 'g13313', 'g17580', 'g25456', 'g19616', 'g25760', 'g12548', 'g9887', 'g1662', 'g10256', 'g2181', 'g17921', 'g7517', 'g11477', 'g31206', 'g23789', 'g14860', 'g10561', 'g25530', 'g25102', 'g15542', 'g15723', 'g23170', 'g22565', 'g24043', 'g15847', 'g13073', 'g23075', 'g6059', 'g28937', 'g30440', 'g31607', 'g29704', 'g2585', 'g1013', 'g22854', 'g27195', 'g24489', 'g7537', 'g11796', 'g26103', 'g16918', 'g30518', 'g20288'])

"""
          
 
