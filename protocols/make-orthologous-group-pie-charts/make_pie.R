# TODO: merge into find-annotation-orthologs pipeline.

# Display C. elegans orthologous groups
t = 26312
o = 17191
slices = c(o, t - o)                                                                                                  
pct = round(slices/sum(slices)*100, 1)                                                                                      
lbls = c('Genes in orthologous groups', 'Genes not in orthologous groups')
lbls = paste(lbls, ' (', pct, '%)', sep='')                                                                                 
svg('orthos_celegans.svg')                                                                                                  
pie(slices, labels=lbls, main='C. elegans genes in orthologous groups')                                                       
dev.off()

# Display C. briggsae orthologous groups
t = 21842
o = 13805
slices = c(o, t - o)                                                                                                  
pct = round(slices/sum(slices)*100, 1)                                                                                      
lbls = c('Genes in orthologous groups', 'Genes not in orthologous groups')
lbls = paste(lbls, ' (', pct, '%)', sep='')                                                                                 
svg('orthos_cbriggsae.svg')                                                                                                  
pie(slices, labels=lbls, main='C. briggsae genes in orthologous groups')                                                       
dev.off()

# Display MHco3 orthologous groups (BLOSUM80, C. elegans outgroup)
t = 21897
o = 8813
slices = c(o, t - o)                                                                                                  
pct = round(slices/sum(slices)*100, 1)                                                                                      
lbls = c('Genes in orthologous groups', 'Genes not in orthologous groups')
lbls = paste(lbls, ' (', pct, '%)', sep='')                                                                                 
svg('orthos_mhco3.svg')                                                                                                  
pie(slices, labels=lbls, main='MHco3 genes in orthologous groups')                                                       
dev.off()

# Display McMaster orthologous groups (BLOSUM80, C. elegans outgroup)
t = 23610
o = 7101
slices = c(o, t - o)                                                                                                  
pct = round(slices/sum(slices)*100, 1)                                                                                      
lbls = c('Genes in orthologous groups', 'Genes not in orthologous groups')
lbls = paste(lbls, ' (', pct, '%)', sep='')                                                                                 
svg('orthos_mcmaster.svg')                                                                                                  
pie(slices, labels=lbls, main='McMaster genes in orthologous groups')                                                       
dev.off()

# Display cardinality relationships of orthologous groups
slices = c(4494, 407, 1951, 181)
pct = round(slices/sum(slices)*100, 1)                                                                                      
lbls = c('1:1', '1:n', 'n:1', 'm:n')
lbls = paste(lbls, ' (', pct, '%)', sep='')                                                                                 
svg('ortho_distribution.svg')
pie(slices, labels=lbls, main='MHco3:McMaster orthologous groups')
dev.off()

# Display tandem distribution of genes in orthologous groups
slices = c(21, 18, 71)
pct = round(slices/sum(slices)*100, 1)                                                                                      
lbls = c('No other genes on scaffold', 'Overlapping gene(s) on same scaffold', 'Non-overlapping gene(s) on same scaffold')
lbls = paste(lbls, ' (', slices, ' scaffolds)', sep='')                                                                                 
svg('ortho_tandem.svg')
pie(slices, labels=lbls, main='Structure of MHco3 orthologous groups')
dev.off()

# Display multi-orthologue scaffolds
slices = c(110, 3)
pct = round(slices/sum(slices)*100, 1)                                                                                      
lbls = c('MHco3', 'McMaster')
lbls = paste(lbls, ' (', slices, ' scaffolds)', sep='')                                                                                 
svg('ortho_multi.svg')
pie(slices, labels=lbls, main='Number of scaffolds with multiple members of single orthologous group')
dev.off()
