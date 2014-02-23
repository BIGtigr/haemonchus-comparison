library(Vennerable)

v <- Venn(
       SetNames = c('MHco3 complete', 'MHco3 partial', 'McMaster complete', 'McMaster partial'),
       Weight = c(
         '0110' = 2,
         '0111' = 0,
         '0000' = 0,
         '0001' = 5,
         '0011' = 0,
         '0010' = 5,
         '0101' = 8,
         '0100' = 6,
         '1111' = 0,
         '1110' = 0,
         '1100' = 0,
         '1101' = 0,
         '1010' = 166,
         '1011' = 0,
         '1001' = 33,
         '1000' = 17))

svg('mhco3_vs_mcmaster_cegs.svg')
plot(v, doWeights = FALSE, type='ellipses', show = list(SetLabels=TRUE))
dev.off()

v <- Venn(SetNames=c('MHco3', 'McMaster'), Weight=c('01'=7, '10'=50, '11'=166))
svg('present_complete.svg')
plot(v, doWeights = FALSE, type='circles', show = list(SetLabels=TRUE))
dev.off()

v <- Venn(SetNames=c('MHco3', 'McMaster'), Weight=c('01'=10, '10'=23, '11'=209))
svg('present_comp+part.svg')
plot(v, doWeights = FALSE, type='circles', show = list(SetLabels=TRUE))
dev.off()

v <- Venn(SetNames=c('MHco3', 'McMaster'), Weight=c('01'=75, '10'=32, '11'=25))
svg('missing_complete.svg')
plot(v, doWeights = FALSE, type='circles', show = list(SetLabels=TRUE))
dev.off()

v <- Venn(SetNames=c('MHco3', 'McMaster'), Weight=c('01'=29, '10'=16, '11'=6))
svg('missing_comp+part.svg')
plot(v, doWeights = FALSE, type='circles', show = list(SetLabels=TRUE))
dev.off()
