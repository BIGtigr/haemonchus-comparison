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
