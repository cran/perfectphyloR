## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----warning=FALSE-------------------------------------------------------
library(perfectphyloR)
 # Haplotype matrix
 haplo_mat <- matrix(c(1,1,1,0,
                       0,0,0,0,
                       1,1,1,1,
                       1,0,0,0), byrow = TRUE, ncol = 4)
 # SNV names
 SNV_names <- c(paste("SNV", 1:4, sep = ""))
 # Haplotype names
 hap_names <- c("h1", "h2", "h3", "h4")
 # SNV positions in base pairs
 SNV_posns <- c(1000, 2000, 3000, 4000)
 ex_hapMat <- createHapMat(hapmat = haplo_mat,
                           snvNames = SNV_names,
                           hapNames = hap_names ,
                           posns = SNV_posns)
 ex_hapMat

## ----message=FALSE, warning=FALSE----------------------------------------
# load the example hapMat data object

data(ex_hapMatSmall_data)

# Reconstruct first dendrogram
rdend <- reconstructPP(hapMat = ex_hapMatSmall_data,
                        focalSNV = 1,
                        minWindow = 1,
                        sep = "-")

# Structure of rdend

str(rdend)

## ------------------------------------------------------------------------
# Lower and upper limits of the window of SNVs in hapMat used to reconstruct 'rdend'. 
ex_hapMatSmall_data$posns[rdend$snvWinIndices]

## ------------------------------------------------------------------------
ex_hapMatSmall_data$hapmat[, rdend$snvWinIndices[1]:rdend$snvWinIndices[2]]

## ----message=FALSE, warning=FALSE----------------------------------------
data(ex_hapMatSmall_data)
                       
# Reconstruct partitions across the region. 
 rdends <- reconsPPregion(hapMat = ex_hapMatSmall_data,
                       minWindow = 1)
 
# Structure of the partition for the first SNV.
 str(rdends[[1]])
 
 

## ----fig1, fig.height=5, fig.width=8, fig.cap = 'Figure 1: Two partition structures reconstructed at the first SNV and the tenth SNV of `ex_hapMatSmall_data`. ', fig.show = 'hold'----

data(ex_hapMatSmall_data)
# Reconstruct first dendrogram
 rdend1 <- reconstructPP(hapMat = ex_hapMatSmall_data,
                        focalSNV = 1,
                        minWindow = 1,
                        sep = "-")

# Reconstruct the second dendrogram
 rdend2 <- reconstructPP(hapMat = ex_hapMatSmall_data,
                        focalSNV = 10,
                        minWindow = 1,
                        sep = "-")
 
 par(mfrow=c(1,2), las=1)
 plotDend(dend = rdend1, direction = "downwards")
 
 plotDend(dend = rdend2, direction = "downwards")
 

## ----eval= FALSE---------------------------------------------------------
#  
#  # Comparator true dendrogram at 975 kbp.
#  data(tdend)
#  
#  # hapMat data object.
#  data(ex_hapMat_data)
#  
#  # Reconstruct dendrograms across the region.
#  allrdends <- reconsPPregion(hapMat = ex_hapMat_data,
#                         minWindow = 1)
#  
#  # Rand index profile based on 6 clusters.
#  RI_profile <- testDendAssoRI(rdend = allrdends, cdend = tdend, k = 6,
#                               hapMat = ex_hapMat_data, nperm = 1000,
#                               xlab = "SNV positions (bp)",
#                               ylab = "Rand indices", main = "Association Profile")
#  
#  # Omnibus P value for overall association.
#  RI_profile$OmPval
#  

## ----echo=FALSE, message=FALSE, warning=FALSE, fig2, fig.height=5.5, fig.width=9, fig.cap = 'Figure 2: Rand index values between the comparator true dendrogram at SNV position 975 kbp and the reconstructed dendrogram at each SNV position across the 2 Mbp genomic region of the example dataset.',fig.show = 'hold'----

pval = 0.000999001

print(pval)
data(ex_hapMat_data)

load(url("https://github.com/cbhagya/perfectphyloR/raw/master/vignette_data/RI_profile.rda"))

plot(ex_hapMat_data$posns, RI_profile$Stats, main = "Association Profile"
     , xlab = "SNV positions (bp)", ylab = "Rand indices")



## ----eval= FALSE---------------------------------------------------------
#  
#  # Comparator true dendrogram at SNV position 975 kbp.
#  data(tdend)
#  
#  # hapMat data object.
#  data(ex_hapMat_data)
#  
#  
#  # Compute rank-based distances between haplotypes based on
#  # the comparator true dendrogram (tdend) using the function 'rdistMatrix()'.
#  tdendDmat = perfectphyloR::rdistMatrix(tdend)
#  
#  # dCor profile, comparing the association between distance matrix of true dendrogram
#  # (comparator dendrogram) and all reconstructed dendrogram across the genomic region.
#  #
#  dCor_profile <- testAssoDist(cdmat = tdendDmat, rdend = allrdends, method = "dCor",
#                               hapMat = ex_hapMat_data, nperm = 1000,
#                               xlab = "SNV positions (bp)",
#                               ylab = "dCor Statistics", main = "Association Profile")
#  
#  # Omnibus P value for overall association.
#  dCor_profile$OmPval

## ----echo=FALSE, fig3, fig.height=5.5, fig.width=9, fig.cap = 'Figure 3: dCor statistics between the comparator distance matrix at SNV position 975 kbp, and the reconstructed dendrograms at each SNV position across the 2 Mbp genomic region of the example dataset. ',fig.show = 'hold'----
print(pval)

load(url("https://github.com/cbhagya/perfectphyloR/raw/master/vignette_data/dCor_profile.rda"))

plot(ex_hapMat_data$posns, dCor_profile$Stats, main = "Association Profile"
     , xlab = "SNV positions (bp)", ylab = "dCor Statistics")

## ----eval = FALSE--------------------------------------------------------
#  
#  # Phenotypic distances
#  data(phenoDist)
#  
#  
#  # After that compute the RV profile.
#  RV_profile <- testAssoDist(cdmat = phenoDist, rdend = allrdends, method = "RV",
#                               hapMat = ex_hapMat_data, nperm = 1000,
#                               xlab = "SNV positions (bp)",
#                               ylab = "RV coeffients", main = "Association Profile")
#  
#  # Omnibus P value for overall association.
#  RV_profile$OmPval

## ----echo=FALSE, fig4, fig.height=5.5, fig.width=9, fig.cap = 'Figure 4: RV coefficients between the phenotypic distance matrix and the reconstructed dendrograms across the 2 Mbp genomic region of the example dataset.', fig.show = 'hold'----
pval2 =  0.4945255

print(pval2)

load(url("https://github.com/cbhagya/perfectphyloR/raw/master/vignette_data/RV_profile.rda"))

plot(ex_hapMat_data$posns, RV_profile$Stats, main = "Association Profile"
     , xlab = "SNV positions (bp)", ylab = "RV coeffients")


## ----eval = TRUE, fig5, fig.height=5.5, fig.width=9, fig.cap = 'Figure 5: RV coefficients between the phenotypic distance matrix and the reconstructed dendrograms across the 2 Mbp genomic region of the example dataset annotated with risk region.', fig.show = 'hold'----

# plot association profile. 
plot(ex_hapMat_data$posns, RV_profile$Stats, main = "Association profile with the risk region"
     , xlab = "SNV positions (bp)", ylab = "RV coefficients")
# Indicate the risk region where the causal SNVs are located.
abline(v=950000); abline(v=1050000)  


