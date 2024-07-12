# Libraries
library(pegas)      # For population genetics analyses
library(adegenet)   # For genetic data analysis
library(vcfR)       # For reading VCF files
library(tidyverse)  # For data manipulation and visualization
library(ape)        # For phylogenetic and evolutionary analyses
library(ade4)       # For multivariate data analysis
library(hierfstat)  # For hierarchical F-statistics
library(poppr)      # For population genetic analyses
library(dartR)      # For SNP data analysis
library(dplyr)      # For data manipulation
library(reshape)    # For data reshaping
library(reshape2)   # For data reshaping

# Read in previously filtered VCF file
no626 <- read.vcfR("VCF6_26.vcf") # Read in VCF file
no626

# Convert to genlight object
mygenlight2 <- vcfR2genlight(no626) # Convert VCF file to genlight object
mygenlight2

# Add population info to genlight object
pop(FF626GL) <- c("Belize", "Belize", "Belize", "Belize", "Belize", "Belize", "Belize", "Belize", "Belize", "Belize",
                  "Belize", "Belize", "Belize", "Belize", "Belize", "Belize", "Belize",
                  "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida",
                  "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida",
                  "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida",
                  "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida",
                  "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida", "Florida",
                  "North Carolina", "North Carolina", "North Carolina", "North Carolina", "North Carolina", "North Carolina", "North Carolina",
                  "Texas", "Texas", "Texas", "Texas", "Texas", "Texas", "Texas")

# Convert into genind objects

# Check population and individual names 
pop(mygenlight2)
indNames(mygenlight2)
pop(mygenlight2)

# Subset Populations into individual genlight objects
Pop1GL <- mygenlight2[pop == "pop1"]
pop(Pop1GL)
Pop2GL <- mygenlight2[pop == "pop2"]
pop(Pop2GL)
Pop3GL <- mygenlight2[pop == "pop3"]
pop(Pop3GL)
Pop4GL <- mygenlight2[pop == "pop4"]
pop(Pop4GL)

# Convert into genind objects
Pop1GI <- dartR::gl2gi(Pop1GL)
Pop1GI
Pop2GI <- dartR::gl2gi(Pop2GL)
Pop2GI
Pop3GI <- dartR::gl2gi(Pop3GL)
Pop3GI
Pop4GI <- dartR::gl2gi(Pop4GL)
Pop4GI

# Run HWE test on each Population
my.hwt1 <- hw.test(Pop1GI, B = 1000) # B = number of draws. Low for testing. Use default for real, B = 1000
summary(my.hwt1)
my.hwt2 <- hw.test(Pop2GI, B = 1000) # B = number of draws. Low for testing. Use default for real, B = 1000
summary(my.hwt2)
my.hwt3 <- hw.test(Pop3GI, B = 1000) # B = number of draws. Low for testing. Use default for real, B = 1000
summary(my.hwt3)
my.hwt4 <- hw.test(Pop4GI, B = 1000) # B = number of draws. Low for testing. Use default for real, B = 1000
summary(my.hwt4)

# Convert genind to hierfstat to view Population 1 basic stats, combine with hw.test and view loci that have too high Ho/Hs ratio
Pop1HS <- genind2hierfstat(Pop1GI)
basicstatPop1 <- basic.stats(Pop1HS)
basicstatPop1
head(basicstatPop1$perloc)
statvaltabPop1 <- data.frame(basicstatPop1$perloc) # Pulls out per locus 
head(statvaltabPop1)
statvaltabPop1$locus <- rownames(statvaltabPop1)
my.hwt1 <- data.frame(my.hwt1)
my.hwt1$locus <- row.names(my.hwt1)
head(my.hwt1)
table(my.hwt1$locus == statvaltabPop1$locus) # Problem is not in table matching
statvaltabPop1 <- left_join(statvaltabPop1, my.hwt1, by = "locus") # Combines with HWE test results
head(statvaltabPop1)
# statvaltabPop1$filter <- ifelse(statvaltabPop1$Pr.exact < 0.01 & statvaltabPop1$Ho/statvaltabPop1$Hs > 1, 0, 1)
# Now using chi-squared probability rather than Pr. exact
statvaltabPop1$filter <- ifelse(statvaltabPop1$Pr.chi.2... < 0.01 & statvaltabPop1$Ho/statvaltabPop1$Hs > 1, 0, 1)
# Note that per ?basic.stats, Ho is observed heterozygosity, Hs is what people mean by expected
ggplot(statvaltabPop1, aes(x = Ho/Hs, fill = factor(filter))) +
  geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop1filtered, aes(x = Ho/Hs, fill = Pr.exact < 0.01)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop1, aes(x = Ho, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop1, aes(x = Hs, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop1, aes(x = Ho/Hs, fill = factor(filter))) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
statvaltabPop1filtered <- filter(statvaltabPop1, filter == 1)
statvaltabPop1filtered[statvaltabPop1filtered$Pr.chi.2... < 0.01,]

# Convert genind to hierfstat to view Population 2 basic stats, combine with hw.test and view loci that have too high Ho/Hs ratio
Pop2HS <- genind2hierfstat(Pop2GI)
basicstatPop2 <- basic.stats(Pop2HS)
basicstatPop2
head(basicstatPop2$perloc)
statvaltabPop2 <- data.frame(basicstatPop2$perloc) # Pulls out per locus 
statvaltabPop2$locus <- rownames(my.hwt2)
statvaltabPop2 <- bind_cols(statvaltabPop2, data.frame(my.hwt2)) # Combines with HWE test results
head(statvaltabPop2)
statvaltabPop2$filter <- ifelse(statvaltabPop2$Pr.chi.2... < 0.01 & statvaltabPop2$Ho/statvaltabPop2$Hs > 1, 0, 1)
# Note that per ?basic.stats, Ho is observed heterozygosity, Hs is what people mean by expected
# ggplot(statvaltabPop2, aes(x = Ho/Hs, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop2, aes(x = Ho, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop2, aes(x = Hs, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop2, aes(x = Ho/Hs, fill = factor(filter))) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
statvaltabPop2filtered <- filter(statvaltabPop2, filter == 0)

# Convert genind to hierfstat to view Population 3 basic stats, combine with hw.test and view loci that have too high Ho/Hs ratio
Pop3HS <- genind2hierfstat(Pop3GI)
basicstatPop3 <- basic.stats(Pop3HS)
basicstatPop3
head(basicstatPop3$perloc)
statvaltabPop3 <- data.frame(basicstatPop3$perloc) # Pulls out per locus 
statvaltabPop3$locus <- rownames(my.hwt3)
statvaltabPop3 <- bind_cols(statvaltabPop3, data.frame(my.hwt3)) # Combines with HWE test results
head(statvaltabPop3)
statvaltabPop3$filter <- ifelse(statvaltabPop3$Pr.chi.2... < 0.01 & statvaltabPop3$Ho/statvaltabPop3$Hs > 1, 0, 1)
# Note that per ?basic.stats, Ho is observed heterozygosity, Hs is what people mean by expected
# ggplot(statvaltabPop3, aes(x = Ho/Hs, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop3, aes(x = Ho, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop3, aes(x = Hs, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop3, aes(x = Ho/Hs, fill = factor(filter))) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")

# Convert genind to hierfstat to view Population 4 basic stats, combine with hw.test and view loci that have too high Ho/Hs ratio
Pop4HS <- genind2hierfstat(Pop4GI)
basicstatPop4 <- basic.stats(Pop4HS)
basicstatPop4
head(basicstatPop4$perloc)
statvaltabPop4 <- data.frame(basicstatPop4$perloc) # Pulls out per locus 
statvaltabPop4$locus <- rownames(my.hwt4)
statvaltabPop4 <- bind_cols(statvaltabPop4, data.frame(my.hwt4)) # Combines with HWE test results
head(statvaltabPop4)
statvaltabPop4$filter <- ifelse(statvaltabPop4$Pr.chi.2... < 0.01 & statvaltabPop4$Ho/statvaltabPop4$Hs > 1, 0, 1)
# Note that per ?basic.stats, Ho is observed heterozygosity, Hs is what people mean by expected
# ggplot(statvaltabPop4, aes(x = Ho/Hs, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop4, aes(x = Ho, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop4, aes(x = Hs, fill = Pr.exact < 0.05/2430)) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
# ggplot(statvaltabPop4, aes(x = Ho/Hs, fill = factor(filter))) +
#   geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")

# Check how many loci pass/fail the filter
table(statvaltabPop1$filter)
table(statvaltabPop2$filter)
table(statvaltabPop3$filter)
table(statvaltabPop4$filter)

# More don't pass under chi squared test

Pop1failedloci <- subset(statvaltabPop1, filter < 1) # Pull pop1 failed loci
dim(Pop1failedloci)

Pop2failedloci <- subset(statvaltabPop2, filter < 1) # Pull pop2 failed loci

Pop1filterlist <- as.list(Pop1failedloci$locus) # Make a list of failed loci for pop1

Pop2filterlist <- as.list(Pop2failedloci$locus) # Make a list of failed loci for pop2

Pop1keeploci <- subset(statvaltabPop1, filter > 0) # Pull pop1 loci to keep, 2053
dim(Pop1keeploci)
dim(Pop1failedloci)

Pop2keeploci <- subset(statvaltabPop2, filter > 0) # Pull pop2 loci to keep, 1964

Pop1keeplist <- as.list(Pop1keeploci$locus) # Make a list pop1 loci to keep
Pop2keeplist <- as.list(Pop2keeploci$locus) # Make a list of pop2 loci to keep

table(Pop1keeploci$locus %in% Pop2keeploci$locus)
table(Pop2keeploci$locus %in% Pop1keeploci$locus)

# Let's isolate the keep loci
Pop1filteredGI <- Pop1GI[loc = Pop1keeploci$locus] # This creates the character vector to use to filter out failed loci
nLoc(Pop1filteredGI) # Check for correct number of loci

# Let's isolate the keep loci
locNames(Pop2GI)
Pop2filteredGI <- Pop2GI[loc = Pop2keeploci$locus]

# This creates a list that only has loci that are in both lists (not what we want)
# completelocilist <- statvaltabPop1filtered$locus[statvaltabPop1filtered$locus %in% statvaltabPop2filtered$locus]

# This keeps loci that passed the test in any population
combinedkeeplist <- base::unique(c(statvaltabPop1filtered$locus, statvaltabPop2filtered$locus))
length(combinedkeeplist) # Check the length of the combined list of loci

# Retest HWE for filtered populations

# Run HWE test on each filtered Population
my.hwt1filtered <- hw.test(Pop1filteredGI, B=1000) # B=number of draws. Low for testing. Use default for real, B=1000
summary(my.hwt1filtered)
head(my.hwt1filtered)

my.hwt2filtered <- hw.test(Pop2filteredGI, B=1000) # B=number of draws. Low for testing. Use default for real, B=1000
summary(my.hwt2filtered)

# Convert filtered genind to hierfstat to view Population 1 basic stats, combine with hw.test and view loci that have too high Ho/Hs ratio
Pop1filteredHS <- genind2hierfstat(Pop1filteredGI)
basicstatPop1filtered <- basic.stats(Pop1filteredHS)
basicstatPop1filtered
head(basicstatPop1filtered$perloc)
statvaltabPop1filtered <- data.frame(basicstatPop1filtered$perloc) # Pulls out per locus statistics
statvaltabPop1filtered$locus <- rownames(statvaltabPop1filtered)
my.hwt1filtered <- data.frame(my.hwt1filtered)
my.hwt1filtered$locus <- row.names(my.hwt1filtered)
table(my.hwt1filtered$locus == statvaltabPop1filtered$locus) # Check for matching loci
statvaltabPop1filtered <- left_join(statvaltabPop1filtered, my.hwt1filtered, by = "locus") # Combines with HWE test results
head(statvaltabPop1filtered)
statvaltabPop1filtered$filter <- ifelse(statvaltabPop1filtered$Pr.chi.2.. < 0.01 & statvaltabPop1filtered$Ho / statvaltabPop1filtered$Hs > 1, 0, 1)
# Note that per ?basic.stats, Ho is observed heterozygosity, Hs is expected heterozygosity
ggplot(statvaltabPop1filtered, aes(x = Ho/Hs, fill = factor(filter))) +
  geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
plot(statvaltabPop1$Pr.chi.2...[statvaltabPop1$filter == 1], statvaltabPop1filtered$Pr.chi.2...)
# Plot shows the probabilities are the same both times

table(statvaltabPop1filtered$filter) # Check if filtering is complete

# Convert filtered genind to hierfstat to view Population 2 basic stats, combine with hw.test and view loci that have too high Ho/Hs ratio
Pop2filteredHS <- genind2hierfstat(Pop2filteredGI)
basicstatPop2filtered <- basic.stats(Pop2filteredHS)
basicstatPop2filtered
head(basicstatPop2filtered$perloc)
statvaltabPop2filtered <- data.frame(basicstatPop2filtered$perloc) # Pulls out per locus statistics
statvaltabPop2filtered$locus <- rownames(statvaltabPop2filtered)
my.hwt2filtered <- data.frame(my.hwt2filtered)
my.hwt2filtered$locus <- row.names(my.hwt2filtered)
table(my.hwt2filtered$locus == statvaltabPop2filtered$locus) # Check for matching loci
statvaltabPop2filtered <- left_join(statvaltabPop2filtered, my.hwt2filtered, by = "locus") # Combines with HWE test results
head(statvaltabPop2filtered)
statvaltabPop2filtered$filter <- ifelse(statvaltabPop2filtered$Pr.chi.2.. < 0.01 & statvaltabPop2filtered$Ho / statvaltabPop2filtered$Hs > 1, 0, 1)
# Note that per ?basic.stats, Ho is observed heterozygosity, Hs is expected heterozygosity
ggplot(statvaltabPop2filtered, aes(x = Ho/Hs, fill = factor(filter))) +
  geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
plot(statvaltabPop2$Pr.chi.2...[statvaltabPop2$filter == 1], statvaltabPop2filtered$Pr.chi.2...)
# Plot shows the probabilities are the same both times

table(statvaltabPop2filtered$filter) # Check if filtering is complete
table(statvaltabPop1filtered$filter) # Check if filtering is complete

# Determine number of loci
nLoc(Pop1filteredGI) # 2053
nLoc(Pop2filteredGI) # 1964
nLoc(Pop3GI) # 2430
nLoc(Pop4GI) # 2430

completelocilist <- statvaltabPop1filtered$locus[statvaltabPop1filtered$locus %in% statvaltabPop2filtered$locus]
length(completelocilist) # Number of common loci

# Both Individual Population filtering and Combined filtering lead to approximately 1,949 SNPs left to work with

# Select complete loci list from the filtered genind object
# Convert into genind objects
mygenind <- dartR::gl2gi(mygenlight2)
mygenind

# Isolate the keep loci
filteredGI <- mygenind[loc = completelocilist] # This creates the character vector to use to filter out failed loci
nLoc(filteredGI) # 1949

# Run HWE test Again
my.hwtfiltered <- hw.test(filteredGI, B = 1000) # B=number of draws. Low for testing. Use default for real, B=1000
summary(my.hwtfiltered)

# Convert filtered genind to hierfstat to view Population 1 basic stats, combine with hw.test and view loci that have too high Ho/Hs ratio
filteredHS <- genind2hierfstat(filteredGI)
basicstatfiltered <- basic.stats(filteredHS)
basicstatfiltered
head(basicstatfiltered$perloc)
statvaltabfiltered <- data.frame(basicstatfiltered$perloc) # Pulls out per locus statistics
statvaltabfiltered$locus <- rownames(statvaltabfiltered)
my.hwtfiltered <- data.frame(my.hwtfiltered)
my.hwtfiltered$locus <- row.names(my.hwtfiltered)
table(my.hwtfiltered$locus == statvaltabfiltered$locus) # Check for matching loci
statvaltabfiltered <- left_join(statvaltabfiltered, my.hwtfiltered, by = "locus") # Combines with HWE test results
head(statvaltabfiltered)
statvaltabfiltered$filter <- ifelse(statvaltabfiltered$Pr.chi.2.. < 0.01 & statvaltabfiltered$Ho / statvaltabfiltered$Hs > 1, 0, 1)
# Note that per ?basic.stats, Ho is observed heterozygosity, Hs is what people mean by expected
ggplot(statvaltabfiltered, aes(x = Ho/Hs, fill = factor(filter))) +
  geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
plot(statvaltabfiltered$Pr.chi.2...[statvaltabfiltered$filter == 1], statvaltabfiltered$Pr.chi.2...)
table(statvaltabfiltered$filter) # Check if filtering is complete

# Run Structure-like test
# Convert filtered GI to GL
FfilteredGL <- gi2gl(filteredGI)
FfilteredGL

# Plot cross-entropy results to assess optimal number of K
# Smaller values of cross-entropy usually mean better runs
# A plateau usually represents the K that best fits the data
plot(project2, col = "blue", cex = 1.5, pch = 19)

# Filter the populations to the 1967 SNPs in Pop2
Pop1FF <- Pop1filteredGI[loc = c()]
Pop3FF <- Pop3GI[loc = c()]
Pop4FF <- Pop4GI[loc = c()]

# Check loci counts again
nLoc(Pop1FF) # 1956
nLoc(Pop2filteredGI) # 1967
nLoc(Pop3FF) # 1967
nLoc(Pop4FF) # 1967

# Filter Again to Remove SNPs missing in Pop1
Pop2FF <- Pop2filteredGI[loc = c()]
Pop3FF <- Pop3FF[loc = c()]
Pop4FF <- Pop4FF[loc = c()]

# Check nLoc AGAIN
nLoc(Pop1FF) # 1956
nLoc(Pop2FF) # 1911
nLoc(Pop3FF) # 1911
nLoc(Pop4FF) # 1911

# Convert Pop1 to GL
Pop1FFGL <- gi2gl(Pop1FF)
Pop1reduced <- gl.keep.loc(Pop1FFGL, first = 1, last = 1911) # Reduce gl object to minimum number of SNPs
nLoc(Pop1reduced)

# Convert to GL
Pop2FFGL <- gi2gl(Pop2FF)
Pop3FFGL <- gi2gl(Pop3FF)
Pop4FFGL <- gi2gl(Pop4FF)

pop(Pop3FFGL)

# Bind
combinedFFGL <- postStrictclean <- rbind.genlight(Pop1reduced, Pop2FFGL, Pop3FFGL, Pop4FFGL)
combinedFFGL # Filtered file

# Convert to genind
combinedFFGI <- gl2gi(combinedFFGL)

# Run HWE test 
my.hwtFFGI <- hw.test(combinedFFGI, B = 1000) # B=number of draws. Low for testing. Use default for real, B=1000
summary(my.hwtFFGI)

# Convert genind to hierfstat, run basic stats, combine with hw.test and view loci that have too high Ho/Hs ratio
HSFFGI <- genind2hierfstat(combinedFFGI)
basicstatFFGI <- basic.stats(HSFFGI)
basicstatFFGI
head(basicstatFFGI$perloc)
statvaltabFFGI <- data.frame(basicstatFFGI$perloc) # Pulls out per locus statistics
statvaltabFFGI$locus <- rownames(my.hwtFFGI)
statvaltabFFGI <- bind_cols(statvaltabFFGI, my.hwtFFGI) # Combines with HWE test results
head(statvaltabFFGI)
statvaltabFFGI$filter <- ifelse(statvaltabFFGI$Pr.exact < 0.01 & statvaltabFFGI$Ho / statvaltabFFGI$Hs > 1, 0, 1)
# Note that per ?basic.stats, Ho is observed heterozygosity, Hs is what people mean by expected
ggplot(statvaltabFFGI, aes(x = Ho/Hs, fill = Pr.exact < 0.05/1911)) +
  geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
ggplot(statvaltabFFGI, aes(x = Ho, fill = Pr.exact < 0.05/1911)) +
  geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
ggplot(statvaltabFFGI, aes(x = Hs, fill = Pr.exact < 0.05/1911)) +
  geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
ggplot(statvaltabFFGI, aes(x = Ho/Hs, fill = factor(filter))) +
  geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
table(statvaltabFFGI$filter) # Check if filtering is complete

# Convert filtered genind object to genlight
FFGLgeno <- gl2geno(combinedFFGL, outfile = "combinedFFGL", outpath = getwd())

# Perform SNMF analysis with different values of K (number of ancestral populations)
project = snmf("combinedFFGL.geno", K = 1:9, repetitions = 10, project = "new",
               entropy = TRUE, ploidy = 2)

# Plot cross-entropy results to assess the optimal number of K
# Smaller values of cross-entropy usually mean better runs
# A plateau usually represents the K that best fits the data
plot(project, col = "blue", cex = 1.5, pch = 19)

# Extract the cross-entropy of all runs where K = 2
ce3 = cross.entropy(project, K = 2)
ce3

# Find the run with the lowest cross-entropy
lowest.ce3 = which.min(ce3)
lowest.ce3

# Extract Q-matrix for the best run (K = 2)
qmatrix3 = as.data.frame(Q(project, K = 2, run = lowest.ce3))
head(qmatrix3)

# Add population labels to the genind object
pop(combinedFFGI) <- c("pop1","pop1","pop1","pop1","pop1",
                       "pop1","pop1","pop1","pop1","pop1",
                       "pop1","pop1","pop1","pop1","pop1",
                       "pop1","pop1", # one through 17
                       "pop2.1","pop2.3","pop2.3","pop2.3","pop2.3", # 18 through 22
                       "pop2.3","pop2.3","pop2.3","pop2.3","pop2.3", # 23 thru 27
                       "pop2.3","pop2.2","pop2.2","pop2.3","pop2.3", # 28 thru 32
                       "pop2.3","pop2.3","pop2.3","pop2.3","pop2.3", # 33 thru 37
                       "pop2.3","pop2.3","pop2.3","pop2.3","pop2.3", # 38 thru 42
                       "pop2.3","pop2.3","pop2.2","pop2.3","pop2.3", # 43 thru 47
                       "pop2.3","pop2.3","pop2.3","pop2.3","pop2.3", # 48 thru 52
                       "pop2.3","pop2.3","pop2.3","pop2.3","pop2.3", # 53 thru 57
                       "pop2.3","pop2.3","pop2.3","pop2.2","pop2.2", # 58 thru 62
                       "pop2.2","pop2.2","pop2.2","pop2.1","pop2.1", # 63 thru 67
                       "pop2.1","pop2.1","pop2.1","pop2.1", # 68 thru 71
                       "pop3","pop3","pop3","pop3","pop3", # 72 thru 76
                       "pop3","pop3", # 77 thru 78
                       "pop4","pop4","pop4","pop4","pop4", # 79 thru 83
                       "pop4","pop4") # 84 and 85

# Plot cross-entropy results to assess the optimal number of K
plot(projectOC, col = "blue", cex = 1.5, pch = 19)

# Extract the cross-entropy of all runs where K = 3
ce = cross.entropy(projectOC, K = 3)
ce

# Find the run with the lowest cross-entropy
lowest.ce = which.min(ce)
lowest.ce

# Extract Q-matrix for the best run (K = 3)
qmatrix = as.data.frame(Q(projectOC, K = 3, run = lowest.ce))
head(qmatrix)

# Label column names of qmatrix
ncol(qmatrix)
cluster_names = c()
for (i in 1:ncol(qmatrix)) {
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) = cluster_names
head(qmatrix)

# Add individual IDs to qmatrix
qmatrix$Ind = indNames(combinedFFGI)
indNames(combinedFFGI)

# Add site IDs to qmatrix
qmatrix$Site = combinedFFGI$pop
head(qmatrix)

# Convert dataframe to long format
qlong = melt(qmatrix, id.vars = c("Ind", "Site"))
head(qlong)

# Subset data to reduce computation time
subsites = sort(c("Belize", "Florida", "North Carolina", "Texas"))
subsites
data_filt = popsub(FF626GL, sublist = subsites)
data_filt

# Change order of sites using the factor function
site.order = c("Belize", "Florida", "North Carolina", "Texas")
qlong$Site_ord = factor(qlong$Site, levels = site.order)

# Adjust facet labels
levels(qlong$Site_ord)
facet.labs = c("Belize", "Florida", "North Carolina", "Texas")
levels(qlong$Site_ord) = facet.labs
levels(qlong$Site_ord)

# Define colour palette
pal = colorRampPalette(c("tomato", "lightblue", "lightgreen"))
cols = pal(length(unique(qlong$variable)))

# Plot admixture barplot
admix.bar = ggplot(data = qlong, aes(x = Ind, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~Site, scales = "free", ncol = 2) +
  scale_fill_manual(values = cols) +
  ylab("Admixture proportion") +
  # xlab("Individual") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour = "black", size = 12),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
admix.bar
# ggsave("4.admixture_barplot.png", width = 6, height = 10, dpi = 300)

# ----------------- #
# Prepare pie charts
# ----------------- #

# Calculate mean admixture proportions for each site
head(qmatrix)
clusters = grep("Cluster", names(qmatrix)) # indexes of cluster columns
clusters
avg_admix1 = aggregate(qmatrix[, clusters], list(qmatrix$Site), mean)

# Order alphabetically by site
avg_admix2 = avg_admix1[order(as.character(avg_admix1$Group.1)), ]
avg_admix2

# Convert dataframe from wide to long format
avg_admix3 = melt(avg_admix2, id.vars = "Group.1")
head(avg_admix3)

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable)) +
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols) +
    theme_void()
}

# Test function on one site
pie_charts(avg_admix3, site = "Texas", cols = cols)

# Apply function to all sites using for loop
pies = list()
for (i in subsites){
  pies[[i]] = pie_charts(admix_df = avg_admix3, site = i, cols = cols) 
}

# ----------------- #
# Prepare basemap
# ----------------- #

# Import csv file containing coordinates
coords = read.csv("coordinates.csv")

# Subset coordinates
coords = coords[coords$Code %in% subsites, ]

# Order alphabetically by site
coords = coords[order(coords$Code), ] 
coords$Code

avg_admix3$Group.1
coords$Code

# Check if the order matches the coords order
as.character(avg_admix3$Group.1) == as.character(coords$Code)

# Load necessary libraries for mapping
library(rworldmap)
library(rworldxtra)
library(sf)

# Get high resolution world map and convert to sf object
world <- getMap(resolution = "high")
world <- st_as_sf(world)
class(world)

# Plot basemap with specified longitude and latitude range
basemap = ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-99, -74), ylim = c(14, 36), expand = FALSE) +
  xlab("Longitude") +
  ylab("Latitude")
basemap

# ----------------- #
# Add pie charts to basemap
# ----------------- #

# Extract coordinates for each site
coord.list = list()
for (i in subsites){
  coord.list[[i]] = c(subset(coords, Code == i)$Lon, subset(coords, Code == i)$Lat)
}
coord.list

# Define pie chart sizes
radius = 1

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(subsites)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

pies.ac

# Add layers to basemap
pie.map = basemap + pies.ac
pie.map
# ggsave("4.pie_charts_map.png", width = 8, height = 10, dpi = 300)

# Load ggpubr for further plotting functionalities
library(ggpubr)

# Combine ggplots
ggarrange(admix.bar + theme(legend.position = "right") + labs(title = "Individual admixture proportions", tag = "A"),
          pie.map + labs(title = "Mean admixture proportions for each site", tag = "B"))
ggsave("3.Admixture_bar_map.png", width = 15, height = 10, dpi = 300)

# Bar plot for admixture coefficients with lines indicating group divisions
OCbp3 <- barplot(t(qmatrix), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions (K=3)", ylab = "Admixture coefficients", border = NA, space = 0)
abline(v = 18, col = "black")  # Adds line right after individual 18
abline(v = 72, col = "black")  # Adds line right after individual 72
abline(v = 79, col = "black")  # Adds line right after individual 79
abline(v = 85, col = "black")  # Adds line right after individual 85

# Analysis without population-specific tests

# Read in vcf file
no626 <- read.vcfR("VCF6_26.vcf")
no626

# Convert into genind objects
mygenind626 <- dartR::gl2gi(mygenlight2)
mygenind626

# Run HWE test
my.hwt <- hw.test(mygenind626, B = 1000)  # B=number of draws. Low for testing. Use default for real, B=1000
summary(my.hwt)

# Convert genind to heirfstat, run basic stats, combine with hw.test, and view loci that have too high Ho/Hs ratio
HS626 <- genind2hierfstat(mygenind626)
basicstat <- basic.stats(HS626)
basicstat
head(basicstat$perloc)
statvaltab <- data.frame(basicstat$perloc)  # Pulls out per locus
statvaltab$locus <- rownames(my.hwt)
statvaltab <- bind_cols(statvaltab, my.hwt)  # Combines with HWE test results
head(statvaltab)
statvaltab$locus <- rownames(statvaltab)
my.hwt <- data.frame(my.hwt)
my.hwt$locus <- row.names(my.hwt)
head(my.hwt)
table(my.hwt$locus == statvaltab$locus)  # Problem is not in table matching
statvaltab <- left_join(statvaltab, my.hwt, by = "locus")  # Combines with HWE test results
head(statvaltab)

# Now using chi squared probability rather than Pr. exact
statvaltab$filter <- ifelse(statvaltab$Pr.chi.2... < 0.01 & statvaltab$Ho / statvaltab$Hs > 1, 0, 1)
table(statvaltab$filter)

Keeploci <- subset(statvaltab, filter > 0)  # Pull loci to keep
dim(Keeploci)

# Filter for loci to keep
Filtered626GI <- mygenind626[loc = c(Keeploci$locus)]

# Retest HWE for filtered populations
nInd(Filtered626GI)
Filtered626GI

indNames(Filtered626GI)

# Run HWE test on filtered populations
my.hwtfiltered <- hw.test(Filtered626GI, B = 1000)  # B=number of draws. Low for testing. Use default for real, B=1000
summary(my.hwtfiltered)

# Convert filtered genind to heirfstst, run basic stats, combine with hw.test, and view loci that have too high Ho/Hs ratio
filteredHS <- genind2hierfstat(Filtered626GI)
basicstatfiltered <- basic.stats(filteredHS)
basicstatfiltered
head(basicstatfiltered$perloc)
statvaltabfiltered <- data.frame(basicstatfiltered$perloc)  # Pulls out per locus
statvaltabfiltered$locus <- rownames(statvaltabfiltered)
my.hwtfiltered <- data.frame(my.hwtfiltered)
my.hwtfiltered$locus <- row.names(my.hwtfiltered)
table(my.hwtfiltered$locus == statvaltabfiltered$locus)
statvaltabfiltered <- left_join(statvaltabfiltered, my.hwtfiltered, by = "locus")  # Combines with HWE test results
head(statvaltabfiltered)
statvaltabfiltered$filter <- ifelse(statvaltabfiltered$Pr.chi.2... < 0.01 & statvaltabfiltered$Ho / statvaltabfiltered$Hs > 1, 0, 1)

# Plot HWE equilibrium deviations
ggplot(statvaltabfiltered, aes(x = Ho / Hs, fill = factor(filter))) +
  geom_histogram() + ggtitle("Significant deviation for HW equilibrium, Bonferroni corrected")
plot(statvaltab$Pr.chi.2...[statvaltab$filter == 1], statvaltabfiltered$Pr.chi.2...)

# Convert from genind to genlight
FF626GL <- gi2gl(Filtered626GI)
FF626GL

# Convert genlight to geno format for SNMF analysis
FF626geno <- gl2geno(FF626GL, outfile = "FF626", outpath = getwd())

# Run SNMF analysis
project626 <- snmf("FF626.geno", K = 1:9, repetitions = 10, project = "new", entropy = TRUE, ploidy = 2)
plot(project626)

# Plot cross-entropy results to assess optimal number of K
# Smaller values of cross-entropy usually mean better runs
# A plateau usually represents the K that best fits the data
plot(project2, col = "blue", cex = 1.5, pch = 19)

# Extract the cross-entropy of all runs where K = 2
ce2 <- cross.entropy(project2, K = 1)
ce2

# Find the run with the lowest cross-entropy
lowest.ce2 <- which.min(ce2)
lowest.ce2

# Extract Q-matrix for the best run
qmatrix2 <- as.data.frame(Q(project2, K = 1, run = lowest.ce2))
head(qmatrix2)

# Calculate Minor Allele Frequency (MAF) for each HWE filter category

# Out ALL
OutAllmfreq <- minorAllele(filteredGI)
mean(OutAllmfreq)
quantile(OutAllmfreq)
range(OutAllmfreq)

# Out Combo
OutCombomfreq <- minorAllele(Filtered626GI)
mean(OutCombomfreq)
quantile(OutCombomfreq)
range(OutCombomfreq)

# No Filter
NoFiltermfreq <- minorAllele(mygenind626)
mean(NoFiltermfreq)
quantile(NoFiltermfreq)
range(NoFiltermfreq)

# Convert genlight to VCF format
gl2vcf(FfilteredGL, outfile = 'OutAll', outpath = getwd())
gl2vcf(FF626GL, outfile = 'OutCombo', outpath = getwd())

# PCAs and DAPCs are in file "PCA comparisons_7_5"

# AMOVAs

# Now doing AMOVA test for differences among populations

# Set strata for FfilteredGL
FfilteredGL@strata <- NULL
strata(FfilteredGL) <- data.frame(population = pop(FfilteredGL))
table(strata(FfilteredGL))

# AMOVA without individual variation
OAamova <- poppr.amova(FfilteredGL, within = FALSE, ~population)
OAamova
OAamova.test <- randtest(OAamova, nrepet = 1000) # Test for significance
plot(OAamova.test)
OAamova.test

# Set strata for FF626GL
FF626GL@strata <- NULL
strata(FF626GL) <- data.frame(population = pop(FF626GL))
table(strata(FF626GL))

# AMOVA without individual variation
OCamova <- poppr.amova(FF626GL, within = FALSE, ~population)
OCamova
OCamova.test <- randtest(OCamova, nrepet = 1000) # Test for significance
plot(OCamova.test)
OCamova.test

# Set strata for mygenlight2
mygenlight2@strata <- NULL
strata(mygenlight2) <- data.frame(population = pop(mygenlight2))
table(strata(mygenlight2))

# AMOVA without individual variation
NFamova <- poppr.amova(mygenlight2, within = FALSE, ~population)
NFamova
NFamova.test <- randtest(NFamova, nrepet = 1000) # Test for significance
plot(NFamova.test)
NFamova.test

# Calculate pairwise Nei's Fst for filtered dataset
OAlane_fst <- genet.dist(filteredGI, method = "Nei87") %>% round(digits = 3)
OAlane_fst

# Bootstrap analysis for confidence intervals
boot.ppfst(dat = filteredGI, nboot = 100, quant = c(0.025, 0.975), diploid = TRUE)

# Calculate pairwise Nei's Fst for "Out Combo" filtered dataset
OClane_fst <- genet.dist(Filtered626GI, method = "Nei87") %>% round(digits = 3)
OClane_fst

# Calculate pairwise Nei's Fst for "No Filter" dataset
NFlane_fst <- genet.dist(mygenind626, method = "Nei87") %>% round(digits = 3)
NFlane_fst

# Gene Flow Metric calculation
(1/4) * ((1 / OAlane_fst) - 1)
(1/4) * ((1 / OClane_fst) - 1)
(1/4) * ((1 / NFlane_fst) - 1)

# Use StAMPP package for pairwise Fst with bootstrapping
OAfst <- stamppFst(FfilteredGL, nboots = 10100, percent = 95, nclusters = 1)
OAfst

OCfst <- stamppFst(FF626GL, nboots = 10100, percent = 95, nclusters = 1)
OCfst

NFfst <- stamppFst(mygenlight2, nboots = 10100, percent = 95, nclusters = 1)
NFfst

# Admixture Analysis using LEA package

# Convert genlight to geno format
FF75geno <- gl2geno(FfilteredGL, outfile = "FF75", outpath = getwd())

# Run SNMF analysis for "Out All" filtered data
projectOA <- snmf("FF75.geno", K = 1:9, repetitions = 100, project = "new", entropy = TRUE, ploidy = 2, alpha = 100)
plot(projectOA, cex.main = 2, cex.sub = 1.5, cex.lab = 1.5, cex.axis = 1.5)

# Run SNMF analysis for "Out Combo" filtered data
projectOC <- snmf("FOC75.geno", K = 1:9, repetitions = 100, project = "new", entropy = TRUE, ploidy = 2, alpha = 100)
plot(projectOC)

# Run SNMF analysis for "No Filter" dataset
NFgeno <- gl2geno(mygenlight2, outfile = "FNF75", outpath = getwd())
projectNF <- snmf("FNF75.geno", K = 1:9, repetitions = 100, project = "new", entropy = TRUE, ploidy = 2, alpha = 100)
plot(projectNF)

# Create admixture bar plot
OCbp <- barplot(t(qmatrix), col = c("green", "blue"), xlab = "Accessions", ylab = "Admixture coefficients", main = "Ancestry Matrix", border = NA, space = 0)
abline(v = 18, col = "black") # adds line right after indv 18
abline(v = 72, col = "black") # adds line right after indv 72
abline(v = 79, col = "black") # adds line right after indv 79
abline(v = 85, col = "black") # adds line right after indv 85

# Extract cross-entropy for K = 2, 3, and 4
OCce2 <- cross.entropy(projectOC, K = 2)
OCce2
OClowest.ce2 <- which.min(OCce2)
OClowest.ce2

OCce3 <- cross.entropy(projectOC, K = 3)
OCce3
OClowest.ce3 <- which.min(OCce3)
OClowest.ce3

OCce4 <- cross.entropy(projectOC, K = 4)
OCce4

# Find the run with the lowest cross-entropy for K=4
OClowest.ce4 <- which.min(OCce4)
OClowest.ce4

# Extract Q-matrices for the best runs at K=2, K=3, and K=4
OCQmatrix2 <- Q(projectOC, K = 2, run = 29)
OCQmatrix3 <- Q(projectOC, K = 3, run = 29)
OCQmatrix4 <- Q(projectOC, K = 4, run = 88)

# Bar plot of Q-matrix
barplot(t(Qmatrix))

# Set up plotting area to have 3 rows
par(mfrow = c(3, 1))

# Admixture bar plots for K=2, K=3, and K=4 for "Out Combo" filtered data
OCbp2 <- barplot(t(OCQmatrix2), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions (K=2)", ylab = "Admixture coefficients",
                 main = "Ancestry Matrix: Out Combo", border = NA, space = 0)
abline(v = c(18, 72, 79, 85), col = "black")

OCbp3 <- barplot(t(OCQmatrix3), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions (K=3)", ylab = "Admixture coefficients", border = NA, space = 0)
abline(v = c(18, 72, 79, 85), col = "black")

OCbp4 <- barplot(t(OCQmatrix4), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions (K=4)", ylab = "Admixture coefficients", border = NA, space = 0)
abline(v = c(18, 72, 79, 85), col = "black")

# Extract the cross-entropy of all runs where K = 2 for "Out All" filtered data
OAce2 <- cross.entropy(projectOA, K = 2)
OAce2

# Find the run with the lowest cross-entropy
OAlowest.ce2 <- which.min(OAce2)
OAlowest.ce2

# Extract the cross-entropy of all runs where K = 3 for "Out All" filtered data
OAce3 <- cross.entropy(projectOA, K = 3)
OAce3

# Find the run with the lowest cross-entropy
OAlowest.ce3 <- which.min(OAce3)
OAlowest.ce3

# Extract the cross-entropy of all runs where K = 4 for "Out All" filtered data
OAce4 <- cross.entropy(projectOA, K = 4)
OAce4

# Find the run with the lowest cross-entropy
OAlowest.ce4 <- which.min(OAce4)
OAlowest.ce4

# Extract Q-matrices for the best runs at K=2, K=3, and K=4 for "Out All" filtered data
OAQmatrix2 <- Q(projectOA, K = 2, run = 70)
OAQmatrix3 <- Q(projectOA, K = 3, run = 53)
OAQmatrix4 <- Q(projectOA, K = 4, run = 53)

# Bar plot of Q-matrix
barplot(t(Qmatrix))

# Set up plotting area to have 3 rows
par(mfrow = c(3, 1))

# Admixture bar plots for K=2, K=3, and K=4 for "Out All" filtered data
OAbp2 <- barplot(t(OAQmatrix2), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions (K=2)", ylab = "Admixture coefficients",
                 main = "Ancestry Matrix: Out All", border = NA, space = 0)
abline(v = c(18, 72, 79, 85), col = "black")

OAbp3 <- barplot(t(OAQmatrix3), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions (K=3)", ylab = "Admixture coefficients", border = NA, space = 0)
abline(v = c(18, 72, 79, 85), col = "black")

OAbp4 <- barplot(t(OAQmatrix4), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions (K=4)", ylab = "Admixture coefficients", border = NA, space = 0)
abline(v = c(18, 72, 79, 85), col = "black")

# Extract the cross-entropy of all runs where K = 2 for "No Filter" dataset
NFce2 <- cross.entropy(projectNF, K = 2)
NFce2

# Find the run with the lowest cross-entropy
NFlowest.ce2 <- which.min(NFce2)
NFlowest.ce2

# Extract the cross-entropy of all runs where K = 3 for "No Filter" dataset
NFce3 <- cross.entropy(projectNF, K = 3)
NFce3

# Find the run with the lowest cross-entropy for K=3
NFlowest.ce3 <- which.min(NFce3)
NFlowest.ce3

# Extract the cross-entropy of all runs where K=4
NFce4 <- cross.entropy(projectNF, K = 4)
NFce4

# Find the run with the lowest cross-entropy for K=4
NFlowest.ce4 <- which.min(NFce4)
NFlowest.ce4

# Extract Q-matrices for the best runs at K=2, K=3, and K=4 for "No Filter" dataset
NFQmatrix2 <- Q(projectNF, K = 2, run = 100)
NFQmatrix3 <- Q(projectNF, K = 3, run = 53)
NFQmatrix4 <- Q(projectNF, K = 4, run = 48)

# Bar plot of Q-matrix
barplot(t(Qmatrix))

# Set up plotting area to have 3 rows
par(mfrow = c(3, 1))

# Admixture bar plots for K=2, K=3, and K=4 for "No Filter" dataset
NFbp2 <- barplot(t(NFQmatrix2), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions(K=2)", ylab = "Admixture coefficients",
                 main = "Ancestry Matrix: No Filter", border = NA, space = 0)
abline(v = c(18, 72, 79, 85), col = "black")

NFbp3 <- barplot(t(NFQmatrix3), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions (K=3)", ylab = "Admixture coefficients", border = NA, space = 0)
abline(v = c(18, 72, 79, 85), col = "black")

NFbp4 <- barplot(t(NFQmatrix4), col = c("tomato", "lightblue", "lightgreen", "gold"),
                 xlab = "Accessions (K=4)", ylab = "Admixture coefficients", border = NA, space = 0)
abline(v = c(18, 72, 79, 85), col = "black")

# Isolation By Distance Analysis

# Datasets:
# - Out All: FfilteredGL (genlight), filteredGI (genind)
# - Out Combo: FF626GL (genlight), Filtered626GI (genind)
# - No Filter: mygenlight2 (genlight), mygenind626 (genind)
