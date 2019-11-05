# Download GenABEL package file from CRAN if not available in our R version  
library(GenABEL)

# loading both genotype and phenotype file.
data = load.gwaa.data(phenofile = "phenotype.dat", genofile = "genotype.raw")

# to access the genotype data at markers 3 through 10 for the first individual.
gtdata(data)[1,3:10]

# to find out on which chromosome markers 3 through 10 are present.
chromosome(gtdata(data)[1,3:10])

# name of a marker.
snpnames(gtdata(data)[,17000])

# to see the genotypes at marker 1177 for the first 4 individuals.
as.character(gtdata(data)[1:4,1])

# to see the reference allele.
refallele(gtdata(data)[,1177])

#summary
summary.snp.data(gtdata(data))

# data for ct(continuous trait) 
ct = phdata(data)$ct

# plot the distribution of ct 
hist(ct, col = "slateblue", main = "Distribution of ct")

# observing this box plot showing the ranges of the variable for both M & F.
# this is to check if the continuos trait depends on one of our other variables.
boxplot(ct~phdata(data)$sex, col = c("blue", "red"), xlab = "sex", names = c("F", "M"), main = "ct by sex")

# a function in GenABEL to check the quality of markers.
?check.marker

# Quality Control check
qc <- check.marker(data, call = 0.95, perid.call = 0.95, maf = 1e-08, p.lev = 1e-08) 

# individuals and SNP that passed the QC check
data.qc <- data[qc$idok, qc$snpok]

# test to see whether or not ct depends on any of the markers
an <- qtscore(ct~1, data=data.qc, trait="gaussian")

# plot the results in a "Manhattan" plot
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")

# plot the genomic inflation corrected plot 
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")

# adding Bonferroni correction(divides each p-value for each test by the number of tests) to avoid false positives.
pval.threshold <- 0.05
bonferroni <- -log10(pval.threshold / nids(data.qc))

# replotting with a line representing Bonferroni score
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")
abline(h=bonferroni, lty=3, color="red")







