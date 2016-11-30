# Background:
# CCA with sparsity http://biorxiv.org/content/biorxiv/early/2016/04/06/047217.full.pdf
# cites https://faculty.washington.edu/dwitten/Papers/pmd.pdf
# R package: https://cran.r-project.org/web/packages/PMA/PMA.pdf
# R package cites:
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2861323/pdf/sagmb1470.pdf


# to prepare data, run the Jupyter notebook prepare_X_m_expression_Y_nmm_expression.ipynb
# Note: need to trim off descriptive columns and transpose.

library("PMA")
library("ggplot2")

m <- read.csv('../../data/m_nmm_expression--sum_by_gene/methanotroph_expression_pooled_on_gene_name.tsv',
        sep = '\t', header = TRUE, stringsAsFactors = FALSE)
dim(m)
# problem loding: the ' in tRNA G18 (ribose-2'-O)-methylase SpoU is giving trouble.
# read in as tRNA G18 (ribose-2-O)-methylase SpoU and the tabs are all part of the first column
# fixed by using read.csv instead of read.table.
        # http://stackoverflow.com/questions/9620155/how-to-read-a-csv-file-containing-apostrophes-into-r
nmm <- read.csv('../../data/m_nmm_expression--sum_by_gene/methylotroph_expression_pooled_on_gene_name.tsv',
                  sep = '\t', header = TRUE, stringsAsFactors = FALSE) #, colClasses=c("product"="character"))

# Note: R is modifying the gene (column) names:
# http://stackoverflow.com/questions/10441437/x-in-my-column-names-of-an-r-data-frame
# (1->4)-alpha-D-glucan 1-alpha-D-glucosylmutase
# becomes "X.1..4..alpha.D.glucan.1.alpha.D.glucosylmutase"

colnames(m)[0:10]
m[0:3, 0:10]
sample_names <- m['X']
head(sample_names)
nmm[0:3, 0:10]

# Cannot have NAs in x or z
# Both should have 0 rows.  (Not true for nmm if read.table is used instaed of read.csv)
dim(m[rowSums(is.na(m)) > 0,])
dim(nmm[rowSums(is.na(nmm)) > 0,])

# get rid of the 'product' column, and transpose.
dim(m)
lapply(m, class)

x = m[, -1]
colnames(x)[0:5]

z = nmm[, -1]
colnames(z)[0:5]

# Can't have columns with standard deviation = 0
dim(x)
dim(Filter(function(q) sd(q) != 0, x))
x = Filter(function(q) sd(q) != 0, x)
dim(z)
z = Filter(function(q) sd(q) != 0, z)
dim(z)


out <- CCA(x, z, typex="standard", typez="standard", K=1)
print(out,verbose=TRUE)

perm.out <- CCA.permute(x,z,typex="standard",typez="standard",nperms=7)

perm.out[0]
names(perm.out)

out07 <- CCA(x, z, typex="standard", typez="standard", K=1, penaltyx = 0.7, penaltyz = 0.7)
head(head(out07))

dim(out07["u"])
dir.create('./results')

getwd()
write.table(out07$u, file = './results/u.csv')
write.table(out07$v, file = './results/v.csv')


