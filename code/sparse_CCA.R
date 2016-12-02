# Background:
# CCA with sparsity http://biorxiv.org/content/biorxiv/early/2016/04/06/047217.full.pdf
# cites https://faculty.washington.edu/dwitten/Papers/pmd.pdf
# R package: https://cran.r-project.org/web/packages/PMA/PMA.pdf
# R package cites:
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2861323/pdf/sagmb1470.pdf
# See /ML_meta-omics/analysis/CCA_methanotroph_expression_methylotroph_expression/sparse_CCA_X_m_expression_Y_nmm_abundance.R
# for work that led to this script.

# Inputs: csv filename for gene expression vectors, penalty_x, penalty_z
# Outputs: u and v vectors from CCA


library("PMA")
library("ggplot2")

args = commandArgs(trailingOnly=TRUE)
x_filename = args[1]
print(paste("x filename:", x_filename))
z_filename = args[2]
print(paste("z filename:", z_filename))
penaltyx = as.numeric(args[3])
penaltyz = as.numeric(args[4])
print(paste('penaltyx:', penaltyx))
print(paste('penaltyz:', penaltyz))

x <- read.csv(x_filename, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
z <- read.csv(z_filename, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Note: R is modifying the gene (column) names:
# http://stackoverflow.com/questions/10441437/x-in-my-column-names-of-an-r-data-frame
# (1->4)-alpha-D-glucan 1-alpha-D-glucosylmutase
# becomes "X.1..4..alpha.D.glucan.1.alpha.D.glucosylmutase"

dim(x) # [1]   83 4593
dim(z) # [1]   83 7355

# Cannot have NAs in x or z
# Both should have 0 rows.  (Not true for nmm if read.table is used instaed of read.csv)
stopifnot(nrow(x[rowSums(is.na(x)) > 0,]) == 0)
stopifnot(nrow(z[rowSums(is.na(z)) > 0,]) == 0)

# get rid of the 'product' column, and transpose.
x = x[, -1]
z = z[, -1]

# train model
print('train model')
model <- CCA(x, z, typex="standard", typez="standard", K=1,
             penaltyx=penaltyx, penaltyz=penaltyz)

# prepare filename to save resulting weights:
u_filename = gsub(x_filename, pattern = '.tsv', replacement = '_u.tsv')
v_filename = gsub(z_filename, pattern = '.tsv', replacement = '_v.tsv')

u <- model$u
v <- model$v

print(paste('save output files to:', u_filename, v_filename))
write.table(u, file = u_filename, row.names = FALSE)
write.table(v, file = v_filename, row.names = FALSE)

