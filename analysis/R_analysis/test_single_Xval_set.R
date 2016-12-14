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

data_path = './../../data/cross_val_data'
list.files(path = data_path, pattern = '.*fold4.*.tsv')

m <- read.csv(paste0(data_path, '/methanotroph_fold4_ss_filtered_train.tsv'),
              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(m)

nmm <- read.csv(paste0(data_path, '/methylotroph_fold4_ss_filtered_train.tsv'),
                sep = '\t', header = FALSE, stringsAsFactors = FALSE) #, colClasses=c("product"="character"))
dim(nmm)

dim(m) # [1]   61 4840
dim(nmm) # [1]    61 10123

colnames(m)[0:10]
m[0:3, 0:10]
nmm[0:3, 0:10]

# Cannot have NAs in x or z
# Both should have 0 rows.  (Not true for nmm if read.table is used instaed of read.csv)
dim(m[rowSums(is.na(m)) > 0,])
dim(nmm[rowSums(is.na(nmm)) > 0,])

x <- m
z <- nmm
# Can't have columns with standard deviation = 0
# I got rid of them in the Python script, so these filters shouldn't be doing anything any longer.
dim(x)
dim(Filter(function(q) sd(q) != 0, x))
x = Filter(function(q) sd(q) != 0, x)
dim(z)
z = Filter(function(q) sd(q) != 0, z)
dim(z)

dir.create('./results')
getwd()

model_stats <- function(CCA_obj){
        x_penalty = CCA_obj$penaltyx
        u = CCA_obj$u[,1]
        u_len = length(u)
        u_zeros = sum(u == 0)
        u_coeffs = sum(u != 0)
        u_frac_zeros = u_zeros/u_len

        z_penalty = CCA_obj$penaltyz
        v = CCA_obj$v[,1]
        v_len = length(v)
        v_zeros = sum(v == 0)
        v_coeffs = sum(v != 0)
        v_frac_zeros = v_zeros/v_len
        return(data.frame(x_penalty=x_penalty,
                          u_len=u_len, u_zeros=u_zeros, u_coeffs=u_coeffs, u_frac_zeros=u_frac_zeros,
                          z_penalty=z_penalty,
                          v_len=v_len, v_zeros=v_zeros, v_coeffs=v_coeffs, v_frac_zeros=v_frac_zeros))
}

penalty_x = 0.0335  # 0.33 --> 7, 0.34 --> 9
penalty_z = 0.022
# remove cols with zero variance (method from Stack Overflow)
x2 = x[,apply(x, 2, var, na.rm=TRUE) != 0]
z2 = z[,apply(z, 2, var, na.rm=TRUE) != 0]
dim(x)
dim(x2)

#final_model = CCA(x2, z2, typex="standard", typez="standard", K=1,
final_model = CCA(x, z, typex="standard", typez="standard", K=1,
                  niter=1000,
                  #standardize = FALSE,  <-- would be great if this magically worked.
                  penaltyx=penalty_x, penaltyz=penalty_z)
final_results = model_stats(final_model)
final_results

# Save these resutls
#write.table(final_model$u, file = './results/u_whole_training_set--8_features.csv', row.names = FALSE)
#write.table(final_model$v, file = './results/v_whole_training_set--8_features.csv', row.names = FALSE)

print(final_model)
