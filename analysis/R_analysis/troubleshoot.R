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
library(reshape2)

standardize = TRUE
num_iter = 1000

# ---- Load methhanotrophy genes and gene names ----
m <- read.csv('../../data/cross_val_data/methanotroph_fold4_ss_filtered_train.tsv',
              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(m)
m_genes <- read.csv('../../data/cross_val_data/methanotroph_fold4_ss_filtered_genes.tsv',
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(m_genes)

# ---- Load methylotrophy genes and gene names ----
nmm <- read.csv('../../data/cross_val_data/methylotroph_fold4_ss_filtered_train.tsv',
                sep = '\t', header = FALSE, stringsAsFactors = FALSE) #, colClasses=c("product"="character"))
dim(nmm)
nmm_genes <- read.csv('../../data/cross_val_data/methylotroph_fold4_ss_filtered_genes.tsv',
                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(nmm_genes)

m[0:3, 0:10]
nmm[0:3, 0:10]

# --- Make sure NAs were removed from X, Z
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

#======  Loop over some different penalty values and find the number of zeros =========

test_pair_of_penalties <- function(penalty_list, x_penalty, z_penalty){
        print(paste("x penalty: " + x_penalty + "z penalty" + z))
        model = CCA(x, z, typex="standard", typez="standard", K=1,
                    standardize = standardize,
                    niter=num_iter,
                    penaltyx=x_penalty, penaltyz=z_penalty_list)
        return(model)}

analyze_penalty_grid <- function(penalty_list, num_iter){
        results_coarse = data.frame()
        print(results_coarse)
        for (x_penalty in penalty_list){
                for (z_penalty in penalty_list){
                        print(sprintf("--- x penalty: %f ---- z_penalty: %f ----",
                                x_penalty, z_penalty))
                        model <- CCA(x, z, typex="standard", typez="standard", K=2,
                                     standardize = standardize,
                                     niter=num_iter,
                                     penaltyx=x_penalty, penaltyz=z_penalty)
                        results_coarse <- rbind(results_coarse, model_stats(model))
                }
        }
        return(results_coarse)
}

#demo_results = analyze_penalty_grid(c(0.1, 0.2), num_iter)
#demo_results
#p <- ggplot(demo_results, aes(x_penalty, z_penalty)) +
#        geom_tile(aes(fill = u_coeffs)) +
#        scale_fill_gradient(low = "white", high = "steelblue", limits=c(0, max(demo_results$u_coeffs))) +
#        ggtitle("# of coefficients for x (u)")
#p

#demo_results_melted <- melt(demo_results, id.vars = c("x_penalty", "z_penalty"), measure.vars = c("u_coeffs", "v_coeffs"))
#demo_results_melted
#
#p <- ggplot(demo_results_melted, aes(x_penalty, z_penalty)) +
#        geom_tile(aes(fill = value)) +
#        scale_fill_gradient(low = "white", high = "steelblue", limits=c(0, max(demo_results_melted$value))) +
#        ggtitle("# of coefficients for each vector") +
#        facet_wrap(~variable)
#p

# -------- Test a real grid ------
grid = analyze_penalty_grid(c(0.001, seq(0.01, 0.1, by=0.01)), num_iter)
grid_melted = melt(grid, id.vars = c("x_penalty", "z_penalty"),
                   measure.vars = c("u_coeffs", "v_coeffs"))
grid_melted
# to make fonts bigger for poster
bold_text = theme(axis.text=element_text(size=16),
                  axis.title=element_text(size=20,face="bold")
                  )

ggplot(grid_melted, aes(x_penalty, z_penalty)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient(low = "white", high = "steelblue",
                            limits=c(0, max(grid_melted$value))) +
        ggtitle("Number of nonzero sparse CCA coefficients") +
        facet_wrap(~variable) +
        theme_bw() + bold_text
        #theme_bw()
dir.create('./plots/')
filename = sprintf('./plots/161207_num_weights_is_independent_of_other_penalty_%d_iter.pdf',
                   num_iter)
filename = sub("pdf", "png", filename)
print(filename)
ggsave(filename, height = 4, width = 8, dpi = 300) # a lot of poster printers are 300 dpi

# -------- Run same penalties for x and z wit higher resolution ------

analyze_penalty_list <- function(penalty_list){
        # do matching x and z penalties, since they don't appear to impact each other.
        results_coarse = data.frame()
        for (penalty in penalty_list){
                model <- CCA(x, z, typex="standard", typez="standard", K=1,
                             standardize = standardize,
                             niter=num_iter,
                             penaltyx=penalty, penaltyz=penalty)
                results_coarse <- rbind(results_coarse, model_stats(model))
        }
        return(results_coarse)
}

penalty_list = c(0.001, seq(0.01, 0.1, by=0.003))
print(penalty_list)
results = analyze_penalty_list(penalty_list)
results$penalty = results$x_penalty
tail(results)

results_melted =  melt(results, id.vars = c('penalty', "u_coeffs", "v_coeffs"),
                       measure.vars = c("u_coeffs", "v_coeffs"),
                       variable.name = 'regularization_penalty',
                       value.name = 'num_nonzero_weights')
tail(results_melted)

bold_text = theme(text = element_text(size=16),
                  title = element_text(size=18),
                  axis.title=element_text(face="bold") )
ggplot(results_melted, aes(penalty, num_nonzero_weights, color=regularization_penalty)) +
        geom_point() + geom_line() +
        ggtitle("tuning regularization for 83 samples") +
        geom_hline(yintercept=8) + theme_bw() + bold_text
filename = sprintf('./plots/161207_num_nonzero_weights_vs_penalty_for_full_training_set_%d_iter.pdf',
                   num_iter)
ggsave(filename)
filename = sub("pdf", "png", filename)
print(filename)
ggsave(filename, height = 4, width = 8)


## with N = 83, we need < ~8 weights
#penalty_list[0:10]
#
#colSums(x)
#rowSums(x)
#
#penalty_x = 0.0335  # 0.33 --> 7, 0.34 --> 9
#penalty_z = 0.022
#x2 = x[,apply(x, 2, var, na.rm=TRUE) != 0]
#z2 = z[,apply(z, 2, var, na.rm=TRUE) != 0]
#dim(x)
#dim(x2)
#
##final_model = CCA(x2, z2, typex="standard", typez="standard", K=1,
#final_model = CCA(x, z, typex="standard", typez="standard", K=1,
#                  niter=num_iter,
#                  standardize = standardize,
#                  penaltyx=penalty_x, penaltyz=penalty_z)
#final_results = model_stats(final_model)
#final_results
#
## Save these resutls
##write.table(final_model$u, file = './results/u_whole_training_set--8_features.csv', row.names = FALSE)
##write.table(final_model$v, file = './results/v_whole_training_set--8_features.csv', row.names = FALSE)
#
#print(final_model)
