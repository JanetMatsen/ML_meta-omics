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

dim(m) # [1]   83 4593
dim(nmm) # [1]   83 7355

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
# I got rid of them in the Python script, so these filters shouldn't be doing anything any longer.
dim(x)
dim(Filter(function(q) sd(q) != 0, x))
x = Filter(function(q) sd(q) != 0, x)
dim(z)
z = Filter(function(q) sd(q) != 0, z)
dim(z)


# first try
out <- CCA(x, z, typex="standard", typez="standard", K=1)
print(out,verbose=TRUE)

# Let the package determine good penalty values.
perm.out <- CCA.permute(x,z,typex="standard",typez="standard",nperms=7)

perm.out[0]
names(perm.out)
perm.out$bestpenaltyx
perm.out$bestpenaltyz

out_best_penalty <- CCA(x, z, typex="standard", typez="standard", K=1,
                        penaltyx=perm.out$bestpenaltyx, penaltyz=perm.out$bestpenaltyz)
dim(out_best_penalty$u)
dir.create('./results')

getwd()
write.table(out_best_penalty$u, file = './results/u_best_penalty.csv', row.names = FALSE)
write.table(out_best_penalty$v, file = './results/v_best_penalty.csv', row.names = FALSE)

# How sparse is the result for the recommended penalty?
qplot(out_best_penalty$u, geom="histogram")
qplot(out_best_penalty$v, geom="histogram")



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

model_stats(out_best_penalty)

#======  Loop over some different penalty values and find the number of zeros =========

test_pair_of_penalties <- function(penalty_list, x_penalty, z_penalty){
        print(paste("x penalty: " + x_penalty + "z penalty" + z))
        model = CCA(x, z, typex="standard", typez="standard", K=1,
                    penaltyx=x_penalty, penaltyz=z_penalty_list)
        return(model)}

analyze_penalty_grid <- function(penalty_list){
        results_coarse = data.frame()
        print(results_coarse)
        for (x_penalty in penalty_list){
                for (z_penalty in penalty_list){
                        model <- CCA(x, z, typex="standard", typez="standard", K=1,
                                     penaltyx=x_penalty, penaltyz=z_penalty)
                        results_coarse <- rbind(results_coarse, model_stats(model))
                }
        }
        return(results_coarse)
}

demo_results = analyze_penalty_grid(c(0.1, 0.2))
demo_results
p <- ggplot(demo_results, aes(x_penalty, z_penalty)) +
        geom_tile(aes(fill = u_coeffs)) +
        scale_fill_gradient(low = "white", high = "steelblue", limits=c(0, max(demo_results$u_coeffs))) +
        ggtitle("# of coefficients for x (u)")
p

library(reshape2)
demo_results_melted <- melt(demo_results, id.vars = c("x_penalty", "z_penalty"), measure.vars = c("u_coeffs", "v_coeffs"))
demo_results_melted

p <- ggplot(demo_results_melted, aes(x_penalty, z_penalty)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient(low = "white", high = "steelblue", limits=c(0, max(demo_results_melted$value))) +
        ggtitle("# of coefficients for each vector") +
        facet_wrap(~variable)
p

# -------- Test a real grid ------
grid = analyze_penalty_grid(c(0.001, seq(0.01, 0.05, by=0.01)))
grid_melted = melt(grid, id.vars = c("x_penalty", "z_penalty"),
                   measure.vars = c("u_coeffs", "v_coeffs"))
grid_melted
ggplot(grid_melted, aes(x_penalty, z_penalty)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient(low = "white", high = "steelblue",
                            limits=c(0, max(grid_melted$value))) +
        ggtitle("# of coefficients for each vector") +
        facet_wrap(~variable) + theme_bw()
dir.create('./plots/')
ggsave('./plots/161201_num_weights_is_independent_of_oter_penalty.pdf')

# -------- Run same penalties for x and z wit higher resolution ------

analyze_penalty_list <- function(penalty_list){
        # do matching x and z penalties, since they don't appear to impact each other.
        results_coarse = data.frame()
        for (penalty in penalty_list){
                model <- CCA(x, z, typex="standard", typez="standard", K=1,
                             penaltyx=penalty, penaltyz=penalty)
                results_coarse <- rbind(results_coarse, model_stats(model))
        }
        return(results_coarse)
}

penalty_list = c(0.001, seq(0.01, 0.05, by=0.003))
print(penalty_list)
results = analyze_penalty_list(penalty_list)
results$penalty = results$x_penalty
tail(results)

results_melted =  melt(results, id.vars = c('penalty', "u_coeffs", "v_coeffs"),
                       measure.vars = c("u_coeffs", "v_coeffs"),
                       variable.name = 'regularization_penalty',
                       value.name = 'num_nonzero_weights')
tail(results_melted)

ggplot(results_melted, aes(penalty, num_nonzero_weights,
                           color=regularization_penalty)) + geom_point() + geom_line() +
        geom_hline(yintercept=8)
ggsave('./plots/161201_num_nonzero_weights_vs_penalty_for_full_training_set.pdf')

# with N = 83, we need < ~8 weights
penalty_list[0:10]

penalty_x = 0.0335  # 0.33 --> 7, 0.34 --> 9
penalty_z = 0.022
final_model = CCA(x, z, typex="standard", typez="standard", K=1,
                  penaltyx=penalty_x, penaltyz=penalty_z)
final_results = model_stats(final_model)
final_results

# Save these resutls
write.table(final_model$u, file = './results/u_whole_training_set--8_features.csv', row.names = FALSE)
write.table(final_model$v, file = './results/v_whole_training_set--8_features.csv', row.names = FALSE)

print(final_model)
