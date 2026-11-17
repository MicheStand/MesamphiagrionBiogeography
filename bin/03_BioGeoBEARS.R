#03_BioGeoBEARS.R
#Create by Miguel Stand, 2026 - mstand20@gmail.com

#Load packages
library(purrr)
library(tidyverse)
library(ape)
library(optimx)   
library(GenSA)    
library(rexpokit)
library(cladoRcpp)
library(snow)   
library(parallel)
library(BioGeoBEARS) #need package Rcpp v1.1.0
library(ggtree)
library(deeptime)
library(paleotree)
library(trans4m8)

#Settings held constant across every model, for valid AICc comparison
trfn    <- "phylogeny.tre"
geogfn  <- "phylip_k4.txt"
distsfn <- "gradualdips_k4.txt"
tipranges       <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)
max_range_size  <- max(rowSums(dfnums_to_numeric(tipranges@df)))

include_null_range <- FALSE 
runslow <- TRUE
outdir  <- "k4"
dir.create(outdir, showWarnings = FALSE)

# ==== Base run-object constructor ====
build_base_run <- function(distsfn, dispersal_multipliers_fn = NULL) {
  ro <- define_BioGeoBEARS_run()
  ro$trfn               <- trfn
  ro$geogfn             <- geogfn
  ro$max_range_size     <- max_range_size
  ro$min_branchlength   <- 1e-6
  ro$include_null_range <- include_null_range
  ro$distsfn            <- distsfn
  if (!is.null(dispersal_multipliers_fn)) ro$dispersal_multipliers_fn <- dispersal_multipliers_fn
  ro$on_NaN_error        <- -1e50
  ro$speedup             <- TRUE
  ro$use_optimx          <- TRUE
  ro$num_cores_to_use    <- 4
  ro$force_sparse        <- FALSE
  ro                     <- readfiles_BioGeoBEARS_run(ro)
  ro$return_condlikes_table                <- TRUE
  ro$calc_TTL_loglike_from_condlikes_table <- TRUE
  ro$calc_ancprobs                         <- TRUE
  ro
}

#Model-family parameter setters
set_DEC_params <- function(ro) ro

set_DIVALIKE_params <- function(ro) {
  pt <- ro$BioGeoBEARS_model_object@params_table
  pt["s","type"]     <- "fixed"; pt["s","init"] <- 0;   pt["s","est"] <- 0
  pt["ysv","type"]   <- "2-j"
  pt["ys","type"]    <- "ysv*1/2"
  pt["y","type"]     <- "ysv*1/2"
  pt["v","type"]     <- "ysv*1/2"
  pt["mx01v","type"] <- "fixed"; pt["mx01v","init"] <- 0.5; pt["mx01v","est"] <- 0.5
  ro$BioGeoBEARS_model_object@params_table <- pt
  ro
}

set_BAYAREALIKE_params <- function(ro) {
  pt <- ro$BioGeoBEARS_model_object@params_table
  pt["s","type"]     <- "fixed"; pt["s","init"] <- 0; pt["s","est"] <- 0
  pt["v","type"]     <- "fixed"; pt["v","init"] <- 0; pt["v","est"] <- 0
  pt["ysv","type"]   <- "1-j"
  pt["ys","type"]    <- "ysv*1/1"
  pt["y","type"]     <- "1-j"
  pt["mx01y","type"] <- "fixed"; pt["mx01y","init"] <- 0.9999; pt["mx01y","est"] <- 0.9999
  ro$BioGeoBEARS_model_object@params_table <- pt
  ro
}

model_setters <- list(DEC         = set_DEC_params,
                      DIVALIKE    = set_DIVALIKE_params,
                      BAYAREALIKE = set_BAYAREALIKE_params)


#Nested "+J", now with family-specific j bounds
j_bounds <- list(DEC         = list(min = 0.00001, max = 2.99999),
                 DIVALIKE    = list(min = 0.00001, max = 1.99999),
                 BAYAREALIKE = list(min = 0.00001, max = 0.99999))


#NULL = leave default d/e bounds untouched (DEC+J, DIVALIKE+J)
de_bounds <- list(DEC = NULL,
                  DIVALIKE = NULL,
                  BAYAREALIKE = list(d = list(min = 0.0000001, max = 4.9999999), 
                                     e = list(min = 0.0000001, max = 4.9999999)))

#NULL = keep default use_optimx = TRUE (DEC+J, DIVALIKE+J)
optimizer_override <- list(DEC = NULL,
                           DIVALIKE = NULL,
                           BAYAREALIKE = "GenSA")

add_founder_event <- function(ro, base_res, family, jstart = 0.0001) {
  dstart <- base_res$outputs@params_table["d","est"]
  estart <- base_res$outputs@params_table["e","est"]
  bounds <- j_bounds[[family]]
  
  pt <- ro$BioGeoBEARS_model_object@params_table
  pt["d","init"] <- dstart; pt["d","est"] <- dstart
  pt["e","init"] <- estart; pt["e","est"] <- estart
  pt["j","type"] <- "free"; pt["j","init"] <- jstart; pt["j","est"] <- jstart
  pt["j","min"]  <- bounds$min; pt["j","max"] <- bounds$max
  
  #BAYAREALIKE+J only: tighten d/e bounds
  de <- de_bounds[[family]]
  if (!is.null(de)) {
    pt["d","min"] <- de$d$min; pt["d","max"] <- de$d$max
    pt["e","min"] <- de$e$min; pt["e","max"] <- de$e$max
  }
  ro$BioGeoBEARS_model_object@params_table <- pt
  
  #BAYAREALIKE+J only: switch optimizer to GenSA
  opt <- optimizer_override[[family]]
  if (!is.null(opt)) ro$use_optimx <- opt
  
  ro
}

#Fit function
fit_or_load <- function(ro, resfn, runslow = TRUE) {
  if (runslow) {
    res <- bears_optim_run(ro)
    save(res, file = resfn)
  } else {
    load(resfn)
  }
  res
}

#Fit and validate the three base models
base_runs <- purrr::imap(model_setters, function(setter, label) {
  ro <- setter(build_base_run(distsfn))
  ro <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = ro)
  check_BioGeoBEARS_run(ro)
  ro
})

resfn_base <- list(DEC         = file.path(outdir, "DEC_k4.RData"),
                   DIVALIKE    = file.path(outdir, "DIVALIKE_k4.RData"),
                   BAYAREALIKE = file.path(outdir, "BAYAREALIKE_k4.RData"))

results_base <- purrr::imap(base_runs, function(ro, label) {
  fit_or_load(ro, resfn_base[[label]], runslow = runslow)
})

#Fit and validate the nested "+J" models
J_runs <- purrr::imap(base_runs, function(ro, label) {
  ro_j <- add_founder_event(ro, results_base[[label]], family = label)
  ro_j <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = ro_j)
  check_BioGeoBEARS_run(ro_j)
  ro_j
})

resfn_J <- list(DEC         = file.path(outdir, "DEC_k4j.RData"),
                DIVALIKE    = file.path(outdir, "DIVALIKE_k4j.RData"),
                BAYAREALIKE = file.path(outdir, "BAYAREALIKE_k4j.RData"))

results_J <- purrr::imap(J_runs, function(ro, label) {
  fit_or_load(ro, resfn_J[[label]], runslow = runslow)
})


#Store the six models together
models <- list(DEC = list(base = results_base$DEC,
                          J    = results_J$DEC),
               DIVALIKE = list(base = results_base$DIVALIKE,
                               J    = results_J$DIVALIKE),
               BAYAREALIKE = list(base = results_base$BAYAREALIKE,
                                  J    = results_J$BAYAREALIKE))

#Pairwise LRT comparison
compare_base_vs_J <- function(base_res, J_res) {
  LnL_null <- get_LnL_from_BioGeoBEARS_results_object(base_res)
  LnL_alt  <- get_LnL_from_BioGeoBEARS_results_object(J_res)
  stats <- AICstats_2models(LnL_alt, LnL_null, numparams1 = 3, numparams2 = 2)
  
  null_row <- extract_params_from_BioGeoBEARS_results_object(
    results_object = base_res, returnwhat = "table", addl_params = "j", paramsstr_digits = 4)
  alt_row <- extract_params_from_BioGeoBEARS_results_object(
    results_object = J_res, returnwhat = "table", addl_params = "j", paramsstr_digits = 4)
  
  list(params = rbind(null_row, alt_row), stats = conditional_format_table(stats))
}

comparisons <- purrr::imap(model_setters, function(setter, label) {
  compare_base_vs_J(results_base[[label]], results_J[[label]])
})

restable   <- purrr::map_dfr(comparisons, "params")
teststable <- purrr::map_dfr(comparisons, "stats")

family_labels          <- names(model_setters)
teststable$alt         <- paste0(family_labels, "+J")
teststable$null        <- family_labels
row.names(restable)    <- as.vector(rbind(family_labels, paste0(family_labels, "+J")))
restable               <- put_jcol_after_ecol(restable)

restable
teststable

#Save intermediate tables 
#save(restable, file = file.path(outdir, "restable_v1.Rdata"))
#save(teststable, file = file.path(outdir, "teststable_v1.Rdata"))
write.table(restable, file.path(outdir, "restable.txt"), 
            quote = FALSE, sep = "\t")
write.table(unlist_df(teststable), file.path(outdir, "teststable.txt"), 
            quote = FALSE, sep = "\t")


#AICc and Akaike weights across all six models
samplesize <- length(ape::read.tree(trfn)$tip.label)
restable <- cbind(restable, calc_AIC_column(LnL_vals = restable$LnL, nparam_vals = restable$numparams))
restable <- cbind(restable, calc_AICc_column(LnL_vals = restable$LnL, nparam_vals = restable$numparams,
                                             samplesize = samplesize))

criteria <- c("AIC", "AICc")
rellike_tables <- purrr::imap(setNames(criteria, criteria), function(colname, label) {
  tbl <- AkaikeWeights_on_summary_table(restable, colname_to_use = colname)
  tbl <- put_jcol_after_ecol(tbl)
  write.table(tbl, file.path(outdir, paste0("restable_", label, "_rellike.txt")),
              quote = FALSE, sep = "\t")
  tbl
})

#Best-supported model under AICc, ranked
rellike_tables$AICc %>% arrange(desc(AICc_wt))

#Plot the best model
#Select best-supported model (highest AICc weight)
best_model <- rownames(rellike_tables$AICc %>%
                         arrange(desc(AICc_wt)) %>%
                         slice(1))
best_model

#Retrieve corresponding BioGeoBEARS results object
best_res <-switch(best_model,
                  "DEC"          = models$DEC$base,
                  "DEC+J"        = models$DEC$J,
                  "DIVALIKE"     = models$DIVALIKE$base,
                  "DIVALIKE+J"   = models$DIVALIKE$J,
                  "BAYAREALIKE"  = models$BAYAREALIKE$base,
                  "BAYAREALIKE+J"= models$BAYAREALIKE$J)

#Plot
#Read tree
tr <- read.tree(trfn)
#Extract ancestral areas
states_list <- best_res$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
state_names <- sapply(states_list, function(x){
  if(length(x)==0) return("NULL")
  paste(c("A","B","C","D")[x+1], collapse="")
})
#Probabilities
probs <- best_res$ML_marginal_prob_each_state_at_branch_top_AT_node
node_probs <- probs[(Ntip(tr)+1):(Ntip(tr)+Nnode(tr)), ]
colnames(node_probs) <- state_names
write.csv(node_probs, "node_probs.csv")
#Now, manually check the node_probs and select the areas with probabilities > 0.15 (just for aesthetic reasons)
cols <- setNames(rep("gray80", length(state_names)), state_names) #Gray for areas with low probabilities
#Color for areas with high probabilities
cols[c("A", "B", "C", "D", "AB", "AD", "BC", "BD", "ABD")] <- c(
  "#C44E52",  # A
  "#DDCC77",  # B
  "#7FA58A",  # C
  "#4C72B0",  # D
  "#556B7A",  # AB
  "#8172B2",  # AD
  "#8C613C",  # BC
  "#64B5CD",  # BD
  "#9E716B")  # ABD

svg("../out/ancestral_areas.svg")
plot.phylo(tr)
nodelabels(pie = node_probs,
           piecol = cols[colnames(node_probs)],
           cex = 0.8,
           frame = "none")
dev.off()

#Save RData 
save.image(file = "03_BioGeoBEARS.RData")
#end :D