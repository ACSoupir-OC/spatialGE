##
# STdiff Helper Functions
# Utility functions for differential expression analysis

##
# count_cores = function(n){
#' @title count_cores: Determine optimal number of cores for parallelization
#' @description Determines the number of cores to use for parallelization
#' @param n an integer representing the number of units to process
#' @return the number of cores to use (max = half of available cores)
#' @keywords internal
#' @importFrom parallel detectCores
count_cores = function(n){
  cores = 1
  if(.Platform$OS.type == 'unix'){
    # Use parallelization (if possible) to read data.
    avail_cores = (parallel::detectCores()) / 2
    if(avail_cores <= n){
      cores = avail_cores
    } else{
      cores = ifelse(n > 0, n, 1)
    }
  }
  return(cores)
}


##
#' @title raise_err: Raise error from exception database
#' @description Stops with an error message from the spatialGE exception database
#' @param err_code error code to look up
#' @param samplename optional sample name to include in error message
#' @keywords internal
#' @importFrom utils read.csv
raise_err = function(err_code=NULL, samplename=NULL){
  pkg_fp = list.files(system.file(package='spatialGE'), full.names=T, pattern='err\\.csv')
  err_db = suppressWarnings(utils::read.csv(pkg_fp, header=F, quote="'"))
  str_to_print = err_db[[2]][ err_db[[1]] == err_code ]
  if(is.null(samplename)){
    stop(paste0('spatialGE exception: ', str_to_print), call.=F)
  } else{
    stop(paste0('spatialGE exception: ', str_to_print, ' Sample: ', samplename, '.'), call.=F)
  }
}


##
#' @title expandSparse: Convert sparse matrix to dense data frame
#' @description Expands a sparse matrix to a dense data frame
#' @param sparsedMatrix a sparse matrix to convert
#' @return a dense data frame
#' @keywords internal
expandSparse = function(sparsedMatrix){
  NonSparse = data.frame(as.matrix(sparsedMatrix), check.names=F)
  return(NonSparse)
}


##
# non_spatial_de
#' @title non_spatial_de: Run non-spatial linear models for gene x cluster combinations
#' @description Runs non-spatial linear models (or t-tests/Wilcoxon) for testing
#' differentially expressed genes between groups of spots/cells
#' @param expr_data a data frame (rows=spots/cells) with normalized expression data
#' as well as coordinates (x,y), and a group variable for random effects.
#' @param combo a data frame with columns (samplename, meta1, meta2, gene) containing
#' all combinations of gene x cluster to test.
#' @param pairwise whether to perform pairwise tests (TRUE) or reference-based (FALSE)
#' @return a list containing non-spatial model results with p-values, log-fold changes,
#' and metadata
#' @keywords internal
#' @importFrom stats as.formula p.adjust
non_spatial_de = function(expr_data=NULL, combo=NULL, pairwise=NULL){
  models_ls = list()
  for(i in 1:nrow(combo)){
    sample_tmp = combo[i, 1]
    meta1_tmp = combo[i, 2]
    meta2_tmp = combo[i, 3]
    gene_tmp = combo[i, 4]

    # Filter data if pairwise or recode annotations if not pairwise
    if(pairwise){
      expr_tmp = expr_data %>%
        dplyr::filter(meta %in% c(meta1_tmp, meta2_tmp))
    } else{
      expr_tmp = expr_data %>%
        dplyr::mutate(meta=factor(dplyr::case_when(meta != meta1_tmp ~ 'other', TRUE ~ meta),
                                  levels=c('other', meta1_tmp)))
    }

    # Get relevant gene column
    expr_tmp = expr_tmp %>%
      dplyr::select(c('group', 'ypos', 'xpos', 'meta'), exprval:=!!gene_tmp)

    # Calculate log-fold change
    avgexpr_cl1 = mean(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select(exprval) %>% unlist())
    avgexpr_cl2 = mean(expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select(exprval) %>% unlist())
    avglogfold_tmp = (avgexpr_cl1 - avgexpr_cl2)

    # Create non-spatial model
    error_message = c()
    warning_message = c()
    res_mod = withCallingHandlers(
      tryCatch(spaMM::fitme(formula=stats::as.formula('exprval~meta'), data=expr_tmp, method="REML"),
               error=function(e){error_message <<- conditionMessage(e)}),
      warning=function(w){warning_message <<- conditionMessage(w)
      invokeRestart("muffleWarning")}
    )

    # Calculate p-value
    res_mod = spaMM::summary.HLfit(res_mod, details=c(p_value="Wald"), verbose=F)[['beta_table']][2, 4]

    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]] = list()
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['nonspmod']] = res_mod
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['samplename']] = sample_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['gene']] = gene_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['avglfc']] = avglogfold_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta1']] = meta1_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta2']] = meta2_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['warning']] = warning_message
  }
  return(models_ls)
}


##
# spatial_de
#' @title spatial_de: Run spatial linear models with Matern covariance
#' @description Fits spatial mixed models with Matern covariance structure to test
#' for spatial autocorrelation in gene expression between clusters
#' @param expr_dat a data frame with expression data and coordinates
#' @param non_sp_mods a list of non-spatial models to update with covariance structures
#' @param annot_dict a data frame mapping coded annotations to original annotation names
#' @param verb verbosity level (0, 1, or 2)
#' @return a list of spatial models with convergence status
#' @keywords internal
#' @importFrom stats as.formula
spatial_de = function(expr_dat=NULL, non_sp_mods=NULL, annot_dict=NULL, verb=NULL){

  # To prevent NOTES in R CMD check
  . = NULL

  res_ls = list()
  for(i in names(non_sp_mods)){
    # Extract sample name
    sample_tmp = non_sp_mods[[i]][['samplename']]
    # Extract gene name
    gene_tmp = non_sp_mods[[i]][['gene']]
    # Extract clusters
    meta1_tmp = non_sp_mods[[i]][['meta1']]
    meta2_tmp = non_sp_mods[[i]][['meta2']]

    # Print information of test
    if(verb == 2L){
      if(grepl('_c[0-9]+_other$', i)){
        stdout_print = i %>%
          stringr::str_extract(., 'c[0-9]+_other$') %>%
          stringr::str_replace(., '_other$', '')
        stdout_print = annot_dict[['orig_annot']][ annot_dict[['coded_annot']] == stdout_print ]
      } else{
        stdout_print = i %>%
          stringr::str_extract(., paste0('c[0-9]+_c[0-9]+$')) %>%
          stringr::str_split(., '_') %>% unlist() %>% unique()
        stdout_print = c(annot_dict[['orig_annot']][ annot_dict[['coded_annot']] == stdout_print[1] ],
                         annot_dict[['orig_annot']][ annot_dict[['coded_annot']] == stdout_print[2] ])
        stdout_print = paste0(stdout_print[1], ' vs. ', stdout_print[2])
      }
      system(sprintf('echo "%s"', paste0('\t\tSample: ', sample_tmp, '; ', gene_tmp, ', ', stdout_print)))
    }

    expr_subset = expr_dat[, c("meta" , "group", "ypos", "xpos", gene_tmp)] %>%
      dplyr::rename(exprval := !!gene_tmp)
    if(meta2_tmp == 'other'){
      expr_subset[['meta']][ expr_subset[['meta']] != meta1_tmp ] = 'other'
    } else{
      expr_subset = expr_subset[ expr_subset[['meta']] %in% c(meta1_tmp, meta2_tmp), ]
    }

    exp_out = tryCatch({
      spaMM::fitme(formula=stats::as.formula(paste0("exprval~meta+Matern(1|xpos+ypos)")),
                   data=expr_subset,
                   fixed=list(nu=0.5), method="REML",
                   control.HLfit=list(algebra="decorr"))
    }, error=function(err){return(err)})

    res_ls[[i]] = list()
    if(any(class(exp_out) == 'TimeoutException')){
      res_ls[[i]][['spmod']] = 'time_out'
    } else if(any(class(exp_out) == 'simpleError')){
      res_ls[[i]][['spmod']] = 'no_conv'
    } else{
      res_ls[[i]][['spmod']] = exp_out
    }
    res_ls[[i]][['sample']] = sample_tmp
    res_ls[[i]][['gene']] = gene_tmp
  }
  return(res_ls)
}


##
# prepare_stdiff_combo
#' @title prepare_stdiff_combo: Create combinations for differential testing
#' @description Creates a data frame with all combinations of genes x cluster pairs
#' to test for differential expression
#' @param to_expand a data frame with columns (samplename, orig_annot, gene, meta)
#' containing all gene x annotation combinations per sample
#' @param user_clusters optional vector of specific clusters to test
#' @param pairwise whether to perform pairwise tests (TRUE) or reference-based (FALSE)
#' @param verbose verbosity level for progress messages
#' @return a data frame with all combinations to test (samplename, meta1, meta2, gene)
#' @keywords internal
prepare_stdiff_combo = function(to_expand=NULL, user_clusters=NULL, pairwise=NULL, verbose=NULL){
  combo_df = tibble::tibble()
  for(sample_name in unique(to_expand[['samplename']])){
    # Extract annotations for a given sample
    to_expand_subset = to_expand %>% dplyr::filter(samplename == sample_name)

    # Define clusters to test if NULL, or subset to those requested by user
    if(is.null(user_clusters)){
      annots_tmp = to_expand_subset %>% dplyr::select('meta') %>% unlist() %>% unique()
    } else{
      annots_tmp = to_expand_subset %>% dplyr::filter(orig_annot %in% user_clusters) %>% dplyr::select('meta') %>% unlist() %>% unique()
    }
    
    if(length(annots_tmp) >= 2){
      if(pairwise){
        # Get unique combinations of each two annotation, without repeating
        combo_meta = tibble::tibble()
        annots_tmp2 = annots_tmp
        for(cl_tmp in annots_tmp){
          combo_meta = dplyr::bind_rows(combo_meta,
                                        expand.grid(cl_tmp,
                                                    grep(cl_tmp, annots_tmp2, value=T, invert=T),
                                                    stringsAsFactors=F))
          annots_tmp2 = grep(cl_tmp, annots_tmp2, value=T, invert=T)
        }
        names(combo_meta) = c('meta1', 'meta2')
        
        # Add genes to each combination
        combo_meta_gene = tibble::tibble()
        for(comb in 1:nrow(combo_meta)){
          combo_meta_gene = dplyr::bind_rows(combo_meta_gene,
                                             tibble::tibble(combo_meta[comb, 1],
                                                            combo_meta[comb, 2],
                                                            gene=to_expand %>% dplyr::filter(samplename == sample_name) %>%
                                                              dplyr::select(gene) %>% unlist() %>% unique()))
        }
        combo_meta = combo_meta_gene
      } else{
        combo_meta = expand.grid(meta1=annots_tmp, meta2='other',
                                 gene=to_expand %>% dplyr::filter(samplename == sample_name) %>%
                                   dplyr::select(gene) %>% unlist() %>% unique(),
                                 stringsAsFactors=F)
      }
      
      combo_df = dplyr::bind_rows(combo_df,
                                  combo_meta %>%
                                    tibble::add_column(samplename=sample_name, .before=1))
    } else{
      if(verbose >= 1L){
        cat(paste0('\t\tSkipping sample ', sample_name, '. Less than two clusters to compare.\n'))
      }
    }
  }
  combo_df = as.data.frame(combo_df) %>%
    dplyr::arrange(samplename, meta1, meta2, gene)

  return(combo_df)
}


##
# stdiff_mean_test
#' @title stdiff_mean_test: Perform t-test or Wilcoxon test for differential expression
#' @description Performs simple mean-based tests (t-test or Wilcoxon) between
#' groups of spots/cells for a given gene
#' @param expr_data a data frame with expression data and metadata
#' @param combo a data frame with columns (samplename, meta1, meta2, gene)
#' @param pairwise whether to perform pairwise tests (TRUE) or reference-based (FALSE)
#' @param test_type type of test: 't_test' or 'wilcoxon'
#' @return a list containing test results with p-values, log-fold changes, and metadata
#' @keywords internal
stdiff_mean_test = function(expr_data=NULL, combo=NULL, pairwise=NULL, test_type=NULL){
  test_ls = list()
  for(i in 1:nrow(combo)){
    sample_tmp = combo[i, 1]
    meta1_tmp = combo[i, 2]
    meta2_tmp = combo[i, 3]
    gene_tmp = combo[i, 4]

    # Filter data if pairwise or recode annotations if not pairwise
    if(pairwise){
      expr_tmp = expr_data %>%
        dplyr::filter(meta %in% c(meta1_tmp, meta2_tmp))
    } else{
      expr_tmp = expr_data %>%
        dplyr::mutate(meta=factor(dplyr::case_when(meta != meta1_tmp ~ 'other', TRUE ~ meta),
                                  levels=c('other', meta1_tmp)))
    }

    # Get relevant gene column
    expr_tmp = expr_tmp %>%
      dplyr::select(c('group', 'ypos', 'xpos', 'meta'), exprval:=!!gene_tmp)

    # Calculate log-fold change
    avgexpr_cl1 = mean(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select(exprval) %>% unlist())
    avgexpr_cl2 = mean(expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select(exprval) %>% unlist())
    avglogfold_tmp = (avgexpr_cl1 - avgexpr_cl2)

    # Perform test for differences in mean
    if(test_type == 't_test'){
      res_test = stats::t.test(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select('exprval') %>% unlist(),
                               expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select('exprval') %>% unlist())
    } else if(test_type == 'wilcoxon'){
      res_test = stats::wilcox.test(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select('exprval') %>% unlist(),
                                    expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select('exprval') %>% unlist())
    }

    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]] = list()
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['mean_test']] = res_test
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['samplename']] = sample_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['gene']] = gene_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['avglfc']] = avglogfold_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta1']] = meta1_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta2']] = meta2_tmp
  }
  return(test_ls)
}