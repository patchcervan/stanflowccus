#' Configure a Workflow for Occupancy Modelling with Stan
#'
#' This function sets up the configuration for an occupancy modelling workflow,
#' specifying directories, species selection, years, covariates, and model details.
#'
#' @param data_dir Character. The directory where data should be looked for.
#' @param out_dir Character. The output directory where results will be saved.
#' @param years Numeric vector. The set of years for the analysis.
#' @param train_years Numeric vector. The subset of years used for model training.
#' @param data_trim Character. Trimming method for data: `"none"`, `"end"`, or `"start"`.
#' @param site_ids Character vector. Variables that serve as unique identifiers for sites.
#' @param visit_ids Character vector. Variables that serve as unique identifiers for visits.
#' @param site_covts Character vector. Site-level covariates.
#' @param visit_covts Character vector. Visit-level covariates.
#' @param stan_model_dir Character. Path to the Stan model directory.
#' @param stan_model_file Character. Name of the Stan model file (optional).
#' @param param_to_plot Character vector. Parameters to plot in diagnostics.
#' @param param_extra Character vector. Additional parameters to monitor.
#' @param dynamic Logical. If `TRUE`, uses a dynamic occupancy model.
#' @param occu_temp_re Character. Temporal random effect structure for occupancy component: `"none"` or `"AR1"`.
#' @param det_temp_re Character. Temporal random effect structure for detection component: `"none"` or `"AR1"`
#' @param prefile Character. String used as a prefix for naming files (e.g. model fits, diagnostics) (optional).
#' @param postfile Character. String used as a postfix for naming files (e.g. model fits, diagnostics) (optional).
#' @param occu_models List of formulas. Occupancy model formula/s.
#' @param phi_models List of formulas. Survival (persistence) model formula/s (for dynamic models).
#' @param gamma_models List of formulas. Colonization model formula/s (for dynamic models).
#' @param det_models List of formulas. Detection model formula/s.
#' @param verbose Logical. If TRUE all configuration parameters are printed. Defaults to FALSE.
#'
#' @details
#' The lists `occu_models`, `phi_models`, `gamma_models` (if they exist) must have the same
#' length. In the case of dynamic models `occu_models` will be used for initial
#' occupancy, `phi_models` will be used for persistance and `gamma_models` will be
#' used for colonization. Models will be paired elementwise, ie., `occu_models[[i]]`,
#' `phi_models[[i]]`, and `gamma_models[[i]]`, for any `i`, will be used for the same dynamic
#' occupancy model. If the model is not dynamic, then `occu_models[[i]]` will be
#' used for occupancy, and the rest will be ignored. Models will be named according
#' to the names of `occu_models` or by a numeric sequence if `occu_models` is unnamed.
#'
#' Temporal random effects defined by `occu_temp_re` affect colonization and survival
#' components of dynamic models and occupancy component of non-dynamic models. Both
#' `occu_temp_re` and `det_temp_re` are ignored if a `stan_model_file` is defined.
#'
#'
#' @return A list containing the workflow configuration with the following components:
#' \itemize{
#'   \item `data_dir`: Data directory.
#'   \item `out_dir`: Output directory.
#'   \item `years`: Years included.
#'   \item `train_years`: Training years.
#'   \item `data_trim`: Data trimming method.
#'   \item `site_ids`: Site identifiers.
#'   \item `visit_ids`: Visit identifiers.
#'   \item `site_covts`: Site covariates.
#'   \item `visit_covts`: Visit covariates.
#'   \item `stan_model`: Stan model file path.
#'   \item `param_to_plot`: Parameters for plotting.
#'   \item `param_extra`: Extra parameters.
#'   \item `dynamic`: Dynamic model flag.
#'   \item `prefile`: File prefix.
#'   \item `postfile`: File postfix.
#'   \item `occu_models`: Occupancy models.
#'   \item `phi_models`: Persistence models.
#'   \item `gamma_models`: Colonization models.
#'   \item `det_models`: Detection models.
#'   \item `raw_data_file`: Raw data file path.
#'   \item `clean_data_file`: Clean and standardized data file path.
#'   \item `occu_data_file`: Occupancy data list file path.
#'   \item `test_data_file`: Occupancy test data list path.
#'   \item `test_fit_file`: Occupancy test-data fit file path.
#'   \item `test_diags_file`: Occupancy test-data fit diagnostics file path.
#'   \item `test_sims_file`: File path to occupancy estimates/predictions from test-data fit.
#' }
#'
#' @examples
#' \dontrun{
#' # Example configuration
#' config <- configWorkflow(
#'   outdir = "results/species1",
#'   years = 2010:2020,
#'   trainyears = 2010:2018,
#'   datatrim = "none",
#'   site_ids = sites$site_id,
#'   visit_ids = visits$visit_id,
#'   site_covts = site_covariates,
#'   visit_covts = visit_covariates,
#'   stanmodel = "occupancy.stan",
#'   paramtoplot = c("psi", "p"),
#'   paramextra = c("gamma", "phi"),
#'   dynamic = TRUE,
#'   prefile = "preprocess.R",
#'   postfile = "postprocess.R",
#'   occu_models = c("~1", "~covariate1"),
#'   phi_models = c("~1", "~covariate2"),
#'   gamma_models = c("~1", "~covariate3"),
#'   det_models = c("~1", "~covariate4")
#' )
#' }
#'
#' @export
configWorkflow <- function(data_dir,
                           out_dir,
                           years,
                           train_years,
                           data_trim = c("none", "end", "start"),
                           site_ids,
                           visit_ids,
                           site_covts,
                           visit_covts,
                           stan_model_dir,
                           stan_model_file = NULL,
                           param_to_plot,
                           param_extra,
                           dynamic = TRUE,
                           occu_temp_re = c("none", "AR1"),
                           det_temp_re = c("none", "AR1"),
                           prefile = NULL,
                           postfile = NULL,
                           occu_models,
                           phi_models = NULL,
                           gamma_models = NULL,
                           det_models,
                           verbose = FALSE){

    data_trim <- match.arg(data_trim)
    occu_temp_re <- match.arg(occu_temp_re)
    det_temp_re <- match.arg(det_temp_re)

    # Join parameters to plot and summarise with other parameters of interest
    paramall <- c(param_to_plot, param_extra)



    # Define file names -------------------------------------------------------

    years_ch <- paste(substr(range(years), 3, 4), collapse = "_")

    if(!is.null(prefile)){
        prefile <- paste0("_", prefile, "_")
    }

    if(!is.null(postfile)){
        postfile <- paste0("_", postfile, "_")
    }

    if(data_trim != "none"){
        trimtext <- paste0("_trim", data_trim)
    } else {
        trimtext <- ""
    }

    # Define model data file name
    datafile <- file.path(data_dir,
                          paste0("model_data", prefile, postfile, ".rds"))

    # Define clean data file
    cleandatafile <- file.path(out_dir,
                               paste0("model_data", prefile, postfile, trimtext, ".rds"))

    # Define occupancy data list file name
    occufile <- file.path(out_dir,
                          paste0("occu_data", prefile, postfile, trimtext, ".rds"))

    # Define occupancy test data list file name
    testfile <- file.path(out_dir,
                          paste0("occu_test_data", prefile, postfile, trimtext, ".rds"))

    # Define test model fit
    testfitfile <- rep(NA, length = length(occu_models)*length(det_models))
    testpredsfile <- rep(NA, length = length(occu_models)*length(det_models))
    testdiagsfile <- rep(NA, length = length(occu_models)*length(det_models))
    testsimsfile <- rep(NA, length = length(occu_models)*length(det_models))

    k <- 1
    for(i in seq_along(occu_models)){
        for(j in seq_along(det_models)){

            if(is.null(names(occu_models))){
                ii <- i
            } else {
                ii <- names(occu_models)[i]
            }

            if(is.null(names(det_models))){
                jj <- j
            } else {
                jj <- names(det_models)[j]
            }

            testfitfile[k] <- file.path(out_dir,
                                        paste0("test_fit", prefile, postfile,
                                               "_occu", ii, "_det", jj, trimtext, ".rds"))

            # Define file for saving diagnostics
            testdiagsfile[k] <- file.path(out_dir,
                                      paste0("diags_stan", prefile, postfile,
                                             "_occu", ii, "_det", jj, trimtext, ".rds"))

            # Define file for saving estimates/predictions
            testsimsfile[k] <- file.path(out_dir,
                                     paste0("sims_stan", prefile, postfile,
                                            "_occu", ii, "_det", jj, trimtext, ".rds"))

            k = k+1
        }
    }


    # Create output list
    setConfig(data_dir = data_dir,
              out_dir = out_dir,
              years = years,
              years_ch = years_ch,
              train_years = train_years,
              data_trim = data_trim,
              site_ids = site_ids,
              visit_ids = visit_ids,
              site_covts = site_covts,
              visit_covts = visit_covts,
              stan_model = file.path(stan_model_dir, stan_model_file),
              param_to_plot = param_to_plot,
              param_all = paramall,
              dynamic = dynamic,
              prefile = prefile,
              postfile = postfile,
              occu_models = occu_models,
              phi_models = phi_models,
              gamma_models = gamma_models,
              det_models = det_models,
              raw_data_file = datafile,
              clean_data_file = cleandatafile,
              occu_data_file = occufile,
              test_data_file = testfile,
              test_fit_file = testfitfile,
              # Add all-data file
              test_diags_file = testdiagsfile,
              test_sims_file = testsimsfile)


    if(verbose){
        message("Analysis configuration set to:")
        print(getConfig())
    } else {
        message("Analysis configuration set. Run getConfig() to retrieve config parameters")
    }

}

