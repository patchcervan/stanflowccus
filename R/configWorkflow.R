#' Configure a Workflow for Occupancy Modelling with Stan
#'
#' This function sets up the configuration for an occupancy modelling workflow,
#' specifying directories, species selection, years, covariates, and model details.
#'
#' @param outdir Character. The output directory where results will be saved.
#' @param sp_sel Character or numeric. The species to be analyzed (name or ID).
#' @param years Numeric vector. The set of years for the analysis.
#' @param trainyears Numeric vector. The subset of years used for model training.
#' @param datatrim Character. Trimming method for data: `"none"`, `"end"`, or `"start"`.
#' @param site_ids Character vector. Variables that serve as unique identifiers for sites.
#' @param visit_ids Character vector. Variables that serve as unique identifiers for visits.
#' @param site_covts Character vector. Site-level covariates.
#' @param visit_covts Character vector. Visit-level covariates.
#' @param stanmodel Character. Path to the Stan model file or model name.
#' @param paramtoplot Character vector. Parameters to plot in diagnostics.
#' @param paramextra Character vector. Additional parameters to monitor.
#' @param dynamic Logical. If `TRUE`, uses a dynamic occupancy model.
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
#'
#' @return A list containing the workflow configuration with the following components:
#' \itemize{
#'   \item `outdir`: Output directory.
#'   \item `sp_sel`: Selected species.
#'   \item `years`: Years included.
#'   \item `trainyears`: Training years.
#'   \item `datatrim`: Data trimming method.
#'   \item `site_ids`: Site identifiers.
#'   \item `visit_ids`: Visit identifiers.
#'   \item `site_covts`: Site covariates.
#'   \item `visit_covts`: Visit covariates.
#'   \item `stanmodel`: Stan model specification.
#'   \item `paramtoplot`: Parameters for plotting.
#'   \item `paramextra`: Extra parameters.
#'   \item `dynamic`: Dynamic model flag.
#'   \item `prefile`: File prefix.
#'   \item `postfile`: File postfix.
#'   \item `occu_models`: Occupancy models.
#'   \item `phi_models`: Persistence models.
#'   \item `gamma_models`: Colonization models.
#'   \item `det_models`: Detection models.
#' }
#'
#' @examples
#' \dontrun{
#' # Example configuration
#' config <- configWorkflow(
#'   outdir = "results/species1",
#'   sp_sel = "species1",
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
configWorkflow <- function(outdir,
                           sp_sel,
                           years,
                           trainyears,
                           datatrim = c("none", "end", "start"),
                           site_ids,
                           visit_ids,
                           site_covts,
                           visit_covts,
                           stanmodel,
                           paramtoplot,
                           paramextra,
                           dynamic = TRUE,
                           prefile = NULL,
                           postfile = NULL,
                           occu_models,
                           phi_models = NULL,
                           gamma_models = NULL,
                           det_models,
                           verbose = FALSE){

    datatrim <- match.arg(datatrim)

    # Join parameters to plot and summarise with other parameters of interest
    paramall <- c(paramtoplot, paramextra)



    # Define file names -------------------------------------------------------

    years_ch <- paste(substr(range(years), 3, 4), collapse = "_")

    if(!is.null(prefile)){
        prefile <- paste0("_", prefile, "_")
    }

    if(!is.null(postfile)){
        postfile <- paste0("_", postfile, "_")
    }

    # Define model data file name
    datafile <- file.path(outdir, sp_sel,
                          paste0("model_data", prefile, sp_sel, postfile, ".rds"))

    # Define occupancy data list file name
    occufile <- file.path(outdir, sp_sel,
                          paste0("occu_data", prefile, sp_sel, postfile, ".rds"))

    # Define occupancy test data list file name
    testfile <- file.path(outdir, sp_sel,
                          paste0("occu_test_data", prefile, sp_sel, postfile, ".rds"))

    # Define test model fit
    testfitfile <- rep(NA, length = length(occu_models)*length(det_models))
    predsfile <- rep(NA, length = length(occu_models)*length(det_models))
    diagsfile <- rep(NA, length = length(occu_models)*length(det_models))
    simsfile <- rep(NA, length = length(occu_models)*length(det_models))

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

            testfitfile[k] <- file.path(outdir, sp_sel,
                                        paste0("test_fit", prefile, sp_sel, postfile,
                                               "_occu", ii, "_det", jj, "_trim", datatrim, ".rds"))
            # Define file for saving predictions
            predsfile[k] <- file.path(outdir, sp_sel,
                                      paste0("preds_stan", prefile, sp_sel, postfile,
                                             "_occu", ii, "_det", jj, "_trim", datatrim, ".rds"))

            # Define file for saving predictions
            diagsfile[k] <- file.path(outdir, sp_sel,
                                      paste0("diags_stan", prefile, sp_sel, postfile,
                                             "_occu", ii, "_det", jj, "_trim", datatrim, ".rds"))

            # Define file for saving simulations
            simsfile[k] <- file.path(outdir, sp_sel,
                                     paste0("sims_stan", prefile, sp_sel, postfile,
                                            "_occu", ii, "_det", jj, "_trim", datatrim, ".rds"))

            k = k+1
        }
    }


    # Create output list
    setConfig(outdir = outdir,
              sp_sel = sp_sel,
              years = years,
              years_ch = years_ch,
              train_years = trainyears,
              datatrim = datatrim,
              site_ids = site_ids,
              visit_ids = visit_ids,
              site_covts = site_covts,
              visit_covts = visit_covts,
              stan_model = stanmodel,
              param_to_plot = paramtoplot,
              param_all = paramall,
              dynamic = dynamic,
              prefile = prefile,
              postfile = postfile,
              datafile = datafile,
              occufile = occufile,
              testfile = testfile,
              testfitfile = testfitfile,
              predsfile = predsfile,
              diagsfile = diagsfile,
              simsfile = simsfile)

    if(verbose){
        message("Analysis configuration set to:")
        message(getConfig())
    } else {
        message("Analysis configuration set. Run getConfig() to retrieve config parameters")
    }

}

