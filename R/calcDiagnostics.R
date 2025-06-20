#' Calculate diagnostics of a Stan occupancy fit
#'
#' @inheritParams makeStanDataBundle
#' @param fit A fitted Stan model.
#' @param test_seasons A numeric vector containing the seasons (`season_id`) used for testing.
#' @param verbose Logical. If TRUE (default), progress of diagnostics is displayed.
#'
#' @return A list of diagnostics
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     calcDiagnostics(occu_data = my_occu_data,
#'                     fit = my_occu_fit,
#'                     test_seasons = 10:15
#'     )
#' }
calcDiagnostics <- function(occu_data, fit, test_seasons = NULL, verbose = TRUE){

    # Start a list
    diags <- list()


    # Calculate AUC -----------------------------------------------------------

    if(verbose){
        message("Calculating out-of-sample AUC")
    }

    # Extract conditional probability of detection from model fit
    preds <- rstan::extract(fit, pars = "p_cond")[[1]]

    # Calculate ROC and AUC for out of sample data

    # This will only work if the response has two levels
    n_levels_out <- occu_data$visit |>
        dplyr::select(site_id, season_id, visit_id, det) |>
        dplyr::filter(!is.na(det)) |>
        dplyr::filter(season_id %in% test_seasons) |>
        dplyr::pull(det) |>
        unique() |>
        length()

    if(n_levels_out == 2){

        pred_roc <- pROC::roc(response = occu_data$visit |>
                                  dplyr::select(site_id, season_id, visit_id, det) |>
                                  dplyr::filter(!is.na(det)) |>
                                  dplyr::filter(season_id %in% test_seasons) |>
                                  dplyr::pull(det),
                              predictor = occu_data$visit |>
                                  dplyr::select(site_id, season_id, visit_id, det) |>
                                  dplyr::mutate(p_cond_mu = apply(preds, 2, mean)) |>
                                  dplyr::filter(!is.na(det)) |>
                                  dplyr::filter(season_id %in% test_seasons) |>
                                  dplyr::pull(p_cond_mu))

        diags$roc_out <- pred_roc

    } else {

        diags$roc_out <- NA

    }

    if(verbose){
        message("Calculating in-sample AUC")
    }

    n_levels_in <- occu_data$visit |>
        dplyr::select(site_id, season_id, visit_id, det) |>
        dplyr::filter(!is.na(det)) |>
        dplyr::filter(!season_id %in% test_seasons) |>
        dplyr::pull(det) |>
        unique() |>
        length()

    if(n_levels_in == 2){

        pred_roc <- pROC::roc(response = occu_data$visit |>
                                  dplyr::select(site_id, season_id, visit_id, det) |>
                                  dplyr::filter(!is.na(det)) |>
                                  dplyr::filter(!season_id %in% test_seasons) |>
                                  dplyr::pull(det),
                              predictor = occu_data$visit |>
                                  dplyr::select(site_id, season_id, visit_id, det) |>
                                  dplyr::mutate(p_cond_mu = apply(preds, 2, mean)) |>
                                  dplyr::filter(!is.na(det)) |>
                                  dplyr::filter(!season_id %in% test_seasons) |>
                                  dplyr::pull(p_cond_mu))

        diags$roc_in <- pred_roc

    } else {

        diags$roc_in <- NA

    }


    # Calculate data likelihood -----------------------------------------------

    if(verbose){
        message("Calculating out-of-sample likelihood")
    }

    # Extract names of all parameters monitored
    all_param <- dimnames(fit)$parameters
    all_param <- unique(gsub("\\[[^]]*\\]", "", all_param))

    if("log_lik_pred" %in% all_param){

        # Extract out-of-sample data likelihood from model fit
        preds <- rstan::extract(fit, pars = "log_lik_pred")[[1]]

        # Calculate data likelihood for out-sample data
        niter <- nrow(preds)

        # Integrate likelihood for each site
        prob_data <- occu_data$visit |>
            dplyr::distinct(site_id) |>
            dplyr::mutate(prob = apply(preds, 2, logSumExp) - log(niter))

        # Take the mean
        diags$lik_out <- prob_data |>
            dplyr::pull(prob) |>
            mean()

    } else {

        diags$lik_out <- NA

    }

    # Calculate model likelihood for in-sample data

    if(verbose){
        message("Calculating in-sample likelihood")
    }

    # Extract in-sample data likelihood from model fit
    preds <- rstan::extract(fit, pars = "log_lik")[[1]]

    # Calculate data likelihood for out-sample data
    niter <- nrow(preds)

    # Integrate likelihood for each site
    prob_data <- occu_data$visit |>
        dplyr::distinct(site_id) |>
        dplyr::mutate(prob = apply(preds, 2, logSumExp) - log(niter))

    diags$lik_in <- prob_data |>
        dplyr::pull(prob) |>
        mean()


    # Proportion of converged parameters --------------------------------------

    rhat <- posterior::rhat
    ess_bulk <- posterior::ess_bulk
    ess_tail <- posterior::ess_tail

    # Some parameters need not converge or we don't want them to compute for conv rate
    no_conv_pars <- c("log_lik", "log_lik_pred", "p_cond", "psi", "mean_psi", "mean_phi", "mean_gamma", "lp__")

    summ <- as.matrix(fit, pars = all_param[!all_param %in% no_conv_pars]) |>
        posterior::as_draws() |>
        posterior::summarise_draws(mean, sd,
                                   rhat = posterior::rhat,
                                   ess_bulk = posterior::ess_bulk,
                                   ess_tail = posterior::ess_tail)

    diags$converge_rate <- sum(abs(summ$rhat - 1) < 0.1) / nrow(summ)


    # Calculate number of divergent transitions -------------------------------

    samp_params <- lapply(rstan::get_sampler_params(fit, inc_warmup = FALSE), as.data.frame)
    div_trans <- sapply(samp_params, function(x) sum(x$divergent_))
    diags$ndivergent <- sum(div_trans)

    return(diags)

}
