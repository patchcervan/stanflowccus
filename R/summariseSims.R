#' Summarise Stan estimates/predictions
#'
#' @inheritParams makeStanDataBundle
#' @param group_vars Character. Variable used for grouping predictions. Either `season_id` or `site_id`.
#' @param pred_var Character. The parameter we want to summarise. One of:
#' "psi", "p", "p_cond", "psi_eq"
#'
#' @return A matrix of summarised estimates/predictions by the variables specified in
#' `group_vars`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     summariseSims(occu_data = my_occu_data,
#'                   fit = my_occu_fit,
#'                   group_var = "season_id",
#'                   pred_var = "psi"
#'     )
#' }
summariseSims <- function(occu_data, fit,
                          group_var =  c("season_id", "site_id"),
                          pred_var = c("psi", "p", "p_cond")){

    group_var <- match.arg(group_var)
    pred_var <- match.arg(pred_var)

    # Extract draws
    if(pred_var == "psi_eq"){
        preds <- rstan::extract(fit, pars = c("mean_phi", "mean_gamma")) # I need to review this
    } else if(pred_var == "p"){
        preds <- rstan::extract(fit, pars = c("psi", "p_cond"))
    } else {
        preds <- rstan::extract(fit, pars = pred_var)
    }

    # if(pred_var %in% c("mean_psi", "mean_phi", "mean_gamma")){
    #
    #     sim_dist <- preds[[pred_var]]
    #     out <- occu_data$visit |>
    #         dplyr::select(dplyr::all_of(group_var)) |>
    #         dplyr::distinct() |>
    #         dplyr::mutate(variable = pred_var,
    #                       mean = apply(sim_dist, 2, mean, na.rm = TRUE),
    #                       ub = apply(sim_dist, 2, quantile, 0.975, na.rm = TRUE),
    #                       lb = apply(sim_dist, 2, quantile, 0.025, na.rm = TRUE))
    # }


    # if(pred_var %in% c("log_lik")){
    #
    #     if(group_var != "site"){
    #         warning(paste("The variable", pred_var, "can only be grouped by site"))
    #     }
    #
    #     sim_dist <- preds[[pred_var]]
    #
    #     out <- occu_data$visit |>
    #         dplyr::select(dplyr::all_of(group_var)) |>
    #         dplyr::distinct() |>
    #         dplyr::mutate(variable = pred_var,
    #                       mean = apply(sim_dist, 2, mean, na.rm = TRUE),
    #                       ub = apply(sim_dist, 2, quantile, 0.975, na.rm = TRUE),
    #                       lb = apply(sim_dist, 2, quantile, 0.025, na.rm = TRUE))
    # }

    if(pred_var %in% c("psi")){

        # We need to find groups
        groups <- occu_data$visit |>
            dplyr::distinct(site_id, season_id) |>
            dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) |>
            dplyr::mutate(group = dplyr::cur_group_id()) |>
            dplyr::ungroup()

        # Mean of simulated detections for each season_id
        ngroups <- max(groups$group)

        sim_dist <- matrix(NA, nrow = ngroups, ncol = nrow(preds[[pred_var]]))

        for(g in 1:ngroups){
            idx <- which(groups$group == g)
            sim_dist[g,] <- apply(preds[[pred_var]][, idx, drop = FALSE], 1, mean, na.rm = TRUE)
        }

        out <- groups |>
            dplyr::distinct(.data[[group_var]]) |>
            dplyr::mutate(variable = pred_var,
                          mean = apply(sim_dist, 1, mean, na.rm = TRUE),
                          ub = apply(sim_dist, 1, quantile, 0.975, na.rm = TRUE),
                          lb = apply(sim_dist, 1, quantile, 0.025, na.rm = TRUE))

    }

    # Need to review this
    # if(pred_var == "psi_eq"){
    #     sim_dist <- preds[c("mean_phi", "mean_gamma")]
    #     sim_dist <- sim_dist[[2]] / (1 - sim_dist[[1]] + sim_dist[[2]])
    #     out <- occu_data$visit |>
    #         dplyr::select(dplyr::all_of(group_var)) |>
    #         dplyr::distinct() |>
    #         dplyr::mutate(variable = pred_var,
    #                       mean = apply(sim_dist, 2, mean, na.rm = TRUE),
    #                       ub = apply(sim_dist, 2, quantile, 0.975, na.rm = TRUE),
    #                       lb = apply(sim_dist, 2, quantile, 0.025, na.rm = TRUE))
    # }

    if(pred_var == "p"){

        # We need to identify years in detection data
        site_years <- occu_data$visit |>
            dplyr::group_by(site_id, season_id) |>
            dplyr::mutate(group = dplyr::cur_group_id()) |>
            dplyr::ungroup() |>
            dplyr::pull(group)

        # Compute marginal probability of detection
        p_marg <- preds$p_cond / preds$psi[ ,site_years]

        # Remove missing data
        p_marg[,is.na(occu_data$visit$det)] <- NA

        # We need to sum simulated detections over each season_id
        groups <- occu_data$visit |>
            dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) |>
            dplyr::mutate(group = dplyr::cur_group_id()) |>
            dplyr::ungroup()

        # Mean of simulated detections for each season_id
        ngroups <- dplyr::n_distinct(groups$group)

        sim_dist <- matrix(NA, nrow = ngroups, ncol = nrow(preds[[1]]))

        for(g in 1:ngroups){
            idx <- which(groups$group == g)
            sim_dist[g,] <- apply(p_marg[, idx, drop = FALSE], 1, mean, na.rm = TRUE)
        }

        out <- occu_data$visit |>
            dplyr::select(dplyr::all_of(group_var)) |>
            dplyr::distinct() |>
            dplyr::mutate(variable = pred_var,
                          mean = apply(sim_dist, 1, mean, na.rm = TRUE),
                          ub = apply(sim_dist, 1, quantile, 0.975, na.rm = TRUE),
                          lb = apply(sim_dist, 1, quantile, 0.025, na.rm = TRUE))
    }

    if(pred_var == "p_cond"){

        # We need to sum simulated detections over each season_id
        groups <- occu_data$visit |>
            dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) |>
            dplyr::mutate(group = dplyr::cur_group_id()) |>
            dplyr::ungroup()

        # Calculate proportion of observed detections for each group
        det_group <- occu_data$visit |>
            dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) |>
            dplyr::mutate(group = dplyr::cur_group_id()) |>
            dplyr::filter(!is.na(det)) |>
            dplyr::summarise(dets = sum(det, na.rm = TRUE),
                             n = dplyr::n(),
                             props = dets/n) |>
            dplyr::ungroup()

        # There might be missing data for some years or site, so they might not
        # be represented
        det_group <- groups |>
            dplyr::distinct(dplyr::across(dplyr::all_of(c(group_var, "group")))) |>
            dplyr::left_join(det_group, by = group_var)

        # Mean of simulated detections for each season_id
        ngroups <- nrow(det_group)

        sim_dist <- matrix(NA, nrow = ngroups, ncol = nrow(preds[[pred_var]]))

        # if(input$mod == "dynamic"){
        #     pred_data$visit_p <- t(pred_data$visit_p)
        # }

        # # If we did this we incorporate some variability each time we sample
        # z <- matrix(rbinom(length(preds[[pred_var]]), 1, preds[[pred_var]]),
        #             ncol = ncol(preds[[pred_var]]))

        # This gives more consistent results but produces expected values
        preds <- preds[[pred_var]]

        for(g in 1:ngroups){
            idx <- which(groups$group == g)
            sim_dist[g,] <- apply(preds[, idx, drop = FALSE], 1, sum, na.rm = TRUE)
            sim_dist[g,] <- sim_dist[g, ] / det_group$n[g]
        }

        out <- det_group |>
            dplyr::mutate(variable = pred_var,
                          mean = apply(sim_dist, 1, mean, na.rm = TRUE),
                          ub = apply(sim_dist, 1, quantile, 0.975, na.rm = TRUE),
                          lb = apply(sim_dist, 1, quantile, 0.025, na.rm = TRUE))

    }

    return(out)

}
