#' Create a Stan data bundle
#'
#' @param occu_data List. A list with two dataframes: one for site data named "site" and one
#' for visit data named "visit". The site data frame must contain at least one column named
#' "site_id" indexing the sites and one named "season_id" indexing the seasons. The
#' visit data frame must at least have one column named "det" with a 1 for detections and
#' 0 for non-detections and a "visit_id" column indexing the visits (overall, not within seasons).
#' @param occ_formula Formula. A formula object specifying the initial occurrence model
#' @param phi_formula Formula. A formula object specifying the persistence dynamics model
#' @param gamma_formula Formula. A formula object specifying the colonization dynamics model
#' @param det_formula Formula. A formula object specifying the detection model
#' @param dynamic Logical. Indicates whether the specified model is a dynamic one (TRUE)
#' or not (FALSE, default)
#'
#' @return A data bundle to pass on to a Stan model
#'
#' @export
#'
#' @examples
#'
#' occu_data <- list(site = data.frame(site_id = rep(1:3, 3),
#'                                     season_id = rep(1:3, each = 3)),
#'                   visit = data.frame(site_id = rep(1:3, 7),
#'                                      season_id = rep(1:3, each = 7),
#'                                      det = rbinom(21, 1, 0.3),
#'                                      effort = rpois(21, 5)))
#'
#' occu_form <- ~ 1
#' det_form <- ~ effort
#'
#' makeStanDataBundle(occu_data,
#'                    occ_formula = occu_form,
#'                    det_formula = det_form,
#'                    dynamic = FALSE)
makeStanDataBundle <- function(occu_data, occ_formula,
                               phi_formula = NULL, gamma_formula = NULL,
                               det_formula, dynamic = FALSE){

    if(dynamic){

        # Create design matrices
        X_psi0 <- formToDesignMatrix(occ_formula,
                                     occu_data$site |>
                                         dplyr::group_by(site_id) |>
                                         dplyr::filter(season_id == min(season_id)) |>
                                         dplyr::ungroup())

        X_psi0 <- X_psi0[[1]] |> as.matrix()

        # if(nrow(X_psi0[[2]]) == 0){
        #   X_psi0 <- X_psi0[[1]]
        # }

        X_phi <- formToDesignMatrix(phi_formula, occu_data$site)

        X_phi <- X_phi[[1]] |> as.matrix()

        # if(nrow(X_phi[[2]]) == 0){
        #   X_phi <- X_phi[[1]]
        # }

        X_gamma <- formToDesignMatrix(gamma_formula, occu_data$site)

        X_gamma <- X_gamma[[1]] |> as.matrix()

        # if(nrow(X_gamma[[2]]) == 0){
        #   X_gamma <- X_gamma[[1]]
        # }

        X_psi <- NULL

    } else {

        # We modify the design matrix to match the dynamic model
        # data generation, where there is an initial occupancy
        # and the following season depend on the covariates of the
        # previous season
        X_psi <- formToDesignMatrix(occ_formula, occu_data$site)

        X_psi <- X_psi[[1]] |> as.matrix()

        X_psi0 <- NULL
        X_phi <- NULL
        X_gamma <- NULL

    }

    X_p <- formToDesignMatrix(det_formula, occu_data$visit |>
                                  dplyr::filter(!is.na(det)))

    X_p <- X_p[[1]] |> as.matrix()
    # if(nrow(X_p[[2]]) == 0){
    #   X_p <- X_p[[1]]
    # }

    # Bundle model data
    K <- c(ncol(X_psi),ncol(X_psi0), ncol(X_phi), ncol(X_gamma), ncol(X_p))
    data.bundle <- list(N_total = nrow(occu_data$visit),
                        N_obs = occu_data$visit |>
                            dplyr::filter(!is.na(det)) |>
                            nrow(),
                        N = dplyr::n_distinct(occu_data$visit$site_id),
                        nseasons = dplyr::n_distinct(occu_data$visit$season_id),
                        N_site = occu_data$visit |>
                            dplyr::count(site_id) |>
                            dplyr::pull(n), # Note that this includes missing data
                        S = occu_data$visit |>
                            dplyr::distinct(site_id, season_id) |>
                            nrow(), # Note that this includes missing data
                        G = occu_data$visit |>
                            dplyr::group_by(site_id) |>
                            dplyr::summarise(n = dplyr::n_distinct(season_id)) |>
                            dplyr::pull(n), # Note that this includes missing data
                        M = occu_data$visit |>
                            dplyr::group_by(site_id, season_id) |>
                            dplyr::summarise(n = dplyr::n()) |>
                            dplyr::pull(n), # Note that this includes missing data
                        K = K,
                        X_psi = X_psi,
                        X_psi0 = X_psi0,
                        X_phi = X_phi,
                        X_gamma = X_gamma,
                        X_p = X_p,
                        ii_obs = which(!is.na(occu_data$visit$det)),
                        ii_mis = which(is.na(occu_data$visit$det)),
                        y = occu_data$visit |>
                            dplyr::filter(!is.na(det)) |>
                            dplyr::pull(det),
                        site = occu_data$visit |>
                            dplyr::filter(!is.na(det)) |>
                            dplyr::pull(site_id),
                        season = occu_data$visit |>
                            dplyr::distinct(site_id, season_id) |>
                            dplyr::pull(season_id),
                        season_visit = occu_data$visit |>
                            dplyr::filter(!is.na(det)) |>
                            dplyr::pull(season_id)
    )

    return(data.bundle)
}


#' Construct a design matrix from a formula object
#'
#' @param form A formula object
#' @param df A data frame with the variables referred to by `form`
#'
#' @return A list with design matrices containing the variables
#'
#' @details
#' The notation "(1|x)" creates a design matrix for random intercepts.
#' At the moment there are no random coefficients other than intercepts.
#'
#' @export
#'
#' @examples
#' df <- data.frame(x = runif(20, 0, 10),
#'                  y = runif(20, -10, 0),
#'                  z = rnorm(20, 0, 1))
#'
#' form1 <- ~ x + y
#'
#' form2 <- ~ x * y # (this does not produce an interaction!)
#'
#' form3 <- ~ X + y + (1|z)
#'
#' formToDesignMatrix(form1, df)
#'
#' formToDesignMatrix(form2, df)
#'
#' formToDesignMatrix(form3, df) # This separates fixed from random effects
formToDesignMatrix <- function(form, df){

    tt <- terms(form)

    vars <- attr(tt, "variables")
    vars <- as.character(vars[-1])

    # Detect functions
    i_funs <- grep("\\(", vars)
    funs <- gsub("\\(", "\\(df$", vars[i_funs])

    # Detect random effects
    i_rms <- grep("\\|", vars)
    rms <- gsub(".*\\| ", "", vars[i_rms])

    # Fill in output variables
    out <- vector("list", length = length(vars))

    for(i in seq_along(i_funs)){
        out[[i_funs[i]]] <- eval(parse(text = funs[i]))
        names(out)[i_funs[i]] <- funs[i]
    }

    for(i in seq_along(i_rms)){
        out[[i_rms[i]]] <- df[,rms[i]]
        names(out)[i_rms[i]] <- rms[i]
    }

    if(length(i_funs) + length(i_rms) == 0){
        for(i in seq_along(out)){
            out[[i]] <-  df[, vars[i]]
            names(out)[i] <- vars[i]
        }
    } else {
        for(i in seq_along(out)[-c(i_funs, i_rms)]){
            out[[i]] <-  df[, vars[i]]
            names(out)[i] <- vars[i]
        }
    }

    # Separate into fixed and random effects
    if(length(i_rms) != 0){
        out_fix <- out[-i_rms]
        out_rm <- out[i_rms]

        # Fix names
        names(out_fix) <- gsub("df\\$", "", names(out_fix))
        names(out_fix) <- gsub("^I\\(|\\)", "", names(out_fix))
        names(out_fix) <- gsub("\\(", "_", names(out_fix))

        names(out_rm) <- gsub("df\\$", "", names(out_rm))
        names(out_rm) <- gsub("^I\\(|\\)", "", names(out_rm))
        names(out_rm) <- gsub("\\(", "_", names(out_rm))

    } else {
        out_fix <- out
        out_rm <- NULL

        # Fix names
        names(out_fix) <- gsub("df\\$", "", names(out_fix))
        names(out_fix) <- gsub("^I\\(|\\)", "", names(out_fix))
        names(out_fix) <- gsub("\\(", "_", names(out_fix))

    }

    if(attr(tt, "intercept")){
        if(length(out_fix) == 0){
            out_fix <- data.frame(intcp = rep(1, nrow(df)))
        } else {
            out_fix <- cbind(data.frame(intcp = rep(1, nrow(df))), out_fix)
        }
    }

    list(X = out_fix |>
             as.data.frame() |>
             as.matrix(),
         X_re = out_rm |>
             as.data.frame() |>
             as.matrix())
}
