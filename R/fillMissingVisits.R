#' Fill in missing visits in detection data
#'
#' @param visit_data A detection/non-detection dataframe.
#' @param site_var Character string. Name of the variable used for identifying sites
#' @param season_var Character string. Name of the variable used for identifying seasons
#' @param seasons An integer vector identifying the seasons the data should cover.
#' @param period Character string. One of "all" (default), "start" and "end". Defining
#' which end of the data should be filled with missing values.
#'
#' @return A dataframe with all possible combinations of `site_var` and `season_var`
#' values. `season_var` values will be taken from `season` if specified or from
#' `range(visit_data)` otherwise. Those cases that were not observed will have NA in detection and any
#' other variables, such as covariates.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' visit_data <- fillMissingVisits(visit_data,
#'                                 site_var = "UTM",
#'                                 season_var = "year",
#'                                 seasons = 1990:2020,
#'                                 period = "all")
#' }

fillMissingVisits <- function(visit_data, site_var, season_var,
                              seasons = NULL,
                              period = c("all", "start", "end")){

    period <- match.arg(period)

    det_data_nest <- visit_data |>
        dplyr::nest_by(.data[[site_var]])

    min_season <- ifelse(is.null(seasons), min(visit_data[season_var]), min(seasons))
    max_season <- ifelse(is.null(seasons), max(visit_data[season_var]), max(seasons))

    if(period == "all"){

        det_periods <- visit_data |>
            dplyr::mutate(min_season = min_season) |>
            dplyr::mutate(max_season = max_season) |>
            dplyr::nest_by(.data[[site_var]])

    } else if(period == "start"){

        det_periods <- visit_data |>
            dplyr::group_by(.data[[site_var]]) |>
            dplyr::summarise(max_season = max(.data[[season_var]])) |>
            dplyr::ungroup() |>
            dplyr::mutate(min_season = min_season) |>
            dplyr::nest_by(.data[[site_var]])

    } else if(period == "end"){

        det_periods <- visit_data |>
            dplyr::group_by(.data[[site_var]]) |>
            dplyr::summarise(min_season = min(.data[[season_var]])) |>
            dplyr::ungroup() |>
            dplyr::mutate(max_season = max_season) |>
            dplyr::nest_by(.data[[site_var]])

    }

    det_data_nest[[season_var]] <- det_periods|>
        dplyr::pull(data) |>
        purrr::map(~unique(.x$min_season):unique(.x$max_season))

    site_seasons <- det_data_nest |>
        dplyr::reframe(.data[[season_var]])

    visit_data <- site_seasons |>
        dplyr::left_join(visit_data, by = c(site_var, season_var))

    return(visit_data)

}
