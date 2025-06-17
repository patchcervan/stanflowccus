#' Create directory structure
#'
#' This function creates the directory structure for an occupancy modelling workflow
#'
#' @param project_dir Character. The path where the structure should be created. Defaults to the current working directory.
#' @param verbose Logical. If TRUE (default) directories created and omitted are printed.
#'
#' @details
#' This function is meant to be run at an empty directory to start working on a workflow.
#' It will check if any of the directories exist and won't create these to avoid losing
#' existing data.
#'
#'
#' @examples
#' \dontrun{
#' # Example configuration
#' createDirs()
#' createDirs("my/preferred/path")
#' }
#'
#' @export
createDirs <- function(project_dir = getwd(), verbose = TRUE){

    dirs <- file.path(project_dir, c("data", "output", "models", "scripts"))

    dirs_create <- dirs[!dir.exists(dirs)]

    if(verbose & any(dir.exists(dirs))){
        warning(paste0("Couldn't create directories: ",
                      paste(dirs[dir.exists(dirs)], collapse = ", "),
                      ", because they already exist"))
    }

    for(i in dirs_create){
        dir.create(i)
    }

    if(verbose){
        message(paste0("Directories created: ", paste(dirs_create, collapse = ", ")))
    }

}
