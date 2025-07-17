#' Create Nearest Neighbours distance matrix
#'
#' @param dist_matrix Matrix. A distance matrix where each entry [i, j] corresponds to the distance between sites i and j.
#' @param nNN Numeric. The number of nearest neighbours to retain. Defaults to 24.
#'
#' @return A distance matrix with the same rows as `dist_matrix` and `nNN` columns
#'
#' @export
#'
#' @examples
#' \dontrun{
#' createNNdistanceMatrix(dist_matrix, 16)
#' }
createNNdistanceMatrix <- function(dist_matrix, nNN = 24){

    # Create an NN distance matrix
    idx_matrix <- matrix(NA, nrow = nrow(dist_matrix), ncol = nNN)
    NN_matrix <- matrix(NA, nrow = nrow(dist_matrix), ncol = nNN)
    sites_matrix <- matrix(NA, nrow = nrow(dist_matrix), ncol = nNN)

    for(i in 1:nrow(dist_matrix)){
        idx_matrix[i, ] <- order(dist_matrix[i,])[1:nNN]
        NN_matrix[i, ] <- dist_matrix[i,idx_matrix[i, ]]
        sites_matrix[i,] <- dimnames(dist_matrix)[[2]][idx_matrix[i,]]
    }

    dimnames(NN_matrix)[[1]] <- dimnames(dist_matrix)[[1]]

    attr(NN_matrix, "sites") <- sites_matrix

    # # Which sites are within the range of other sites
    # sites_keep <- unique(as.vector(sites_matrix)) |>
    #     as.numeric() |>
    #     sort()
    #
    # # Find the indices of these sites
    # sites_keep_idx <- which(dimnames(dist_matrix)[[1]] %in% sites_keep)
    #
    # # Subset NN matrix to those sites whitin the range of other sites
    # NN_matrix <- NN_matrix[sites_keep_idx,]

    return(list(idx_matrix = idx_matrix,
                NN_matrix = NN_matrix))

}
