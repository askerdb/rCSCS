#' Calculate CSCS distance
#'
#' @param features Matrix of feature intensities in each sample (pxn)
#' @param css Square matrix of cosine similarities of features (pxp)
#' @param dissimilarity Output dissimilarity matrix
#' @param cosine_threshold Only include features below this threshold
#' @param weighted Weight features by intensity (TRUE) or absence/presence (FALSE)
#'
#' @details The value of cosine_threshold is 0.6 in Sedio et al. but for certain applications other values might be better.
#' @references Sedio et al, 2016
#' @return A pxp matrix of CSCS distances
#' @export
#'
#' @examples
#'
#' #Get arbitrary GNPS project
#' gnps <- prepare_GNPS("0310e20491314ddbbf12d56b592548b4")
#' dist <- cscs(gnps$features, gnps$css)
#'
#' #GlobalEuphorbia data
#' cscs(GEfeatures, GEcss)
#' @importFrom foreach %dopar%
cscs <- function(features, css, dissimilarity = T, cosine_threshold = 0.6, weighted = T,  normalize = T){

  stopifnot(nrow(features) == nrow(css))
  if (normalize==T){
    features <- apply(features, 2, function(x) x/sum(as.numeric(x)))
  }

  diag(css) <- 1
  css[css < cosine_threshold] <- 0
  css <- as.matrix(css)
  sample_names <- colnames(features)
  distm <- matrix(0, nrow = length(sample_names), ncol = length(sample_names),
                  dimnames = list(sample_names, sample_names))
  spn <- combn(sample_names, 2, simplify=FALSE)
  distlist <- foreach::foreach( pair = spn) %dopar% {
    i <- pair[1]
    j <- pair[2]
    feature_union <- rownames(features)[which(rowSums(features[,c(i,j)]) > 0)]
    css_tmp <- css[feature_union,feature_union]
    a <-  features[feature_union,i]
    b <- features[feature_union,j]
    if (weighted == FALSE){
      a[a > 0] <- 1
      b[b > 0] <- 1
    }
    abt <- a %*% t(b)
    aat <- a %*% t(a)
    bbt <- b %*% t(b)
    cssAB <- css_tmp * abt
    d <- sum(cssAB)/max(sum(css_tmp * aat), sum(css_tmp * bbt))
  }
  distm[lower.tri(distm)] <- as.numeric(distlist)
  distmt <- t(distm)
  distmt[lower.tri(distmt)] <- as.numeric(distlist)
  colnames(distmt) <- sample_names
  row.names(distmt) <- sample_names

  # Make dissimilarity
  if (dissimilarity == T){
    distmt <- 1 - distmt
    diag(distmt) <- 0
  return(as.dist(distmt))
  }
  diag(distmt) <- 1
  return(distmt)
}

