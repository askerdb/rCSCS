#' download_GNPS
#'
#' @param ID The GNPS id to retrieve
#' @param dir Location of the output
#'
#' @return Downloads and extracts a file of the GNPS content, returns a list with paths to the buckettable and network files
#' @export
#'
#' @examples  download_GNPS("0310e20491314ddbbf12d56b592548b4", ".")

download_GNPS <- function(id, dir, overwrite = F, mgf=F){
  path = file.path(paste(dir,"/",id,sep=""))
  if (!dir.exists(path)) {dir.create(dir); dir.create(path); }
  if (read_GNPS_dir(path)$buckettable != -1 & overwrite == F ) stop("Buckettable found in dir and overwrite = F")
  download_buckettable(id = id,file.path(path, "buckettable.tsv"))
  download_edges(id = id,file.path(path, "edges_file.txt"))
  if (mgf == T) download_mgf(id = id,file.path(path, "Ions.mgf"))
}

#' Prepare the CSS matrix from an edges table
#'
#' @param edge_path Path to edges file from GNPS
#'
#' @return CSS matrix
#' @export
#'
#' @examples
#' @import igraph

prepare_css <- function(edge_path){
  csslong <- read.table(edge_path, header=T, sep = "\t")
  g <- igraph::graph.data.frame(csslong, directed=FALSE)
  css <- igraph::get.adjacency(g, attr="Cosine", sparse=FALSE)
  diag(css) <- 1
  return(css)
}

#' prepare_GNPS
#'
#' @param ID GNPS ID to prepare
#' @param dir Path to put GNPS files
#' @param select Only output features in the feature table found in the CSS matrix
#' @param overwrite Overwrite the directory?
#' @param mgf Download mgf file?
#' @return A list of a feature matrix and a chemical similarity matrix (CSS)
#' @export
#' @import RCurl
#'
#' @examples prepare_GNPS("0310e20491314ddbbf12d56b592548b4")

prepare_GNPS <- function(id=NULL, dir = ".", select = T, overwrite=F, mgf=F){
  path = file.path(paste(dir,"/",id,sep=""))
  if (!is.null(id)) download_GNPS(id, dir, overwrite, mgf)
  gnps <- read_GNPS_dir(path)
  features <- read.table(gnps$buckettable, sep = "\t", header=T, row.names=1, comment.char="")
  css <- prepare_css(gnps$edges)
  return(list(features = features, css = css))
}

#' Read a previously downloaded GNPS dir
#'
#' @param dir
#'
#' @return return a list of GNPS files, or -1 if no buckettable file was found
#' @export
#'
#' @examples
#'
read_GNPS_dir <- function(ID, dir = "."){
  files <- Sys.glob(file.path(paste(dir,"/",ID,sep=""), "*"))
  netattr <- list.files(grep("clusterinfosummarygroup_attributes_withIDs_withcomponentID", files, value = T), full.names = T)
  buckettable <- grep("buckettable.tsv", files, value = T)
  edges <- grep("edges_file.txt", files, value = T)
  if (length(buckettable) == 0) return(list(buckettable = -1))
  #ids <- list.files(paste(dir,"/",ID,"/","clusterinfosummarygroup_attributes_withIDs",sep=""), full.names = T)
  return(list(buckettable = buckettable, edges = edges, attr = netattr))#, ids = ids))
}

#' Download buckettable from GNPS
#'
#' @param id GNPS task ID
#'
#' @return The downloaded file is written to disk
#' @export
#'
#' @examples
#'
download_buckettable <- function(id, path = "buckettable.tsv"){
  file = paste0("http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=",id,"&block=main&file=cluster_buckets/")
  download_file(file, path)
}

#' Download edges from GNPS
#'
#' @param id GNPS task ID
#'
#' @return The downloaded file is written to disk
#' @export
#'
#' @examples
#'

download_edges <- function(id, path = "edges_file.txt"){
  file = paste0("http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=",id,"&block=main&file=networkedges_selfloop/")
  download_file(file, path)
}
#' Download MGF file from GNPS
#'
#' @param id GNPS task ID
#'
#' @return The downloaded file is written to disk
#' @export
#'
#' @examples
#'

download_mgf <- function(id, path = "ions.mgf"){
  file = paste0("http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=",id,"&block=main&file=spectra/specs_ms.mgf")
  download_file(file, path)
}

#' Download and save a file
#'
#' @param file URL to the file
#' @param path Path to save the file
#'
#' @return The downloaded file is written to disk
#' @export
#'
#' @examples
#'
download_file <- function(file, path = "NULL"){
  if (is.null(path)){
    error("Must enter a path")
  }
  f = CFILE(path, mode="wb")
  curlPerform(url = file, writedata = f@ref)
  close(f)
}

