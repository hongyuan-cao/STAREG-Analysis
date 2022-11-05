## A function to calculate the Moran's I and Geary's C

corrValue <- function(counts, location){
  require(spdep)
  coords = as.matrix(location)
  k1 = knn2nb(knearneigh(coords))
  all.linked <- max(unlist(nbdists(k1, coords)))
  col.nb <- dnearneigh(coords, 0, all.linked, 
                       row.names=row.names(coords), 
                       longlat = FALSE)
  col.S <- nb2listw(col.nb, style = "S")
  MI.S <- apply(counts, 1, moran, listw = col.S, n = length(col.nb), S0 = Szero(col.S))
  f <- function(x){x = x[[1]]}
  MI.S <- lapply(MI.S, f)
  MI <- NULL
  for(i in 1:length(MI.S)){
    MI <- c(MI, MI.S[[i]])
  }
  names(MI) <- row.names(counts)
  
  GC.S <- apply(counts, 1, geary, listw = col.S, 
                      n = length(col.nb), n1 = length(col.nb) - 1, S0 = Szero(col.S))
  GC.S <- lapply(GC.S, f)
  GC <- NULL
  for(i in 1:length(GC.S)){
    GC <- c(GC, GC.S[[i]])
  }
  names(GC) <- row.names(counts)
  
  return(list(MI = MI, GC = GC))
}

gearyValur <- function(counts, location){
  
}