#' Filter out rows with NA or low value
#'
#' @param m A matrix-like object.
#' @param design A vector, indicating the group of columns and specifying to
#' filter out low-quality rows within each condition.
#' @param minS A numeric, specifying the minimum proportion/number of values should
#' be quantified for each row.
#' @param out A vector, specifying low-quality values (NA and 0).
#'
#' @return A matrix with the same columns as input matrix.
#' @author Wubing Zhang
#' @export
#'
filterN <- function(m, design=NULL, minS = 3, out = c(NA, 0),
                    imputeNA = TRUE, ...){
  m = as.matrix(m)

  if(is.null(design)){
    if(minS<1) minS = minS*ncol(m)
    idx = is.na(m)
    if(length(out[!is.na(out)])>0)
      idx = idx | (m<=max(out[!is.na(out)]))
    sel = rowSums(!idx)>=minS
    return(m[sel,])
  }else{
    if(length(design)!=ncol(m)) stop("parameter design error ...")
    nonsel = c()
    for(con in unique(design)){
      tmp = which(design==con)
      idx = is.na(m[, tmp])
      if(length(out[!is.na(out)])>0)
        idx = idx | (m[,tmp]<=max(out[!is.na(out)]))
      if(minS<1)
        nonsel = c(nonsel, which(rowSums(!idx)<minS*length(tmp)))
      else
        nonsel = c(nonsel, which(rowSums(!idx)<minS))
    }
    if(imputeNA) m = imputeNA(m[-nonsel, ], design, ...)
    return(m)
  }
}

imputeNA <- function(m, design = NULL, k = 3, maxIteration = 10){
  getDist <- function(m){
    tmp = filterN(m, minS = 0.8, out = NA)
    tmp = t(apply(tmp, 1, function(x){
      x[is.na(x)] = mean(x, na.rm = TRUE)
      return(x)
    }))
    distance = cor(tmp)
    diag(distance) = 0
    return(distance)
  }
  imputed_m = m
  iteration = 0
  while(sum(is.na(imputed_m))>0 & iteration<=maxIteration){
    if(is.null(design)){# impute NA as mean of k-nearest neighbor samples.
      if(iteration==maxIteration) k = ncol(m)-2
      distance = getDist(m)
      imputed_m <- sapply(1:ncol(m), function(j, m, distance){
        values = m[,j]
        kmean = rowMeans(m[, order(-distance[j,])[1:k], drop = FALSE], na.rm = TRUE)
        values[is.na(values)] = kmean[is.na(values)]
        return(values)
      }, m, distance)
      # if(sum(is.na(imputed_m))>0){ # Impute NA as mean of k-nearest neighbor genes.
      #   distance = getDist(t(imputed_m))
      #   imputed_m <- t(sapply(1:nrow(m), function(i, m, distance){
      #     values = m[i,]
      #     kmean = colMeans(m[order(-distance[i,])[1:k], ], na.rm = TRUE)
      #     values[is.na(values)] = kmean[is.na(values)]
      #     return(values)
      #   }, imputed_m, distance))
      # }
    }else{ # impute NA within each group
      if(length(design)!=ncol(m)) stop("parameter design error ...")
      for(con in unique(design)){
        idx = design==con
        if(iteration==maxIteration) k = sum(idx)-2
        distance = getDist(m[, idx])
        imputed_m[, idx] <- sapply(1:sum(idx), function(j, m, distance){
          values = m[,j]
          kmean = rowMeans(m[, order(-distance[j,])[1:k], drop = FALSE], na.rm = TRUE)
          values[is.na(values)] = kmean[is.na(values)]
          return(values)
        }, m[, idx], distance)
      }
    }
    iteration = iteration + 1
  }
  return(imputed_m)
}

