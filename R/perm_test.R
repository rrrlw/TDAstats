
# Script for Statistical Inference using Persistent Hom
# based on functions in TDAstats source code

# IMPORT LIBRARIES ----
library(clue)
library(rdist)
library(TDA)        # for circleUnif
library(data.table) # for as.data.table

# FUNCTIONS ----

#FUNCTION CHECKED (RRRLW)
error_check <- function(obj1, obj2, iterations,
                        p, q, dim, format, standardize,
                        input_type){
  # function which performs all necessary error checks for the permutation_test_two_groups function
  # each parameter is exactly the same as in the permutation_test_two_groups function
  
  if (missing(obj1) | missing(obj2)) stop("both parameters must be supplied.")
  
  if (input_type == "list")
  {
    
    if ((!is.list(obj1) & !is.vector(obj1)) |
        (!is.list(obj2) & !is.vector(obj2))) stop("l1 and l2 must either be lists or vectors.")
    
    df_mat_check1 <- vapply(X = 1:length(obj1),
                            FUN.VALUE = logical(1),
                            FUN = function(i) {
                              curr_obj <- obj1[[i]]
                              
                              if (!is.data.frame(curr_obj) &
                                  !is.matrix(curr_obj)) return(TRUE)
                              
                              if (nrow(curr_obj) <= 1 |
                                  ncol(curr_obj) == 0) return(TRUE)
                              
                              return(FALSE)
                            })
    
    df_mat_check2 <- vapply(X = 1:length(obj2),
                            FUN.VALUE = logical(1),
                            FUN = function(i) {
                              curr_obj <- obj2[[i]]
                              
                              if (!is.data.frame(curr_obj) &
                                  !is.matrix(curr_obj)) return(TRUE)
                              
                              if (nrow(curr_obj) <= 1 |
                                  ncol(curr_obj) == 0) return(TRUE)
                              
                              return(FALSE)
                            })
    
    if (any(df_mat_check1) |
        any(df_mat_check2)) stop("all list elements must either be multi-row dataframes or multi-row matrices of at least 1 column.")
    
  } else if (input_type == "sample") {
    
    if(!is.matrix(obj1) & !is.data.frame(obj1)) stop("df1 must either be a matrix or a data frame.")
    
    if(nrow(obj1) <= 1 | ncol(obj1) == 0) stop("df1 must have at least 2 rows and 1 column.")
    
    if(!is.matrix(obj2) & !is.data.frame(obj2)) stop("df2 must be either a matrix or a data frame.")
    
    if(nrow(obj2) <= 1 | ncol(obj2) == 0) stop("df2 must have at least 2 rows and 1 column.")
    
  } else stop("Invalid input type")
  
  if(missing(iterations)) stop("iterations must be supplied.")
  
  if(!is.integer(iterations) &
     !is.numeric(iterations)) stop("iterations must be an integer or numeric")
  
  if(iterations < 1) stop("iterations must be at least 1.")
  
  if(!is.numeric(p) |
     p < 1 |
     is.infinite(p)) stop("p must be a finite number (at minimum equal to 1)")
  
  if(!is.numeric(q) |
     q < 1 |
     is.infinite(q)) stop("q must be a finite number (at minimum equal to 1)")
  
  if(missing(dim)) stop("dim must be supplied.")
  
  if(!is.integer(dim) |
     dim < 0 |
     is.infinite(dim)) stop("dim must be a finite and positive integer.")
  
  if(!is.character(format)) stop("format must be a character string.")
  
  if(format %in% c("distmat","cloud") == FALSE) stop("format must either be distmat or cloud.")
  
  if(!is.logical(standardize)) stop("standardize must be either TRUE or FALSE.")
}

enclosing_radius <- function(point_cloud){
  # function which finds radius beyond which no homology changes
  # point_cloud is a point cloud data frame
  
  pc_dist <- dist(point_cloud)
  pc_nrow <- nrow(point_cloud)
  
  #I DONT UNDERSTAND HOW/WHY THIS WORKS (RRRLW)
  dist_calc <- vapply(X = 1:(pc_nrow - 1),
                      FUN.VALUE = numeric(1),
                      FUN = function(i) {
                        max(
                          pc_dist[(1 + (i - 1) * n - (i - 1) * i / 2):(i * n - i * (i + 1) / 2)]
                        )
                      })
  
  return(min(dist_calc))
  
}

loss <- function(barcodes1,barcodes2,dim,p,q){
  # function to compute the F_{p,q} loss between two sets of barcodes
  # barcodes1 and barcodes2 are two lists of barcodes as outputted from the calculate homology function
  # dim is the maximum dimension of barcodes to consider
  # p and q are finite numbers >=1
  
  # get all possible pairs of distinct members of each list of barcodes
  comb1 = combn(1:length(barcodes1),m = 2,simplify = FALSE)
  comb2 = combn(1:length(barcodes2),m = 2,simplify = FALSE)
  
  # compute the pairwise barcode distances for both lists
  d1_tot = do.call(rbind,lapply(X = comb1,FUN = function(X){ return(d(B1 = as.data.frame(barcodes1[[X[[1]]]]),B2 = as.data.frame(barcodes1[[X[[2]]]]),dim = dim,p = p)^q) }))
  d2_tot = do.call(rbind,lapply(X = comb2,FUN = function(X){ return(d(B1 = as.data.frame(barcodes2[[X[[1]]]]),B2 = as.data.frame(barcodes2[[X[[2]]]]),dim = dim,p = p)^q) }))
  
  # for each dimension we compute the total loss function as described in Robinson & Turner
  res = lapply(X = 0:dim,FUN = function(X){
    
    l = (1/(2*length(barcodes1)*(length(barcodes1) - 1))) * 2 * sum(d1_tot[,(X+1)])
    l = l + (1/(2*length(barcodes2)*(length(barcodes2) - 1))) * 2 * sum(d2_tot[,(X+1)])
    return(l)
    
  })
  
  # return list of distances, one for each dimension
  return(res)
  
}

d <- function(B1,B2,dim,p){
  # function to compute the wasserstein metric between two barcodes
  # B1 and B2 are barcodes as outputted from the function calculate_homology
  # dim is the maximum dimension to consider
  # p is the finite power of the wasserstein distance, p >= 1
  
  # calculate the wasserstein metric in each dimension
  metrics = lapply(X = 0:dim,FUN = function(X){
    
    # subset both barcodes by dimension X
    B1_subset = B1[which(B1$dimension == X),]
    B2_subset = B2[which(B2$dimension == X),]
    B1_subset = B1_subset[,2:3]
    B2_subset = B2_subset[,2:3]
    
    # create empty diagonals for the persistence landscapes
    diag1 = B1_subset[0,]
    diag2 = B2_subset[0,]
    
    # if both subsets are empty then set their distance to 0
    if(nrow(B1_subset) == 0 & nrow(B2_subset) == 0)
    {
      return(0)
    }
    
    if(nrow(B1_subset) > 0)
    {
      for(i in 1:nrow(B1_subset))
      {
        if(B1_subset[i,1] != B1_subset[i,2])
        {
          # for each non-trivial element in B1_subset we add its projection onto the diagonal in diag1
          proj_diag = mean(as.numeric(B1_subset[i,]))
          diag1 = rbind(diag1,data.frame(birth = proj_diag,death = proj_diag))
        }
      }
    }
    
    if(nrow(B2_subset) > 0)
    {
      for(i in 1:nrow(B2_subset))
      {
        if(B2_subset[i,1] != B2_subset[i,2])
        {
          # for each non-trivial element in B2_subset we add its projection onto the diagonal in diag2
          proj_diag = mean(as.numeric(B2_subset[i,]))
          diag2 = rbind(diag2,data.frame(birth = proj_diag,death = proj_diag))
        }
      }
    }
   
    # since an element b of B1_subset is either matched to an element of B2 or to the projection of b onto the diagonal
    # we form the two sets to be matched by row binding B1_subset with diag2 and B2_subset with diag1
    B1_subset = rbind(B1_subset,diag2)
    B2_subset = rbind(B2_subset,diag1)
    
    # use cdist function from the rdist package to find the wasserstein distance matrix between rows of the updated B1_subset and B2_subset
    dist_mat = as.matrix(cdist(B1_subset,B2_subset,metric = "minkowski",p = p))
    
    # use the Hungarian algorithm from the clue package to find the minimal weight matching
    best_match = solve_LSAP(x = dist_mat,maximum = FALSE)
    
    # return the distance for dimension X
    return(sum(dist_mat[cbind(seq_along(best_match), best_match)]^(p))^(1/p))
    
  })
  
  # merge and organize 
  ret_dt = as.data.table(t(as.matrix(unlist(metrics))))
  setnames(ret_dt,old = colnames(ret_dt),new = paste("dimension_",c(0:dim),sep = ""))
  
  return(ret_dt)
  
}

permutation_test_two_groups <- function(l1,l2,iterations,p,q,dim,format,standardize){
  
  # function to test whether or not two sets of labelled data come from the same geometric process
  # l1 is a list of datasets in the first group and likewise for l2 and the second group
  # iterations is the number of permutations we will calculate for class labels
  # q is the finite exponent used in distance calculations between barcodes, q >= 1
  # dim is the maximum dimension in which to calculate homology for barcodes
  # p is the finite wasserstein distance parameter, p >= 1
  
  # perform all error checks, stop() called from function if error found
  error_check(obj1 = l1, obj2 = l2, iterations = iterations,
              p = p, q = q, dim = dim, format = format,
              standardize = standardize, input_type = "list")
  
  # calculate barcodes for each group
  barcodes1 = lapply(X = l1,FUN = function(X){return(calculate_homology(as.matrix(X), format = format, standardize = standardize, dim = dim,threshold = enclosing_radius(X)))})
  barcodes2 = lapply(X = l2,FUN = function(X){return(calculate_homology(as.matrix(X), format = format, standardize = standardize, dim = dim,threshold = enclosing_radius(X)))})
  
  # get group sizes
  n1 = length(l1)
  n2 = length(l2)
  
  # compute loss function on observed data
  test_loss = loss(barcodes1 = barcodes1,barcodes2 = barcodes2,dim = dim,p = p,q = q)
  
  # get permutation values
  perm_values = lapply(X = 1:iterations,FUN = function(X){
    
    # sample two groups of size n1 and n2 from l1 union l2
    ind = sample(1:(n1+n2),size = n1,replace = FALSE)
    ind1 = ind[which(ind <= n1)]
    ind2 = setdiff(ind,ind1) - n1
    barcodes1_sample = c(barcodes1[ind1],barcodes2[ind2])
    barcodes2_sample = c(barcodes1[setdiff(1:n1,ind1)],barcodes2[setdiff(1:n2,ind2)])
    
    # return loss function
    ret_loss = t(as.matrix(unlist(loss(barcodes1 = barcodes1_sample,barcodes2 = barcodes2_sample,dim = dim,p = p,q = q))))
    return(ret_loss)
    
  })
  perm_values = do.call(rbind,perm_values)
  
  # organize results by dimension
  answer <- lapply(X = 0:dim,
                   FUN = function(curr.dim) {
                     curr.ans <- list()
                     curr.ans$dimension <- curr.dim
                     curr.ans$permvals <- perm_values[,curr.dim + 1] # this is the distance between the barcodes in the current dimension across each iteration
                     curr.ans$wasserstein <- test_loss[[curr.dim + 1]] # test statistic distance between barcodes in this dimension
                     curr.ans$pvalue <- (sum(curr.ans$permvals > curr.ans$wasserstein)+1)/(length(curr.ans$permvals)+1)
                     return(curr.ans)
                   })
  
  return(answer)
  
}

permutation_test_two_samples <- function(df1,df2,iterations,p,q,dim,format,standardize){
  
  # function to test whether or not two sets of labelled data come from the same geometric process
  # df1 and df2 are the two datasets
  # iterations is the number of permutations we will calculate for class labels
  # q is the finite exponent used in distance calculations between barcodes, q >= 1
  # dim is the maximum dimension in which to calculate homology for barcodes
  # p is the finite wasserstein distance parameter, p >= 1
  
  # perform all error checks, return from function if error found
  error_check(obj1 = df1, obj2 = df2, iterations = iterations,
              p = p, q = q, dim = dim, format = format,
              standardize = standardize, input_type = "sample")
  
  # calculate barcode for each dataset
  barcode1 = calculate_homology(as.matrix(df1), format = format, standardize = standardize, dim = dim,threshold = enclosing_radius(df1))
  barcode2 = calculate_homology(as.matrix(df2), format = format, standardize = standardize, dim = dim,threshold = enclosing_radius(df2))
  
  # get group sizes
  n1 = nrow(df1)
  n2 = nrow(df2)
  
  # compute loss function on observed data
  test_loss = d(B1 = as.data.frame(barcode1),B2 = as.data.frame(barcode2),dim = dim,p = p)^q
  
  # get permutation values
  perm_values = lapply(X = 1:iterations,FUN = function(X){
    
    # sample two datasets of size n1 and n2 from df1 union df2
    ind = sample(1:(n1+n2),size = n1,replace = FALSE)
    ind1 = ind[which(ind <= n1)]
    ind2 = setdiff(ind,ind1) - n1
    sample1 = rbind(df1[ind1,],df2[ind2,])
    sample2 = rbind(df1[setdiff(1:n1,ind1),],df2[setdiff(1:n2,ind2),])
    barcode_sample1 = calculate_homology(as.matrix(sample1),format = format,standardize = standardize,dim = dim,threshold = enclosing_radius(X = sample1))
    barcode_sample2 = calculate_homology(as.matrix(sample2),format = format,standardize = standardize,dim = dim,threshold = enclosing_radius(X = sample2))
    
    # return loss function
    ret_loss = t(as.matrix(unlist(d(B1 = as.data.frame(barcode_sample1),B2 = as.data.frame(barcode_sample2),dim = dim,p = p)^q)))
    return(ret_loss)
    
  })
  perm_values = do.call(rbind,perm_values)
  
  # organize results by dimension
  answer <- lapply(X = 0:dim,
                   FUN = function(curr.dim) {
                     curr.ans <- list()
                     curr.ans$dimension <- curr.dim
                     curr.ans$permvals <- perm_values[,curr.dim + 1] # this is the distance between the barcodes in the current dimension across each iteration
                     curr.ans$wasserstein <- test_loss[[curr.dim + 1]] # test statistic distance between barcodes in this dimension
                     curr.ans$pvalue <- (sum(curr.ans$permvals > curr.ans$wasserstein)+1)/(length(curr.ans$permvals)+1)
                     return(curr.ans)
                   })
  
  return(answer)
  
}

#NEED TO ADD PARAM ERROR CHECKING OR DEFAULT VALUES
permutation_test <- function(param1, param2, iterations, p = 2, q = 2, dim,
                             format = "cloud", standardize = FALSE, type){
  # wrapper function for the two kinds of tests
  if (type == "samples") {
    permutation_test_two_samples(df1 = param1, df2 = param2, iterations = iterations,
                                 p = p, q = q, dim = dim, format = format,
                                 standardize = standardize)
  } else if (type == "groups") {
    permutation_test_two_groups(l1 = param1, l2 = param2, iterations = iterations,
                                p = p, q = q, dim = dim, format = format,
                                standardize = standardize)
  }
}

# TEST ----
set.seed(42)
l1 = list(data.frame(circleUnif(n = 20)), data.frame(2*circleUnif(n = 15)))
l2 = list(data.frame(circleUnif(n = 20)), data.frame(1.5*circleUnif(n = 20)))
iterations = 2L
dim = 1L
test1 = permutation_test(l1, l2, iterations = iterations,
                         dim = dim, type = "groups")

df1 = l1[[1]]
colnames(df1) <- c("x","y")
df2 = data.frame(x = 1:20,
                 y = runif(n = 20, min = 0, max = 1))
test2 = permutation_test(df1, df2, iterations = iterations,
                         dim = dim, type = "samples")
