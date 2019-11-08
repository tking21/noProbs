#' This function simulates trait data under BM model
#'
#' simBM simulates trait evolution along a tree using fitted parameters for rate and root and calculates the phylogenetic contrast along the tree for each simulation.
#'
#' @param means matrix of mean trait values with corresponding species of interest
#' @param tree of class 'phylo'
#' @param n interger for number of simulations
#' @param species  vector of class 'character' which lists the species in the tree
#'
#' @return data matrix of phylogenetic independent contrasts calculated for each node in the tree. Number of elements in matrix is equvilent to number of simulations performed.
#'
#' @example
#' trait_data <- c(40,57,29,31,31,32,53,40,41,31)
#' my_tree <- sample_tree.RData
#' n = 10
#' tip_names <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8","s9", "s10")
#'simBM(trait_data, my_tree, n, tip_names)
#'
#'@export


simBM <- function(means,tree, n, species) {

  #if data is not continuous
  if (class(means) != 'numeric'){
    warning("Argument 'means' must be a vector of class 'numeric'")
    stop()
  }


  #if # of species and # of values in means do not match
  if( length(means) !=  length(species)){
    warning("Number of trait values in argument 'means' does not match length of augument 'species'")
    stop()
  }

  #if 'n' is not an integer
  if (n %% 1 != 0){
    warning("Argument 'n' must be an interger of class 'numeric'")
    stop()
  }

  if(class(species) != 'character'){
    warning("Argument 'species' must be a vector of class 'character'")
    stop()
  }

  if (class(tree) != 'phylo'){
    warning("Argument 'tree' must be a phylogenetic tree of class 'phylo'")
    stop()
  }

  results <- vector("list", n)

  d <- matrix(means)
  colnames(d) <- c("value")
  rownames(d) <- species

  #estimates the evolutionary/phylogenetic variance-covariance matrix
  Q <- ratematrix(tree, d)

  #simulate data along tree; root=mean trait value across all species
  sim_results<- sim.char(tree, Q, n, model = "BM", root=mean(means))


  for (x in 1:n){
    results[[x]] <- sim_results[,,x]
  }

  #return(results)

  for (x in 1:n){
    pic <- pic(results[[x]], tree )
    if (x == 1){
      trait_pic <- pic
    }
    else{
      trait_pic = rbind(trait_pic, pic)
    }
  }
  return(trait_pic)


}
