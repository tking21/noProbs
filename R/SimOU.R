#' This function simulates trait data under OU model and calculated phylogenetic independent contrast for each simulation
#'
#' SimOU simulates trait evolution along a tree using fitted parameters for rate, alpha parameter, and root value. It then calculated pylogenetic independent contrast for each node of the phylogeny across all simulations.
#'
#'
#' @param means vector of class 'numeric' with continuous trait values for each species
#' @param tree of class 'phylo'
#' @param n interger for number of simulations to perform
#' @param species vector of class 'character' with corresponding species names
#'
#' @return data matrix of phylogenetic independent contrasts calculated for each node in the tree. Number of elements in matrix is equvilent to number of simulations performed.
#'
#' @example
#'#' trait_data <- c(40,57,29,31,31,32,53,40,41,31)
#' my_tree <- sample_tree.RData
#' n = 10
#' tip_names <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8","s9", "s10")
#'simOU(trait_data, my_tree, n, tip_names)
#'
#'
#'@export

simOU <- function(means,tree, n, species) {


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


  #colnames(means) <- c("value")
  d <- matrix(means)
  colnames(d) <- c("value")

  rownames(d) <- species
  fitParams <- fitContinuous(tree, d, model="OU", bounds = list(alpha = c(0, 1000)))

  #make empty vector to hold results in, what is returned
  results <- vector("list", n)

  for (x in 1:n){
    results[[x]] <- rTraitCont(tree, model ="OU", sigma = fitParams$opt[2], alpha = fitParams$opt[1] , theta = fitParams$opt[3], root.value = mean(means))
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





