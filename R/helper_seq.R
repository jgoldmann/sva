####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
  weights <- pos_res <- list()
  for(i in 1:nrow(dat)){
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    #LH <- sapply(1:G_sub, function(j){sum(log2(dnbinom(x, mu=m[j,], size=1/phi_sub[j])+1))})  
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
    LH[is.nan(LH)]=0; 
    if(sum(LH)==0 | is.na(sum(LH))){
      pos_res[[i]] <- c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i]))
    }else{
      pos_res[[i]] <- c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
    }
    
    weights[[i]] <- as.matrix(LH/sum(LH))
  }
  pos_res <- do.call(rbind, pos_res)
  weights <- do.call(cbind, weights)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"], weights=weights)	
  return(res)
} 


####  Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
  new_counts_sub <- 
    bplapply(1:length(counts_sub),
             match_quantiles_inner,
             counts_sub, old_mu, old_phi, new_mu, new_phi)
  new_counts_sub_mx <- 
    matrix(as.numeric(new_counts_sub), nrow = nrow(counts_sub))
  return(new_counts_sub_mx)
}

match_quantiles_inner <- 
  function(n, counts_sub, old_mu, old_phi, new_mu, new_phi) {
  if(counts_sub[n] <= 1){
    res <- counts_sub[n]
  }else{
    index <- ifelse(n%%length(old_phi)==0, length(old_phi), n%%length(old_phi))
    tmp_p <- pnbinom(counts_sub[n]-1, mu=old_mu[n], size=1/old_phi[index])
    if(abs(tmp_p-1)<1e-4){
      res <- counts_sub[n]  
      # for outlier count, if p==1, will return Inf values -> use original count instead
    }else{
      index <- ifelse(n%%length(new_phi)==0, length(new_phi), n%%length(new_phi))
      res <- 1+qnbinom(tmp_p, mu=new_mu[n], size=1/new_phi[index])
    }
  }
  return(res)
}




mapDisp <- function(old_mu, new_mu, old_phi, divider){
  new_phi <- matrix(NA, nrow=nrow(old_mu), ncol=ncol(old_mu))
  for(a in 1:nrow(old_mu)){
    for(b in 1:ncol(old_mu)){
      old_var <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
      new_var <- old_var / (divider[a, b]^2)
      new_phi[a, b] <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
    }
  }
  return(new_phi)
}
