iimpute = function (raw_count, out_file_path, n = NULL, lam = 0) 
{
  "Hi! I am really sorry that I do not have a mic here. So I will do 
  the demo and I have written what you may need to know about the demo.
  
  My project is more about algorithm, so there is a little boring to watch 
  this demo.
  
  First, the user should have installed R and two packages, penalized and 
  kernlab. The method that I want to show is call iimpute.
  
  I wrap it in a single function so that it is easier to use. There are 
  five parameters, count_path and out_dir must be input.
  
  Once the file path and output path are decided, you can run the function 
  by loading it into your R enviornment. 
  
  It may takes some time for running, and the program will show what step 
  it is doing.
  
  Detailed description can be found in my final paper. 
  
  Thank you!"
  
  library(penalized)
  library(kernlab)
  
  raw_count = as.matrix(raw_count)
  print(paste("number of genes in raw count matrix", nrow(raw_count)))
  print(paste("number of cells in raw count matrix", ncol(raw_count)))
  
  totalCounts_by_cell = colSums(raw_count)
  totalCounts_by_cell[totalCounts_by_cell == 0] = 1
  raw_count = sweep(raw_count, MARGIN = 2, 10^6/totalCounts_by_cell, FUN = "*")
  if (min(raw_count) < 0) {
    stop("smallest read count cannot be negative!")
  }
  raw_count=raw_count[which(rowSums(raw_count)!=0),]
  raw_count=raw_count[,which(colSums(raw_count)!=0)]
  count_lnorm = log10(raw_count + 1.01)
  print("reading finished!")
  if(is.null(n)){
    n = ncol(raw_count)/2
  }
  genenames = rownames(count_lnorm)
  cellnames = colnames(count_lnorm)
  count = count_lnorm
  point = log10(1.01)
  drop_thre = 0.5
  count = as.matrix(count)
  I = nrow(count)
  J = ncol(count)
  count_imp = count
  
  count_nzero = lapply(1:I, function(i) setdiff(count[i, ], log10(1.01)))
  mu = sapply(count_nzero, mean)
  mu[is.na(mu)] = 0
  sd = sapply(count_nzero, sd)
  sd[is.na(sd)] = 0
  cv = sd/mu
  cv[is.na(cv)] = 0
  # sum(mu >= 1 & cv >= quantile(cv, 0.25), na.rm = TRUE)
  high_var_genes = which(mu >= 1 & cv >= quantile(cv, 0.25))
  if(length(high_var_genes) < 500){ 
    high_var_genes = 1:I}
  count_hv = count[high_var_genes, ]
  
  

  print("inferring cell similarities ...")
  
  
  set.seed(n)
  npc=100
  ## dimeansion reduction
  print("dimension reduction ...")
  if(J < 5000){
    var_thre = 0.4
    pca = prcomp(t(count_hv))
    eigs = (pca$sdev)^2
    var_cum = cumsum(eigs)/sum(eigs)
    if(max(var_cum) <= var_thre){
      npc = length(var_cum)
    }else{
      npc = which.max(var_cum > var_thre)
    }
  }else{
    var_thre = 0.6
    pca = rpca(t(count_hv), k = 1000, center = TRUE, scale = FALSE) 
    eigs = (pca$sdev)^2
    var_cum = cumsum(eigs)/sum(eigs)
    if(max(var_cum) <= var_thre){
      npc = length(var_cum)
    }else{
      npc = which.max(var_cum > var_thre)
    }
  }
  if (npc < 3){ npc = 3 }
  mat_pcs = t(pca$x[, 1:npc]) # columns are cells
  
  ## detect outliers
  print("calculating cell distances ...")
  dist_cells_list = lapply(1:J, function(id1){
    d = sapply(1:id1, function(id2){
      sse = sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
      sqrt(sse)
    })
    return(c(d, rep(0, J-id1)))
  })
  dist_cells = matrix(0, nrow = J, ncol = J)
  for(cellid in 1:J){dist_cells[cellid, ] = dist_cells_list[[cellid]]}
  dist_cells = dist_cells + t(dist_cells)
  
  min_dist = sapply(1:J, function(i){
    min(dist_cells[i, -i])
  })
  iqr = quantile(min_dist, 0.75) - quantile(min_dist, 0.25)
  outliers = which(min_dist > 1.5 * iqr + quantile(min_dist, 0.75))
  
  non_out = setdiff(1:J, outliers)
  
  width = n
  x=t(mat_pcs[, non_out])
  km = as.matrix(dist(x, method = "euclidean"))
  
  for(i in 1:nrow(km)){
    
    wid=km[i,order(km[i,])[width]]
    for(k in 1:ncol(km)){
      
      if(km[i,k]<=wid){
        km[i,k]=exp(-(km[i,k]^2)/(wid^2))
      }
      else{
        km[i,k]=0
      }
    }
    km[i,i]=mean(km[i,-c(i,which(km[i,-i]==0))])
  }
  
  
  count1 = as.matrix(count)
  parslists = lapply(1:nrow(count1), function(ii) {
    
    xdata = count1[ii, ]
    inits = rep(0, 5)
    inits[1] = sum(xdata == point)/length(xdata)
    if (inits[1] == 0) {inits[1] = 0.01}
    inits[2:3] = c(0.5, 1)
    xdata_rm = xdata[xdata > point]
    inits[4:5] = c(mean(xdata_rm), sd(xdata_rm))
    if (is.na(inits[5])) {inits[5] = 0}
    paramt = inits
    eps = 10
    iter = 0
    loglik_old = 0
    try(while(eps > 0.5) {
      print(params[0])
      pz1 = paramt[1] * dgamma(xdata, shape = paramt[2], rate = paramt[3])
      pz2 = (1 - paramt[1]) * dnorm(xdata, mean = paramt[4], sd = paramt[5])
      pz = pz1/(pz1 + pz2)
      pz[pz1 == 0] = 0
      wt=cbind(pz, 1 - pz)

      paramt[1] = sum(wt[, 1])/nrow(wt)
      paramt[4] = sum(wt[, 2] * xdata)/sum(wt[, 2])
      paramt[5] = sqrt(sum(wt[, 2] * (xdata - paramt[4])^2)/sum(wt[, 2]))
      
      
      
      wt1=wt[,1]
      tp_s = sum(wt)
      tp_t = sum(wt * xdata)
      tp_u = sum(wt * log(xdata))
      tp_v = log(tp_t / tp_s) - tp_u / tp_s
      
      if (tp_v <= 0){
        alpha = 20
      }else{
        alpha = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
        if (alpha >= 20){alpha = 20
        }
      }
      
      beta = tp_s / tp_t * alpha
      paramt[2:3]=c(alpha, beta)
      
      
      
      loglik = sum(log10(paramt[1] * dgamma(xdata, shape = paramt[2], rate = paramt[3]) + (1 - paramt[1]) * dnorm(xdata, mean = paramt[4], sd = paramt[5])))
      
      
      eps = (loglik - loglik_old)^2
      loglik_old = loglik
      iter = iter + 1
      if (iter > 100) 
        break
    }, silent = TRUE)
    
    
    if (class(paramt) == "try-error"){
      paramt = rep(NA, 5)
    }
    return(paramt)
  })
  
  
  
  
  
  parslist = Reduce(rbind, parslists)
  colnames(parslist) = c("rate", "alpha", "beta", "mu", "sigma")
  point = log10(1.01)
  valid_genes = which(complete.cases(parslist) )
  if(length(valid_genes) <= 10){ next }
  # find out genes that violate assumption
  mu = parslist[, "mu"]
  sgene1 = which(mu <= log10(1+1.01))
  # sgene2 = which(mu <= log10(10+1.01) & mu - parslist[,5] > log10(1.01))
  
  dcheck1 = dgamma(mu+1, shape = parslist[, "alpha"], rate = parslist[, "beta"])
  dcheck2 = dnorm(mu+1, mean = parslist[, "mu"], sd = parslist[, "sigma"])
  sgene3 = which(dcheck1 >= dcheck2 & mu <= 1)
  sgene = union(sgene1, sgene3)
  valid_genes = setdiff(valid_genes, sgene)
  
  
  
  for(cc in c(1:length(non_out))){
    gc()
    print(paste("Impute dropout values for cell", cc, "..."))
    round=cc
    weigth=km[cc,]
    
    subcount = count[valid_genes, non_out[cc] , drop = FALSE]
    
    Ic = length(valid_genes)
    parslists = parslist[valid_genes, , drop = FALSE]
    
    
    
    droprates = t(sapply(1:Ic, function(i) {
      paramt=parslists[i, ]
      x2=subcount[i, ]
      pz1 = paramt[1] * dgamma(x2, shape = paramt[2], rate = paramt[3])
      pz2 = (1 - paramt[1]) * dnorm(x2, mean = paramt[4], sd = paramt[5])
      pz = pz1/(pz1 + pz2)
      pz[pz1 == 0] = 0
      wt=cbind(pz, 1 - pz)
      return(wt[, 1])
    }))
    print(droprates)
    mucheck = sweep(subcount, MARGIN = 1, parslists[, "mu"], FUN = ">")
    droprates[intersect(which(mucheck),which(droprates>drop_thre))] = 0
    droprate=1-droprates
    nbs = setdiff(c(1:length(non_out)), cc)
    yimpute = rep(0, Ic)
    count_xy=count[valid_genes,]
    count_xy=sweep(count_xy, 1, droprate, "*")
    xx = count_xy[,non_out[nbs]]
    yy = count_xy[,non_out[cc]]
    
    xx=sweep(xx, 2, weigth[-cc], "*")
    ximpute = count[valid_genes, non_out[nbs]]
    ximpute=sweep(ximpute, 2, weigth[-cc], "*")
    set.seed(cc)
    nnls = penalized(yy, penalized = xx, unpenalized = ~0,
                     positive = TRUE, lambda1 = lam, lambda2 = 0, 
                     maxiter = 3000, trace = FALSE)
    ynew = penalized::predict(nnls, penalized = ximpute, unpenalized = ~0)[,1]
    y = ynew
    maxobs = apply(count[valid_genes,], 1, max)
    y[y > maxobs] = maxobs[y > maxobs]
    
    
    if (class(y) == "try-error") {
      # print(y)
      y = subcount[, non_out[cc], drop = FALSE]
    }
    rewrite_gene=valid_genes[which(droprates>=drop_thre)]
    nonrewrite_gene=setdiff(c(1:nrow(count_imp)),rewrite_gene)
    
    count_imp[valid_genes, non_out[cc]] = y
    count_imp[nonrewrite_gene, non_out[cc]] = count[nonrewrite_gene, non_out[cc]]
    
  }
  
  
  
  count_imp[count_imp < point] = point
  
  count_imp = 10^count_imp - 1.01
  count_imp = sweep(count_imp, MARGIN = 2, totalCounts_by_cell/10^6, 
                    FUN = "*")
  write.csv(count_imp,file = out_file_path)
  
}
