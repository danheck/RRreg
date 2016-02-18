

checkPerfectSeparation <- function(y,x,n.response){
  for(i in 1:ncol(x)){
    if(length(unique(x[,i])) != 1){
      xorder <- order(x[,i])
      checkSeparation <- all(y[xorder] == sort(y[xorder])) | 
        all(rev(y[xorder]) == sort(y[xorder]))
      if(checkSeparation)
        warning("A covariate (i.e., ",colnames(x)[i], ") perfectly separates",
                "\n  responsess of the RR-variable, which will probably lead to unstable estimates.")
      
    }
  }
  
}