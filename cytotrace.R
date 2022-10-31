library(CytoTRACE)
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){

  int k = z.size() ;

  IntegerMatrix  mat(nrows, ncols);

  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }

  return mat;
}
' )


as_matrix <- function(mat){
  
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
results <- CytoTRACE(as_matrix(a[discard_gene2(rownames(a)),mac.metatable$Well_ID[mac.metatable$time != "Adult"]]))

mac_prog_metatable<-metatable[metatable$major%in%c("macrophage","progenitor"),]
results_prog_mf <- CytoTRACE(as_matrix(a[discard_gene2(rownames(a)),mac_prog_metatable$Well_ID[mac_prog_metatable$time != "Adult"]]))
