## PCA exploration plots from Shannon Ellis taken from
## http://research.libd.org/recount-brain/example_multistudy/recount_brain_multistudy.html#run-pca

## Check http://genomicsclass.github.io/book/pages/pca_svd.html
## for differences/similarities between svd() and prcomp()

## Run PCA
pc_function <- function(exp){
  svd(exp - rowMeans(exp))
}

## Plot PCA
pc_plot <- function(pca,legend=NULL, color=NULL, main=NULL,
                    type=class(color),
                    ptsize=1.5,position="bottomright",
                    legend.pan = 4, ...){
  var <-  (pca$d^2/sum(pca$d^2))*100
  var_lab <- function(i) {
      paste0('PC', i, ': ', round(var[i], 2), '%')
  }
  par(mfrow=c(2,2))

  if(type=='character'){
    color = colors
    c2 = names(table(color))
  }else{
    color = as.factor(color)
    c2 = 1:length(unique(as.factor(color)))
  }

  make_len <- function() {
       legend(position, col= c2,
         names(summary(as.factor(legend))),bty="n", ...)
  }

  par(font.lab = 2, cex.lab = 1.2, font.axis = 2, cex.axis = 1.2)
  plot(pca$v[, 1], pca$v[, 2], col= color,
       pch = 19, cex = ptsize,
       xlab = var_lab(1), ylab = var_lab(2),
       main = main)
  if(legend.pan == 1) make_len()

  par(font.lab = 2, cex.lab = 1.2, font.axis = 2, cex.axis = 1.2)
  plot(pca$v[, 3], pca$v[, 4], col= color,
       pch = 19, cex = ptsize,
       xlab = var_lab(3), ylab = var_lab(4),
       main = main)
  if(legend.pan == 2) make_len()

  par(font.lab = 2, cex.lab = 1.2, font.axis = 2, cex.axis = 1.2)
  plot(pca$v[, 5], pca$v[, 6], col= color ,
       pch = 19, cex = ptsize,
       xlab = var_lab(5), ylab = var_lab(6),
       main = main)
  if(legend.pan == 3) make_len()

  par(font.lab = 2, cex.lab = 1.2, font.axis = 2, cex.axis = 1.2)
  plot(pca$v[, 7], pca$v[, 8], col= color,
       pch = 19, cex = ptsize,
       xlab = var_lab(7), ylab = var_lab(8),
       main = main)
  if(legend.pan == 4) make_len()
}

## Plot Variance Explained
var_plot <- function(pca){
  par(mfrow=c(1,1))
  var <-  (pca$d^2/sum(pca$d^2))*100
  print(round(var[1:10], 2))
  plot(var, xlim = c(0, 15), type = "b",
      pch = 16, xlab = "principal components",
      ylab = "variance explained (%)", ylim = c(0, round(max(var) / 10 + .5, 0) * 10))
}
