## Run PCA
pc_function <- function(exp){
  svd(exp - rowMeans(exp))
}

## Plot PCA
pc_plot <- function(pca,legend=NULL, color=NULL, main=NULL,
                    type=NULL,
                    ptsize=1.5,position="bottomright"){
  par(mfrow=c(2,2))

  if(type=='character'){
    color = colors
    c2 = names(table(color))
  }else{
    color = as.factor(color)
    c2 = 1:length(unique(as.factor(color)))
  }

  par(font.lab = 2, cex.lab = 1.2, font.axis = 2, cex.axis = 1.2)
  plot(pca$v[, 1], pca$v[, 2], col= color,
       pch = 19, cex = ptsize,
       xlab = 'PC1', ylab = 'PC2',
       main = main)

  par(font.lab = 2, cex.lab = 1.2, font.axis = 2, cex.axis = 1.2)
  plot(pca$v[, 3], pca$v[, 4], col= color,
       pch = 19, cex = ptsize,
       xlab = 'PC3', ylab = 'PC4',
       main = main)

  par(font.lab = 2, cex.lab = 1.2, font.axis = 2, cex.axis = 1.2)
  plot(pca$v[, 5], pca$v[, 6], col= color ,
       pch = 19, cex = ptsize,
       xlab = 'PC5', ylab = 'PC6',
       main = main)

  par(font.lab = 2, cex.lab = 1.2, font.axis = 2, cex.axis = 1.2)
  plot(pca$v[, 7], pca$v[, 8], col= color,
       pch = 19, cex = ptsize,
       xlab = 'PC7', ylab = 'PC8',
       main = main)
  legend(position, pch = 19, col= c2,
         names(summary(as.factor(legend))),bty="n")
}

## Plot Variance Explained
var_plot <- function(pca){
  par(mfrow=c(1,1))
  plot((pca$d^2/sum(pca$d^2))*100, xlim = c(0, 15), type = "b",
      pch = 16, xlab = "principal components",
      ylab = "variance explained (%)")
}
