cor.plot <- function(x, y, xlab, ylab, lim, main) {
    smoothScatter(x, y, colramp=white.blue, xlim=lim, ylim=lim, main=main, xlab=xlab, ylab=ylab)
    fit <- lm(y~x)
    abline(fit,col='red',lwd=2)
    r2 <- summary(fit)$r.squared
    text(lim[1],lim[2]-0.5,pos=4,paste("Linear Regression R^2", round(r2, 3), sep="="),col='red')
    cor <- cor.test(x, y, method="spearman")$estimate
    text(lim[1],lim[2]-1,pos=4,paste("Spearman Rho", round(cor, 3), sep="="),col='red')
    N <- length(x)
    text(lim[1],lim[2]-1.5,pos=4,paste("N", N, sep="="),col='red')
    list(r2=r2, cor=cor, N=N)
}