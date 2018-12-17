#' @title Multivariate Shapiro-Wilk Test for Multivariate Normality
#'
#' @description
#' \code{mvShapiro} function tests whether the data set has multivariate 
#' normality distribution or not.
#'
#' @details
#' This function computes the test statistic and p-value of the Shapiro-Wilk
#' test for multivariate normality proposed by Villasenor-Alva and
#' GonzalezEstrada (2009). Moreover, it can perform normality test for all
#' groups in grouped datasets. The hypotheses are defined as H0:The dataset comes
#' from multivariate normal distribution and H1:The dataset does not
#' come from multivariate normal distribution.
#'
#' @importFrom stats cor cov pchisq pf pnorm qchisq qf shapiro.test var
#' @param data a data frame.
#' @param group a logical argument. If \code{group="TRUE"}, the normality 
#' tests are performed for the groups in data.
#' When \code{group=TRUE}, \code{G=NULL} cannot be. default 
#' \code{group="FALSE"} 
#' @param G a group vector. default \code{G=NULL} 
#'
#' @export
#' @return a list with 2 elements:
#' \item{Stat}{The value of Test Statistic(s)}
#' \item{p.value}{p value(s)}
#' @references Villasenor Alva, J. A., & Estrada, E. G. (2009). 
#' A generalization of Shapiro Wilk's test for multivariate normality. 
#' Communications in Statistics Theory and Methods, 38(11), 1870-1883.
#' @references Elizabeth Gonzalez Estrada and Jose A. Villasenor-Alva (2013).
#' mvShapiroTest: Generalized Shapiro Wilk test for multivariate   normality.
#' R package version 1.0.
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#' @examples
#' 
#' data(iris) 
#' # For raw data
#' result <- mvShapiro(data=iris[,1:4])
#' summary(result)
#' #For grouped data
#' result.group <- mvShapiro(data=iris[,1:4],group=TRUE,G=iris[,5])
#' summary(result.group)

 mvShapiro<-function(data,group=FALSE,G=NULL){
 Name<-"mvShapiro"
 data<-as.matrix(data)
 n<-nrow(data)
 p<-ncol(data)
 
if (n < 12 || n > 5000) 
        stop("Sample size must be between 12 and 5000.")
    if (n <= p) 
        stop("Sample size must be larger than vector dimension.")
    
        if (group=="FALSE") {
        z <- scale(data, scale = FALSE)
        eigenv <- eigen(var(data), symmetric = TRUE)
        eigen.val<-eigenv$values
        e.vec <- as.matrix(eigenv$vectors)
        sqrS <- e.vec %*% diag(1/sqrt(eigen.val), ncol = p) %*%t(e.vec)
        Z <- t(sqrS %*% t(z))
        Wi <-NULL
        for (i in 1:p) {
        Wi[i]<-shapiro.test(Z[,i])$statistic
        }
        W <- mean(Wi)
        y <- log(n)
        w1 <- log(1 - W)
        m <- -1.5861 - 0.31082 * y - 0.083751 * y^2 + 0.0038915 *y^3
        s <- exp(-0.4803 - 0.082676 * y + 0.0030302 * y^2)
        s2 <- s^2
        sigma2 <- log((p - 1 + exp(s2))/p)
        mu1 <- m + s2/2 - sigma2/2
        pval <- pnorm(w1, mean = mu1, sd = sqrt(sigma2), lower.tail = FALSE)
        results <- list(Stat= W, p.value = pval,Test = Name,group="FALSE")
        
        } 
        if (group=="TRUE")  {
        G<-as.factor(G)
        Levels<-levels(G) ; g<-length(Levels); pval=Wast=NULL
       
        for (i in 1:g) {
        data.group<-data[which(G==Levels[i]),]
	  z <- scale(data.group, scale = FALSE)
        eigenv <- eigen(var(data.group), symmetric = TRUE)
        eigen.val<-eigenv$values
        e.vec <- as.matrix(eigenv$vectors)
        sqrS <- e.vec %*% diag(1/sqrt(eigen.val), ncol = p) %*%t(e.vec)
        Z <- t(sqrS %*% t(z))
        Wi <-NULL
        for (j in 1:p) {
        Wi[j]<-shapiro.test(Z[,i])$statistic
        }
        W <- mean(Wi)
        y <- log(n)
        w1 <- log(1 - W)
        m <- -1.5861 - 0.31082 * y - 0.083751 * y^2 + 0.0038915 *y^3
        s <- exp(-0.4803 - 0.082676 * y + 0.0030302 * y^2)
        s2 <- s^2
        sigma2 <- log((p - 1 + exp(s2))/p)
        mu1 <- m + s2/2 - sigma2/2
        pval[i] <- pnorm(w1, mean = mu1, sd = sqrt(sigma2), lower.tail = FALSE)
        Wast[i]<-W
        }
        
        PVAL<-data.frame(GROUPS=Levels, P.Values=pval)
        WAST<-data.frame(GROUPS=Levels, Statistic=Wast)

        results <- list(Stat= WAST, p.value = PVAL,Test = Name,group="TRUE")        
        }


        class(results)<-c("MVTests","list")
 	  return(results) 
}

