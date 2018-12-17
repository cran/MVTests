#' @title One Way Multivariate Analysis of Variance (MANOVA)
#'
#' @description
#' \code{Manova} function computes one-way MANOVA test and gives confidence 
#' intervals
#'
#' @details
#' This function computes the MANOVA test for more than two independent 
#' samples with and without the assumption that covariance matrices are 
#' homogeneity. The hypotheses are \code{H0:Mu_1=Mu_2=...=Mu_g} and
#' \code{H1:At least a Mu_j is different from others (j=1,2,...,g)}.
#' When \code{H0} is rejected, this function computes confidence intervals 
#' for all variables to determine variable(s) affecting on rejection decision.
#' 
#' @importFrom stats cor cov pchisq pf pnorm qchisq qf shapiro.test var
#' @param data a data frame.
#' @param group a group vector
#' @param method The method that will be used for MANOVA
#' @param CI a logical argument. If \code{CI="TRUE"}, then the decision 
#' matrices based on confidence intervals obtained for all variables and 
#' groups are calculated.
#' CI values are calculated only in Roy method. \code{default CI=FALSE}. 
#' @param alpha Significance Level that will be used for confidence intervals. \code{default alpha=0.05}
#' @param Homogeneity a logical argument. If \code{Homogeneity=TRUE}, 
#' then classical MANOVA methods which are under the assumption that 
#' covariance matrices of groups are homogeneity are used.
#' Otherwise, the MANOVA approach (James, 1954) without homogeneous 
#' covariance matrices assumption is used. Default \code{Homogeneity=TRUE}. 
#' The homogeneity test for covariance matrices is performed \code{BoxM} function.
#' @export
#'
#' @return a list with 9 elements:
#' \item{Method}{The method used for MANOVA}
#' \item{Test.Stat}{The value of the test statistic for the selected method}
#' \item{App.Stat}{The approximate value of F or Chi-Squared statistic. 
#' The \code{Manova} function uses F statistic when the covariance matrices 
#' are homogeneity; otherwise, it uses Chi-squared Statistic.}
#' \item{df}{The F or Chi-Squared statistic's degree of freedom}
#' \item{p.value}{p value}
#' \item{Decision.Matrices}{When \code{method="Roy"} and \code{CI="TRUE"},  it gives the 
#' decision matrices based on confidence intervals obtained for all variables 
#' and groups.
#' If between groups decision is TRUE for any variable, it means that 
#' the mean vectors of these groups are statistically different.
#' When the covariance matrices are not homogeneity, confidence intervals 
#' are not calculated.}
#' \item{alpha}{The alpha value using in confidence intervals}
#' \item{Descriptive}{Descriptive Statistics for all groups}
#' \item{Homogeneity}{It gives the result of the assumption that covariance
#' matrices are homogeneity.}
#'
#' @references Rencher, A. C. (2003). Methods of multivariate analysis 
#' (Vol. 492). John Wiley & Sons.
#' @references Tatlidil, H. (1996). Uygulamali Cok Degiskenli Istatistiksel 
#' Yontemler. Cem Web.
#' @references James, G. S. (1954). Tests of linear hypotheses in univariate 
#' and multivariate analysis when the ratios of the population variances are unknown. Biometrika, 41(1/2), 19-43.
#' @references Tsagris M. T. (2014). Multivariate statistical functions 
#' in R. Athens, Nottingham and Abu Halifa (Kuwait).
#'
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#' @examples
#' 
#' data(iris)
#' 
#' # Wilk's Test
#' results.Wilks <- Manova(data=iris[,1:4],group=iris[,5],alpha=0.01)
#' summary(results.Wilks)
#' 
#' # Pillai's Test
#' results.Pillai <- Manova(data=iris[,1:4],group=iris[,5],method="Pillai")
#' summary(results.Pillai)
#'  
#' # Hotelling and Lawley's Test
#' results.HL <- Manova(data=iris[,1:4],group=iris[,5], method="Hotelling-Lawley")
#' summary(results.HL)
#' 
#' # Roy's Test
#' results.Roy <- Manova(data=iris[,1:4],group=iris[,5], method="Roy")
#' summary(results.Roy)
#' results.RoyCI <- Manova(data=iris[,1:4],group=iris[,5], method="Roy",CI=TRUE)
#' summary(results.RoyCI)
#'
#' # James's Test
#' results.James <- Manova(data=iris[,1:4],group=iris[,5],Homogeneity=FALSE)
#' summary(results.James)
 
 Manova<-function(data,group,method="Wilks",CI="FALSE",alpha=0.05,Homogeneity=TRUE) {
 Name<-"MANOVA"
 if (CI=="TRUE" && method!="Roy"){
	stop("Confidence Intervals might be only calculated by using the Roy Method.")
 }
 
 if (Homogeneity==TRUE) {
 group<-as.factor(group)
 Levels<-levels(group)
 X<-cbind(data,group)
 p<-ncol(data);n<-nrow(data);g=length(levels(group));VE<-n-g; VH<-g-1
 N.star<-(VE-p-1)/2; s<-min(VH,p); r<-max(VH,p); m<-(abs(VH-p)-1)/2
 
 ## Calculation B and W
 ns<-NULL;X.mean.group<-matrix(c(rep(0,p*g)),nrow=p,ncol=g)
 W<-B<-matrix(rep(0,p*p),p,p)
 for (i in 1:g) {
 ns[i]<-length(which(group==Levels[i]))
 X.mean<-colMeans(data)
 X.mean.group[,i]<-colMeans(data[which(group==Levels[i]),])
 W<-W+(ns[i]-1)*cov(data[which(group==Levels[i]),])
 B<-B+ns[i]*(X.mean.group[,i]-X.mean)%*%t((X.mean.group[,i]-X.mean))
 }
 WB<-solve(W)%*%B
 eigens<-eigen(WB)$values


 ##Descriptive Statistics
 Variables<-colnames(data)
 
 Descriptive<-NULL
 for (A in 1:g) { 
 MU<-colMeans(data[which(group==Levels[A]),]) 
 S<-cov(data[which(group==Levels[A]),]) 
 Desc<-rbind(MU,sqrt(diag(S)))
 colnames(Desc)<-Variables
 rownames(Desc)<-c("Means","Sd")
 Descriptive[[Levels[A]]]<-Desc
 }
 
 ##Wilk's Method
 if (method=="Wilks") {
 lamda<-det(W)/det(W+B)
 t<-sqrt((p^2*VH^2-4)/(p^2+VH^2-5))
 w<-VE+VH-0.5*(p+VH+1)
 df1<-p*VH; df2<-w*t-0.5*(p*VH-2)
 F<-((1-lamda^(1/t))/lamda^(1/t))*(df2/df1)
 df<-c(df1,df2)
 pval<-pf(q=F, df1=df1,df2=df2,lower.tail=FALSE)
 Result<-list(Method="Wilks",Test.Stat=lamda,App.Stat=F,df=df,p.value=pval,Descriptive=Descriptive,Homogeneity=TRUE,Test=Name,CI=FALSE)
 } 

 ## Hotelling-Lawley Method
 if (method=="Hotelling-Lawley") {
 T0.sq<-sum(eigens)
 F<-(2*(s*N.star+1)*T0.sq)/(s^2*(2*m+s+1))
 df<-c((s*(2*m+s+1)),(2*(s*N.star+1)))
 pval<-pf(q=F, df1=(s*(2*m+s+1)),df2=(2*(s*N.star+1)), lower.tail=FALSE)
 Result<-list(Method="Hotelling-Lawley",Test.Stat=T0.sq, App.Stat=F,df=df, p.value=pval,Descriptive=Descriptive,Homogeneity=TRUE,Test=Name,CI=FALSE)
 }

 # Pillai Method
 if (method=="Pillai") {
 T=0
 for (i in 1:p){
 T<-T+eigens[i]/(eigens[i]+1)
 }
 F<-((2*N.star+s+1)*T)/((2*m+s+1)*(s-T))
 df<-c((s*(2*m+s+1)),(s*(2*N.star+s+1)))
 pval<-pf(q=F, df1=(s*(2*m+s+1)), df2=(s*(2*N.star+s+1)), lower.tail=FALSE)
 Result<-list(Method="Pillai",Test.Stat=T,App.Stat=F,df=df,p.value=pval,
 Descriptive=Descriptive,Homogeneity=TRUE,Test=Name,CI=FALSE)
 }


 # Roy Method
 if (method=="Roy") {
 Teta<-eigens[1]/(eigens[1]+1)
 F<-(eigens[1]*(VE-r+VH))/r
 df<-c(r,(VE-r+VH))
 pval<-pf(q=F, df1=r,df2=(VE-r+VH),lower.tail=FALSE)
 
 if (CI=="FALSE") {   # None Confidence Intervals
 Result<-list(Method="Roy",Test.Stat=Teta,App.Stat=F,df=df,p.value=pval,
 Descriptive=Descriptive,Homogeneity=TRUE,Test=Name,CI=FALSE)
 }
 
 else if (CI=="TRUE") {   # with Confidence Intervals
 
 Ft<-qf(p=alpha,df1=r,df2=(VE-r+VH),lower.tail=FALSE)
 lamda.star<-(VE-r+VH)/(r*Ft) 

 
 C<-choose(g, 2)
 CI.list<-NULL
 for (d in 1:p) {  
 mcp=d   # for dth variable confidence intervals
 
 Karar<-matrix(rep('.',C*C),g,g)

 a<-rep(0,p);a[mcp]<-1;
 for (i in 1:(g-1)) {  # 1. Group selection
 for (k in (i+1):g) {  # 2. Group selection

 c<-rep(0,g)
 c[i]=1;c[k]=-1
 D=0  # The total term with C^2 
 for (z in 1:g) {
 D<-D+c[z]^2/ns[z]}
 mean.diff<-c%*%t(a%*%X.mean.group)
 aWa<-t(a)%*%W%*%a
 Lower<-mean.diff-sqrt(lamda.star*aWa*D)
 Upper<-mean.diff+sqrt(lamda.star*aWa*D)
 if (Lower<0 & Upper>0) {
 Karar[i,k]<-Karar[k,i]<-"TRUE"
 }else{Karar[i,k]<-Karar[k,i]<-"FALSE"
 }}}
 colnames(Karar)<-rownames(Karar)<-Levels
 CI.list[[Variables[d]]]<-Karar
 }
 
 Result<-list(Method="Roy",Test.Stat=Teta,App.Stat=F,df=df,p.value=pval,
 Decision.Matrices=CI.list,Descriptive=Descriptive,Homogeneity=TRUE,Test=Name,CI=TRUE)
 }}
 
 } 

 if (Homogeneity==FALSE) {
 group<-as.factor(group)
 Levels<-levels(group)
 X<-cbind(data,group)
 p<-ncol(data);n<-nrow(data);g=length(levels(group));VE<-n-g; VH<-g-1
 ns<-as.vector(table(group))


 #Calculation W, Wi, Xi.mean, X.mean
 Means<-matrix(0,ncol=p,nrow=g)
 Wi<-array(dim=c(p,p,g))
 sumwiXi<-0
 for (i in 1:g) {  
 Means[i,]<-colMeans(data[which(group==Levels[i]),])
 Wi[,,i]<-ns[i]*solve(cov(data[which(group==Levels[i]),]))
 sumwiXi<-sumwiXi+Wi[,,i]%*%Means[i,]
 }
 W.sum<-apply(Wi,1:2,sum)
 X.mean<-solve(W.sum)%*%sumwiXi
 
 ##Calculation J
 J=0 
 for (i in 1:g) {
 J=J+t(Means[i,]-X.mean)%*%Wi[,,i]%*%(Means[i,]-X.mean)
 } 


 # Calculation A and B
 Ip<-diag(p); r=p*(g-1)
 A=1;B=0
 for (i in 1:g) {
 ki<-(sum(diag(Ip-solve(W.sum)%*%Wi[,,i]))^2)/(ns[i]-1)
 A=A+(1/2*r)*ki
 B<-B+(1/(r*(r+2)))*(ki+ki/2)
 }
 
 Chisq.table<-qchisq(alpha,df=r,lower.tail=FALSE)
 Chisq.adjucted.table<-Chisq.table*(A+B*Chisq.table)
 pval=pchisq(q=J/(A+B*Chisq.table),df=r,lower.tail=FALSE)

##Descriptive Statistics
 Variables<-colnames(data)
 
 Descriptive<-NULL
 for (A in 1:g) { 
 MU<-colMeans(data[which(group==Levels[A]),]) 
 S<-cov(data[which(group==Levels[A]),]) 
 Desc<-rbind(MU,sqrt(diag(S)))
 colnames(Desc)<-Variables
 rownames(Desc)<-c("Means","Sd")
 Descriptive[[Levels[A]]]<-Desc
 }


 
 Result<-list(Method="James (1954)",Test.Stat=J,App.Stat=Chisq.adjucted.table,df=r,
 p.value=pval,Descriptive=Descriptive,Homogeneity=FALSE,Test=Name)
 }




 class(Result)<-c("MVTests","list")
 return(Result)
 }

