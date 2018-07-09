# code for R (http://www.R-project.org)
# Copyright (C) 2018 Natalia Lombardi de Oliveira,
#                    Marcio Alves Diniz,
#                    Carlos Alberto de Bragança Pereira,
#                    Adriano Polpo (polpo@ufscar.br).
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    For the GNU General Public License see <http://www.gnu.org/licenses/>.
#
#    Please see the paper:
#    DOI: 10.1371/journal.pone.0199102 (To cite this code, please cite the paper).

require("plot3D")

#to have arial font in plot, it is may necessary to run the 3 lines below.
#install.packages("extrafont")
library(extrafont)
#font_import()
loadfonts()
loadfonts(device = "postscript")


## this function returns the current amount of spent time by
## the current R process (user time + system time) in minutes.
.gettime <- function() {
  t <- proc.time()
  ## note that t[1]+t[2] are the process time, not the real computer time,
  ## so we still take proper values even if the machine is running other stuff.
  out <- (t[1]+t[2])/60
  names(out) <- NULL
  return(out)
}

###################################
## Density of a Dirichlet Distribution
## x: 
## a: parameters (a >= 2)
ddirichlet <- function(x,a,log=TRUE) {
  if (!is.matrix(x)) {
    x <- matrix(x,1,length(a))
  }
  
  aux <- matrix(NA,length(x[,1]),length(a))
  for (i in 1:length(a)) {
    aux[,i] <- (a[i]-1)*log(x[,i])
  }
  out <- lgamma(sum(a))-sum(lgamma(a))+apply(aux,1,sum)
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}

###################################
## Sampling from a Dirichlet Distribution
## n: sample size
## a: parameters (a >= 2)
rdirichlet <- function(n=1,a=c(1,1)) {
  out <- matrix(NA,n,length(a))
  for (i in 1:length(a)) {
    out[,i] <- rgamma(n,shape=a[i],scale=1)
  }
  out <- out/apply(out,1,sum)
  return(out)
}

###################################
## Probability of a Multinomial Distribution
## x: 
## a: parameters (a >= 2)
dmult <- function(x,a,log=FALSE) {
  if (!is.matrix(a)) {
    a <- matrix(a,1,length(x))
  }
  
  aux <- matrix(NA,length(a[,1]),length(x))
  for (i in 1:length(x)) {
    aux[,i] <- x[i]*log(a[,i])
  }
  out <- lgamma(sum(x)+1)-sum(lgamma(x+1))+apply(aux,1,sum)
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}

###################################
## Probability of a Multinomial Distribution
## x: 
## a: parameters (a >= 2)
dmult2 <- function(x,a,log=FALSE) {
  if (!is.matrix(x)) {
    x <- matrix(x,1,length(x))
  }
  
  aux <- matrix(NA,length(x[,1]),length(a))
  for (i in 1:length(a)) {
    aux[,i] <- x[,i]*log(a[i])
  }
  n <- sum(x[1,])
  out <- lgamma(n+1)-apply(lgamma(x+1),1,sum)+apply(aux,1,sum)
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}


###################################
# Function to build all possible tables 2x2 for a homogeneity test.
# n1:sample size of row 1.
# n2:sample size of row 2.
homog.2x2 <- function(n1=10,n2=10,n.mcmc=1000,verbose=TRUE,verbose.i=5000) {
  ini <- .gettime()

  if (n1 %% 1 != 0) {
    stop("'n1' must be integer number!")
  }
  if (n2 %% 1 != 0) {
    stop("'n2' must be integer number!")
  }

  n.tab <- (n1+1)*(n2+1)
  out <- matrix(NA,n.tab,17)
  h <- 100
  aux.barnard <- matrix(NA,n.tab,h)
  theta.barnard <- seq(0.0001,0.9999,(0.9999-0.0001)/(h-1))
  i <- 1
  n1. <- n1
  n2. <- n2
  n.. <- n1. + n2.
  for (x11 in 0:n1.) {
    x12 <- n1. - x11
    for (x21 in 0:n2.) {
      x22 <- n2. - x21
      n.1 <- x11+x21
      n.2 <- x12+x22

      #P-value (LRT: exact computation)
      log.lambda <- ( ifelse(n.1==0,0,n.1*(log(n.1)-log(n..)))
                     +ifelse(n.2==0,0,n.2*(log(n.2)-log(n..)))
                     -ifelse(x11==0,0,x11*(log(x11)-log(n1.)))
                     -ifelse(x12==0,0,x12*(log(x12)-log(n1.)))
                     -ifelse(x21==0,0,x21*(log(x21)-log(n2.)))
                     -ifelse(x22==0,0,x22*(log(x22)-log(n2.))))
      h <- exp( lgamma(n1.+1)-lgamma(x11+1)-lgamma(x12+1)   #choose(n1.,x11)
               +lgamma(n2.+1)-lgamma(x21+1)-lgamma(x22+1)   #choose(n2.,x21)
               -lgamma(n..+2)+lgamma(n.1+1)+lgamma(n.2+1))  #choose(n..+1,n.1)

      out[i,1:8] <- c(n1.,n2.,
                      x11,x12,
                      x21,x22,
                      exp(log.lambda),h)

      #p-value (LRT: asymptotic computation)
      out[i,10] <- pchisq(-2*log.lambda,df=1,lower.tail=FALSE)

      #p-value (Chi-Square: asymptotic computation)
      stat.chi <- ( exp(log((x11-(n1.*n.1/n..))^2)-log(n1.)-log(n.1)+log(n..))
                   +exp(log((x12-(n1.*n.2/n..))^2)-log(n1.)-log(n.2)+log(n..))
                   +exp(log((x21-(n2.*n.1/n..))^2)-log(n2.)-log(n.1)+log(n..))
                   +exp(log((x22-(n2.*n.2/n..))^2)-log(n2.)-log(n.2)+log(n..)))
      out[i,11:12] <- c(stat.chi,pchisq(stat.chi,df=1,lower.tail=FALSE))

      #e-value (FBST: asymptotic computation)
      out[i,13] <- pchisq(-2*log.lambda,df=2,lower.tail=FALSE)

      #e-value (FBST: approximated computation)
      theta1 <- rbeta(n.mcmc,x11+1,x12+1)
      theta2 <- rbeta(n.mcmc,x21+1,x22+1)

      stat.fbst <- ( lgamma(n1.+2)+lgamma(n2.+2)
                    -lgamma(x11+1)-lgamma(x12+1)
                    -lgamma(x21+1)-lgamma(x22+1)
                    +ifelse(n.1 == 0,0,n.1*(log(n.1)-log(n..)))
                    +ifelse(n.2 == 0,0,n.2*(log(n.2)-log(n..))))
      out[i,14] <- exp(stat.fbst)
      out[i,15] <- 1-sum(( dbeta(theta1,x11+1,x12+1,log=TRUE)
                          +dbeta(theta2,x21+1,x22+1,log=TRUE)) >= stat.fbst)/n.mcmc

      #p-value (Barnard test)
      aux.barnard[i,] <- exp( dbinom(x11,n1.,theta.barnard,log=TRUE)
                             +dbinom(x21,n2.,theta.barnard,log=TRUE))

      #p-value (Fisher)
      out[i,17] <- fisher.test(x=matrix(out[i,3:6],2,2,byrow=T),
                               alternative="two.sided")$p.value

      if ((verbose) && ((i %% verbose.i) == 0)) {
        spent <- .gettime() - ini
        cat("\n",sep="")
        cat(date(),"\n",sep="")
        cat("i: ",i," of ",n.tab,
            "   Time spent: ",sprintf("%.2f",spent),"min",
            "   Expected to finish in: ",sprintf("%.2f",(spent/i)*(n.tab-i)),"min",
            "\n",sep="")
      }
      i <- i+1
    }
  }
  colnames(out) <- c("n1.","n2.","x11","x12","x21","x22",
                     "7 LRT.stat","8 h(x)","9 LRT.P","10 LRT.p",
                     "11 Chi.stat","12 Chi.p",
                     "13 FBST.easymp","14 FBST.sup","15 FBST.e",
                     "16 Barnard.p","17 Fisher.p")
  rownames(out) <- NULL

  if (verbose) {
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Just a few more computations, it should be fast.\n",sep="")
  }

  test.order <- order(out[,7])

  aux.barnard <- aux.barnard[test.order,]
  aux.barnard <- apply(aux.barnard,2,cumsum)
  aux.barnard <- apply(aux.barnard,1,max)

  out <- out[test.order,]
  aux <- cumsum(out[,8])/sum(out[,8])
  out[length(out[,1]),9]  <- aux[length(out[,1])]
  out[length(out[,1]),16] <- aux.barnard[length(out[,1])]
  for (i in (n.tab-1):1) {
    if (out[i,7] == out[(i+1),7]) {
      out[i,9]  <- out[(i+1),9]
      out[i,16] <- out[(i+1),16]
    } else {
      out[i,9]  <- aux[i]
      out[i,16] <- aux.barnard[i]
    }
  }

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Total Time spent: ",sprintf("%.2f",spent),"min\n",sep="")
  }

  return(out)
}

###################################
# Power function of a Hard-Weinberg test.
# p1: number of points in the grid for theta1.
# p2: number of points in the grid for theta2.
# alpha: significance level to evaluate the power function.
# n: sample size.
# homog: output of homog.2x2 function.
# (n or homog must be specified, if both are defined, n will be used)
# n.rep: number of tables to sample.
power.homog.2x2 <- function(p1=100,p2=100,alpha=0.05,n=NULL,homog=NULL,n.rep=1000,
                            n.mcmc=1000,verbose=TRUE,verbose.i=100) {
  ini <- .gettime()

  if (!is.null(n)) {
    if (length(n) != 2) {
      stop("n must be a vector of size 2!")
    }
    if (verbose) {
      spent <- .gettime() - ini
      cat("\n",sep="")
      cat(date(),"\n",sep="")
      cat("Evaluating the p.values of the Homogeneity test...\n",sep="")
    }
    homog <- homog.2x2(n1=n[1],n2=n[2],n.mcmc=n.mcmc,verbose=verbose,verbose.i=verbose.i)
  } else if (is.null(homog)) {
    stop("n or homog must be defined.")
  }

  out <- NULL
  out$n <- homog[1,1:2]
  out$p <- homog

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Evaluating the power function of the Homogeneity test...\n",sep="")
  }

  ini2 <- .gettime()

  P.exa <- stepfun(x=homog[,7],y=c(0,homog[,9]),f=0)
  p.LRT <- stepfun(x=homog[,7],y=c(0,homog[,10]),f=0)

  n.size <- p1*p2
  out$power <- matrix(NA,n.size,7)
  theta1 <- seq(0.0001,0.9999,(0.9999-0.0001)/(p1-1))
  theta2 <- seq(0.0001,0.9999,(0.9999-0.0001)/(p2-1))
  aux.bar1 <- homog[(homog[,16] <= alpha),7]
  aux.bar1 <- aux.bar1[length(aux.bar1)]
  aux.bar <- homog[(homog[,7] <= aux.bar1),c(3,5)]

  i <- 1
  for (i1 in 1:p1) {
    for (i2 in 1:p2) {
      x1 <- rbinom(n=n.rep,size=out$n[1],prob=theta1[i1])
      x2 <- rbinom(n=n.rep,size=out$n[2],prob=theta2[i2])
      x11 <- x1
      x12 <- n[1]-x1
      x21 <- x2
      x22 <- n[2]-x2
      n1. <- n[1]
      n2. <- n[2]
      n.. <- n1. + n2.
      n.1 <- x11+x21
      n.2 <- x12+x22

      lambda <- exp( ifelse(n.1==0,0,n.1*(log(n.1)-log(n..)))
                    +ifelse(n.2==0,0,n.2*(log(n.2)-log(n..)))
                    -ifelse(x11==0,0,x11*(log(x11)-log(n1.)))
                    -ifelse(x12==0,0,x12*(log(x12)-log(n1.)))
                    -ifelse(x21==0,0,x21*(log(x21)-log(n2.)))
                    -ifelse(x22==0,0,x22*(log(x22)-log(n2.))))

      p.bar  <- sum(dbinom(aux.bar[,1],size=out$n[1],prob=theta1[i1])*
                    dbinom(aux.bar[,2],size=out$n[2],prob=theta2[i2]))

      stat.chi <- ( exp(log((x11-(n1.*n.1/n..))^2)-log(n1.)-log(n.1)+log(n..))
                   +exp(log((x12-(n1.*n.2/n..))^2)-log(n1.)-log(n.2)+log(n..))
                   +exp(log((x21-(n2.*n.1/n..))^2)-log(n2.)-log(n.1)+log(n..))
                   +exp(log((x22-(n2.*n.2/n..))^2)-log(n2.)-log(n.2)+log(n..)))
      stat.chi <- stat.chi[!is.nan(stat.chi)]
      p.chi    <- pchisq(stat.chi,df=1,lower.tail=FALSE)

      aux.fisher <- 0
      for (j in 1:n.rep) {
        aux.p <- fisher.test(x=matrix(c(x11[j],x12[j],x21[j],x22[j]),2,2,byrow=T),
                             alternative="two.sided")$p.value
        if (aux.p < alpha) {
          aux.fisher <- aux.fisher+1
        }
      }

      out$power[i,] <- c(theta1[i1],theta2[i2],
                         mean(P.exa(lambda) < alpha),mean(p.LRT(lambda) < alpha),
                         p.bar,sum(p.chi < alpha)/n.rep,
                         aux.fisher/n.rep)

      if ((verbose) && ((i %% verbose.i) == 0)) {
        spent <- .gettime() - ini2
        cat("\n",sep="")
        cat(date(),"\n",sep="")
        cat("i: ",i," of ",n.size,
            "   Time spent: ",sprintf("%.2f",spent),"min",
            "   Expected to finish in: ",sprintf("%.2f",(spent/i)*(n.size-i)),"min",
            "\n",sep="")
      }
      i <- i+1
    }
  }
  colnames(out$power) <- c("theta1","theta2","3 P","4 LRT",
                           "5 Barnard","6 Chi","7 Fisher")
  rownames(out$power) <- NULL

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Total Time spent: ",sprintf("%.2f",spent),"min\n",sep="")
  }
  return(out)
}

##############################
fig.power.homog.2x2 <- function(homog,x="P",y="LRT",tp=TRUE) {
  fn.aux <- function(x) {
    if (x == "P") {
      out <- 3
    } else if (x == "LRT") {
      out <- 4
    } else if (x == "Chi") {
      out <- 6
    } else if (x == "Bar") {
      out <- 5
    } else if (x == "Fish") {
      out <- 7
    }
    return(out)
  }
  x <- fn.aux(x)
  y <- fn.aux(y)
  textcolor <- rgb(0.4,0.4,0.4,1)

  plot(homog$power[,x],homog$power[,y],xlab="",ylab="",
       xlim=c(0,1),ylim=c(0,1),pch=18,cex=0.5)
  lines(c(-1,2),c(-1,2),col=rgb(0.7,0.7,0.7,0.7),lwd=4)
  op <- par(font = 2)
  legend("topleft",paste(round(100*mean(homog$power[,x] < homog$power[,y]),2),
                         "% above\nthe gray line",sep=""),
         bty="n",text.col=textcolor)
  legend("bottomright",paste(round(100*mean(homog$power[,x] > homog$power[,y]),2),
                         "% bellow\nthe gray line",sep=""),
         bty="n",text.col=textcolor)
  if (tp) {
    legend("topright",paste(round(100*mean(homog$power[,x] == homog$power[,y]),2),
                           "% equal to                  \n",
                           "the gray line               ",sep=""),
           bty="n",text.col=textcolor)
  } else {
    legend("bottomleft",paste("               ",
                              round(100*mean(homog$power[,x] == homog$power[,y]),2),
                              "% equal to\n",
                              "               the gray line",sep=""),
           bty="n",text.col=textcolor)
  }
  par(op)
}

###################################
# Function to build all possible tables 2x3 for a homogeneity test.
# n1, n2: sample size of each row.
homog.2x3 <- function(n1=10,n2=10,n.mcmc=1000,verbose=TRUE,verbose.i=1000) {
  ini <- .gettime()

  if (n1 %% 1 != 0) {
    stop("'n1' must be integer number!")
  }
  if (n2 %% 1 != 0) {
    stop("'n2' must be integer number!")
  }

  n.tab <- (n1+2)*(n1+1)*(n2+2)*(n2+1)/4
  out <- matrix(NA,n.tab,17)
  i <- 1
  n1. <- n1
  n2. <- n2
  n.. <- n1. + n2.
  for (x11 in 0:n1.) {
    for (x12 in 0:(n1.-x11)) {
      x13 <- n1.-x11-x12
      for (x21 in 0:n2.) {
        for (x22 in 0:(n2.-x21)) {
          x23 <- n2.-x21-x22
          n.1 <- x11+x21
          n.2 <- x12+x22
          n.3 <- x13+x23

          #P-value (LRT: exact computation)
          log.lambda <- ( ifelse(n.1==0,0,n.1*(log(n.1)-log(n..)))
                         +ifelse(n.2==0,0,n.2*(log(n.2)-log(n..)))
                         +ifelse(n.3==0,0,n.3*(log(n.3)-log(n..)))
                         -ifelse(x11==0,0,x11*(log(x11)-log(n1.)))
                         -ifelse(x12==0,0,x12*(log(x12)-log(n1.)))
                         -ifelse(x13==0,0,x13*(log(x13)-log(n1.)))
                         -ifelse(x21==0,0,x21*(log(x21)-log(n2.)))
                         -ifelse(x22==0,0,x22*(log(x22)-log(n2.)))
                         -ifelse(x23==0,0,x23*(log(x23)-log(n2.))))
          h <- exp( lgamma(n1.+1)+lgamma(n2.+1)
                   +lgamma(n.1+1)+lgamma(n.2+1)+lgamma(n.3+1)
                   -lgamma(x11+1)-lgamma(x12+1)-lgamma(x13+1)
                   -lgamma(x21+1)-lgamma(x22+1)-lgamma(x23+1)
                   -lgamma(n..+3))
    
          out[i,1:10] <- c(n1.,n2.,
                          x11,x12,x13,
                          x21,x22,x23,
                          exp(log.lambda),h)
    
          #p-value (LRT: asymptotic computation)
          out[i,12] <- pchisq(-2*log.lambda,df=2,lower.tail=FALSE)
    
          #p-value (Chi-Square: asymptotic computation)
          stat.chi <- ( exp(log((x11-(n1.*n.1/n..))^2)-log(n1.)-log(n.1)+log(n..))
                       +exp(log((x12-(n1.*n.2/n..))^2)-log(n1.)-log(n.2)+log(n..))
                       +exp(log((x13-(n1.*n.3/n..))^2)-log(n1.)-log(n.3)+log(n..))
                       +exp(log((x21-(n2.*n.1/n..))^2)-log(n2.)-log(n.1)+log(n..))
                       +exp(log((x22-(n2.*n.2/n..))^2)-log(n2.)-log(n.2)+log(n..))
                       +exp(log((x23-(n2.*n.3/n..))^2)-log(n2.)-log(n.3)+log(n..)))
          out[i,13:14] <- c(stat.chi,pchisq(stat.chi,df=2,lower.tail=FALSE))
        
          #e-value (FBST: asymptotic computation)
          out[i,15] <- pchisq(-2*log.lambda,df=4,lower.tail=FALSE)
    
          #e-value (FBST: approximated computation)
          theta1 <- rdirichlet(n.mcmc,c(x11+1,x12+1,x13+1))
          theta2 <- rdirichlet(n.mcmc,c(x21+1,x22+1,x23+1))
    
          stat.fbst <- ( lgamma(n1.+3)+lgamma(n2.+3)
                        -lgamma(x11+1)-lgamma(x12+1)-lgamma(x13+1)
                        -lgamma(x21+1)-lgamma(x22+1)-lgamma(x23+1)
                        +ifelse(n.1 == 0,0,n.1*(log(n.1)-log(n..)))
                        +ifelse(n.2 == 0,0,n.2*(log(n.2)-log(n..)))
                        +ifelse(n.3 == 0,0,n.3*(log(n.3)-log(n..))))
          out[i,16] <- exp(stat.fbst)
          out[i,17] <- 1-sum( ddirichlet(theta1,c(x11+1,x12+1,x13+1))
                             +ddirichlet(theta2,c(x21+1,x22+1,x23+1))
                               >= stat.fbst)/n.mcmc
        
          if ((verbose) && ((i %% verbose.i) == 0)) {
            spent <- .gettime() - ini
            cat("\n",sep="")
            cat(date(),"\n",sep="")
            cat("i: ",i," of ",n.tab,
                "   Time spent: ",sprintf("%.2f",spent),"min",
                "   Expected to finish in: ",sprintf("%.2f",(spent/i)*(n.tab-i)),"min",
                "\n",sep="")
          }
          i <- i+1
        }
      }
    }
  }
  colnames(out) <- c("n1.","n2.","x11","x12","x13","x21","x22","x23",
                     "9 LRT.stat","10 h(x)","11 LRT.P","12 LRT.p",
                     "13 Chi.stat","14 Chi.p",
                     "15 FBST.easymp","16 FBST.sup","17 FBST.e")
  rownames(out) <- NULL

  if (verbose) {
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Just a few more computations, it should be fast.\n",sep="")
  }

  test.order <- order(out[,9])
  print(c(i,n.tab))

  out <- out[test.order,]
  aux <- cumsum(out[,10])/sum(out[,10])
  out[n.tab,11] <- aux[n.tab]
  for (i in (n.tab-1):1) {
    if (out[i,9] == out[(i+1),9]) {
      out[i,11] <- out[(i+1),11]
    } else {
      out[i,11] <- aux[i]
    }
  }

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Total Time spent: ",sprintf("%.2f",spent),"min\n",sep="")
  }

  return(out)
}


###################################
# Function to build all possible tables 3x3 for a homogeneity test.
# n1, n2, n3: sample size of each row.
homog.3x3 <- function(n1=10,n2=10,n3=10,n.mcmc=1000,verbose=TRUE,verbose.i=1000) {
  ini <- .gettime()

  if (n1 %% 1 != 0) {
    stop("'n1' must be integer number!")
  }
  if (n2 %% 1 != 0) {
    stop("'n2' must be integer number!")
  }
  if (n3 %% 1 != 0) {
    stop("'n3' must be integer number!")
  }


  n.tab <- round( exp(lgamma(n1+3)-lgamma(n1+1)-lgamma(3))
                 *exp(lgamma(n2+3)-lgamma(n2+1)-lgamma(3))
                 *exp(lgamma(n3+3)-lgamma(n3+1)-lgamma(3)),0)
  out <- matrix(NA,n.tab,20)
  i <- 1
  n1. <- n1
  n2. <- n2
  n3. <- n3
  n.. <- n1. + n2. + n3.
  for (x11 in 0:n1.) {
    for (x12 in 0:(n1.-x11)) {
      x13 <- n1.-x11-x12
      for (x21 in 0:n2.) {
        for (x22 in 0:(n2.-x21)) {
          x23 <- n2.-x21-x22
          for (x31 in 0:n3.) {
            for (x32 in 0:(n3.-x31)) {
              x33 <- n2.-x31-x32

              n.1 <- x11+x21+x31
              n.2 <- x12+x22+x32
              n.3 <- x13+x23+x33

              #P-value (LRT: exact computation)
              log.lambda <- ( ifelse(n.1==0,0,n.1*(log(n.1)-log(n..)))
                             +ifelse(n.2==0,0,n.2*(log(n.2)-log(n..)))
                             +ifelse(n.3==0,0,n.3*(log(n.3)-log(n..)))
                             -ifelse(x11==0,0,x11*(log(x11)-log(n1.)))
                             -ifelse(x12==0,0,x12*(log(x12)-log(n1.)))
                             -ifelse(x13==0,0,x13*(log(x13)-log(n1.)))
                             -ifelse(x21==0,0,x21*(log(x21)-log(n2.)))
                             -ifelse(x22==0,0,x22*(log(x22)-log(n2.)))
                             -ifelse(x23==0,0,x23*(log(x23)-log(n2.)))
                             -ifelse(x31==0,0,x31*(log(x31)-log(n3.)))
                             -ifelse(x32==0,0,x32*(log(x32)-log(n3.)))
                             -ifelse(x33==0,0,x33*(log(x33)-log(n3.))))
              h <- exp( lgamma(n1.+1)+lgamma(n2.+1)+lgamma(n3.+1)
                       +lgamma(n.1+1)+lgamma(n.2+1)+lgamma(n.3+1)
                       -lgamma(x11+1)-lgamma(x12+1)-lgamma(x13+1)
                       -lgamma(x21+1)-lgamma(x22+1)-lgamma(x23+1)
                       -lgamma(x31+1)-lgamma(x32+1)-lgamma(x33+1)
                       -lgamma(n..+3))
        
              out[i,1:13] <- c(n1.,n2.,
                              x11,x12,x13,
                              x21,x22,x23,
                              x31,x32,x33,
                              exp(log.lambda),h)
        
              #p-value (LRT: asymptotic computation)
              out[i,15] <- pchisq(-2*log.lambda,df=4,lower.tail=FALSE)
        
              #p-value (Chi-Square: asymptotic computation)
              stat.chi <- ( exp(log((x11-(n1.*n.1/n..))^2)-log(n1.)-log(n.1)+log(n..))
                           +exp(log((x12-(n1.*n.2/n..))^2)-log(n1.)-log(n.2)+log(n..))
                           +exp(log((x13-(n1.*n.3/n..))^2)-log(n1.)-log(n.3)+log(n..))
                           +exp(log((x21-(n2.*n.1/n..))^2)-log(n2.)-log(n.1)+log(n..))
                           +exp(log((x22-(n2.*n.2/n..))^2)-log(n2.)-log(n.2)+log(n..))
                           +exp(log((x23-(n2.*n.3/n..))^2)-log(n2.)-log(n.3)+log(n..))
                           +exp(log((x31-(n3.*n.1/n..))^2)-log(n3.)-log(n.1)+log(n..))
                           +exp(log((x32-(n3.*n.2/n..))^2)-log(n3.)-log(n.2)+log(n..))
                           +exp(log((x33-(n3.*n.3/n..))^2)-log(n3.)-log(n.3)+log(n..)))
              out[i,16:17] <- c(stat.chi,pchisq(stat.chi,df=4,lower.tail=FALSE))
            
              #e-value (FBST: asymptotic computation)
              out[i,18] <- pchisq(-2*log.lambda,df=6,lower.tail=FALSE)
        
              #e-value (FBST: approximated computation)
              theta1 <- rdirichlet(n.mcmc,c(x11+1,x12+1,x13+1))
              theta2 <- rdirichlet(n.mcmc,c(x21+1,x22+1,x23+1))
              theta3 <- rdirichlet(n.mcmc,c(x31+1,x32+1,x33+1))
        
              stat.fbst <- ( lgamma(n1.+3)+lgamma(n2.+3)+lgamma(n3.+3)
                            -lgamma(x11+1)-lgamma(x12+1)-lgamma(x13+1)
                            -lgamma(x21+1)-lgamma(x22+1)-lgamma(x23+1)
                            -lgamma(x31+1)-lgamma(x32+1)-lgamma(x33+1)
                            +ifelse(n.1 == 0,0,n.1*(log(n.1)-log(n..)))
                            +ifelse(n.2 == 0,0,n.2*(log(n.2)-log(n..)))
                            +ifelse(n.3 == 0,0,n.3*(log(n.3)-log(n..))))
              out[i,19] <- exp(stat.fbst)
              out[i,20] <- 1-sum( ddirichlet(theta1,c(x11+1,x12+1,x13+1))
                                 +ddirichlet(theta2,c(x21+1,x22+1,x23+1))
                                 +ddirichlet(theta2,c(x31+1,x32+1,x33+1))
                                   >= stat.fbst)/n.mcmc

              if ((verbose) && ((i %% verbose.i) == 0)) {
                spent <- .gettime() - ini
                cat("\n",sep="")
                cat(date(),"\n",sep="")
                cat("i: ",i," of ",n.tab,
                    "   Time spent: ",sprintf("%.2f",spent),"min",
                    "   Expected to finish in: ",sprintf("%.2f",(spent/i)*(n.tab-i)),
                    "min","\n",sep="")
              }
              i <- i+1
            }
          }
        }
      }
    }
  }
  colnames(out) <- c("n1.","n2.","x11","x12","x13","x21","x22","x23","x31","x32","x33",
                     "12 LRT.stat","13 h(x)","14 LRT.P","15 LRT.p",
                     "16 Chi.stat","17 Chi.p",
                     "18 FBST.easymp","19 FBST.sup","20 FBST.e")
  rownames(out) <- NULL

  if (verbose) {
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Just a few more computations, it should be fast.\n",sep="")
  }

  test.order <- order(out[,12])
  print(c(i,n.tab))

  out <- out[test.order,]
  aux <- cumsum(out[,13])/sum(out[,13])
  out[n.tab,14] <- aux[n.tab]
  for (i in (n.tab-1):1) {
    if (out[i,12] == out[(i+1),12]) {
      out[i,14] <- out[(i+1),14]
    } else {
      out[i,14] <- aux[i]
    }
  }

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Total Time spent: ",sprintf("%.2f",spent),"min\n",sep="")
  }

  return(out)
}


###################################
# Function to build all possible tables 2x2 for a independence test.
# n: sample size.
indep.2x2 <- function(n=10,n.mcmc=1000,verbose=TRUE,verbose.i=1000) {
  n.. <- n
  ini <- .gettime()

  if (n %% 1 != 0) {
    stop("'n' must be an integer number!")
  }

  n.tab <- (n..+3)*(n..+2)*(n..+1)/6

  out <- matrix(NA,n.tab,14)
  i <- 1
  for (x11 in 0:n..) {
    for (x12 in 0:(n..-x11)) {
      for (x21 in 0:(n..-x11-x12)) {
        x22 <- (n..-x11-x12-x21)
        n1. <- x11+x12
        n2. <- x21+x22
        n.1 <- x11+x21
        n.2 <- x12+x22

        #P-value (LRT: exact computation)
        log.lambda <- ( ifelse(n.1==0,0,n.1*(log(n.1)-log(n..)))
                       +ifelse(n.2==0,0,n.2*(log(n.2)-log(n..)))
                       +ifelse(n1.==0,0,n1.*(log(n1.)-log(n..)))
                       +ifelse(n2.==0,0,n2.*(log(n2.)-log(n..)))
                       -ifelse(x11==0,0,x11*(log(x11)-log(n..)))
                       -ifelse(x12==0,0,x12*(log(x12)-log(n..)))
                       -ifelse(x21==0,0,x21*(log(x21)-log(n..)))
                       -ifelse(x22==0,0,x22*(log(x22)-log(n..))))
        h <- exp( lgamma(n..+1)+2*lgamma(n..+3)
                 -lgamma(x11+1)-lgamma(x12+1)
                 -lgamma(x21+1)-lgamma(x22+1)
                 -lgamma(n.1+1)-lgamma(n.2+1)
                 -lgamma(n1.+1)-lgamma(n2.+1))
    
        out[i,1:7] <- c(n..,
                        x11,x12,
                        x21,x22,
                        exp(log.lambda),h)
    
        #p-value (LRT: asymptotic computation)
        out[i,9] <- pchisq(-2*log.lambda,df=1,lower.tail=FALSE)
    
        #p-value (Chi-Square: asymptotic computation)
        stat.chi <- ( exp(log((x11-(n1.*n.1/n..))^2)-log(n1.)-log(n.1)+log(n..))
                     +exp(log((x12-(n1.*n.2/n..))^2)-log(n1.)-log(n.2)+log(n..))
                     +exp(log((x21-(n2.*n.1/n..))^2)-log(n2.)-log(n.1)+log(n..))
                     +exp(log((x22-(n2.*n.2/n..))^2)-log(n2.)-log(n.2)+log(n..)))
        out[i,10:11] <- c(stat.chi,pchisq(stat.chi,df=1,lower.tail=FALSE))
      
        #e-value (FBST: asymptotic computation)
        out[i,12] <- pchisq(-2*log.lambda,df=3,lower.tail=FALSE)
  
        #e-value (FBST: approximated computation)
        theta <- rdirichlet(n.mcmc,c(x11+1,x12+1,x21+1,x22+1))
    
        stat.fbst <- (-lgamma(x11+1)-lgamma(x12+1)
                      -lgamma(x21+1)-lgamma(x22+1)
                      +lgamma(n..+4)
                      +ifelse(n.1 == 0,0,n.1*(log(n.1)-log(n..)))
                      +ifelse(n.2 == 0,0,n.2*(log(n.2)-log(n..)))
                      +ifelse(n1. == 0,0,n1.*(log(n1.)-log(n..)))
                      +ifelse(n2. == 0,0,n2.*(log(n2.)-log(n..))))
        out[i,13] <- exp(stat.fbst)
        out[i,14] <- 1-sum(ddirichlet(theta,c(x11+1,x12+1,x21+1,x22+1))
                            >= stat.fbst)/n.mcmc
        
        if ((verbose) && ((i %% verbose.i) == 0)) {
          spent <- .gettime() - ini
          cat("\n",sep="")
          cat(date(),"\n",sep="")
          cat("i: ",i," of ",n.tab,
              "   Time spent: ",sprintf("%.2f",spent),"min",
              "   Expected to finish in: ",sprintf("%.2f",(spent/i)*(n.tab-i)),"min",
              "\n",sep="")
        }
        i <- i+1
      }
    }
  }
  colnames(out) <- c("n..","x11","x12","x21","x22",
                     "6 LRT.stat","7 h(x)","8 LRT.P","9 LRT.p",
                     "10 Chi.stat","11 Chi.p",
                     "12 FBST.easymp","13 FBST.sup","14 FBST.e")
  rownames(out) <- NULL

  if (verbose) {
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Just a few more computations, it should be fast.\n",sep="")
  }

  test.order <- order(out[,6])
  out <- out[test.order,]
  aux <- cumsum(out[,7])/sum(out[,7])
  out[n.tab,8] <- aux[n.tab]
  for (i in (n.tab-1):1) {
    if (out[i,6] == out[(i+1),6]) {
      out[i,8] <- out[(i+1),8]
    } else {
      out[i,8] <- aux[i]
    }
  }

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Total Time spent: ",sprintf("%.2f",spent),"min\n",sep="")
  }

  return(out)
}

###################################
# Function to build all possible tables 2x3 for a independence test.
# n: sample size.
indep.2x3 <- function(n=10,n.mcmc=1000,verbose=TRUE,verbose.i=1000) {
  n.. <- n
  ini <- .gettime()

  if (n %% 1 != 0) {
    stop("'n' must be an integer number!")
  }

  n.tab <- round(exp(lgamma(n..+2*3)-lgamma(n..+1)-lgamma(2*3)),0)

  out <- matrix(NA,n.tab,16)
  i <- 1
  for (x11 in 0:n..) {
    for (x12 in 0:(n..-x11)) {
      for (x13 in 0:(n..-x11-x12)) {
        for (x21 in 0:(n..-x11-x12-x13)) {
          for (x22 in 0:(n..-x11-x12-x13-x21)) {
            x23 <- (n..-x11-x12-x13-x21-x22)
            n1. <- x11+x12+x13
            n2. <- x21+x22+x23
            n.1 <- x11+x21
            n.2 <- x12+x22
            n.3 <- x13+x23

            #P-value (LRT: exact computation)
            log.lambda <- ( ifelse(n.1==0,0,n.1*(log(n.1)-log(n..)))
                           +ifelse(n.2==0,0,n.2*(log(n.2)-log(n..)))
                           +ifelse(n.3==0,0,n.3*(log(n.3)-log(n..)))
                           +ifelse(n1.==0,0,n1.*(log(n1.)-log(n..)))
                           +ifelse(n2.==0,0,n2.*(log(n2.)-log(n..)))
                           -ifelse(x11==0,0,x11*(log(x11)-log(n..)))
                           -ifelse(x12==0,0,x12*(log(x12)-log(n..)))
                           -ifelse(x13==0,0,x13*(log(x13)-log(n..)))
                           -ifelse(x21==0,0,x21*(log(x21)-log(n..)))
                           -ifelse(x22==0,0,x22*(log(x22)-log(n..)))
                           -ifelse(x23==0,0,x23*(log(x23)-log(n..))))
            h <- exp( lgamma(n..+1)+lgamma(n..+3)+lgamma(n..+4)
                     -lgamma(x11+1)-lgamma(x12+1)-lgamma(x13+1)
                     -lgamma(x21+1)-lgamma(x22+1)-lgamma(x23+1)
                     -lgamma(n.1+1)-lgamma(n.2+1)-lgamma(n.3+1)
                     -lgamma(n1.+1)-lgamma(n2.+1))
    
            out[i,1:9] <- c(n..,
                            x11,x12,x13,
                            x21,x22,x23,
                            exp(log.lambda),h)
    
            #p-value (LRT: asymptotic computation)
            out[i,11] <- pchisq(-2*log.lambda,df=2,lower.tail=FALSE)
    
            #p-value (Chi-Square: asymptotic computation)
            stat.chi <- ( exp(log((x11-(n1.*n.1/n..))^2)-log(n1.)-log(n.1)+log(n..))
                         +exp(log((x12-(n1.*n.2/n..))^2)-log(n1.)-log(n.2)+log(n..))
                         +exp(log((x13-(n1.*n.3/n..))^2)-log(n1.)-log(n.3)+log(n..))
                         +exp(log((x21-(n2.*n.1/n..))^2)-log(n2.)-log(n.1)+log(n..))
                         +exp(log((x22-(n2.*n.2/n..))^2)-log(n2.)-log(n.2)+log(n..))
                         +exp(log((x23-(n2.*n.3/n..))^2)-log(n2.)-log(n.3)+log(n..)))
            out[i,12:13] <- c(stat.chi,pchisq(stat.chi,df=2,lower.tail=FALSE))
      
            #e-value (FBST: asymptotic computation)
            out[i,14] <- pchisq(-2*log.lambda,df=5,lower.tail=FALSE)
  
            #e-value (FBST: approximated computation)
            theta <- rdirichlet(n.mcmc,c(x11+1,x12+1,x13+1,x21+1,x22+1,x23+1))
    
            stat.fbst <- (-lgamma(x11+1)-lgamma(x12+1)-lgamma(x13+1)
                          -lgamma(x21+1)-lgamma(x22+1)-lgamma(x23+1)
                          +lgamma(n..+5)
                          +ifelse(n.1 == 0,0,n.1*(log(n.1)-log(n..)))
                          +ifelse(n.2 == 0,0,n.2*(log(n.2)-log(n..)))
                          +ifelse(n.3 == 0,0,n.3*(log(n.3)-log(n..)))
                          +ifelse(n1. == 0,0,n1.*(log(n1.)-log(n..)))
                          +ifelse(n2. == 0,0,n2.*(log(n2.)-log(n..))))
            out[i,15] <- exp(stat.fbst)
            out[i,16] <- 1-sum(ddirichlet(theta,c(x11+1,x12+1,x13+1,x21+1,x22+1,x23+1))
                                >= stat.fbst)/n.mcmc
        
            if ((verbose) && ((i %% verbose.i) == 0)) {
              spent <- .gettime() - ini
              cat("\n",sep="")
              cat(date(),"\n",sep="")
              cat("i: ",i," of ",n.tab,
                  "   Time spent: ",sprintf("%.2f",spent),"min",
                  "   Expected to finish in: ",sprintf("%.2f",(spent/i)*(n.tab-i)),"min",
                  "\n",sep="")
            }
            i <- i+1
          }
        }
      }
    }
  }
  colnames(out) <- c("n..","x11","x12","x13","x21","x22","x23",
                     "8 LRT.stat","9 h(x)","10 LRT.P","11 LRT.p",
                     "12 Chi.stat","13 Chi.p",
                     "14 FBST.easymp","15 FBST.sup","16 FBST.e")
  rownames(out) <- NULL

  if (verbose) {
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Just a few more computations, it should be fast.\n",sep="")
  }

  test.order <- order(out[,8])
  out <- out[test.order,]
  aux <- cumsum(out[,9])/sum(out[,9])
  out[n.tab,10] <- aux[n.tab]
  for (i in (n.tab-1):1) {
    if (out[i,8] == out[(i+1),8]) {
      out[i,10] <- out[(i+1),10]
    } else {
      out[i,10] <- aux[i]
    }
  }

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Total Time spent: ",sprintf("%.2f",spent),"min\n",sep="")
  }

  return(out)
}


###################################
# Function to build all possible tables 3x3 for a independence test.
# n: sample size.
indep.3x3 <- function(n=10,n.mcmc=1000,verbose=TRUE,verbose.i=1000) {
  n.. <- n
  ini <- .gettime()

  if (n %% 1 != 0) {
    stop("'n' must be an integer number!")
  }

  n.tab <- round(exp(lgamma(n..+3*3)-lgamma(n..+1)-lgamma(3*3)),0)
  out <- matrix(NA,nrow=n.tab,ncol=19)

  i <- 1
  for (x11 in 0:n..) {
    for (x12 in 0:(n..-x11)) {
      for (x13 in 0:(n..-x11-x12)) {
        for (x21 in 0:(n..-x11-x12-x13)) {
          for (x22 in 0:(n..-x11-x12-x13-x21)) {
            for (x23 in 0:(n..-x11-x12-x13-x21-x22)) {
              for (x31 in 0:(n..-x11-x12-x13-x21-x22-x23)) {
                for (x32 in 0:(n..-x11-x12-x13-x21-x22-x23-x31)) {
                  x33 <- (n..-x11-x12-x13-x21-x22-x23-x31-x32)
                  n1. <- x11+x12+x13
                  n2. <- x21+x22+x23
                  n3. <- x31+x32+x33
                  n.1 <- x11+x21+x31
                  n.2 <- x12+x22+x32
                  n.3 <- x13+x23+x33

                  #P-value (LRT: exact computation)
                  log.lambda <- ( ifelse(n.1==0,0,n.1*(log(n.1)-log(n..)))
                                 +ifelse(n.2==0,0,n.2*(log(n.2)-log(n..)))
                                 +ifelse(n.3==0,0,n.3*(log(n.3)-log(n..)))
                                 +ifelse(n1.==0,0,n1.*(log(n1.)-log(n..)))
                                 +ifelse(n2.==0,0,n2.*(log(n2.)-log(n..)))
                                 +ifelse(n3.==0,0,n3.*(log(n3.)-log(n..)))
                                 -ifelse(x11==0,0,x11*(log(x11)-log(n..)))
                                 -ifelse(x12==0,0,x12*(log(x12)-log(n..)))
                                 -ifelse(x13==0,0,x13*(log(x13)-log(n..)))
                                 -ifelse(x21==0,0,x21*(log(x21)-log(n..)))
                                 -ifelse(x22==0,0,x22*(log(x22)-log(n..)))
                                 -ifelse(x23==0,0,x23*(log(x23)-log(n..)))
                                 -ifelse(x31==0,0,x31*(log(x31)-log(n..)))
                                 -ifelse(x32==0,0,x32*(log(x32)-log(n..)))
                                 -ifelse(x33==0,0,x33*(log(x33)-log(n..))))
                  h <- exp( lgamma(n..+1)+2*lgamma(n..+4)
                           -lgamma(x11+1)-lgamma(x12+1)-lgamma(x13+1)
                           -lgamma(x21+1)-lgamma(x22+1)-lgamma(x23+1)
                           -lgamma(x31+1)-lgamma(x32+1)-lgamma(x33+1)
                           -lgamma(n.1+1)-lgamma(n.2+1)-lgamma(n.3+1)
                           -lgamma(n1.+1)-lgamma(n2.+1)-lgamma(n3.+1))

                  out[i,1:12] <- c(n..,
                                   x11,x12,x13,
                                   x21,x22,x23,
                                   x31,x32,x33,
                                   exp(log.lambda),h)
    
                  #p-value (LRT: asymptotic computation)
                  out[i,14] <- pchisq(-2*log.lambda,df=4,lower.tail=FALSE)
    
                  #p-value (Chi-Square: asymptotic computation)
                  stat.chi <- ( exp(log((x11-(n1.*n.1/n..))^2)-log(n1.)-log(n.1)+log(n..))
                               +exp(log((x12-(n1.*n.2/n..))^2)-log(n1.)-log(n.2)+log(n..))
                               +exp(log((x13-(n1.*n.3/n..))^2)-log(n1.)-log(n.3)+log(n..))
                               +exp(log((x21-(n2.*n.1/n..))^2)-log(n2.)-log(n.1)+log(n..))
                               +exp(log((x22-(n2.*n.2/n..))^2)-log(n2.)-log(n.2)+log(n..))
                               +exp(log((x23-(n2.*n.3/n..))^2)-log(n2.)-log(n.3)+log(n..))
                               +exp(log((x31-(n3.*n.1/n..))^2)-log(n3.)-log(n.1)+log(n..))
                               +exp(log((x32-(n3.*n.2/n..))^2)-log(n3.)-log(n.2)+log(n..))
                               +exp(log((x33-(n3.*n.3/n..))^2)-log(n3.)-log(n.3)+log(n..)))
                  out[i,15:16] <- c(stat.chi,pchisq(stat.chi,df=4,lower.tail=FALSE))
      
                  #e-value (FBST: asymptotic computation)
                  out[i,17] <- pchisq(-2*log.lambda,df=8,lower.tail=FALSE)
  
                  #e-value (FBST: approximated computation)
                  theta <- rdirichlet(n.mcmc,c(x11+1,x12+1,x13+1,
                                               x21+1,x22+1,x23+1,
                                               x31+1,x32+1,x33+1))
    
                  stat.fbst <- (-lgamma(x11+1)-lgamma(x12+1)-lgamma(x13+1)
                                -lgamma(x21+1)-lgamma(x22+1)-lgamma(x23+1)
                                -lgamma(x31+1)-lgamma(x32+1)-lgamma(x33+1)
                                +lgamma(n..+8)
                                +ifelse(n.1 == 0,0,n.1*(log(n.1)-log(n..)))
                                +ifelse(n.2 == 0,0,n.2*(log(n.2)-log(n..)))
                                +ifelse(n.3 == 0,0,n.3*(log(n.3)-log(n..)))
                                +ifelse(n1. == 0,0,n1.*(log(n1.)-log(n..)))
                                +ifelse(n2. == 0,0,n2.*(log(n2.)-log(n..)))
                                +ifelse(n3. == 0,0,n3.*(log(n3.)-log(n..))))
                  out[i,18] <- exp(stat.fbst)
                  out[i,19] <- 1-sum(ddirichlet(theta,c(x11+1,x12+1,x13+1,
                                                        x21+1,x22+1,x23+1,
                                                        x31+1,x32+1,x33+1))
                                      >= stat.fbst)/n.mcmc
        
                  if ((verbose) && ((i %% verbose.i) == 0)) {
                    spent <- .gettime() - ini
                    cat("\n",sep="")
                    cat(date(),"\n",sep="")
                    cat("i: ",i," of ",n.tab,
                        "   Time spent: ",sprintf("%.2f",spent),"min",
                        "   Expected to finish in: ",sprintf("%.2f",(spent/i)*(n.tab-i)),
                        "min","\n",sep="")
                  }
                  i <- i+1
                }
              }
            }
          }
        }
      }
    }
  }
  colnames(out) <- c("n..","x11","x12","x13","x21","x22","x23","x31","x32","x33",
                     "11 LRT.stat","12 h(x)","13 LRT.P","14 LRT.p",
                     "15 Chi.stat","16 Chi.p",
                     "17 FBST.easymp","18 FBST.sup","19 FBST.e")
  rownames(out) <- NULL

  if (verbose) {
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Just a few more computations, it should be fast.\n",sep="")
  }

  test.order <- order(out[,11])
  out <- out[test.order,]
  aux <- cumsum(out[,12])/sum(out[,11])
  out[n.tab,13] <- aux[n.tab]
  for (i in (n.tab-1):1) {
    if (out[i,11] == out[(i+1),11]) {
      out[i,13] <- out[(i+1),13]
    } else {
      out[i,13] <- aux[i]
    }
  }

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Total Time spent: ",sprintf("%.2f",spent),"min\n",sep="")
  }

  return(out)
}


###################################
# Function to build all possible tables for a Hard-Weinberg test.
# n: sample size.
HW <- function(n=10,n.mcmc=1000,verbose=TRUE,verbose.i=1000) {
  ini <- .gettime()

  if (n %% 1 != 0) {
    stop("'n' must be an integer number!")
  }

  n.tab <- (n+2)*(n+1)/2
  out <- matrix(NA,n.tab,14)
  i <- 1
  h <- 100
  aux.barnard <- matrix(NA,n.tab,h)
  theta.barnard <- seq(0.0001,0.9999,(0.9999-0.0001)/(h-1))
  for (x1 in 0:n) {
    for (x2 in 0:(n-x1)) {
      x3 <- (n-x1-x2)

      log.lambda <- ( x2*log(2)
                     +ifelse((2*x1+x2)==0,0,(2*x1+x2)*(log(2*x1+x2)-log(2)-log(n)))
                     +ifelse((2*x3+x2)==0,0,(2*x3+x2)*(log(2*x3+x2)-log(2)-log(n)))
                     -ifelse(x1==0,0,x1*(log(x1)-log(n)))
                     -ifelse(x2==0,0,x2*(log(x2)-log(n)))
                     -ifelse(x3==0,0,x3*(log(x3)-log(n))))
      h <- exp( lgamma(n+1)+x2*log(2)
               +lgamma(2*x1+x2+1)+lgamma(2*x3+x2+1)
               -lgamma(x1+1)-lgamma(x2+1)-lgamma(x3+1)
               -lgamma(2*n+2))

      out[i,1:6] <- c(n,x1,x2,x3,exp(log.lambda),h)
    
      #p-value (LRT: asymptotic computation)
      out[i,8] <- pchisq(-2*log.lambda,df=1,lower.tail=FALSE)
    
      #p-value (Chi-Square: asymptotic computation)
      est.theta <- (2*x1+x2)/(2*n)
      est1 <- n*(est.theta^2)
      est2 <- n*2*est.theta*(1-est.theta)
      est3 <- n*((1-est.theta)^2)
      stat.chi <- ( ((x1-est1)^2)/est1
                   +((x2-est2)^2)/est2
                   +((x3-est3)^2)/est3)
      out[i,9:10] <- c(stat.chi,pchisq(stat.chi,df=1,lower.tail=FALSE))
      
      #e-value (FBST: asymptotic computation)
      out[i,11] <- pchisq(-2*log.lambda,df=2,lower.tail=FALSE)
  
      #e-value (FBST: approximated computation)
      theta <- rdirichlet(n.mcmc,c(x1+1,x2+1,x3+1))
    
      stat.fbst <- (lgamma(n+3)-lgamma(x1+1)-lgamma(x2+1)-lgamma(x3+1)
                    +x2*log(2)
                    +ifelse((2*x1+x2) == 0,0,(2*x1+x2)*(log(2*x1+x2)-log(2)-log(n)))
                    +ifelse((2*x3+x2) == 0,0,(2*x3+x2)*(log(2*x3+x2)-log(2)-log(n))))
 
      out[i,12] <- stat.fbst
      out[i,13] <- 1-sum(ddirichlet(theta,c(x1+1,x2+1,x3+1)) >= stat.fbst)/n.mcmc

      #p-value (Barnard test)
      aux.barnard[i,] <- dmult(x=c(x1,x2,x3),a=cbind(theta.barnard^2,
                                                     2*theta.barnard*(1-theta.barnard),
                                                     (1-theta.barnard)^2))
 

      if ((verbose) && ((i %% verbose.i) == 0)) {
        spent <- .gettime() - ini
        cat("\n",sep="")
        cat(date(),"\n",sep="")
        cat("i: ",i," of ",n.tab,
            "   Time spent: ",sprintf("%.2f",spent),"min",
            "   Expected to finish in: ",sprintf("%.2f",(spent/i)*(n.tab-i)),"min",
            "\n",sep="")
      }
      i <- i+1
    }
  }
  colnames(out) <- c("n","x1","x2","x3",
                     "5 LRT.stat","6 h(x)","7 LRT.P","8 LRT.p",
                     "9 Chi.stat","10 Chi.p",
                     "11 FBST.easymp","12 FBST.sup","13 FBST.e","14 Barnard.p")
  rownames(out) <- NULL

  if (verbose) {
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Just a few more computations, it should be fast.\n",sep="")
  }

  test.order <- order(out[,5])
  aux.barnard <- aux.barnard[test.order,]
  aux.barnard <- apply(aux.barnard,2,cumsum)
  aux.barnard <- apply(aux.barnard,1,max)

  out <- out[test.order,]
  aux <- cumsum(out[,6])/sum(out[,6])
  aux <- cumsum(out[,6])/sum(out[,6])
  out[n.tab,7] <- aux[n.tab]
  out[n.tab,14] <- aux.barnard[n.tab]
  for (i in (n.tab-1):1) {
    if (out[i,5] == out[(i+1),5]) {
      out[i,7]  <- out[(i+1),7]
      out[i,14] <- out[(i+1),14]
    } else {
      out[i,7]  <- aux[i]
      out[i,14] <- aux.barnard[i]
    }
  }

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Total Time spent: ",sprintf("%.2f",spent),"min\n",sep="")
  }
  return(out)
}


###################################
# Power function of a Hard-Weinberg test.
# p1: number of points in the grid for theta1.
# p2: number of points in the grid for theta2.
# alpha: significance level to evaluate the power function.
# n: sample size.
# hw: output of HW function.
# (n or hw must be specified, if both are defined, n will be used)
# n.rep: number of tables to sample.
power.hw <- function(p1=100,p2=100,alpha=0.05,n=NULL,hw=NULL,n.rep=1000,
                     verbose=TRUE,verbose.i=1000) {
  ini <- .gettime()

  if (!is.null(n)) {
    if (verbose) {
      spent <- .gettime() - ini
      cat("\n",sep="")
      cat(date(),"\n",sep="")
      cat("Evaluating the p.values of the HW test...\n",sep="")
    }
    hw <- HW(n=n,verbose=verbose,verbose.i=verbose.i)
  } else if (is.null(hw)) {
    stop("n or hw must be defined.")
  }

  out <- NULL
  out$n <- hw[1,1]
  out$p <- hw

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Evaluating the power function of the HW test...\n",sep="")
  }

  ini2 <- .gettime()

  P.exa <- stepfun(x=hw[,5],y=c(0,hw[,7]),f=0)
  p.LRT <- stepfun(x=hw[,5],y=c(0,hw[,8]),f=0)

  theta1 <- seq(0.0001,0.9999,(0.9999-0.0001)/(p1-1))
  theta2 <- seq(0.0001,0.9999,(0.9999-0.0001)/(p2-1))
  aux.bar1 <- hw[(hw[,14] <= alpha),5]
  aux.bar1 <- aux.bar1[length(aux.bar1)]
  aux.bar <- hw[(hw[,5] <= aux.bar1),c(2,3,4)]

  aux <- matrix(NA,p1*p2,3)
  i <- 1
  for (i1 in 1:p1) {
    for (i2 in 1:p2) {
      if (theta1[i1]+theta2[i2] <= 1) {
        aux[i,] <- c(theta1[i1],theta2[i2],1-(theta1[i1]+theta2[i2]))
        i <- i+1
      }
    }
  }
  aux <- aux[!is.na(aux[,1]),]
  n.size <- length(aux[,1])

  out$power <- matrix(NA,n.size,6)
  for (i in 1:n.size) {
    x <- rmultinom(n=n.rep,size=out$n,prob=aux[i,])

    lambda <- exp( x[2,]*log(2)
                  +ifelse((2*x[1,]+x[2,])==0,0,
                          (2*x[1,]+x[2,])*(log(2*x[1,]+x[2,])-log(2)-log(out$n)))
                  +ifelse((2*x[3,]+x[2,])==0,0,
                          (2*x[3,]+x[2,])*(log(2*x[3,]+x[2,])-log(2)-log(out$n)))
                  -ifelse(x[1,]==0,0,x[1,]*(log(x[1,])-log(out$n)))
                  -ifelse(x[2,]==0,0,x[2,]*(log(x[2,])-log(out$n)))
                  -ifelse(x[3,]==0,0,x[3,]*(log(x[3,])-log(out$n))))

    est.theta <- (2*x[1,]+x[2,])/(2*out$n)
    est1 <- out$n*(est.theta^2)
    est2 <- out$n*2*est.theta*(1-est.theta)
    est3 <- out$n*((1-est.theta)^2)
    stat.chi <- ( ((x[1,]-est1)^2)/est1
                 +((x[2,]-est2)^2)/est2
                 +((x[3,]-est3)^2)/est3)
    stat.chi <- stat.chi[!is.nan(stat.chi)]
    p.chi    <- pchisq(stat.chi,df=1,lower.tail=FALSE)

    if (aux[i,3] == 0) {
      aux2 <- c(aux[i,1],aux[i,2]-10^(-10),10^(-10))
    } else {
      aux2 <- aux[i,]
    }
    p.bar <- sum(dmult2(x=aux.bar,a=aux2))

    out$power[i,] <- c(aux[i,1],aux[i,2],
                       mean(P.exa(lambda) < alpha),mean(p.LRT(lambda) < alpha),
                       p.bar,sum(p.chi < alpha)/n.rep)

    if ((verbose) && ((i %% verbose.i) == 0)) {
      spent <- .gettime() - ini2
      cat("\n",sep="")
      cat(date(),"\n",sep="")
      cat("i: ",i," of ",n.size,
          "   Time spent: ",sprintf("%.2f",spent),"min",
          "   Expected to finish in: ",sprintf("%.2f",(spent/i)*(n.size-i)),"min",
          "\n",sep="")
    }
  }
  colnames(out$power) <- c("theta1","theta2","3 P","4 LRT",
                           "5 Barnard","6 Chi")
  rownames(out$power) <- NULL

  if (verbose) {
    spent <- .gettime() - ini
    cat("\n",sep="")
    cat(date(),"\n",sep="")
    cat("Total Time spent: ",sprintf("%.2f",spent),"min\n",sep="")
  }
  return(out)
}

#############################
fig.power.hw <- function(hw,x="P",y="LRT",tp=TRUE) {
  fn.aux <- function(x) {
    if (x == "P") {
      out <- 3
    } else if (x == "LRT") {
      out <- 4
    } else if (x == "Chi") {
      out <- 6
    } else if (x == "Bar") {
      out <- 5
    }
    return(out)
  }
  x <- fn.aux(x)
  y <- fn.aux(y)
  textcolor <- rgb(0.4,0.4,0.4,1)
  
  plot(hw$power[,x],hw$power[,y],xlab="",ylab="",
       xlim=c(0,1),ylim=c(0,1),pch=18,cex=0.5)
  lines(c(-1,2),c(-1,2),col=rgb(0.7,0.7,0.7,0.7),lwd=4)
  op <- par(font = 2)
  legend("topleft",paste(round(100*mean(hw$power[,x] < hw$power[,y]),2),
                         "% above\nthe gray line",sep=""),
         bty="n",text.col=textcolor)
  legend("bottomright",paste(round(100*mean(hw$power[,x] > hw$power[,y]),2),
                         "% bellow\nthe gray line",sep=""),
         bty="n",text.col=textcolor)
  if (tp) {
    legend("topright",paste(round(100*mean(hw$power[,x] == hw$power[,y]),2),
                           "% equal to                  \n",
                           "the gray line               ",sep=""),
           bty="n",text.col=textcolor)
  } else {
    legend("bottomleft",paste("               ",
                              round(100*mean(hw$power[,x] == hw$power[,y]),2),
                              "% equal to\n",
                              "               the gray line",sep=""),
           bty="n",text.col=textcolor)
  }
  par(op)
}

##############################
##############################
##############################

pairs2 <- 
  function (x, labels, panel = points, ..., lower.panel = panel, 
            upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
            label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
            row1attop = TRUE, gap = 1) 
  {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex,
                                                                 font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                          oma, ...) {
      if (side%%2 == 1) 
        Axis(x, side = side, xpd = NA, at=c(0,0.2,0.4,0.6,0.8,1),
             labels=c(0,NA,NA,NA,NA,1), cex.axis=0.8, ...)
      else Axis(y, side = side, xpd = NA, at=c(0,0.2,0.4,0.6,0.8,1),
                labels=c(0,NA,NA,NA,NA,1), cex.axis=0.8, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
      x <- as.data.frame(x)
      for (i in seq_along(names(x))) {
        if (is.factor(x[[i]]) || is.logical(x[[i]])) 
          x[[i]] <- as.numeric(x[[i]])
        if (!is.numeric(unclass(x[[i]]))) 
          stop("non-numeric argument to 'pairs'")
      }
    }
    else if (!is.numeric(x)) 
      stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
      lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
      upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
      diag.panel <- match.fun(diag.panel)
    if (row1attop) {
      tmp <- lower.panel
      lower.panel <- upper.panel
      upper.panel <- tmp
      tmp <- has.lower
      has.lower <- has.upper
      has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
      stop("only one column in the argument to 'pairs'")
    has.labs <- TRUE
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      has.labs <- FALSE
    oma <- if ("oma" %in% nmdots) 
      dots$oma
    else NULL
    main <- if ("main" %in% nmdots) 
      dots$main
    else NULL
    if (is.null(oma)) {
      oma <- c(4, 4, 4, 4)
      if (!is.null(main)) 
        oma[3L] <- 6
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    for (i in if (row1attop) 
      1L:nc
         else nc:1L) for (j in 1L:nc) {
           localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                     type = "n", xlim=c(0,1), ylim=c(0,1), ...)
           if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
             box()
             # edited here...
             #           if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
             #           localAxis(1 + 2 * row1attop, x[, j], x[, i], 
             #                       ...)
             # draw x-axis
             if (i == nc & j != nc) 
               localAxis(1, x[, j], x[, i], 
                         ...)
             # draw y-axis
             if (j == 1 & i != 1) 
               localAxis(2, x[, j], x[, i], ...)
             #           if (j == nc && (i%%2 || !has.upper || !has.lower)) 
             #             localAxis(4, x[, j], x[, i], ...)
             mfg <- par("mfg")
             if (i == j) {
               if (has.diag) 
                 localDiagPanel(as.vector(x[, i]), ...)
               if (has.labs) {
                 par(usr = c(0, 1, 0, 1))
                 if (is.null(cex.labels)) {
                   l.wid <- strwidth(labels, "user")
                   cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                 }
                 text.panel(0.5, label.pos, labels[2,i], cex = cex.labels, 
                            font = font.labels)
                 text.panel(0.5, label.pos+0.17, labels[1,i], cex = cex.labels, 
                            font = font.labels)
                 text.panel(0.5, label.pos-0.15, labels[3,i], cex = cex.labels, 
                            font = font.labels)
               }
             }
             else if (i < j) {
               localLowerPanel(as.vector(x[, j]), as.vector(x[,i]), ...)
             }
             else {
               localUpperPanel(as.vector(x[, j]), as.vector(x[,i]), ...)
             #  lines(c(-1,2),c(-1,2),type="l",col=rgb(0.5,0.5,0.5,alpha=0.7),lwd=2)
              lines(c(-1,2),c(-1,2),type="l",col=rgb(0.7,0.7,0.7,0.7),lwd=1)
             }
             if (any(par("mfg") != mfg)) 
               stop("the 'panel' function made a new plot")
           }
           else par(new = FALSE)
         }
    if (!is.null(main)) {
      font.main <- if ("font.main" %in% nmdots) 
        dots$font.main
      else par("font.main")
      cex.main <- if ("cex.main" %in% nmdots) 
        dots$cex.main
      else par("cex.main")
      mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
  }

##############################
##############################
##############################
plot.power3d <- function(power,x="P",y="Bar",name=NULL,flag=TRUE,a1=60,a2=10) {
  fn.aux <- function(x) {
    if (x == "P") {
      out <- list(i=3,col=4)
    } else if (x == "LRT") {
      out <- list(i=4,col=3)
    } else if (x == "Chi") {
      out <- list(i=6,col=2)
    } else if (x == "Bar") {
      out <- list(i=5,col=1)
    } else if (x == "Fish") {
      out <- list(i=7,col="magenta")
    }
    return(out)
  }
  x <- fn.aux(x)
  y <- fn.aux(y)

  smooth.x <- loess(power[,x$i] ~ power[,1] + power[,2],span=0.01)
  smooth.y <- loess(power[,y$i] ~ power[,1] + power[,2],span=0.01)
  theta1 <- sort(unique(power[,1]))
  theta2 <- sort(unique(power[,2]))

  z3d.x <- matrix(NA,length(theta1),length(theta2))
  z3d.y <- matrix(NA,length(theta1),length(theta2))
  k <- 1
  for (i in 1:length(theta1)) {
    for (j in 1:length(theta2)) {
      if (flag) {
        if ((power[k,1] == theta1[i]) && (power[k,2] == theta2[j])) {
          if (x$i != 5) {
            z3d.x[i,j] <- smooth.x$fitted[k]
          }
          if (y$i != 5) {
            z3d.y[i,j] <- smooth.y$fitted[k]
          }
          k <- ifelse(k < length(power[,1]),k+1,k)
        }
      }
      if ((!flag) || (x$i == 5)) {
        aux1 <- power[(power[,1] == theta1[i]),c(2,x$i)]
        if (!is.matrix(aux1)) {
          aux1 <- matrix(aux1,1,2)
        }
        aux1 <- aux1[(aux1[,1] == theta2[j]),2]
        if (length(aux1) != 0) {
          z3d.x[i,j] <- aux1
        }
      }
      if ((!flag) || (y$i == 5)) {
        aux2 <- power[(power[,1] == theta1[i]),c(2,y$i)]
        if (!is.matrix(aux2)) {
          aux2 <- matrix(aux2,1,2)
        }
        aux2 <- aux2[(aux2[,1] == theta2[j]),2]
        if (length(aux2) != 0) {
          z3d.y[i,j] <- aux2
        }
      }

      if ((!is.na(z3d.x[i,j])) && (z3d.x[i,j] < 0)) {
        z3d.x[i,j] <- 0
      }
      if ((!is.na(z3d.x[i,j])) && (z3d.x[i,j] > 1)) {
        z3d.x[i,j] <- 1
      }
      if ((!is.na(z3d.y[i,j])) && (z3d.y[i,j] < 0)) {
        z3d.y[i,j] <- 0
      }
      if ((!is.na(z3d.y[i,j])) && (z3d.y[i,j] > 1)) {
        z3d.y[i,j] <- 1
      }
    }
  }

#  print(mean(z3d.x > z3d.y, na.rm=TRUE))
#  print(mean(z3d.x == z3d.y, na.rm=TRUE))
#  print(mean(z3d.x < z3d.y, na.rm=TRUE))

  if (!is.null(name)) {
#    bitmap(paste(name,".tiff",sep=""),height=10,width=10,units ='cm',type="tiff24nc",
#           res=600,method="pdf")
#    png(paste(name,".png",sep=""),width=4000,height=4000,res=570)
  }
  persp3D(z=z3d.y,facets=F,theta=a1,phi=a2,col=y$col,ticktype="detailed",family="Arial",
          cex.axis=0.7,shade=0.2,add=F,xlab="θ1",ylab="θ2",zlab="power(θ1,θ2)")
  persp3D(z=z3d.x,facets=F,col=x$col,shade=0.2,add=T)
  if (!is.null(name)) {
    dev2bitmap(paste(name,".tiff",sep=""),height=10,width=10,units ='cm',
               type="tiff24nc",res=300,method="pdf")
    dev.off()
  }
}




