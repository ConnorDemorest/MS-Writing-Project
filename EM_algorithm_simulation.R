#Simulate EM algorithm
set.seed(10)

## Sum over only finite values because log yields infinite values sometimes
sum.finite <- function(x) {
  sum(x[is.finite(x)])
}

#Simulation scenario
pi <- 0.35; n <- 1000
mu0 <- 0; sig0 <- 1;
mu1 <- 5; sig1 <- 2
theta <- rbinom(n, 1, pi)
x <- theta*rnorm(n, mu1, sig1) + (1-theta)*rnorm(n, mu0, sig0)

plot(density(x))


#Implement EM Algorithm#

#initial values
Q <- 0
#-------Clustering and finding initial values--------#
mem <- kmeans(x,2)$cluster
# Think these are backwards! How does this affect my data?
# xtabs(~ x + mem)
mu1.new <- mean(x[mem==1])
mu0.new <- mean(x[mem==2])
sig1.new <- var(x[mem==1])
sig0.new <- var(x[mem==2])
pi1 <- p <- sum(mem==1)/length(mem)
pi2 <- p.inv <- sum(mem==2)/length(mem)
#----------------------------------------------------#
#First iteration
p.inv.numerator.update <- p.inv * dnorm(x, mu0.new, sqrt(sig0.new)) #p.inv = 1-pi
p.numerator.update <- p * dnorm(x, mu1.new, sqrt(sig1.new))

p.tot <- p.inv.numerator.update + p.numerator.update
p.new.x <- p.numerator.update /p.tot
p.inv.new.x <- p.inv.numerator.update/p.tot


Q <- sum.finite((p.new.x * log(dnorm(x, mu1.new, sqrt(sig1.new)) * p))) + 
  sum.finite(p.inv.new.x * log(dnorm(x, mu0.new, sqrt(sig0.new)) * p.inv))

k <- 2; Q[2] <- Q+1 #Q[2] is a dummy so loop is entered.

while(abs(Q[k]-Q[k-1])> 1e-6){
  
  p.numerator.update <- p * dnorm(x, mu1.new, sqrt(sig1.new))
  p.inv.numerator.update <- p.inv * dnorm(x, mu0.new, sqrt(sig0.new)) #p.inv = 1-pi
  
  
  p.tot <- p.inv.numerator.update + p.numerator.update
  p.new.x <- p.numerator.update /p.tot
  p.inv.new.x <- p.inv.numerator.update/p.tot
  
  p <- sum.finite(p.new.x)/n
  p.inv <- sum.finite(p.inv.new.x)/n
  
  mu0.new <- sum.finite(x*p.inv.new.x)/sum.finite(p.inv.new.x)
  mu1.new <- sum.finite(x*p.new.x)/sum.finite(p.new.x)
  
  sig0.new <- sum.finite( (x-mu0.new)^2*p.inv.new.x ) / sum.finite(p.inv.new.x)
  sig1.new <- sum.finite( (x-mu1.new)^2*p.new.x ) / sum.finite(p.new.x)
  
   k <- k+1
  
  Q[k] <- sum.finite((p.new.x * log(dnorm(x, mu1.new, sqrt(sig1.new)) * p))) + 
                sum.finite(p.inv.new.x * log(dnorm(x, mu0.new, sqrt(sig0.new)) * p.inv))
  
  if(k%%100 ==0) 
    print(k)
  
}
#Check performance
k
p
p.inv
mu0.new; sqrt(sig0.new);
mu1.new; sqrt(sig1.new);

abs(mu0-mu0.new);
abs(sig0-sig0.new)

abs(mu1-mu1.new);
abs(sig1-sig1.new)

#_______________________________________________________________________________#
#Another way of implementing EM using built-in  function
library(mixtools)
gm<-normalmixEM(x,k=2,lambda=c(pi1, pi2),mu=c(mean(x[mem==1], mean(x[mem==2]))),sigma=c(sd(x[mem==1]), sd(x[mem==2])))
gm$lambda #pi and 1-pi
gm$mu
gm$sigma
# not always a good result