library(tictoc)
library(mclust)
library(datasets)
library(gridExtra)
library(grid)
euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))

kmeanspp <- function(X,K){
  len <- dim(X)[1]
  res <- c(sample(1:len,1))
  max <- K-1
  for (k in seq(1:max)){# Pour avoir K points
    dist_tab <- c() # Notre vecteur de proba distance
    idx <- seq(1:len)
    idx <- idx[!(seq(1:len) %in% res)]
    for (i in idx){
      dist <- .Machine$integer.max
      for(j in res){ # Pour tous les points c_i déjà sélectionner
        aux <- euclidean_dist(X[i,],X[j,])
        if (aux < dist){
          dist <- aux
        } # On prend le min de la distance
      }
      dist_tab <- append(dist_tab,dist)
    }
    dist_tab <- dist_tab^2/sum(dist_tab^2)
    # On doit choisir un point avec une proba pondéré par la distance
    c_k <- sample(idx,1,prob=dist_tab)
    res <- append(res,c_k)
  }
  res
}



kmean <- function(X,K){
  X <- as.matrix(X)
  len <- dim(X)[1]
  size <- dim(X)[2]
  C <- kmeanspp(X,K)
  centerMass <- matrix(0,nrow = K, ncol = size)
  meanMass <- matrix(0,nrow = K, ncol = size)
  for(mass in seq(1:K)){
    l <- C[mass]
    centerMass[mass,] <- X[l,]
  }
  C <- seq(1:K)
  while(euclidean_dist(centerMass,meanMass) > 0.1){
    cluster <- seq(1:len)
    for (i in seq(1:K)){ # Pour chaque i
      for (x in seq(1:len)){ # On va tester tous les x (oui c'est pas ouf)
        b <- TRUE
        loop <- seq(1:K)
        loop <- loop[!(seq(1:K) %in% i)]
        for (j in loop){ # On va faire la comparaison de la distance avec tous les autres
          if (euclidean_dist(centerMass[i,],X[x,])>euclidean_dist(centerMass[j,],X[x,])){
            b <- FALSE
          }
        }
        if(b){
          cluster[x] <- C[i]
        }
      }
    }
    #cat(cluster)
    #cat("\n")
    meanMass <- matrix(0,nrow = K, ncol = size)
    for(i in seq(1:len)){
      for (c in C)
        if(cluster[i] == c){
          meanMass[c,] <- meanMass[c,] + X[i,]
        }
    }
    summ <- table(cluster)
    
    for(c in C){
      meanMass[c,] <- meanMass[c,]/summ[c]
    }
    
    aux <- centerMass
    centerMass <- meanMass
    meanMass <- aux
    
  }
  centerMass
}

# --------------------------------------------------
kmean2 <- function(X,K){
  len <- dim(X)[1]
  size <- dim(X)[2]
  C <- c(sample(1:len,K))
  centerMass <- matrix(nrow = K, ncol = size)
  meanMass <- matrix(0,nrow = K, ncol = size)
  for(mass in seq(1:K)){
    l <- C[mass]
    centerMass[mass,] <- X[l,]
  }
  C <- seq(1:K)
  while(euclidean_dist(centerMass,meanMass) > 0.1){
    cluster <- seq(1:len)
    for (i in seq(1:K)){ # Pour chaque i
      for (x in seq(1:len)){ # On va tester tous les x (oui c'est pas ouf)
        b <- TRUE
        loop <- seq(1:K)
        loop <- loop[!(seq(1:K) %in% i)]
        for (j in loop){ # On va faire la comparaison de la distance avec tous les autres
          if (euclidean_dist(centerMass[i,],X[x,])>euclidean_dist(centerMass[j,],X[x,])){
            b <- FALSE
          }
        }
        if(b){
          cluster[x] <- C[i]
        }
      }
    }
    #cat(cluster)
    #cat("\n")
    meanMass <- matrix(0,nrow = K, ncol = size)
    for(i in seq(1:len)){
      for (c in C)
        if(cluster[i] == c){
          meanMass[c,] <- meanMass[c,] + X[i,]
        }
    }
    summ <- table(cluster)
    
    for(c in C){
      meanMass[c,] <- meanMass[c,]/summ[c]
    }
    
    aux <- centerMass
    centerMass <- meanMass
    meanMass <- aux
    
  }
  centerMass
}

# n doit être un multiple de x
normx <- function(x,d,n){
  if (n%%x != 0){
    stop("n n est pas un multiple de x")
  }
  len <- n/x
  center <- matrix(nrow = x, ncol = d)
  res <- matrix(nrow = n, ncol = d)
  for (i in seq(1:x)){
    center[i,] <- runif(d,0,500)
    #cat(center[i,])
    #cat("\n")
    for (j in seq(1:d)){
      vec <- rnorm(len,mean=center[i,j],sd=1)
      debut <- ((i-1)*len)+1
      fin <- i*len
      res[debut:fin,j] <- vec
    }
    res[i*len,] <- center[i,]
  }
  res
}


phi <- function(X,C){
  res <- 0
  for (x in X){
    minimum <- .Machine$integer.max
    for (c in C){
      d <- euclidean_dist(x,c)
      d <- d^2
      if (d < minimum){
        minimum <- d
      }
    }
    res <- res + minimum
  }
  res
}

classif <- function(X,C){
  X <- as.matrix(X)
  lx <- dim(X)[1]
  res <- c()
  l = dim(C)[1]
  for (x in seq(1:lx)){
    minimum <- .Machine$integer.max
    idx <- 1
    for (i in seq(1:l)){
      d <- euclidean_dist(X[x,],C[i,])
      if (d < minimum){
        minimum <- d
        idx <- i
      }
    }
    res <- append(res,idx)
  }
  res
}

# Variable importante 

n = 1000
X <- normx(25,5,n)
# -------------------------------------#
# Pour k = 10
cat("pour k = 10")
cat("\n")
cat("Boucle pour kmean")

mean_kclas10 <- c()
mean_timek10 <- c()
for (loop in seq(1:20)){
  tic()
  cat("Boucle ",loop)
  cat("\n")
  #kpp <- kmean(X,10)
  #mean_kpp <- append(mean_kpp,phi(X,kpp)/n)
  kclas <- kmean2(X,10)
  mean_kclas10 <- append(mean_kclas10,phi(X,kclas)/n)
  time <- toc(quiet = TRUE)
  mean_timek10 <- append(mean_timek10, time$toc - time$tic)
  
}

#kclas <- kmeans(X,10,iter.max= 15, nstart = 1)

cat("\n")
cat("Boucle pour kmean++")

mean_kpp10 <- c()
mean_timekpp10 <- c()
for (loop in seq(1:20)){
  tic()
  cat("Boucle ",loop)
  cat("\n")
  kpp <- kmean(X,10)
  mean_kpp10 <- append(mean_kpp10,phi(X,kpp)/n)
  time <- toc(quiet = TRUE)
  mean_timekpp10 <- append(mean_timekpp10, time$toc - time$tic)
}

# -------------------------------------#
# Pour k = 25

cat("pour k = 25")
cat("\n")
cat("Boucle pour kmean")

mean_kclas25 <- c()
mean_timek25 <- c()
for (loop in seq(1:20)){
  tic()
  cat("Boucle ",loop)
  cat("\n")
  kclas <- kmean2(X,25)
  mean_kclas25 <- append(mean_kclas25,phi(X,kclas)/n)
  time <- toc(quiet = TRUE)
  mean_timek25 <- append(mean_timek25, time$toc - time$tic)
  
}


cat("\n")
cat("Boucle pour kmean++")

mean_kpp25 <- c()
mean_timekpp25 <- c()
for (loop in seq(1:20)){
  tic()
  cat("Boucle ",loop)
  cat("\n")
  kpp <- kmean(X,25)
  mean_kpp25 <- append(mean_kpp25,phi(X,kpp)/n)
  time <- toc(quiet = TRUE)
  mean_timekpp25 <- append(mean_timekpp25, time$toc - time$tic)
}

# -------------------------------------#
# Pour k = 50

cat("pour k = 50")
cat("\n")
cat("Boucle pour kmean")

mean_kclas50 <- c()
mean_timek50 <- c()
for (loop in seq(1:20)){
  tic()
  cat("Boucle ",loop)
  cat("\n")
  kclas <- kmean2(X,50)
  mean_kclas50 <- append(mean_kclas50,phi(X,kclas)/n)
  time <- toc(quiet = TRUE)
  mean_timek50 <- append(mean_timek50, time$toc - time$tic)
  
}

cat("\n")
cat("Boucle pour kmean++")


mean_kpp50 <- c()
mean_timekpp50 <- c()
for (loop in seq(1:20)){
  tic()
  cat("Boucle ",loop)
  cat("\n")
  kpp <- kmean(X,50)
  mean_kpp50 <- append(mean_kpp50,phi(X,kpp)/n)
  time <- toc(quiet = TRUE)
  mean_timekpp50 <- append(mean_timekpp50, time$toc - time$tic)
}

k <- c(10,25,50)
mean_phi_k <- c(mean(mean_kclas10),mean(mean_kclas25),mean(mean_kclas50))
mean_phi_kpp <- c(mean(mean_kpp10),mean(mean_kpp25),mean(mean_kpp50))
min_phi_k <- c(min(mean_kclas10),min(mean_kclas25),min(mean_kclas50))
min_phi_kpp <- c(min(mean_kpp10),min(mean_kpp25),min(mean_kpp50))
time_k <- c(mean(mean_timek10),mean(mean_timek25),mean(mean_timek50))
time_kpp <- c(mean(mean_timekpp10),mean(mean_timekpp25),mean(mean_timekpp50))

# mean_phi_k <- c(mean(c(1,2,3,4,5,6)),12,23)
# mean_phi_kpp <- c(mean(c(1,2,3,4,5,6)),12,42)
# min_phi_k <- c(mean(c(1,2,3,4,5,6)),2,31)
# min_phi_kpp <- c(mean(c(1,2,3,4,5,6,1,1,1,1)),2,3)
# time_k <- c(mean(c(1,2,3,4,5,6)),2,23)
# time_kpp <- c(mean(c(1,2,3,4,5,6)),82,3)
final = data.frame(k,mean_phi_k,mean_phi_kpp,min_phi_k,min_phi_kpp,time_k,time_kpp)
g <- tableGrob(final, rows = NULL)
grid.draw(g)

##### Data IRIS

data(iris)

kpp <- kmean(iris[,1:4],3)
km <- kmeans(iris[,1:4],3)
kclust = Mclust(iris[-5],G=3,modelNames="VEV")
t(kclust$parameters$mean)

clustkpp <- classif(iris[-5],kpp)

table(clustkpp)
table(km$cluster)
table(kclust$classification)

pcaIris <- prcomp(iris[-5])

plot(pcaIris$x[,1],pcaIris$x[,2],col=clustkpp)