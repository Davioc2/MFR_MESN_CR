# Estudo de simulação 03 - Classificação

##### ----- Simulação 01 - Skew-t ----- #####

CCR <- matrix(NA, nrow = 100, ncol = 2)
colnames(CCR) <- c('SN', 'ST')

i <- 69
while(i < 101){
  # Simulação
  
  n <- 400
  g <- 2
  p <- 2
  X <- cbind(rep(1, n), runif(n, -2, 2))
  
  beta01 <- c(0, -2)
  sigma2_01 <- 4
  lambda01 <- 3
  mu01 <- as.numeric(X%*%beta01)
  
  beta02 <- c(0, 4)
  sigma2_02 <- 1
  lambda02 <- 4
  mu02 <- as.numeric(X%*%beta02)
  
  nu <- 4
  
  grupo <- numeric(n)
  y <- numeric(n)
  
  for(ii in 1:n){
    arg1 <- list(mu = mu01[ii], sigma2 = sigma2_01, shape = lambda01, nu = nu)
    arg2 <- list(mu = mu02[ii], sigma2 = sigma2_02, shape = lambda02, nu = nu)
    
    obs <- rmix(1, pii = c(0.6, 0.4), family = 'Skew.t', 
                arg = list(arg1, arg2), cluster = T)
    y[ii] <- obs$y
    grupo[ii] <- obs$cluster
  }
  
  nivelC <- 0.10
  ki <- as.numeric(quantile(y, probs = nivelC))
  y[y <= ki] <- ki
  phi <- as.numeric(y == ki)
  
  rm(beta01, sigma2_01, lambda01, mu01, beta02, sigma2_02, lambda02, 
     mu02, arg1, arg2, obs, nivelC, ii, nu)
  
  
  # Skew-normal
  
  modelo01 <- MFMRCSN(y, X, ki, g = 2, maxIter = 1000, tol = 1e-05, showEP = F,
                      ordem = c(1,2))
  
  tabela01 <- as.numeric(table(modelo01$Classification))
  
  if(tabela01[1] > tabela01[2]){
    classSN <- sum(diag(table(grupo, modelo01$Classification)))/400
  }else{
    classSN <- 1 - sum(diag(table(grupo, modelo01$Classification)))/400
  }
  
  
  # Skew-t
  
  modelo02 <- MFMRCST03(y, X, ki, g = 2, maxIter = 1000, tol = 1e-05, showEP = F,
                        ordem = c(1,2), nu = 4)
  
  tabela02 <- as.numeric(table(modelo02$Classification))
  
  if(tabela02[1] > tabela02[2]){
    classST <- sum(diag(table(grupo, modelo02$Classification)))/400
  }else{
    classST <- 1 - sum(diag(table(grupo, modelo02$Classification)))/400
  }
  
  
  # Final
  
  CCR[i,] <- c(classSN, classST)
  
  print(i)
  
  print(CCR[i,])
  
  i <- i + 1
}

write.csv2(CCR, file = 'SkewT_C10_Parte01.csv')


##### ----- Simulação 02 - Skew-slash ----- #####

CCR <- matrix(NA, nrow = 100, ncol = 2)
colnames(CCR) <- c('SN', 'ST')

i <- 91
while(i < 101){
  # Simulação
  
  n <- 400
  g <- 2
  p <- 2
  X <- cbind(rep(1, n), runif(n, -2, 2))
  
  beta01 <- c(0, -2)
  sigma2_01 <- 4
  lambda01 <- 3
  mu01 <- as.numeric(X%*%beta01)
  
  beta02 <- c(0, 4)
  sigma2_02 <- 1
  lambda02 <- 4
  mu02 <- as.numeric(X%*%beta02)
  
  nu <- 2
  
  grupo <- numeric(n)
  y <- numeric(n)
  
  for(ii in 1:n){
    arg1 <- list(mu = mu01[ii], sigma2 = sigma2_01, shape = lambda01, nu = nu)
    arg2 <- list(mu = mu02[ii], sigma2 = sigma2_02, shape = lambda02, nu = nu)
    
    obs <- rmix(1, pii = c(0.6, 0.4), family = 'Skew.slash', 
                arg = list(arg1, arg2), cluster = T)
    y[ii] <- obs$y
    grupo[ii] <- obs$cluster
  }
  
  nivelC <- 0.10
  ki <- as.numeric(quantile(y, probs = nivelC))
  y[y <= ki] <- ki
  phi <- as.numeric(y == ki)
  
  rm(beta01, sigma2_01, lambda01, mu01, beta02, sigma2_02, lambda02, 
     mu02, arg1, arg2, obs, nivelC, ii, nu)
  
  
  # Skew-normal
  
  modelo01 <- MFMRCSN(y, X, ki, g = 2, maxIter = 1000, tol = 1e-05, showEP = F,
                      ordem = c(1,2))
  
  tabela01 <- as.numeric(table(modelo01$Classification))
  
  if(tabela01[1] > tabela01[2]){
    classSN <- sum(diag(table(grupo, modelo01$Classification)))/400
  }else{
    classSN <- 1 - sum(diag(table(grupo, modelo01$Classification)))/400
  }
  
  
  # Skew-t
  
  modelo02 <- MFMRCST03(y, X, ki, g = 2, maxIter = 1000, tol = 1e-05, showEP = F,
                        ordem = c(1,2), nu = 15)
  
  tabela02 <- as.numeric(table(modelo02$Classification))
  
  if(tabela02[1] > tabela02[2]){
    classST <- sum(diag(table(grupo, modelo02$Classification)))/400
  }else{
    classST <- 1 - sum(diag(table(grupo, modelo02$Classification)))/400
  }
  
  
  # Final
  
  CCR[i,] <- c(classSN, classST)
  
  print(i)
  
  print(CCR[i,])
  
  i <- i + 1
}

write.csv2(CCR, file = 'SkewSlash_C10_Parte01.csv')


