# Estudo de simulação 02

library(tictoc)

##### ----- Opção 01 ----- #####

MatrizAIC <- MatrizBIC <- matrix(NA, nrow = 28, ncol = 8)
colnames(MatrizAIC) <- colnames(MatrizBIC) <- c('N1', 'N2', 'N3', 'N4', 'T1', 
                                                'T2', 'T3', 'T4')

i <- 1
while(i < 29){
  
  # Simulando um modelo skew-t com 3 grupos e nu = 4
  
  n <- 1000
  g <- 3
  p <- 2
  X <- cbind(rep(1, n), runif(n, -2, 2))
  
  beta01 <- c(-4, 4)
  sigma2_01 <- 1
  lambda01 <- -2
  mu01 <- as.numeric(X%*%beta01)
  
  beta02 <- c(0, -2)
  sigma2_02 <- 4
  lambda02 <- 3
  mu02 <- as.numeric(X%*%beta02)
  
  beta03 <- c(0, 4)
  sigma2_03 <- 1
  lambda03 <- 4
  mu03 <- as.numeric(X%*%beta03)
  
  nu <- 4
  
  grupo <- numeric(n)
  y <- numeric(n)
  
  for(ii in 1:n){
    arg1 <- list(mu = mu01[ii], sigma2 = sigma2_01, shape = lambda01, nu = nu)
    arg2 <- list(mu = mu02[ii], sigma2 = sigma2_02, shape = lambda02, nu = nu)
    arg3 <- list(mu = mu03[ii], sigma2 = sigma2_03, shape = lambda03, nu = nu)
    
    obs <- rmix(1, pii = c(0.5, 0.2, 0.3), family = 'Skew.t', 
                arg = list(arg1, arg2, arg3), cluster = T)
    y[ii] <- obs$y
    grupo[ii] <- obs$cluster
  }
  
  nivelC <- 0.30
  ki <- as.numeric(quantile(y, probs = nivelC))
  y[y <= ki] <- ki
  phi <- as.numeric(y == ki)
  
  rm(beta01, sigma2_01, lambda01, mu01, beta02, sigma2_02, lambda02, 
     mu02, beta03, sigma2_03, lambda03, mu03, arg1, 
     arg2, arg3, obs, nivelC, ii, nu)
  
  
  # Ajustando skew-normal
  
  resultSN_1 <- MFMRCSN(y, X, ki, g = 1, maxIter = 2000, tol = 1e-04, showEP = F,
                        ordem = c(1))
  
  resultSN_2 <- MFMRCSN(y, X, ki, g = 2, maxIter = 2000, tol = 1e-04, showEP = F,
                        ordem = c(1, 2))
  
  resultSN_3 <- MFMRCSN(y, X, ki, g = 3, maxIter = 2000, tol = 1e-04, showEP = F,
                        ordem = c(2, 1, 3))
  
  resultSN_4 <- MFMRCSN(y, X, ki, g = 4, maxIter = 2000, tol = 1e-04, showEP = F,
                        ordem = c(1, 2, 3, 4))
  
  
  # Ajustando skew-t
  
  resultST_1 <- MFMRCST03(y, X, ki, g = 1, maxIter = 2000, tol = 1e-04, 
                          showEP = F, ordem = c(1), nu = 4)
  
  resultST_2 <- MFMRCST03(y, X, ki, g = 2, maxIter = 2000, tol = 1e-04, 
                          showEP = F, ordem = c(1,2), nu = 4)
  
  resultST_3 <- MFMRCST03(y, X, ki, g = 3, maxIter = 2000, tol = 1e-04, 
                          showEP = F, ordem = c(2, 1, 3), nu = 4)
  
  resultST_4 <- MFMRCST03(y, X, ki, g = 4, maxIter = 2000, tol = 1e-04, 
                          showEP = F, ordem = c(1, 2, 3, 4), nu = 4)
  
  
  # Alocando resultados
  
  MatrizAIC[i,] <- c(resultSN_1$Performance[2], resultSN_2$Performance[2],
                     resultSN_3$Performance[2], resultSN_4$Performance[2],
                     resultST_1$Performance[2], resultST_2$Performance[2],
                     resultST_3$Performance[2], resultST_4$Performance[2])
  
  MatrizBIC[i,] <- c(resultSN_1$Performance[3], resultSN_2$Performance[3],
                     resultSN_3$Performance[3], resultSN_4$Performance[3],
                     resultST_1$Performance[3], resultST_2$Performance[3],
                     resultST_3$Performance[3], resultST_4$Performance[3])
  
  print(i)
  
  i <- i + 1
}

write.csv2(MatrizAIC, file = 'AIC_C30_Parte03.csv')
write.csv2(MatrizBIC, file = 'BIC_C30_Parte03.csv')


##### ----- Opção 02 ----- #####

# Simulando um modelo skew-t com 2 grupos, nu = 3 e censura de 5%

n <- 1000
g <- 2
p <- 4
X <- cbind(rep(1, n), runif(n, 1, 5), runif(n, -2, 2), runif(n, 1, 4))

beta01 <- c(0, -1, -2, -3)
sigma2_01 <- 1
lambda01 <- -2
mu01 <- as.numeric(X%*%beta01)

beta02 <- c(-1, 1, 2, 3)
sigma2_02 <- 2
lambda02 <- 3
mu02 <- as.numeric(X%*%beta02)

nu <- 4

grupo <- numeric(n)
y <- numeric(n)

for(i in 1:n){
  arg1 <- list(mu = mu01[i], sigma2 = sigma2_01, shape = lambda01, nu = nu)
  arg2 <- list(mu = mu02[i], sigma2 = sigma2_02, shape = lambda02, nu = nu)
  
  obs <- rmix(1, pii = c(0.7, 0.3), family = 'Skew.t', 
              arg = list(arg1, arg2), cluster = T)
  y[i] <- obs$y
  grupo[i] <- obs$cluster
}

nivelC <- 0.05
ki <- as.numeric(quantile(y, probs = nivelC))
y[y <= ki] <- ki
phi <- as.numeric(y == ki)

rm(beta01, sigma2_01, lambda01, mu01, beta02, sigma2_02, lambda02, mu02, arg1, 
   arg2, obs, nivelC, i, nu)

tic()
# Ajustando skew-normal

resultSN_1 <- MFMRCSN(y, X, ki, g = 1, maxIter = 2000, tol = 1e-04, showEP = F,
                      ordem = c(1))

resultSN_2 <- MFMRCSN(y, X, ki, g = 2, maxIter = 2000, tol = 1e-04, showEP = F,
                      ordem = c(2, 1))

resultSN_3 <- MFMRCSN(y, X, ki, g = 3, maxIter = 2000, tol = 1e-04, showEP = F,
                      ordem = c(1, 2, 3))

# Ajustando skew-t

resultST_1 <- MFMRCST03(y, X, ki, g = 1, maxIter = 2000, tol = 1e-04, showEP = F,
                        ordem = c(1), nu = 4)

resultST_2 <- MFMRCST03(y, X, ki, g = 2, maxIter = 2000, tol = 1e-04, showEP = F,
                        ordem = c(2,1), nu = 4)

resultST_3 <- MFMRCST03(y, X, ki, g = 3, maxIter = 2000, tol = 1e-04, showEP = F,
                        ordem = c(1, 2, 3), nu = 4)
toc()