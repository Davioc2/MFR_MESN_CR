# Simulações

MonteCarloDavi <- matrix(NA, nrow = 500, ncol = 15)
colnames(MonteCarloDavi) <- c('Beta10', 'Beta11', 'Beta12', 'Beta13',
                              'Beta20', 'Beta21', 'Beta22', 'Beta23',
                              'Sigma2_1', 'Sigma2_2', 
                              'Lambda1', 'Lambda2', 
                              'Prop1', 'Prop2', 'Nu')

i <- 1
while(i < 501){
  ## Simulando os dados
  
  ## TAMANHO DA AMOSTRA ##
  n <- 1000
  ## ------------------ ##
  
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
  
  nu <- 5
  
  y <- numeric(n)
  for(ii in 1:n){
    arg1 <- list(mu = mu01[ii], sigma2 = sigma2_01, shape = lambda01, nu = nu)
    arg2 <- list(mu = mu02[ii], sigma2 = sigma2_02, shape = lambda02, nu = nu)
    
    obs <- rmix(1, pii = c(0.7, 0.3), family = 'Skew.t', 
                arg = list(arg1, arg2))
    y[ii] <- obs
  }
  
  ## Nível de censura ##
  nivelC <- 0.08
  ## ---------------- ##
  
  ki <- as.numeric(quantile(y, probs = nivelC))
  y[y <= ki] <- ki
  phi <- as.numeric(y == ki)
  
  rm(beta01, sigma2_01, lambda01, mu01, beta02, sigma2_02, lambda02, mu02, arg1, 
     arg2, obs, nivelC, ii, nu)
  
  # Resultados
  
  result <- MFMRCST03(y, X, ki, 2, 5, maxIter = 1000, ordem = c(1,2), 
                      showEP = F)
  
  explod <- sum(sum(result$Lambda > 10), sum(result$Lambda < -10))
  
  if(explod == 0){
    MonteCarloDavi[i,] <- c(as.vector(t(result$Beta)), result$Sigma2, 
                            result$Lambda, result$Prop, result$Nu)
    
    print(i)
    
    i <- i + 1
  }else{
    print('Explodiu')
  }
}


colMeans(MonteCarloDavi, na.rm = T)
write.csv2(MonteCarloDavi, paste0('N1000C08_ST_Parte01.csv'))

rm(list = ls())