# Simulação MFMRCSN

MonteCarloDavi <- matrix(NA, nrow = 500, ncol = 12)
colnames(MonteCarloDavi) <- c('Beta10', 'Beta11', 'Beta12',
                              'Beta20', 'Beta21', 'Beta22',
                              'Sigma2_1', 'Sigma2_2', 
                              'Lambda1', 'Lambda2', 
                              'Prop1', 'Prop2')

i <- 1
while(i < 501){
  ## Simulando os dados
  
  ## TAMANHO DA AMOSTRA ##
  n <- 500
  ## ------------------ ##
  
  g <- 2
  p <- 3
  X <- cbind(rep(1, n), rnorm(n, 5, 3), rnorm(n, 20, 4))
  
  beta01 <- c(-1, -4, -3)
  sigma2_01 <- 3
  lambda01 <- 4
  mu01 <- X%*%beta01
  
  beta02 <- c(3, 7, 9)
  sigma2_02 <- 1
  lambda02 <- 3
  mu02 <- X%*%beta02
  
  y <- numeric(n)
  for(k in 1:n){
    arg1 <- list(mu = mu01[k,], sigma2 = sigma2_01, shape = lambda01)
    arg2 <- list(mu = mu02[k,], sigma2 = sigma2_02, shape = lambda02)
    
    obs <- rmix(1, pii = c(0.7, 0.3), family = 'Skew.normal', 
                arg = list(arg1, arg2))
    y[k] <- obs
  }
  
  ## Nível de censura ##
  nivelC <- 0.08
  ## ---------------- ##
  
  ki <- as.numeric(quantile(y, probs = nivelC))
  y[y <= ki] <- ki
  phi <- as.numeric(y == ki)
  
  rm(beta01, sigma2_01, lambda01, mu01, beta02, sigma2_02, lambda02, mu02, arg1, 
     arg2, obs, nivelC, k)
  
  # Resultados
  
  result <- MFMRCSN(y, X, cen = ki, g = 2, maxIter = 2000, ordem = c(2,1),
                    showEP = F, printOutput = F)
  
  explod <- sum(sum(result$Lambda > 10), sum(result$Lambda < -10))
  
  if(explod == 0){
    MonteCarloDavi[i,] <- c(as.vector(t(result$Beta)), result$Sigma2, 
                            result$Lambda, result$Prop)
    
    print(i)
    
    i <- i + 1
  }else{
    print('Explodiu')
  }
}


colMeans(MonteCarloDavi, na.rm = T)
write.csv2(MonteCarloDavi, 'N500C8_SN.csv')