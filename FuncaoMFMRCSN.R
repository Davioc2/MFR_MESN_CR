# Função MFMRCSN

library(MomTrunc) # Sempre passar a sigma² como parâmetro nas funções
library(sn)
library(moments)
library(mixsmsn)
library(numDeriv)

MFMRCSN <- function(y, X, cen, g, maxIter = 100, tol = 1e-6, ordem, 
                    showEP = T, printOutput = T){
  ki <- cen
  n <- length(y)
  phi <- as.numeric(y == ki)
  
  p <- ncol(X)
  
  ## Função auxiliar
  
  dUnc <- function(yUnc, XUnc, beta, sigma, lambda, prop){
    auxMatriz <- matrix(NA, nrow = n-m, ncol = g)
    for(j in 1:g){
      mediaUnc <- as.numeric(XUnc%*%beta[j,])
      auxMatriz[,j] <- prop[j]*dsn(yUnc, mediaUnc, sigma[j], lambda[j])
    }
    return(rowSums(auxMatriz))
  }
  
  dCen <- function(yCen, XCen, beta, sigma, lambda, prop){
    auxMatriz <- matrix(NA, nrow = m, ncol = g)
    for(j in 1:g){
      mediaCen <- as.numeric(XCen%*%beta[j,])
      auxMatriz[,j] <- prop[j]*psn(yCen, mediaCen, sigma[j], lambda[j])
    }
    return(rowSums(auxMatriz))
  }
  
  logVero <- function(y, X, beta, sigma, lambda, prop){
    yUnc <- y[phi == 0]
    yCen <- y[phi == 1]
    
    XUnc <- X[phi == 0,]
    XCen <- X[phi == 1,]
    
    lv <- sum(log(dUnc(yUnc, XUnc, beta, sigma, lambda, prop))) + 
      sum(log(dCen(yCen, XCen, beta, sigma, lambda, prop)))
    
    return(lv)
  }
  
  EP_MFMRCSN <- function(X, y, phi, beta, sigma2, lambda, prop){
    
    # Funções auxiliares
    
    dUncEP <- function(yUnc, XUnc, beta, sigma, lambda, prop){
      g <- nrow(beta)
      aux <- numeric(g)
      
      for(j in 1:g){
        mediaUnc <- as.numeric(t(XUnc)%*%beta[j,])
        aux[j] <- prop[j]*dsn(yUnc, mediaUnc, sigma[j], lambda[j])
      }
      return(sum(aux))
    }
    
    dCenEP <- function(yCen, XCen, beta, sigma, lambda, prop){
      g <- nrow(beta)
      aux <- numeric(g)
      
      for(j in 1:g){
        mediaCen <- as.numeric(t(XCen)%*%beta[j,])
        aux[j] <- prop[j]*psn(yCen, mediaCen, sigma[j], lambda[j])
      }
      return(sum(aux))
    }
    
    logVeroEP <- function(yi, xi, g, phi, theta){
      p <- length(xi)
      
      beta <- matrix(theta[(1:(p*g))], nrow = g, byrow = T)
      sigma <- sqrt(theta[(p*g + 1):(p*g + g)])
      lambda <- theta[(p*g + g + 1):(p*g + 2*g)]
      pAux <- theta[(p*g + 2*g + 1):(p*g + 3*g - 1)]
      
      prop <- c(pAux, 1-sum(pAux))
      
      if(phi == 1){
        lv <- log(dCenEP(yi, xi, beta, sigma, lambda, prop))
      }else{
        lv <- log(dUncEP(yi, xi, beta, sigma, lambda, prop))
      }
      
      return(lv)
    }
    
    
    # Encontrando o erro-padrão
    
    n <- length(y)
    g <- nrow(beta)
    p <- ncol(X)
    
    theta <- c(as.vector(t(beta)), sigma2, lambda, prop[-g])
    i <- 1
    
    while(i < n+1){
      score <- grad(logVeroEP, theta, yi = y[i], xi = X[i,], g = g, phi = phi[i])
      
      w <- score%*%t(score)
      
      if(i == 1){
        InfObs <- w
      }else{
        InfObs <- InfObs + w
      }
      
      i <- i + 1
    }
    
    EPGrad <- round(sqrt(diag(solve(InfObs))), 4)
    
    names(EPGrad)[(1:(p*g))] <- paste0('beta', rep(1:g, rep(p, g)), rep(0:(p-1), g))
    names(EPGrad)[(p*g + 1):(p*g + g)] <- paste0('sigma2_', 1:g)
    names(EPGrad)[(p*g + g + 1):(p*g + 2*g)] <- paste0('lambda', 1:g)
    names(EPGrad)[(p*g + 2*g + 1):(p*g + 3*g - 1)] <- paste0('p', 1:(g-1))
    
    
    # Output
    
    valCritico <- theta/EPGrad
    valorP <- round(1 - pnorm(abs(valCritico)), 4)
    names(valorP) <- names(EPGrad)
    
    output <- rbind(EPGrad, valorP)
    rownames(output) <- c('Erro-padrão', 'Valor-p')
    
    return(output)
  }
  
  printMFMRCSN <- function(theta, ep, lv, aic, bic){
    sig <- ifelse(ep[2,] > 0.05, ' ', 
                  ifelse(ep[2,] > 0.01, '*', ifelse(ep[2,] > 0.001, '**', '***')))
    tabela <- data.frame(Estimate = theta, `Std error` = ep[1,], 
                         `p-value` = paste0(ep[2,], sig))
    
    cat('------------------------------------------------------------------------ \n')
    cat('  \n')
    print(tabela)
    cat('  \n')
    cat('Loglikelihood =', lv, '   AIC =', aic, '   BIC =', bic, '\n')
    cat('  \n')
    cat('  \n')
    cat('Ps: the results take into account the normality of the ML estimators for large n.  \n')
    cat('------------------------------------------------------------------------')
  }
  
  
  ## Chute inicial
  
  dados <- as.data.frame(cbind(Y = y, X[,-1]))
  km <- kmeans(dados, centers = g, nstart = 50, iter.max = 100)
  aux <- as.numeric(names(sort(table(km$cluster))))
  
  betaNova <- matrix(NA, nrow = g, ncol = p)
  dpNova <- pNova <- assNova <- numeric(g)
  
  for(i in 1:g){
    dados0 <- dados[km$cluster == aux[ordem[i]],]
    modelo <- lm(Y ~ ., data = dados0)
    betaNova[i,] <- as.numeric(modelo$coefficients)
    dpNova[i] <- sqrt(sum((dados0$Y - modelo$fitted.values)^2)/(nrow(dados0) - p))
    assNova[i] <- skewness(modelo$residuals)
    pNova[i] <- nrow(dados0)/n
  }
  
  rm(dados, km, aux, ordem, dados0, modelo, i)
  
  
  ## Processo iterativo
  
  yUnc <- y[phi == 0]
  yCen <- y[phi == 1]
  
  XUnc <- X[phi == 0,]
  XCen <- X[phi == 1,]
  
  m <- sum(phi == 1)
  
  crit <- 1
  c <- 0
  
  while(crit > tol & c < maxIter){
    ## Etapa E
    
    betaAtual <- betaNova
    dpAtual <- dpNova
    assAtual <- assNova
    pAtual <- pNova
    
    Z <- ZY <- ZY2 <- ZT <- ZT2 <- ZTY <- matrix(NA, nrow = n, ncol = g)
    
    delta <- assAtual/sqrt(1 + assAtual^2)
    D <- dpAtual*delta
    G <- (dpAtual^2)*(1 - delta^2)
    
    for(j in 1:g){
      M <- sqrt(G[j]/(G[j] + D[j]^2))
      
      ### Dados sem censura
      
      mediaUnc <- as.numeric(XUnc%*%betaAtual[j,])
      aUnc <- assAtual[j]*(yUnc - mediaUnc)/dpAtual[j]
      aux <- ifelse(pnorm(aUnc) == 0, .Machine$double.xmin, pnorm(aUnc))
      
      #### E[Zi | Yi]
      
      ZUnc <- pAtual[j]*dsn(yUnc, mediaUnc, dpAtual[j], assAtual[j])/
        dUnc(yUnc, XUnc, betaAtual, dpAtual, assAtual, pAtual)
      
      #### E[ZiYi | Yi]
      
      YUnc <- yUnc
      
      #### E[ZiYi² | Yi]
      
      Y2Unc <- yUnc^2
      
      #### E[ZiTi | Yi]
      
      p01_T1 <- (M^2)*D[j]*(yUnc - mediaUnc)/G[j]
      p02_T1 <- M*dnorm(aUnc)/aux
      
      TUnc <- p01_T1 + p02_T1
      
      #### E[ZiTi² | Yi]
      
      T2Unc <- p01_T1^2 + p01_T1*p02_T1 + M^2
      
      rm(p01_T1, p02_T1, mediaUnc, aUnc, aux)
      
      #### E[ZiYiTi | Yi]
      
      TYUnc <- YUnc*TUnc
      
      
      ### Dados censurados
      
      mediaCen <- as.numeric(XCen%*%betaAtual[j,])
      
      auxMean <- auxEYY <- w0 <- numeric(m)
      for(i in 1:m){
        aux02 <- meanvarTMD(lower = -Inf, upper = yCen[i], mu = mediaCen[i], 
                            Sigma = dpAtual[j]^2, lambda = assAtual[j], dist = 'SN')
        auxMean[i] <- aux02$mean
        auxEYY[i] <- aux02$EYY
        
        w0[i] <- meanvarTMD(lower = -Inf, upper = yCen[i], mu = mediaCen[i], 
                            Sigma = G[j], dist = 'normal')$mean
      }
      
      P0 <- pnorm(yCen, mediaCen, sqrt(G[j]))
      R0 <- ifelse(psn(yCen, mediaCen, dpAtual[j], assAtual[j]) == 0, 
                   .Machine$double.xmin,
                   psn(yCen, mediaCen, dpAtual[j], assAtual[j]))
      gi <- P0/(sqrt(pi*(1 + assAtual[j]^2)/2)*R0)
      
      #### E[Zi | Yi < ki]
      
      ZCen <- pAtual[j]*psn(yCen, mediaCen, dpAtual[j], assAtual[j])/
        dCen(yCen, XCen, betaAtual, dpAtual, assAtual, pAtual)
      
      #### E[ZiYi | Yi < ki]
      
      YCen <- auxMean
      
      #### E[ZiYi² | Yi < ki]
      
      Y2Cen <- auxEYY
      
      #### E[ZiTi | Yi < ki]
      
      TCen <- (M^2)*D[j]*(YCen - mediaCen)/G[j] + M*gi
      
      #### E[ZiTi² | Yi < ki]
      
      p01_T2 <- (M^4)*(D[j]^2)*(Y2Cen - 2*YCen*mediaCen + mediaCen^2)/(G[j]^2)
      p02_T2 <- (M^3)*D[j]*(w0 - mediaCen)*gi/G[j] + M^2
      
      T2Cen <- p01_T2 + p02_T2
      
      rm(p01_T2, p02_T2)
      
      #### E[ZiTiYi | Yi < ki]
      
      TYCen <- (M^2)*D[j]*(Y2Cen - YCen*mediaCen)/G[j] + M*w0*gi
      
      
      ### Vetores finais 
      
      Z[phi == 0, j] <- ZUnc; Z[phi == 1, j] <- ZCen
      ZY[phi == 0, j] <- ZUnc*YUnc; ZY[phi == 1, j] <- ZCen*YCen
      ZY2[phi == 0, j] <- ZUnc*Y2Unc; ZY2[phi == 1, j] <- ZCen*Y2Cen
      ZT[phi == 0, j] <- ZUnc*TUnc; ZT[phi == 1, j] <- ZCen*TCen
      ZT2[phi == 0, j] <- ZUnc*T2Unc; ZT2[phi == 1, j] <- ZCen*T2Cen
      ZTY[phi == 0, j] <- ZUnc*TYUnc; ZTY[phi == 1, j] <- ZCen*TYCen
      
      rm(ZUnc, YUnc, Y2Unc, TUnc, T2Unc, TYUnc, ZCen, YCen, Y2Cen, TCen, 
         T2Cen, TYCen, auxEYY, auxMean, gi, i, j, P0, R0, w0, mediaCen, aux02, M)
    }
    
    
    ## Etapa M
    
    betaNova <- matrix(NA, nrow = g, ncol = p)
    DNova <- GNova <- pNova <- numeric(g)
    
    for(j in 1:g){
      pNova[j] <- sum(Z[,j])/n
      
      betaNova[j,] <- solve(t(X)%*%diag(Z[,j])%*%X)%*%(t(X)%*%(ZY[,j] - ZT[,j]*D[j]))
      
      DNova[j] <- sum(ZTY[,j] - ZT[,j]*(X%*%betaNova[j,]))/sum(ZT2[,j])
      
      GNova[j] <- sum(ZY2[,j] - 2*ZY[,j]*(X%*%betaNova[j,]) - 2*ZTY[,j]*DNova[j] + 
                        Z[,j]*(X%*%betaNova[j,])^2 + 
                        2*ZT[,j]*(X%*%betaNova[j,])*DNova[j] + ZT2[,j]*(DNova[j]^2))/
        sum(Z[,j])
    }
    
    dpNova <- sqrt(GNova + DNova^2)
    
    assNova <- DNova/sqrt(GNova)
    
    
    ## Critério de parada
    
    crit <- abs(logVero(y, X, betaNova, dpNova, assNova, pNova)/
                  logVero(y, X, betaAtual, dpAtual, assAtual, pAtual) - 1)
    c <- c + 1
  }
  
  
  ## Resultados
  
  theta <- round(c(as.vector(t(betaNova)), dpNova^2, assNova, pNova[g-1]), 4)
  
  
  ### Desempenho
  
  nPar <- length(theta)
  lv <- logVero(y, X, betaNova, dpNova, assNova, pNova)
  aic <- -2*lv + 2*nPar
  bic <- -2*lv + log(n)*nPar
  
  desempenho <- c(lv, aic, bic)
  names(desempenho) <- c('Loglikelihood', 'AIC', 'BIC')
  
  
  ### Classficação das observações
  
  class <- apply(Z, 1, FUN = function(x) which.max(x))
  
  
  ### Erro-padrão
  
  if(showEP){
    erroPadrao <- EP_MFMRCSN(X, y, phi, betaNova, dpNova^2, assNova, pNova)
    
    resultado <- list(Iterations = c, Prop = pNova, Beta = betaNova, 
                      Sigma2 = dpNova^2, Lambda = assNova,
                      StdError = erroPadrao, Performance = desempenho,
                      Classification = class)
  }else{
    resultado <- list(Iterations = c, Prop = pNova, Beta = betaNova, 
                      Sigma2 = dpNova^2, Lambda = assNova, 
                      Performance = desempenho, Classification = class)
  }
  
  
  ### Output
  
  if(showEP & printOutput){
    printMFMRCSN(theta, erroPadrao, lv, aic, bic)
  }
  
  
  return(resultado)
}