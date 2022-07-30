# Mistura finita de modelos de regressão ST com resposta censurada à esquerda

library(MomTrunc) # Sempre passar a sigma² como parâmetro nas funções
library(sn)
library(moments)
library(mixsmsn)
library(numDeriv)

# Função com apenas 1 nu

MFMRCST03 <- function(y, X, cen, g, nu, maxIter = 100, tol = 1e-6, showEP = T, 
                      printOutput = T, ordem){
  ki <- cen
  n <- length(y)
  phi <- as.numeric(y == ki)
  
  p <- ncol(X)
  
  
  ## Função auxiliar
  
  dUnc <- function(yUnc, XUnc, beta, sigma, lambda, nu, prop){
    auxMatriz <- matrix(NA, nrow = n-m, ncol = g)
    for(j in 1:g){
      mediaUnc <- as.numeric(XUnc%*%beta[j,])
      auxMatriz[,j] <- prop[j]*dst(yUnc, mediaUnc, sigma[j], lambda[j], nu)
    }
    return(rowSums(auxMatriz))
  }
  
  dCen <- function(yCen, XCen, beta, sigma, lambda, nu, prop){
    auxMatriz <- matrix(NA, nrow = m, ncol = g)
    for(j in 1:g){
      mediaCen <- as.numeric(XCen%*%beta[j,])
      auxMatriz[,j] <- prop[j]*pst(yCen, mediaCen, sigma[j], lambda[j], nu)
    }
    return(rowSums(auxMatriz))
  }
  
  logVero <- function(y, X, beta, sigma, lambda, nu, prop){
    yUnc <- y[phi == 0]
    yCen <- y[phi == 1]
    
    XUnc <- X[phi == 0,]
    XCen <- X[phi == 1,]
    
    lv <- sum(log(dUnc(yUnc, XUnc, beta, sigma, lambda, nu, prop))) + 
      sum(log(dCen(yCen, XCen, beta, sigma, lambda, nu, prop)))
    
    return(lv)
  }
  
  EP_MFMRCST <- function(X, y, phi, beta, sigma2, lambda, prop, nu){
    
    # Funções auxiliares
    
    dUncEP <- function(yUnc, XUnc, beta, sigma, lambda, prop, nu){
      g <- nrow(beta)
      aux <- numeric(g)
      
      for(j in 1:g){
        mediaUnc <- as.numeric(t(XUnc)%*%beta[j,])
        aux[j] <- prop[j]*dst(yUnc, mediaUnc, sigma[j], lambda[j], nu)
      }
      return(sum(aux))
    }
    
    dCenEP <- function(yCen, XCen, beta, sigma, lambda, prop, nu){
      g <- nrow(beta)
      aux <- numeric(g)
      
      for(j in 1:g){
        mediaCen <- as.numeric(t(XCen)%*%beta[j,])
        aux[j] <- prop[j]*pst(yCen, mediaCen, sigma[j], lambda[j], nu)
      }
      return(sum(aux))
    }
    
    logVeroEP <- function(yi, xi, g, phi, nu, theta){
      p <- length(xi)
      
      beta <- matrix(theta[(1:(p*g))], nrow = g, byrow = T)
      sigma <- sqrt(theta[(p*g + 1):(p*g + g)])
      lambda <- theta[(p*g + g + 1):(p*g + 2*g)]
      pAux <- theta[(p*g + 2*g + 1):(p*g + 3*g - 1)]
      
      prop <- c(pAux, 1-sum(pAux))
      
      if(phi == 1){
        lv <- log(dCenEP(yi, xi, beta, sigma, lambda, prop, nu))
      }else{
        lv <- log(dUncEP(yi, xi, beta, sigma, lambda, prop, nu))
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
      score <- grad(logVeroEP, theta, yi = y[i], xi = X[i,], g = g, phi = phi[i], 
                    nu = nu)
      
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
  
  printMFMRCST <- function(theta, nu, ep, lv, aic, bic){
    sig <- ifelse(ep[2,] > 0.05, ' ', 
                  ifelse(ep[2,] > 0.01, '*', ifelse(ep[2,] > 0.001, '**', '***')))
    tabela <- data.frame(Estimate = theta, `Std error` = ep[1,], 
                         `p-value` = paste0(ep[2,], sig))
    
    cat('------------------------------------------------------------------------ \n')
    cat('  \n')
    print(tabela)
    cat('  \n')
    cat('nu =', nu, '\n')
    cat('  \n')
    cat('Loglikelihood =', lv, '   AIC =', aic, '   BIC =', bic, '\n')
    cat('  \n')
    cat('  \n')
    cat('Ps: the results take into account the normality of the ML estimators for large n. \n')
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
  
  nuNova <- nu
  
  rm(dados, km, aux, ordem, dados0, modelo, i, nu)
  
  
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
    nuAtual <- nuNova
    pAtual <- pNova
    
    Z <- ZU <- ZUY <- ZUY2 <- ZUT <- ZUT2 <- ZUTY <- matrix(NA, nrow = n, ncol = g)
    
    delta <- assAtual/sqrt(1 + assAtual^2)
    D <- dpAtual*delta
    G <- (dpAtual^2)*(1 - delta^2)
    
    for(j in 1:g){
      M <- sqrt(G[j]/(G[j] + D[j]^2))
      
      ### Dados sem censura
      
      mediaUnc <- as.numeric(XUnc%*%betaAtual[j,])
      aUnc <- assAtual[j]*(yUnc - mediaUnc)/dpAtual[j]
      d2Unc <- ((yUnc - mediaUnc)^2)/(dpAtual[j]^2)
      
      #### E[Zi | Yi]
      
      ZUnc <- pAtual[j]*dst(yUnc, mediaUnc, dpAtual[j], assAtual[j], nuAtual)/
        dUnc(yUnc, XUnc, betaAtual, dpAtual, assAtual, nuAtual, pAtual)
      
      #### E[Ui | Yi]
      
      p01_U <- 4*(nuAtual^(nuAtual/2))*gamma((nuAtual + 3)/2)/
        (sqrt(pi)*gamma(nuAtual/2)*dpAtual[j])
      p02_U <- (nuAtual + d2Unc)^(-(nuAtual + 3)/2)*
        pt(aUnc*sqrt((nuAtual + 3)/(nuAtual + d2Unc)), nuAtual + 3)
      p03_U <- ifelse(dst(yUnc, mediaUnc, dpAtual[j], assAtual[j], nuAtual) == 0,
                      .Machine$double.xmin,
                      dst(yUnc, mediaUnc, dpAtual[j], assAtual[j], nuAtual))
      
      UUnc <- p01_U*p02_U/p03_U
      
      rm(p01_U, p02_U)
      
      #### E[UiYi | Yi]
      
      UYUnc <- UUnc*yUnc
      
      #### E[UiYi² | Yi]
      
      UY2Unc <- UUnc*(yUnc^2)
      
      #### E[UiTi | Yi]
      
      p01_UT <- (M^2)*D[j]*(yUnc - mediaUnc)*UUnc/G[j]
      p02_UT <- 2*(nuAtual^(nuAtual/2))*gamma((nuAtual + 2)/2)/
        (pi*gamma(nuAtual/2)*dpAtual[j])
      p03_UT <- (nuAtual + d2Unc + aUnc^2)^(-(nuAtual + 2)/2)
      
      UTUnc <- p01_UT + M*p02_UT*p03_UT/p03_U
      
      rm(p01_UT)
      
      #### E[UiTi² | Yi]
      
      p01_UT2 <- (((M^2)*D[j]*(yUnc - mediaUnc)/G[j])^2)*UUnc
      p02_UT2 <- M^2 + (M^3)*D[j]*(yUnc - mediaUnc)*p02_UT*p03_UT/(p03_U*G[j])
      
      UT2Unc <- p01_UT2 + p02_UT2
      
      rm(p03_U, p02_UT, p03_UT, p01_UT2, p02_UT2)
      
      #### E[UiYiTi | Yi]
      
      UTYUnc <- yUnc*UTUnc
      
      rm(d2Unc, aUnc)
      
      
      ### Dados censurados
      
      mediaCen <- as.numeric(XCen%*%betaAtual[j,])
      
      dpAtual02 <- sqrt((nuAtual*(dpAtual[j]^2))/(nuAtual + 2))
      
      dpAtual03 <- sqrt((nuAtual*(dpAtual[j]^2))/((nuAtual + 1)*(1 + assAtual[j]^2)))
      
      cNu <- 2*gamma((nuAtual + 1)/2)/
        (gamma(nuAtual/2)*sqrt(nuAtual*(1 + assAtual[j]^2)*pi))
      
      aux <- ifelse(pst(yCen, mediaCen, dpAtual[j], assAtual[j], nuAtual) == 0,
                    .Machine$double.xmin,
                    pst(yCen, mediaCen, dpAtual[j], assAtual[j], nuAtual))
      
      wPhi <- pst(yCen, mediaCen, dpAtual03, 0, nuAtual + 1)/aux
      
      auxMean <- auxEYY <- w0 <- numeric(m)
      for(i in 1:m){
        aux02 <- meanvarTMD(lower = -Inf, upper = yCen[i], mu = mediaCen[i], 
                            Sigma = dpAtual02^2, lambda = assAtual[j], 
                            nu = nuAtual + 2, dist = 'ST')
        auxMean[i] <- aux02$mean
        auxEYY[i] <- aux02$EYY
        
        w0[i] <- meanvarTMD(lower = -Inf, upper = yCen[i], mu = mediaCen[i], 
                            Sigma = dpAtual03^2, lambda = 0, nu = nuAtual + 1, 
                            dist = 'ST')$mean
      }
      
      #### E[Zi | Yi < ki]
      
      ZCen <- pAtual[j]*pst(yCen, mediaCen, dpAtual[j], assAtual[j], nuAtual)/
        dCen(yCen, XCen, betaAtual, dpAtual, assAtual, nuAtual, pAtual)
      
      #### E[Ui | Yi < ki]
      
      UCen <- pst(yCen, mediaCen, dpAtual02, assAtual[j], nuAtual + 2)/aux
      
      #### E[UiYi | Yi < ki]
      
      UYCen <- auxMean*UCen
      
      #### E[UiYi² | Yi < ki]
      
      UY2Cen <- auxEYY*UCen
      
      #### E[UiTi | Yi < ki]
      
      p01_UT <- D[j]*(UYCen - UCen*mediaCen)/(G[j] + D[j]^2)
      p02_UT <- M*cNu*wPhi
      
      UTCen <- p01_UT + p02_UT
      
      rm(p01_UT)
      
      #### E[UiTi² | Yi < ki]
      
      p01_UT2 <- ((D[j]/(G[j] + D[j]^2))^2)*
        (UY2Cen - 2*UYCen*mediaCen + UCen*(mediaCen^2))
      p02_UT2 <- p02_UT*D[j]*(w0 - mediaCen)/(G[j] + D[j]^2) + M^2
      
      UT2Cen <- p01_UT2 + p02_UT2
      
      rm(p01_UT2, p02_UT2)
      
      #### E[UiTiYi | Yi < ki]
      
      p01_UTY <- D[j]*(UY2Cen - UYCen*mediaCen)/(G[j] + D[j]^2)
      p02_UTY <- p02_UT*w0
      
      UTYCen <- p01_UTY + p02_UTY
      
      rm(p02_UT, p01_UTY, p02_UTY, dpAtual02, dpAtual03, cNu, aux, wPhi, aux02, 
         auxMean, auxEYY, w0)
      
      
      ### Vetores finais 
      
      Z[phi == 0, j] <- ZUnc; Z[phi == 1, j] <- ZCen
      ZU[phi == 0, j] <- ZUnc*UUnc; ZU[phi == 1, j] <- ZCen*UCen 
      ZUY[phi == 0, j] <- ZUnc*UYUnc; ZUY[phi == 1, j] <- ZCen*UYCen
      ZUY2[phi == 0, j] <- ZUnc*UY2Unc; ZUY2[phi == 1, j] <- ZCen*UY2Cen
      ZUT[phi == 0, j] <- ZUnc*UTUnc; ZUT[phi == 1, j] <- ZCen*UTCen
      ZUT2[phi == 0, j] <- ZUnc*UT2Unc; ZUT2[phi == 1, j] <- ZCen*UT2Cen
      ZUTY[phi == 0, j] <- ZUnc*UTYUnc; ZUTY[phi == 1, j] <- ZCen*UTYCen
      
      rm(ZUnc, UUnc, UYUnc, UY2Unc, UTUnc, UT2Unc, UTYUnc, ZCen, UCen, UYCen, UY2Cen,
         UTCen, UT2Cen, UTYCen, i, j, mediaCen)
    }
    
    
    ## Etapa M
    
    betaNova <- matrix(NA, nrow = g, ncol = p)
    DNova <- GNova <- pNova <- numeric(g)
    
    for(j in 1:g){
      pNova[j] <- sum(Z[,j])/n
      
      betaNova[j,] <- solve(t(X)%*%diag(ZU[,j])%*%X)%*%(t(X)%*%(ZUY[,j] - ZUT[,j]*D[j]))
      
      DNova[j] <- sum(ZUTY[,j] - ZUT[,j]*(X%*%betaNova[j,]))/sum(ZUT2[,j])
      
      GNova[j] <- sum(ZUY2[,j] - 2*ZUY[,j]*(X%*%betaNova[j,]) - 2*ZUTY[,j]*DNova[j] + 
                        ZU[,j]*(X%*%betaNova[j,])^2 + 
                        2*ZUT[,j]*(X%*%betaNova[j,])*DNova[j] + ZUT2[,j]*(DNova[j]^2))/
        sum(Z[,j])
    }
    
    dpNova <- sqrt(GNova + DNova^2)
    
    assNova <- DNova/sqrt(GNova)
    
    nuNova <- nuAtual
    
    
    ## Critério de parada
    
    crit <- abs(logVero(y, X, betaNova, dpNova, assNova, nuNova, pNova)/
                  logVero(y, X, betaAtual, dpAtual, assAtual, nuAtual, pAtual) - 1)
    c <- c + 1
  }
  
  ## Resultados
  
  theta <- round(c(as.vector(t(betaNova)), dpNova^2, assNova, pNova[-g]), 4)
  
  
  ### Desempenho
  
  nPar <- length(theta)
  lv <- logVero(y, X, betaNova, dpNova, assNova, nuNova, pNova)
  aic <- -2*lv + 2*nPar
  bic <- -2*lv + log(n)*nPar
  
  desempenho <- c(lv, aic, bic)
  names(desempenho) <- c('Loglikelihood', 'AIC', 'BIC')
  
  
  ### Classficação das observações
  
  class <- apply(Z, 1, FUN = function(x) which.max(x))
  
  
  ### Erro-padrão
  
  if(showEP){
    erroPadrao <- EP_MFMRCST(X, y, phi, betaNova, dpNova^2, assNova, pNova, 
                             nuNova)
    
    resultado <- list(Iterations = c, Prop = pNova, Beta = betaNova, 
                      Sigma2 = dpNova^2, Lambda = assNova, Nu = nuNova,
                      StdError = erroPadrao, Performance = desempenho,
                      Classification = class)
  }else{
    resultado <- list(Iterations = c, Prop = pNova, Beta = betaNova, 
                      Sigma2 = dpNova^2, Lambda = assNova, Nu = nuNova,
                      Performance = desempenho, Classification = class)
  }
  
  
  ### Output
  
  if(showEP & printOutput){
    printMFMRCST(theta, nuNova, erroPadrao, lv, aic, bic)
  }
  
  
  return(resultado)
}
