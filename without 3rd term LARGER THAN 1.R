#Revised Markov Chain model-----------------------------------------------------------------------
LFC <- function(Yt2, Yt1, Yt,phi1, phi2, mu0, mu1){
  
  gamma1 <- log((1-mu1)/(1-mu0))
  
  gamma2 <- log((mu1*(1-mu0))/(mu0*(1-mu1)))
  
  
  ZT1 <- (Yt - mu1) / (sqrt(mu1 * (1 - mu1)))
  
  ZT11 <- (Yt1 - mu1) / (sqrt(mu1 * (1 - mu1)))
  
  ZT21 <- (Yt2 - mu1) / (sqrt(mu1 * (1 - mu1)))
  
  ZT0 <- (Yt - mu0) / (sqrt(mu0 * (1 - mu0)))
  
  ZT10 <- (Yt1 - mu0) / (sqrt(mu0 * (1 - mu0)))
  
  ZT20 <- (Yt2 - mu0) / (sqrt(mu0 * (1 - mu0)))
  
  rho_term1 <- (1 + phi1 * ZT1 * ZT11 + phi1 * ZT11 * ZT21 + phi2 * ZT1 * ZT21) / (1 + phi1 * ZT11 * ZT21) 
  rho_term2 <- (1 + phi1 * ZT0 * ZT10 + phi1 * ZT10 * ZT20 + phi2 * ZT0 * ZT20) / (1 + phi1 * ZT10 * ZT20) 
  
  rho_term <- log(rho_term1 / rho_term2)
  
  L <- Yt + (rho_term / gamma2) + gamma1 / gamma2
  
  return(L)
}

L_second <- function(mu0, mu1, phi1, phi2){
  
  l111 <- LFC(1, 1, 1, phi1, phi2,mu0, mu1)
  
  l011 <- LFC(0, 1, 1, phi1, phi2, mu0, mu1)
  
  l001 <- LFC(0, 0, 1, phi1, phi2, mu0, mu1)
  
  l101 <- LFC(1, 0, 1, phi1, phi2, mu0, mu1)
  
  l110 <- LFC(1, 1, 0, phi1, phi2, mu0, mu1)
  
  l000 <- LFC(0, 0, 0, phi1, phi2, mu0, mu1) 
  
  l010 <- LFC(0, 1, 0, phi1, phi2, mu0, mu1)
  
  l100 <- LFC(1, 0, 0, phi1, phi2, mu0, mu1)
  
  m <- round(1/abs(l000),digits=0)
  
  l111_star <- round(l111 * m,digits=0) / m
  
  l011_star <- round(l011 * m,digits=0) / m
  
  l001_star <- round(l001 * m,digits=0) / m
  
  l101_star <- round(l101 * m,digits=0) / m
  
  l110_star <- round(l110 * m,digits=0) / m
  
  l000_star <- round(l000 * m,digits=0) / m
  
  l010_star <- round(l010 * m,digits=0) / m
  
  l100_star <- round(l100 * m,digits=0) / m
  
  
  return(list(l111=l111, l011=l011, l001=l001, l101=l101, l110=l110, l000=l000, l010=l010, l100=l100, m=m, 
              l111_star = l111_star, l011_star=l011_star, l001_star=l001_star, l101_star=l101_star, 
              l110_star=l110_star, l000_star=l000_star, l010_star=l010_star, l100_star=l100_star))
  
}



#Function to calculate the proability, given in-control parameter mu0 
prob <- function(mu, Yt2, Yt1, Yt, phi1, phi2 ,phi3){
  
  term1 <- (mu^Yt) * ((1 - mu)^(1-Yt))
  
  Zt1 <- (Yt1 - mu) / sqrt(mu * (1 - mu))
  
  Zt2 <- (Yt2 - mu) / sqrt(mu * (1 - mu))
  
  Zt <- (Yt- mu) / sqrt(mu * (1 - mu))
  
  p <- (term1 * (1 + phi1 * Zt1 * Zt2 + phi1 * Zt * Zt1 + phi2 * Zt * Zt2 + phi3*Zt*Zt1*Zt2)) / (1 + phi1 * Zt1 * Zt2)
  
  return(list(p=p))
  
}



N_vector <- function(h, mu0, mu1, p_true, phi1, phi2,phi3){
  
  m = L_second(mu0,mu1,phi1,phi2)$m
  l111_star <- abs(L_second(mu0,mu1,phi1,phi2)$l111_star)  
  l011_star <- abs(L_second(mu0,mu1,phi1,phi2)$l011_star)
  l001_star <- abs(L_second(mu0,mu1,phi1,phi2)$l001_star)
  l101_star <- abs(L_second(mu0,mu1,phi1,phi2)$l101_star)
  l110_star <- abs(L_second(mu0,mu1,phi1,phi2)$l110_star)
  l000_star <- abs(L_second(mu0,mu1,phi1,phi2)$l000_star)
  l010_star <- abs(L_second(mu0,mu1,phi1,phi2)$l010_star)
  l100_star <- abs(L_second(mu0,mu1,phi1,phi2)$l100_star)
  
  p000 <- prob(p_true, 0, 0, 0, phi1,phi2,phi3)$p
  p010 <- prob(p_true, 0, 1, 0, phi1,phi2,phi3)$p
  p100 <- prob(p_true, 1, 0, 0, phi1,phi2,phi3)$p
  p110 <- prob(p_true, 1, 1, 0, phi1,phi2,phi3)$p
  p001 <- prob(p_true, 0, 0, 1, phi1,phi2,phi3)$p
  p011 <- prob(p_true, 0, 1, 1, phi1,phi2,phi3)$p
  p101 <- prob(p_true, 1, 0, 1, phi1,phi2,phi3)$p
  p111 <- prob(p_true, 1, 1, 1, phi1,phi2,phi3)$p
  
  H <- m * h
  
  t <- rep(NA, (4*H))
  t <- as.matrix(t)
  
  for(i in 1:(H)){
    
    t[(4*i-3):(4*i),1] <- as.matrix(rep(i, 4), 4, 1)
    
  }
  
  Q <- matrix(rep(0, (4*H * 4*H)), nrow = 4*H, ncol=4*H)
  
  for(i in 1:(4*H)){
    
    #l000
    if((4*t[i]-3) > 5|(4*t[i]-3) == 5){
      Q[(4*t[i]-3),((4*t[i]-3) - 4)] = p000
    }else{Q[(4*t[i]-3),1] = p000}
    
    #l001
    if((4*t[i]-3) < (4*H-4*m*l001_star-3)|(4*t[i]-3) == (4*H-4*m*l001_star-3)){
      Q[(4*t[i]-3),((4*t[i]-3) + 4*m*l001_star + 1)] = p001
    }
    
    #l010
    if((4*t[i]-2) < (4*H - 4*m*l010_star - 2)|(4*t[i]-2) == (4*H - 4*m*l010_star - 2)){
      Q[(4*t[i]-2),((4*t[i]-2) + 4*m*l010_star + 1)] = p010
    }
    
    #l011
    if((4*t[i]-2) > (2 + 4*m*l011_star)|(4*t[i]-2) == (2 + 4*m*l011_star)){
      Q[(4*t[i]-2),((4*t[i]-2) - 4*m*l011_star + 2)] = p011
    }else{Q[(4*t[i]-2),1]=p011} 
    
    #l100
    if((4*t[i]-1) < (4*H - 4*m*l100_star -1)|(4*t[i]-1) == (4*H - 4*m*l100_star -1)){
      Q[(4*t[i]-1),((4*t[i]-1) + 4*m*l100_star - 2)] = p100
    }
    
    #l101
    if((4*t[i]-1) > (4*m*l101_star + 3)|(4*t[i]-1) == (4*m*l101_star + 3)){
      Q[(4*t[i]-1),((4*t[i]-1) - 4*m*l101_star - 1)] = p101
    } else{Q[(4*t[i]-1),1]=p101} 
    
    #l110
    if((4*t[i]) > (4 + 4*m*l110_star)|(4*t[i]) == (4 + 4*m*l110_star)){
      Q[(4*t[i]),((4*t[i]) - 4*m*l110_star - 1)] = p110
    }else{Q[(4*t[i]),1] = p110}
    
    #l111
    if((4*t[i]) < (4*H - 4*m*l111_star )|(4*t[i]) == (4*H - 4*m*l111_star)){
      Q[(4*t[i]),((4*t[i]) + 4*m*l111_star )] = p111
    } 
    
  }
  
  I <- diag(4*H)
  
  I1 <- matrix(rep(1, 4*H), nrow=4*H, ncol=1)
  
  N = solve(I - Q) %*% I1
  
  return(list(N=N, m=m,l111_star = l111_star, l011_star=l011_star, l001_star=l001_star, l101_star=l101_star, 
              l110_star=l110_star, l000_star=l000_star, l010_star=l010_star, l100_star=l100_star,p000=p000, p001=p001,
              p010=p010, p011=p011, p101=p101, p110=p110, p111=p111, p100=p100))
  
}

#Generating AR2---------------------------------------------------------------

#>+++++++++++++++++++++++++++++++++++++++++++
RL_CUSUM_rho_second <- function(mu0, mu1, h, X, phi1, phi2){
  
  m <- L_second(mu0, mu1, phi1, phi2)$m
  l111_star_real <- (L_second(mu0,mu1,phi1,phi2)$l111_star)  
  l011_star_real <- (L_second(mu0,mu1,phi1,phi2)$l011_star)
  l001_star_real <- (L_second(mu0,mu1,phi1,phi2)$l001_star)
  l101_star_real <- (L_second(mu0,mu1,phi1,phi2)$l101_star)
  l110_star_real <- (L_second(mu0,mu1,phi1,phi2)$l110_star)
  l000_star_real <- (L_second(mu0,mu1,phi1,phi2)$l000_star)
  l010_star_real <- (L_second(mu0,mu1,phi1,phi2)$l010_star)
  l100_star_real <- (L_second(mu0,mu1,phi1,phi2)$l100_star)
  
  RL <- 3
  
  S.old <- 0
  
  n1 <- length(X)
  
  repeat{
    
    if(X[RL-2]==1&X[RL-1]==1&X[RL]==1){LP=l111_star_real}else if(X[RL-2]==0&X[RL-1]==1&X[RL]==1){LP=l011_star_real}else if(X[RL-2]==0&X[RL-1]==0&X[RL]==1){LP=l001_star_real}else if(X[RL-2]==1&X[RL-1]==0&X[RL]==1){LP=l101_star_real}else if(X[RL-2]==1&X[RL-1]==1&X[RL]==0){LP=l110_star_real}else if(X[RL-2]==0&X[RL-1]==0&X[RL]==0){LP=l000_star_real}else if(X[RL-2]==0&X[RL-1]==1&X[RL]==0){LP=l010_star_real}else if(X[RL-2]==1&X[RL-1]==0&X[RL]==0){LP=l100_star_real}
    
    S.new <- max(0, S.old  + LP)
    
    S.old <- S.new
    
    RL <- RL + 1
    
    if(S.new > h | RL > n1){break}
    
  }
  
  RL <- RL-1
  
  return(list(RL=RL,S.new = S.new))
}

#Monitoring CUSUM++++++++++++++++++++++++++++++++++++++++++++++++++++++++

M_CUSUMrho_AR2_second_OUT <- function(n1, mu0, mu1, p_true, h, phi1, phi2, phi3, n){
  
  RL <- rep(NA, n)
  
  for(i in 1:n){
    
    X <-  Basim(phi1, phi2, phi3, p_true, n1)$Y
    
    RL[i] <- RL_CUSUM_rho_second(mu0, mu1, h, X, phi1, phi2)$RL
    
  }
  return(list(ARL = mean(RL), sd = sd(RL)/sqrt(n)))
}




mu0=0.050
mu1=0.075

phi1=0.4132
phi2=0.2454
phi3=0.8095
p_true = 0.4
m=L_second(mu0, mu1, phi1, phi2 )$m
h=260/m
N=N_vector(h, mu0, mu1,  p_true,  phi1 ,  phi2, phi3 )$N
head(N)
M_CUSUMrho_AR2_second_OUT(3000, mu0, mu1, p_true, h, phi1, phi2, phi3, 10000)

