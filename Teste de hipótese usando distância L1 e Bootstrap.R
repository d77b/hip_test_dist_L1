### Distância L1 ###
rm(list=ls())
############### inicio da funcao que calcula dist L1 ###############
f <- function(x,y,k){
  cx<-(k^(1/5))/(2*sd(x)) # Largura de banda do alisamento
  cy<-(k^(1/5))/(2*sd(y))
  kF <- density(x,cx ,kernel="r", n=k) # kernel density F^
  kG <- density(y,cy, kernel="r", n=k) # kernel density G^
    #plot(kF)
    #plot(kG)
  F <- kF$y # os valores da funçao estimada
  G <- kG$y
  S <- kF$x # pontos onde a densidade F é estimada
  T <- kG$x
  NF <- kF$n # a quant de pontos amostrados; tamanho do array
  NG <- kG$n
  S[1] <- min(S[1], T[1])
  T[1] <- S[1]
  S[NF] <- max(S[NF-1], T[NG-1])
  T[NG] <- S[NF]
  End <- S[NF]
  Fcur <- F[1] # valor corrente do F^ no loop
  Gcur <- G[1]
  InextF <- 2 # índice do prox valor de F^
  InextG <- 2
  Low <- S[1]
  L1 <- 0
  while (Low < End) {
    if(T[InextG] < S[InextF]) {
      High <- T[InextG] #atualiza pontos da abicissa
      L1 <- L1 + (abs(Fcur - Gcur))*(High-Low) #Distância norma 1
      Gcur <- G[InextG] #atualiza pontos da ordenada
      InextG <- InextG + 1
    }
    else{
      High <- S[InextF] #atualiza pontos da abicissa
      L1 <- L1 + (abs(Fcur - Gcur))*(High-Low) #Distância norma 1
      Fcur <- F[InextF] #atualiza pontos da ordenada
      InextF <- InextF + 1
    } # fim do if
    Low <- High
  } # fim do while
  return(L1)
} 
############### fim da funcao que calcula dist L1 ###############

### inicio da simulação ###
repli <- 250
rej <- 0
quantil <- rep(0,c(repli))
for (r in 1:repli){
  n <- 80
  x <- rnorm(80,0,1)
  #y <- rnorm(80,0.25,1)
  y <- rcauchy(80,location = 0, scale = 1)
#  y <- rnorm(80, 0.5, 0.01)
#  y <- runif(1,0,1)
#  if (z <= 0.05) y <- y2
  cont <- 0
  nboot <- 500 # B #
  dist <- rep(0,c(nboot))
  # bootstrap: 
  bsL1x <- matrix(sample(x, size = n * nboot, replace = TRUE), nrow = nboot)
  bsL1y <- matrix(sample(y, size = n * nboot, replace = TRUE), nrow = nboot)
  # calculo das estatisticas do bootstrap:
  for (i in 1:nboot){
    bsL1x[i,] <- sort(bsL1x[i,])
    bsL1y[i,] <- sort(bsL1y[i,])
    dist[i] <- f(bsL1x[i,],bsL1y[i,],n)
  }
  quantil[r] <- summary(dist)[5]
  # teste de hipotese
  # ponto de corte empirico definido apartir do intervalo de confiança 
  # para distro t de 95% para cada rodada de bootstrap.
  
  for (j in 1:nboot){
	if (dist[j] >= 0.2960) cont<-cont+1
	}
  if (cont/nboot >= 0.05) rej <- rej + 1
  print(r)
} 
### fim da simulação ###
hist(quantil)
summary(quantil)


rej_r <- rej/repli
print(rej_r)


#________________________________________
### Teste KS ###

rm(list=ls())
repli <- 250
ac <- 0
pks <- 0
act <- 0
nboot <- 499
sqx <-array(0,c(nboot))
sqy <-array(0,c(nboot))
den <-array(0,c(nboot))
n <- 80
############### Replicação (for principal) ###############
for (i in 1:repli){
  x <- rnorm(80)
  #y <- rnorm(80, .25, 1)
  y <- rcauchy(80,location = 0, scale = 1)
  ### Bootstrap
  bootsamx <- matrix(sample(x, size = n * nboot, replace = TRUE), nrow = nboot)
  bootsamy <- matrix(sample(y, size = n * nboot, replace = TRUE), nrow = nboot)
  ### Estatística t
  xbar <- apply(bootsamx, 1, mean)
  for(j in 1:nboot)
    ybar <- apply(bootsamy, 1, mean)
  for(j in 1:nboot)
  {sqy[j] <- sum((bootsamy[j,]-ybar[j])^2)}
  {sqx[j] <- sum((bootsamx[j,]-xbar[j])^2)}
  t <-array(0,c(nboot)) 
  den <- sqrt((sqx + sqy)/(2*n - 2))*sqrt(2/n)
  t <- (xbar-ybar)/den
  t1 <- t.test(x,y)
  tchap <- t1$s
  qi <- quantile(t, 0.025)
  qs <- quantile(t, 0.975)
  print(i)
  #print(qi)
  #print(qs)
  if (qi < 0 & qs > 0) ac <- ac + 1 
  sdteta <- sd(t) # fazer conta com teta
  vet <- (t-tchap)/sdteta
  qi <- quantile(t, 0.025)
  qs <- quantile(t, 0.975)
  if (qi < 0 & qs > 0) act <- act + 1
  
  ### Kolmogorov - Smirnov
  ks <- rep(0,repli)
  res <- ks.test(x,y) # guarda a lista de resultados
  ks[i]<- res$p.value # guarda os p-valores de cada teste
  if (ks[i]<=0.05) pks<- pks+1
}
############### Fim da Replicação (for principal) ###############

prop <- ac/repli #
rejt <- 1-prop # rejeição da t
print(rejt)
rejteta <- 1 - act/repli
print(rejteta)
rejks <- pks/repli
print(rejks)

rm(list=ls())
