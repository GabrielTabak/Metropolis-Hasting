# Vamos fazer um algoritmo de Metropolis Hasting 
# para estimar o traço latente dos alunos
# usando TRI, aqui vamos assumir que os 
# parâmetros de itens são conhecidos


# Bibliotecas usadas
library(mirt)
library(matrixStats)


# CCI -> TRI
Pfunction <- function(X){
  return(c + (1 - c)/(1 + exp(-a*(X - b))))
}


#Metropolis Hasting, com burning e lag
# P -> número de pessoas que fizeram a prova
# I -> número de itens da prova
# a, b, c -> parâmetros dos itens (vetores em que a posição X 
#            representa o X-ésimo item)
# It -> número de iterações
# U -> Matriz com a resposta dos alunos(linhas são os alunos e
#                                       colunas os itens) 
# Br -> Burning
# lag -> thinning


MH <- function(P,I,a,b,c, It, U, Br = 1, lag = 1){

  U <- t(U) 
  
  # Matriz para armazenar as cadeias, colunas são as iterações, linhas
  # são os alunos
  thetas2 <- matrix(0,nrow = P,ncol = It)
  acc <- rep(0,P) #Aceitação de cada linha 

  #Valores iniciais (normal 0,1)
  thetas2[,1] <- rnorm(P)

  
  #Gerar as cadeias
  for (k in 2:It){
    #Geramos os candidatos
    candidate <- as.array(rnorm(P,thetas2[,k-1],1))
    
    P1 <- apply(X = candidate,MARGIN = 1, FUN = Pfunction)
    Q1 <- 1 - P1
    
    #Pegamos os anteriores
    atual <- as.array(thetas2[,k-1])
    
    P2 <- apply(X = atual,MARGIN = 1, FUN = Pfunction)
    Q2 <- 1 - P2
    
    #Calculamos R, e se o candidato é aceito ou não
    Ac <- colProds(as.matrix(P1^U*Q1^(1-U)/(P2^U*Q2^(1-U))))
    Ac <- Ac*(dnorm(candidate,0,1)/dnorm(atual,0,1))*(dnorm(atual,candidate,1)/dnorm(candidate,atual,1))
    Ac <- ifelse(runif(P) < Ac, 1, 0)
    
    thetas2[,k] <- ifelse(Ac==1,candidate,thetas2[,k-1])
    acc <- acc + Ac
  }

  
  #Vamos retornar também o mesmo modelo usando o MIRT
  U <- t(U)
  model <- mirt(as.data.frame(U),model = 1, itemtype = "3PL")
  
  #Burning
  cadeiaFinal <- thetas2[,Br:It]
  
  #Thinning
  randomX <- seq(1,length(cadeiaFinal[1,]),lag)
  cadeiaFinal <- cadeiaFinal[,randomX]

  lista <- list("thetas" = cadeiaFinal, "model" = model, "acc" = acc)
  
  return(lista)        
}




################
#Gerando dados simulados, antes de ir para o MCMC
#Vamos estimar apenas os thetas

set.seed(3110)
P <- 400 #alunos
I <- 100   #Items

a <- runif(I,0.8,2)
b <- runif(I,-2,2)
c <- rep(0.1, I)

thetas <- rnorm(P)

  
  #Matriz de Probabilidades (todos os parâmetros conhecidos)
  Probs <- matrix(0,P,I)
  for (i in 1:P) {
    Probs[i,] <- c  + (1-c)/(1 + exp(-a*(thetas[i] - b)))  
  }
  
  #Respostas dos alunos simuladas, geramos a matriz U
  U <- ifelse(runif(I*P)<Probs,1,0)
  U <- as.matrix(U) 

It <- 5000
thetas2 <- MH(P,I,a,b,c,It,U)


#################


#Vamos analisar todas as combinações de traços latentes
#para valores determinados de a,b,c

a <- rep(1.8,5)
b <- c(-3, -1.5, 0, 1.5, 3)
c <- rep(0.2,5)

# c <- rep(0,5)
# a <- rep(1,5)

It <- 10000 
P <- 32
I <- 5

U <- expand.grid(0:1, 0:1, 0:1, 0:1, 0:1)

U <- as.matrix(U)

thetas2 <- MH(P,I,a,b,c,It,U, Br = 500, lag = 4)


#########################


#Analisar convergência

library(ggplot2)


ggplot() + geom_point(mapping = aes(1:P,thetas2$acc/It)) + 
  xlab("Alunos") + ylab("Taxa de Aceitação")

# Nesse caso, ficou muito alto o nível de aceitação, poderiamos aumentar 
# a variância da Normal para os candidatos

#Grafico da cadeia
ts.plot(thetas2$thetas[4,])


plot(cumsum(thetas2$thetas[4,]) / (1:length(thetas2$thetas[4,])), type = "l", xlab = "Iteration",
     ylab = "Ergodic mean of y")

acf(thetas2$thetas[5,], lag.max = 10, main = "x", xlab = "Lag",
    ylab = "ACF")



#Burning e Thinning fora da função
Br <- 1000
lag <- 10

cadeiaFinal <- thetas2$thetas[,Br:5000]

randomX <- seq(1,length(cadeiaFinal[1,]),lag)
cadeiaFinal <- cadeiaFinal[,randomX]

# Vamos comparar os resultados, os thetas reais (usados na simulação),
# com o resultado da cadeia(usando a média) e com o resultado do MIRT
# Podemos usar rowMedians() também

resultado <- data.frame("cadeia" = rowMeans(cadeiaFinal),
                        "cadeiaMediana" = rowMedians(cadeiaFinal),
                        "reais" = thetas, 
                        "fscores" = fscores(thetas2$model),
                        "F - C" = abs(fscores(thetas2$model) - rowMeans(cadeiaFinal)),
                        "Real - cadeia" = abs(thetas - rowMeans(cadeiaFinal)))


library(ggplot2)

ggplot(resultado) + geom_line(aes(x = 1:length(cadeia), y = cadeia),col=3)+
  geom_line(aes(x = 1:length(cadeia), y = reais),col=2)

ggplot(resultado) + geom_line(aes(x = reais, y = cadeia),col=3)

ggplot(resultado) + geom_line(aes(x = cadeiaMediana, y = cadeia),col=3)

ggplot(resultado) + geom_line(aes(x = F1, y = cadeia),col=3)

ggplot(resultado) + geom_line(aes(x = F1, y = reais),col=3)

ggplot(resultado) + geom_line(aes(x = 1:length(cadeia), y = Real...cadeia),col=3)

ggplot(resultado) + geom_line(aes(x = 1:length(cadeia), y = F1.1),col=3)

ggplot(resultado) + geom_density(aes(x = cadeia))


View(rowMeans(thetas2$thetas))

Fim <- cbind(U, rowMeans(thetas2$thetas))




