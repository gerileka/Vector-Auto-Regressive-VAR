rm(list = ls())
setwd("C:/Users/geril/OneDrive/Bureau/MoSEF/Ielpo/Projet/Final")
data_entry= read.csv("Data.csv", head= TRUE, sep= ";", dec=",")

# Transformation data type
sapply(data_entry, typeof)
data_entry$X <- as.Date(data_entry$X, format = "%d.%m.%Y")  
indx <- sapply(data_entry, is.factor)
data_entry[indx] <- lapply(data_entry[indx], function(x) as.numeric(as.character(x)))
sapply(data_entry, typeof)
rm(indx)


# Ploting initial data
dev.off()
plot(data_entry$X,data_entry$China, xlab="Date", ylab="GDP Growth",type='l', main="China GDP", col="blue")
plot(data_entry$X,data_entry$USA, xlab="Date", ylab="GDP Growth",type='l', main="USA GDP", col="red")
plot(data_entry$X,data_entry$Eurozone, xlab="Date", ylab="GDP Growth",type='l', main="Eurozone", col="darkgreen")
# We are certenly working with non stationary data so we better take the difference


dates = as.Date(data_entry[, 1], format= "%d.%m.%Y")
data=as.matrix(data_entry[,2:4])
data = data/100
print("Correlation Matrix")
cor(data)


select_lag_Var = function(y,maxLag)
{
  K = ncol(y) # Number of series (in our case 3)
  lag = abs(as.integer(maxLag + 1)) # Numbers of lags (maxlag + 1)
  ylagged = embed(y, lag)[, -c(1:K)] # Transforms the matrix, into a low dim euclidien space (the new table with the lagged columns)
  yendog = y[-c(1:maxLag), ] # y but with the lags applied (meaning that the first dates are deleted)
  sample = nrow(ylagged) # Just the numbers of row
  idx = seq(K, K * maxLag, K) # Where our data is situated (index of column of the real data not lagged)
  
  # Taking the three criterions and creating empty lists
  criteria <- matrix(NA, nrow = 3, ncol = maxLag)
  rownames(criteria) = c("AIC(n)", "HQ(n)", "BIC(n)")
  colnames(criteria) = paste(seq(1:maxLag))
  
  for (i in 1:maxLag) {
    ys.lagged = cbind(ylagged[, c(1:idx[i])], NULL) # Creates the dataframe of lagged plus original data
    sampletot = nrow(y)
    nstar = ncol(ys.lagged)
    resids = lm.fit(x=ys.lagged, y=yendog)$residuals #Takes the residual of every regression for each lag
    sigma.det = det(crossprod(resids)/sample) # Residu/number of rows
    # Bare in mind that K is number columns (3 in our equation)
    criteria[1, i] = log(sigma.det) + (2/sample) * (i * K^2) # AIC
    criteria[2, i] = log(sigma.det) + (2 * log(log(sample))/sample) * (i * K^2) # HQ
    criteria[3, i] = log(sigma.det) + (log(sample)/sample) * (i * K^2) # BIC
    
  }
  order = apply(criteria, 1, which.min)
  final_list = (list(selection = order, criteria = criteria))
  return(order)
  
}


best_lag_list = select_lag_Var(data,8)
best_lag_list
lag_ch = min(best_lag_list)


mvnorm<-function(X,mu,sigma)
{
  # We calculate the density of a normal distribution 
  A=(2*pi)^(ncol(sigma)/2)
  B=det(sigma)^(1/2)
  C=-1/2*t(X-mu)%*%solve(sigma)%*%(X-mu)
  D=exp(C)
  return(1/(A*B)*D)
}


VAR_loglik<-function(para,X)
{
  phi_0=para[1:ncol(X)] # it will get the first paramaters of phi zero
  list_matrix <- list() # Empty list that will stock the other matrix of parameters. 
  E = matrix(phi_0,nrow(X)-lag_ch,3,byrow=TRUE) # Starting point of expectation, this will not change, only updated through interations.
  for (i in 1:lag_ch){
    list_matrix[[i]]=matrix(para[ncol(X)*(1+(i-1)*ncol(X))+1:length(para)],ncol(X),ncol(X)) # We create phi_1 ..... ph_p
    
    E = E+X[1:(nrow(X)-lag_ch),]%*%list_matrix[[i]] # we calculate the expectation (n-p dimension) 
  }
  
  residus=X[1:(nrow(X)-lag_ch),]-E # We get the errors from real X - Expected
  sigma=var(residus) # Compute the simga of the residuals
  
  loglik=0
  for (i in 2:(nrow(data)-(lag_ch+1)))
  {
    temp=mvnorm(data[i,], E[i-1,],sigma)
    temp=log(temp)
    loglik=loglik-temp
    # Calculate the Log-L which should be minimized
  }
  #print(loglik)
  return(loglik)
}



para = numeric(ncol(data)*(1+lag_ch*ncol(data)))
para

VAR_loglik(para,data)


estimation = optim(para,VAR_loglik,,data,method = "BFGS")
estimation



para = estimation$par

get_sigma<-function(para,X)
  ## the new para that we found and we apply the same function as before
{
  phi_0=para[1:ncol(X)]
  list_matrix <- list()
  E = matrix(phi_0,nrow(X)-lag_ch,3,byrow=TRUE) 
  for (i in 1:lag_ch){
    list_matrix[[i]]=matrix(para[ncol(X)*(1+(i-1)*ncol(X))+1:length(para)],ncol(X),ncol(X))
    
    E = E+X[1:(nrow(X)-lag_ch),]%*%list_matrix[[i]] 
  }
  
  residus=X[1:(nrow(X)-lag_ch),]-E
  sigma=var(residus) 
  results <- list(sigma = sigma, phi_1 = list_matrix[[1]])
  return(results)
}



sigma = get_sigma(para,data)$sigma
sigma
phi_1 = get_sigma(para,data)$phi_1
phi_1



P = t(chol(sigma))
P


E = numeric(3) # Create a dim 3 of zeros list
E[1] = -0.08 # Input a shock of negative -0.8
horizon = 10 # We'll comment the first 4 timeframes (1 year) but it looks better shown like this
IRF = c()
for(i in 1:horizon)
{
  #print(i)
  phi = phi_1^i
  #print(phi)
  temp = phi%*%P%*%E
  #print(temp)
  IRF=cbind(IRF,temp)
}
# Construction des trois impact sur les trois series
#dev.off()
for (i in 1:3)
{
  plot(IRF[i,])
}



E = numeric(3)
E[3] = -0.05
horizon = 10
IRF = c()
for(i in 1:horizon)
{
  #print(i)
  phi = phi_1^i
  #print(phi)
  temp = phi%*%P%*%E
  #print(temp)
  IRF=cbind(IRF,temp)
}

# Construction des trois impact sur les trois series
#dev.off()
for (i in 1:3)
{
  plot(IRF[i,])
}



E = numeric(3)
E[2] = -0.05
horizon = 10
IRF = c()
for(i in 1:horizon)
{
  #print(i)
  phi = phi_1^i
  #print(phi)
  temp = phi%*%P%*%E
  #print(temp)
  IRF=cbind(IRF,temp)
}
for (i in 1:3)
{
  plot(IRF[i,])
}

