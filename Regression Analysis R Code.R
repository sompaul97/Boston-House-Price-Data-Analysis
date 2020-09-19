#################################### Loading and Overviewing Data ####################################

data = read.table("Boston_Housing.txt",header=T)
attach(data)
y = MEDV
n = length(y)
X = cbind(CRIM,INDUS,NOX,RM,AGE,DIS,TAX,PTRATIO,B,LSTAT)
X = cbind(rep(1,n),X)
colnames(X)[1] = 'INTERCEPT'
p = ncol(X)
detach(data)

for(i in 2:p)	#Standardizing regressors
	X[,i] = (X[,i]-mean(X[,i])) / sd(X[,i])



#################################### Dividing into Training & Test Data ####################################

set.seed(153)
n.test = 100
n.train = n - n.test
ind.test = sort(sample(1:n,n.test,replace=FALSE))
y.test = y[ind.test]
X.test = X[ind.test,]
y.train = y[-ind.test]
X.train = X[-ind.test,]



#################################### Implementing OLS ####################################

OLS = function(X,y)
{
	beta.hat = solve(crossprod(X)) %*% crossprod(X,y)
	y.hat = X %*% beta.hat
	err = y - y.hat
	RSS = sum(err^2)
	MSRes = RSS / (nrow(X)-ncol(X))
	r = {}
	r$beta.hat = beta.hat
	r$y.hat = y.hat
	r$err = err
	r$RSS = RSS
	r$MSRes = MSRes
	r
}

#################################### Implementing Box-Cox Transformation ####################################

BoxCox = function(x,L)
{
	if(L!=0) u=(x^L-1)/L else u=log(x)
	u
}


#################################### Computing RMSE ####################################

OLS.train = OLS(X.train,y.train)
beta.hat = OLS.train$beta.hat
y.test.hat = X.test %*% beta.hat
RMSE = sqrt(sum((y.test-y.test.hat)^2) / n.test)
message("RMSE before dropping observations: ",round(RMSE,4))



#################################### Dealing with Leverage & Influential Points ####################################

# Leverage Points

#.....................................................................#
leverage.cutoff = 2*p / n.train + 2*0.02

# Note: Cutoff is subject to modification as per requirement!
#.....................................................................#

H = X.train %*% solve(crossprod(X.train)) %*% t(X.train)
plot(1:n.train,diag(H),pch=16,col='blue', xlab='Observation No.',ylab='Leverage Value')
leverage.ind = which(diag(H)>leverage.cutoff)
message("No. of Leverage Points (at Cutoff ",round(leverage.cutoff,4),"): ",length(leverage.ind))
X.train = X.train[-leverage.ind,]
y.train = y.train[-leverage.ind]
n.train = nrow(X.train)
OLS.train = OLS(X.train,y.train)
beta.hat = OLS.train$beta.hat
y.test.hat = X.test %*% beta.hat
RMSE = sqrt(sum((y.test-y.test.hat)^2) / n.test)
message("RMSE after dropping leverage points: ",round(RMSE,4))


# Influential Points (after dropping leverage points)

#.....................................................................#
CooksD.cutoff = 1
DFFITS.cutoff = 2 * sqrt(p/n.train) + 2*0.25
COVRATIO.cutoff = 3*p/n.train + 1*0.1
DFBETAS.cutoff = 2 / sqrt(n.train)

# Note: Cutoffs are subject to modifications as per requirement!
#.....................................................................#

OLS.train = OLS(X.train,y.train)
beta.hat = OLS.train $ beta.hat
MSRes = OLS.train $ MSRes
y.train.hat = OLS.train $ y.hat
C = solve(crossprod(X.train))
H = X.train %*% C %*% t(X.train)

CooksD = 0
DFFITS = 0
COVRATIO = 0
DFBETAS = matrix(0,n.train,p)
for(i in 1:n.train)
{
	OLS.i = OLS(X.train[-i,],y.train[-i])
	beta.hat.i = OLS.i $ beta.hat
	y.i.hat = X.train %*% beta.hat.i
	Si2 = OLS.i $ MSRes
	CooksD[i] = sum((y.i.hat - y.train.hat)^2) / (p*MSRes)
	DFFITS[i] = (y.train.hat[i] - y.i.hat[i]) / sqrt(H[i,i] * Si2)
	COVRATIO[i] = (Si2/MSRes)^p / (1-H[i,i])
	for(j in 1:p)
		DFBETAS[i,j] = (beta.hat[j]-beta.hat.i[j]) / sqrt(Si2*C[j,j])
}

Inf.Analysis_CooksD = data.frame(i=1:n.train,im=CooksD,d=rep(0,n.train))
Inf.Analysis_DFFITS = data.frame(i=1:n.train,im=DFFITS,d=rep(0,n.train))
Inf.Analysis_COVRATIO = data.frame(i=1:n.train,im=COVRATIO,d=rep(0,n.train))
Inf.Analysis_DFBETAS = cbind(1:n.train, matrix(0,n.train,p+1))
for(i in 1:n.train)
{
	count.1 = 0
	for(j in 1:p)
	{
		if (abs(DFBETAS[i,j]) > DFBETAS.cutoff)
		{
			Inf.Analysis_DFBETAS[i,j+1] = 1
			count.1 = count.1 + 1
		}
	}
	if(count.1>0.6*p) Inf.Analysis_DFBETAS[i,p+2] = 1
}

Inf.Analysis_CooksD $ d[which(CooksD>CooksD.cutoff)] = 1
message("No. of Points Suspected to be Influential by Cooks Distance (at Cutoff ",round(CooksD.cutoff,4),"): ",length(which(Inf.Analysis_CooksD$d == 1)))
Inf.Analysis_DFFITS $ d[which(abs(DFFITS)>DFFITS.cutoff)] = 1
message("No. of Points Suspected to be Influential by DFFITS (at Cutoff ",round(DFFITS.cutoff,4),"): ",length(which(Inf.Analysis_DFFITS$d == 1)))
Inf.Analysis_COVRATIO $ d[which(abs(COVRATIO-1)>COVRATIO.cutoff)] = 1
message("No. of Points Suspected to be Influential by COVRATIO (at Cutoffs ",round(1-COVRATIO.cutoff,4)," and ",round(1+COVRATIO.cutoff,4),"): ",length(which(Inf.Analysis_COVRATIO$d == 1)))
message("No. of Points Suspected to be Influential by DFBETAS (at Cutoff ",round(DFBETAS.cutoff,4),"): ",length(which(Inf.Analysis_DFBETAS[,p+2] == 1)))

par(mfrow=c(2,2))
	plot(Inf.Analysis_CooksD$i,Inf.Analysis_CooksD$im, pch=21, bg=c('blue')[unclass(as.factor(Inf.Analysis_CooksD$d))], xlab='Observation No.',ylab='Cooks Distance', xlim=c(1,n.train))
	plot(Inf.Analysis_DFFITS$i,Inf.Analysis_DFFITS$im, pch=21, bg=c('blue','red')[unclass(as.factor(Inf.Analysis_DFFITS$d))], xlab='Observation No.',ylab='DFFITS', xlim=c(1,n.train))
	plot(Inf.Analysis_COVRATIO$i,Inf.Analysis_COVRATIO$im, pch=21, bg=c('blue','red')[unclass(as.factor(Inf.Analysis_COVRATIO$d))], xlab='Observation No.',ylab='COVRATIO', xlim=c(1,n.train))
mtext("Red Colour indicates Influential Point          ",
	side=1,line=-10,adj=1,outer=TRUE,col='red',cex=0.8)
Inf.Analysis = cbind(1:n.train,Inf.Analysis_CooksD$d,Inf.Analysis_DFFITS$d,Inf.Analysis_COVRATIO$d,Inf.Analysis_DFBETAS[,p+2])
count.inf = array(0,n.train)
for(i in 1:n.train)
{
	for(j in 2:5)
	{
		if(Inf.Analysis[i,j]==1)
			count.inf[i] = count.inf[i] + 1
	}
}
influential.ind = which(count.inf>=3 & Inf.Analysis[,4]==1)
X.train.temp = X.train
y.train.temp = y.train
if(length(influential.ind)!=0)
{
	X.train.temp = X.train[-influential.ind,]
	y.train.temp = y.train[-influential.ind]
	n.train.temp = nrow(X.train.temp)
}
message("No. of Influential Points: ",length(influential.ind))

OLS.temp = OLS(X.train.temp,y.train.temp)
beta.hat = OLS.temp$beta.hat
y.test.hat = X.test %*% beta.hat
RMSE.temp = sqrt(sum((y.test-y.test.hat)^2) / n.test)
message("RMSE after dropping influential observations: ",round(RMSE.temp,4))

X.train = X.train.temp
y.train = y.train.temp
n.train = nrow(X.train)
RMSE = RMSE.temp



#################################### Dealing with Curvature ####################################

OLS.train = OLS(X.train,y.train)
beta.hat = OLS.train $ beta.hat
err = OLS.train $ err

par(mfrow=c(2,5))	#Plotting Residuals vs Regressors
for(j in 2:p)
	plot(X.train[,j],err,pch=16,col='blue', xlab=colnames(X.train)[j],ylab='Residual')

par(mfrow=c(4,5))	#APR and CPR Plots
for(j in 2:p)
{
	CPR.j = err + beta.hat[j]*X.train[,j]
	X.star = cbind(X.train,X.train[,j]^2)
	OLS.star = OLS(X.star,y.train)
	beta.star = OLS.star $ beta.hat
	err.star = OLS.star $ err
	APR.j = err + beta.star[j]*X.star[,j] + beta.star[length(beta.star)]*X.star[,ncol(X.star)]
	
	plot(X.train[,j],CPR.j,pch=16,col='blue', xlab=colnames(X.train)[j],ylab='CPR')
	plot(X.train[,j],APR.j,pch=16,col='blue', xlab=colnames(X.train)[j],ylab='APR')
}


# Non-linearity suspected for RM and LSTAT

# Removing Non-Linearity

# For RM
reg.RM = which(colnames(X.train)=='RM')	#identifying the regressor which seems to be non-linearly to E(y)
X_RM = X.train[,-reg.RM]
x.RM = X.train[,reg.RM]
err.y_RM = OLS(X_RM,y.train) $ err
err.RM.X_RM = OLS(X_RM,x.RM) $ err
par(mfrow=c(1,1))
plot(err.y_RM,err.RM.X_RM,pch=16,col='blue', xlab=colnames(X.train)[reg.RM],ylab='Partial Residual')

# For LSTAT
reg.LSTAT = which(colnames(X.train)=='LSTAT')
X_LSTAT = X.train[,-reg.LSTAT]
x.LSTAT = X.train[,reg.LSTAT]
err.y_LSTAT = OLS(X_LSTAT,y.train) $ err
err.LSTAT.X_LSTAT = OLS(X_LSTAT,x.LSTAT) $ err
par(mfrow=c(1,1))
plot(err.y_LSTAT,err.LSTAT.X_LSTAT,pch=16,col='blue', xlab=colnames(X.train)[reg.LSTAT],ylab='Partial Residual')

OLS.LSTAT = OLS(x.LSTAT,log(y.train))	#Transformation used: y = b^x
parameter.nonlinear.LSTAT = as.numeric(exp(OLS.LSTAT $ beta.hat))
LSTAT.new = parameter.nonlinear.LSTAT ^ x.LSTAT
cor.old.LSTAT = cor(y.train,x.LSTAT)
cor.new.LSTAT = cor(y.train,LSTAT.new)
message("Correlation:\tAfter transformation: ",cor.new.LSTAT,"\tBefore transformation: ",cor.old.LSTAT)


X.new = X.train
X.new[,reg.LSTAT] = LSTAT.new
beta.hat.new = OLS(X.new,y.train) $ beta.hat
y.test.hat.new = X.test %*% beta.hat.new
RMSE.new = sqrt(sum((y.test-y.test.hat.new)^2) / n.test)
message("RMSE:\tBefore transformation: ",round(RMSE,4),"\tAfter transformation: ",round(RMSE.new,4))

if(RMSE.new<=RMSE)
{
	X.train = X.new
	RMSE = RMSE.new
}



#################################### Dealing with Heteroscedasticity ####################################

OLS.train = OLS(X.train,y.train)
beta.hat = OLS.train $ beta.hat
err = OLS.train $ err
err2 = err * err
plot(y.train,err,pch=16,col='blue', xlab='MEDV',ylab='Residual')
par(mfrow=c(2,5))
for(j in 2:p)
	plot(X.train[,j],err,pch=16,col='blue', xlab=colnames(X.train)[j],ylab='Residual')
par(mfrow=c(2,5))
for(j in 2:p)
	plot(X.train[,j],err2,pch=16,col='blue', xlab=colnames(X.train)[j],ylab='Estimate of Error Variance (e^2)')


# Breusch-Pagan Test

di.star = n.train*err2 / sum(err2)
ind = which(colnames(X.train)=='CRIM' | colnames(X.train)=='AGE' | 
		colnames(X.train)=='DIS' | colnames(X.train)=='B')
Z = X.train[,c(1,ind)]	#obtained from the plots of ei^2 vs xj
di.star.hat = OLS(Z,di.star) $ y.hat
Q.BP = n.train * cor(di.star,di.star.hat)^2
message("Value of test statistic for Breusch Pagan Test: ",round(Q.BP,4))
if(Q.BP>qchisq(0.95,ncol(Z))) message("Reject Homoscedasticity Assumption") else message("Accept Homoscedasticity Assumption")


# Estimating Regression Coefficients in the presence of Heteroscedasticity

beta.hat.0 = OLS(X.train,y.train) $ beta.hat
beta.hat.1 = beta.hat.0
alpha.hat.0 = OLS(Z,log(err2)) $ beta.hat
alpha.hat.1 = alpha.hat.0
Sigma.hat = diag(0,n.train)
while(sum(beta.hat.1-beta.hat.0)^2 <= 0.1 & sum(alpha.hat.1-alpha.hat.0)^2 <= 0.1)
{
	beta.hat.0 = beta.hat.1
	alpha.hat.0 = alpha.hat.1
	err2.I = (y.train - X.train %*% beta.hat.0)^2
	OLS.I = OLS(Z,log(err2.I))
	alpha.hat.1 = OLS.I $ beta.hat
	Sigma.hat = diag(as.vector(OLS.I $ y.hat))
	beta.hat.1 = solve(t(X.train)%*%solve(Sigma.hat)%*%X.train) %*% t(X.train)%*%solve(Sigma.hat)%*%y.train
}

y.test.hat.new = X.test %*% beta.hat.1
RMSE.new = sqrt(sum((y.test-y.test.hat.new)^2) / n.test)
message("RMSE:\tBased on New Estimate: ",round(RMSE.new,4),"\tBased on Old Estimate: ",round(RMSE,4))

if(RMSE>RMSE.new)
{
	RMSE = RMSE.new
	X.train = sqrt(solve(Sigma.hat)) %*% X.train
	y.train = sqrt(solve(Sigma.hat)) %*% y.train
}


#################################### Dealing with Non-Normality ####################################

# Q-Q Plot (before transformation)

H = X.train %*% solve(crossprod(X.train)) %*% t(X.train)
OLS.train = OLS(X.train,y.train)
err = OLS.train $ err
R.student = array(0)
for(i in 1:n.train)
{
	Si2 = OLS(X.train[-i,],y.train[-i]) $ MSRes
	R.student[i] = err[i] / sqrt(Si2*(1-H[i,i]))
}
par(mfrow=c(1,1))
plot(qt(seq(0.01,0.99,0.01),n.train-p),quantile(R.student,seq(0.01,0.99,0.01)),
	xlab='Population Quantile',ylab='Sample Quantile',main='Q-Q Plot')
abline(a=0,b=1,col='red')


# Box-Cox Transformation & Profile Likelihood Method

lambda.val = seq(-2,2,0.01)
g.lambda = array(0)
for(i in 1:length(lambda.val))
{
	y.lambda = BoxCox(y.train,lambda.val[i])
	OLS.lambda = OLS(X.train,y.lambda)
	beta.hat.lambda = OLS.lambda $ beta.hat
	RSS.lambda = OLS.lambda $ RSS
	sigma2.lambda = RSS.lambda / n.train
	g.lambda[i] = -0.5*n.train*(log(2*pi*sigma2.lambda) + 1) + n.train*(lambda.val[i]-1)*mean(log(y.train))
}
plot(lambda.val,g.lambda,type='l',col='red', xaxt='n',xlab=expression(lambda),ylab=expression(paste('g(',lambda,')')))
axis(side=1,at=seq(-2,2,0.1),tck=1)

lambda = 0.2	#obtained from the plot


# Q-Q Plot (after transformation)

y.BoxCox = BoxCox(y.train,lambda)
OLS.BoxCox = OLS(X.train,y.BoxCox)
err = OLS.BoxCox $ err
R.student = array(0)
for(i in 1:n.train)
{
	Si2 = OLS(X.train[-i,],y.BoxCox[-i]) $ MSRes
	R.student[i] = err[i] / sqrt(Si2*(1-H[i,i]))
}
plot(qt(seq(0.01,0.99,0.01),n.train-p),quantile(R.student,seq(0.01,0.99,0.01)),
	xlab='Population Quantile',ylab='Sample Quantile',main='Q-Q Plot')
abline(a=0,b=1,col='red')

#######################################################################################################