

#Flexible Quantile Contours for Multivariate Functional Data: Beyond Convexity#


library(splines)
library(MASS)
library(matrixcalc) ## for vec
library(Matrix) ## for kronecker product
library(gurobi)
library(mvtnorm)
x.points <- seq(1,9,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
mean <- c(5,5)
sigma <- matrix(c(0.1,0,0,0.1),nrow=2)


### regression coeffcients
## since we have generated b-spline matrix we provide gammas instead of betas
gamma1 = c(1, 0.7, 0.8, 0.9, 0.7, 0.8, 0.9)
gamma2 = c(1, 0.8, 0.7, 0.8, 0.9, 0.7, 0.6)


##### flexible contours for bivriate normal distribution
##############################################

t_length=10
t =seq(0.1,1,length.out = t_length)

library(splines)
## intercept is not included by default.
# Number of knots = df-degree(3)-1(if intercept)
basis <- bs(x = t, df=7, degree = 3,
            Boundary.knots = c(0,1), intercept = TRUE)
gamma1 = c(1, 0.7, 0.8, 0.9, 0.7, 0.8, 0.9)
gamma2 = c(1, 0.8, 0.7, 0.8, 0.9, 0.7, 0.6)



n_each = 300
Y1=matrix(0, nrow=length(t), ncol = n_each); mean1=matrix(0, nrow=length(t), ncol = n_each)
Y2=matrix(0, nrow=length(t), ncol = n_each); mean2=matrix(0, nrow=length(t), ncol = n_each)
for(i in 1:length(t)){
  
  biv.normal = rmvnorm(n=n_each, mean = mean, sigma)
  
  Y1[i, ] = t(basis[i,])%*%gamma1+ biv.normal[,1]
  Y2[i, ] = t(basis[i,])%*%gamma2+ biv.normal[,2]
  
  mean1[i,] = t(basis[i,])%*%gamma1
  mean2[i,] = t(basis[i,])%*%gamma2
  
}


library(scatterplot3d)
##plotting the simulated launches
s=scatterplot3d(Y1[,1],Y2[,1],t,color="grey60", xlim = c(5.5,6.5), ylim = c(5.5,6),
                xlab = "x", ylab = "y",zlab = "t",angle=10, cex.symbols = 0.2)

for(i in 2:n_each){
  s$points3d(Y1[,i],Y2[,i],t,col="grey60", cex=0.2)
}
## observed spline function
s$points3d(mean1[,1]+5,mean2[,1]+5,t,col="red", cex=0.2, type = "l")


##############################################
grid_size=21
m = grid_size^2

Y1_alltime = vec(t(Y1)) 
Y2_alltime = vec(t(Y2)) 
Y= cbind(Y1_alltime, Y2_alltime)



#plot(Y1[5,], Y2[5,])
#points(Y_cap, col="red")
#min(Y_cap)
n = n_each*length(t)

df = 7
basis_model = bs(x = t, df= df, degree = 3,
                 Boundary.knots = c(0,1), intercept = TRUE)

vec_n= vector()
for(i in 1:t_length){
  vec_i = c(rep(basis_model[i,], times= n_each))
  vec_n = c(vec_n, vec_i)
}


X =  matrix(vec_n, nrow = n, ncol =df, byrow=T)

#X = rep(1,n)  ## without slope
nu = rep(1/n, n)
mu = rep(1/m, m)
# simulating U from U[0,1]^2
#U = cbind(runif(m, 0,1), runif(m,0,1))


## simulating U from a regular grid
x <- seq(0, 1, length.out = grid_size)
y <- seq(0, 1, length.out = grid_size)
U <- as.matrix(expand.grid(x = x, y = y))




######################## DUAL ############################
obj1 = vec(as.matrix(nu))
obj2 = vec(mu%*%t(nu)%*%X)

model_dual = list()

##objective function
model_dual$obj = c(obj1, obj2)

## A matrix of constraints (lhs)
A1_dual = Matrix(kronecker(X= diag(n), Y= rep(1,m)), sparse = T)
A2_dual = Matrix(kronecker(X= X , Y= diag(m)), sparse = T)


model_dual$A = cbind(A1_dual,A2_dual)

## b vector of constraints (rhs)
rhs_dual = vec(U%*%t(Y))
model_dual$rhs <- rhs_dual

model_dual$vtype <- 'C'  ## Continuous variable type (default)
model_dual$sense <- '>'
model_dual$modelsense <- 'min'

params <- list(OutputFlag=0)
result_dual <- gurobi(model_dual, params)


print(result_dual$objval)
#print(result_dual$x)
psi = result_dual$x[1:n]
vec_b = result_dual$x[-(1:n)]
b = matrix(vec_b, nrow = m, ncol = df) ## b constructed according to x

## obtaining the beta's using numererical differentiation
library(numDeriv)

#covariate = c(1, t[15], t[15]^2, t[15]^3)
## predicting at  t = t[15]
covariate = predict(basis_model, newx = t[5])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)

########## function for numeical differentiation #######
finite.differences <- function(x, y) {
  
  stepsize = abs(x[1,1]-x[2,1])
  n <- dim(x)[1]
  
  # Initialize a vector of length n to enter the derivative approximations
  fdx_u <- vector(length = grid_size)
  fdx_f = rep(0,grid_size); fdx_b = rep(0,grid_size)
  # Iterate through the values using the forward differencing method
  fdx = matrix(0, nrow=grid_size, ncol=grid_size)
  
  for(u in 1:grid_size){
    for (i in 1:(grid_size-1)) {
      fdx_f[i] <- (y[i+1 +grid_size*(u-1)] - y[i +grid_size*(u-1)]) / stepsize
    }
    
    # Iterate through the values using the backward differencing method
    for (i in 2:grid_size) {
      fdx_b[i] <- (y[i +grid_size*(u-1)] - y[i-1 +grid_size*(u-1)]) / stepsize
    }
    
    fdx_u <- (fdx_f+fdx_b)/2
    fdx_u[1] <- fdx_f[1]
    fdx_u[grid_size] <- fdx_b[grid_size]
    
    fdx[,u] = fdx_u
  }
  
  return(fdx)
}

beta1Tx = finite.differences(U, bTx)
fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))


# Here U=(u2,u1) so we take the transpose
beta2Tx = t(finite.differences(U, b2Tx))
#par(mfrow=c(1,2))
#plot(Y1[15,], Y2[15,], xlab="Y1_cap", ylab="Y2_cap",type="n")
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

norm_vec <- function(x) sqrt(sum(x^2))
norm_max = function(x) max(abs(x[1]), abs(x[2]))

#plot(U_centered)
#points(U_centered[which(norm_rec<0.125),], col="red")
U_centered = U-0.5
#norm_U_c = apply(U_centered, 1, norm_vec)
#summary(norm_U_c)
#norm = norm_U_c/norm_U_c[m]
tau = seq(0.1,1,length.out=7)


## transform uniform reference to spherical uniform
radius = seq(0.00001, 1, length.out = grid_size)

num <- grid_size # number of points you want on the unit circle

pts_all <- t(sapply(1:num,function(p) c(radius[1]*cos(2*p*pi/num),radius[1]*sin(2*p*pi/num))))
for(i in 2:grid_size){
  pts.circle <- t(sapply(1:num,function(p) c(radius[i]*cos(2*p*pi/num),radius[i]*sin(2*p*pi/num))))
  pts_all = rbind(pts_all, pts.circle)
}

U_n = pts_all

#plot(U_n)
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

library(pdist)
## U_n as rows, Y_n as columns
dist_mat = as.matrix(pdist(X=U_n, Y=U_centered))
ranks = rank(c(dist_mat))

##defining cost matrix
cmat = matrix(ranks, nrow=grid_size^2, ncol=grid_size^2)


library(adagio)
matching = assignment(cmat)

perm = matching$perm

#tau = 0.75
#U_number = which(apply(U_n, 1, norm_vec) < tau)
#U_centered_corresponding = perm[U_number]

library(alphahull)
plot(Y1[5,], Y2[5,], xlab=expression(paste(Y[1])), ylab=expression(paste(Y[2])), pch=20, cex=0.5)
#points(Y_cap[U_centered_corresponding,1],Y_cap[U_centered_corresponding,2], xlab="Y1_cap", ylab= "Y2_cap", col="red", pch=20, cex=0.5)

tau=seq(0.1,0.9, length.out = 6)
tau = c(0.25,0.5,0.75)
for (i in 1:length(tau)) {
  
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 800)
  plot(ahull.obj,lwd=c(1,1,2), add=T, wpoints=F)
  
  
}


##### Overylay true bivariate normal density contours 
mean_t = vector()
mean_t[1] = mean[1] + t(basis[5,])%*%gamma1
mean_t[2] = mean[2] + t(basis[5,])%*%gamma2
sn = function(y1, y2){
  sndensity = mvtnorm::dmvnorm(c(y1, y2), mean_t, sigma)
  return(sndensity)}


# making the data-frame at which the contour-plot should be made
z1 = seq(1, 10, length.out = 100)
z2 = seq(1, 10, length.out = 100)
n = length(z1)
z3 = numeric(0)
for(i in 1:n){
  for(j in 1:n){
    z = sn(z1[i], z2[j])
    z3 = append(z3, z)
  }
}
z3 = matrix(z3, nrow = n, byrow=T)


# contour-plot
#require(grDevices)
#graphics::filled.contour(z1, z2, z3, nlevels = 20, las = 1, color.palette = rainbow)
#contour(z1, z2, z3, las = 1, nlevels = 10, drawlabels = F, method = "flattest", vfont = c("sans serif", "bold"), col = "red", xlim = c(-5, 5), ylim = c(-5, 5), lwd = 1.2, las = 1)

tmpz <- sort(as.vector(z3), decreasing = TRUE)
tmpz2 <- cumsum(tmpz)

color_density = c("red","blue", "green")
for (i in 1:length(tau)) {
  w <- which(tmpz2 > tau[i]*tail(tmpz2,1))[1]
  contour(z1,z2,z3, levels = tmpz[w],
          drawlabels = FALSE, col = color_density[i], xlab = "", ylab = "", main = "", add=T)
}




## connect contours over time i.e. plot estimated spline function for each u
Y1_cap = matrix(0, nrow = m, ncol = length(t)); Y2_cap = matrix(0, nrow = m, ncol = length(t))
for (i in 1:length(t)) {
  covariate = predict(basis_model, newx = t[i])[1,]
  bTx= apply(b, 1, function(x) t(x)%*%covariate)
  beta1Tx = finite.differences(U, bTx)
  
  fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
  b2Tx = vec(t(fmatrix))
  beta2Tx = t(finite.differences(U, b2Tx))
  
  Y1_cap[,i] = c(beta1Tx)
  Y2_cap[,i] = c(beta2Tx)
  
}

### plotting the estimated spline function
library(scatterplot3d)
s=scatterplot3d(Y1_cap[1,],Y2_cap[1,],t,type="l" , xlim=c(5,6.5), ylim = c(5,7),
                cex.symbols =0.4, pch = 20,
                xlab = expression(paste(hat(Y)[1])),ylab = "", zlab = "t",
                color ="grey60")

for(i in 2:dim(U)[1]){
  s$points3d(Y1_cap[i,],Y2_cap[i,],t, type = "l", cex=0.4,pch=20,col="grey60")
}
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])-1
y <- dims[3]+ 0.08*diff(dims[3:4])-0.3
text(x,y,expression(paste(hat(Y)[2])),srt=10)

## observed spline function
s$points3d(mean1[,1]+5,mean2[,1]+5,t,col="red", cex=0.2, type = "l")

##median
#s$points3d(Y1_cap[145,],Y2_cap[145,],t, type = "l", cex=0.4,pch=20,col="blue")

legend("topright", legend=c("Predicted", "True"),
       col=c("black", "red"), lty=1, cex=0.8)




#### Repeating Simulations ####

set.seed(151)

library(matrixcalc) ## for vec
library(Matrix) ## for kronecker product
library(gurobi)

true_mean_ls = list()
predicted_ls = list()
euc_dis_curve = vector()
Y_cap_list =list()
map_list = list()
simulationsize = 10
for (k in 1:simulationsize) {
  
  
  
  t_length=10
  t =seq(0.1,1,length.out = t_length)
  library(splines)
  ## intercept is not included by default.
  # Number of knots = df-degree(3)-1(if intercept)
  basis <- bs(x = t, df=7, degree = 3,
              Boundary.knots = c(0,1), intercept = TRUE)
  
  
  ### regression coeffcients
  ## since we have generated b-spline matrix we provide gammas instead of betas
  gamma1 = c(1, 0.7, 0.8, 0.9, 0.7, 0.8, 0.9)
  gamma2 = c(1, 0.8, 0.7, 0.8, 0.9, 0.7, 0.6)
  
  
  n_each = 300
  Y1=matrix(0, nrow=length(t), ncol = n_each); mean1=matrix(0, nrow=length(t), ncol = n_each)
  Y2=matrix(0, nrow=length(t), ncol = n_each); mean2=matrix(0, nrow=length(t), ncol = n_each)
  #for(i in 1:length(t)){
  
  #  X = runif(n_each, -0.5, 0.5)
  #  Z = runif(n_each, 0, 0.5)
  #  phi = runif(n_each, 0, 2*pi)
  #  R = 0.2*Z*(0.5+(0.5-abs(X))/2)
  
  #  epsilon1_t = X+ R*cos(phi)
  #  epsilon2_t = X^2 + R*sin(phi)
  
  #  mean1[i,] = t(basis[i,])%*%gamma1
  #  mean2[i,] = t(basis[i,])%*%gamma2
  
  #  Y1[i, ] = t(basis[i,])%*%gamma1 + epsilon1_t 
  #  Y2[i, ] = t(basis[i,])%*%gamma2 + epsilon2_t 
  
  
  #}
  
  mean <- c(5,5)
  sigma <- matrix(c(0.1,0,0,0.1),nrow=2)
  for(i in 1:length(t)){
    
    biv.normal = rmvnorm(n=n_each, mean = mean, sigma)
    
    Y1[i, ] = t(basis[i,])%*%gamma1+ biv.normal[,1]
    Y2[i, ] = t(basis[i,])%*%gamma2+ biv.normal[,2]
    
    mean1[i,] = t(basis[i,])%*%gamma1
    mean2[i,] = t(basis[i,])%*%gamma2
    
  }
  
  ##############################################
  grid_size=17
  m = grid_size^2
  
  Y1_alltime = vec(t(Y1))
  Y2_alltime = vec(t(Y2))
  Y= cbind(Y1_alltime, Y2_alltime)
  
  #plot(Y)
  
  
  
  n = n_each*length(t)
  
  df = 7
  basis_model = bs(x = t, df= df, degree = 3,
                   Boundary.knots = c(0,1), intercept = TRUE)
  
  vec_n= vector()
  for(i in 1:t_length){
    vec_i = c(rep(basis_model[i,], times= n_each))
    vec_n = c(vec_n, vec_i)
  }
  
  
  X =  matrix(vec_n, nrow = n, ncol =df, byrow=T)
  
  #X = rep(1,n)  ## without slope
  nu = rep(1/n, n)
  mu = rep(1/m, m)
  # simulating U from U[0,1]^2
  #U = cbind(runif(m, 0,1), runif(m,0,1))
  
  
  ## simulating U from a regular grid
  x <- seq(0, 1, length.out = grid_size)
  y <- seq(0, 1, length.out = grid_size)
  U <- as.matrix(expand.grid(x = x, y = y))
  
  
  ######################## DUAL ############################
  obj1 = vec(as.matrix(nu))
  obj2 = vec(mu%*%t(nu)%*%X)
  
  model_dual = list()
  
  ##objective function
  model_dual$obj = c(obj1, obj2)
  
  ## A matrix of constraints (lhs)
  A1_dual = Matrix(kronecker(X= diag(n), Y= rep(1,m)), sparse = T)
  A2_dual = Matrix(kronecker(X= X , Y= diag(m)), sparse = T)
  
  
  model_dual$A = cbind(A1_dual,A2_dual)
  
  ## b vector of constraints (rhs)
  rhs_dual = vec(U%*%t(Y))
  model_dual$rhs <- rhs_dual
  
  model_dual$vtype <- 'C'  ## Continuous variable type (default)
  model_dual$sense <- '>'
  model_dual$modelsense <- 'min'
  
  params <- list(OutputFlag=0)
  result_dual <- gurobi(model_dual, params)
  
  
  print(result_dual$objval)
  #print(result_dual$x)
  psi = result_dual$x[1:n]
  vec_b = result_dual$x[-(1:n)]
  b = matrix(vec_b, nrow = m, ncol = df) ## b constructed according to x
  
  ## obtaining the beta's using numererical differentiation
  library(numDeriv)
  
  #covariate = c(1, t[15], t[15]^2, t[15]^3)
  ## predicting at  t = t[15]
  covariate = predict(basis_model, newx = t[5])[1,]
  bTx= apply(b, 1, function(x) t(x)%*%covariate)
  beta1Tx = finite.differences(U, bTx)
  fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
  b2Tx = vec(t(fmatrix))
  
  # Here U=(u2,u1) so we take the transpose
  beta2Tx = t(finite.differences(U, b2Tx))
  #par(mfrow=c(1,2))
  #plot(Y1[15,], Y2[15,], xlab="Y1_cap", ylab="Y2_cap",type="n")
  Y_cap =cbind(c(beta1Tx), c(beta2Tx))
  
  Y_cap_list[[k]] = Y_cap
  norm_vec <- function(x) sqrt(sum(x^2))
  
  U_centered = U-0.5
  tau= seq(0.1,1,length.out=7)
  
  
  
  ## transform uniform reference to spherical uniform
  radius = seq(0.00001, 1, length.out = grid_size)
  
  num <- grid_size # number of points you want on the unit circle
  
  pts_all <- t(sapply(1:num,function(p) c(radius[1]*cos(2*p*pi/num),radius[1]*sin(2*p*pi/num))))
  for(i in 2:grid_size){
    pts.circle <- t(sapply(1:num,function(p) c(radius[i]*cos(2*p*pi/num),radius[i]*sin(2*p*pi/num))))
    pts_all = rbind(pts_all, pts.circle)
  }
  
  U_n = pts_all
  
  #plot(U_n)
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  
  library(pdist)
  ## U_n as rows, Y_n as columns
  dist_mat = as.matrix(pdist(X=U_n, Y=U_centered))
  ranks = rank(c(dist_mat))
  
  ##defining cost matrix
  cmat = matrix(ranks, nrow=grid_size^2, ncol=grid_size^2)
  
  library(adagio)
  matching = assignment(cmat)
  
  perm = matching$perm
  
  map_list[[k]] = perm 
  
  
  ## connect contours over time i.e. plot estimated spline function for each u
  Y1_cap = matrix(0, nrow = m, ncol = length(t)); Y2_cap = matrix(0, nrow = m, ncol = length(t))
  for (i in 1:length(t)) {
    covariate = predict(basis_model, newx = t[i])[1,]
    bTx= apply(b, 1, function(x) t(x)%*%covariate)
    beta1Tx = finite.differences(U, bTx)
    
    fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
    b2Tx = vec(t(fmatrix))
    beta2Tx = t(finite.differences(U, b2Tx))
    
    Y1_cap[,i] = c(beta1Tx)
    Y2_cap[,i] = c(beta2Tx)
    
  }
  
  
  ## observed spline function
  true_mean = cbind(mean1[,1],mean2[,1])+5
  predicted_ls[[k]] = cbind(Y1_cap[145,],Y2_cap[145,])
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  
  euc=0
  for (j in 1:length(t)) {
    euc= euc + euc.dist(true_mean[j,],predicted_ls[[k]][j,])
  }
  
  euc_dis_curve[k] = euc
  
}



mean(euc_dis_curve)/t_length


library(scatterplot3d)
s=scatterplot3d(mean1[,1]+5,mean1[,1]+5,t,type="l" , xlim=c(5, 6.5), ylim = c(5,7),
                cex.symbols =0.4, pch = 20,
                xlab = expression(paste(hat(Y)[1])),ylab = "", zlab = "t",
                color ="red")

for(i in 2:simulationsize){
  s$points3d(predicted_ls[[i]][,1],predicted_ls[[i]][,2],t, type = "l", cex=0.4,pch=20,col="grey60")
}
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])-0.2
y <- dims[3]+ 0.08*diff(dims[3:4])-0.2
text(x,y,expression(paste(hat(Y)[2])),srt=10)

## observed spline function
s$points3d(mean1[,1],mean2[,1],t,col="red", cex=0.2, type = "l")



plot(Y1[5,],Y2[5,], xlab=expression(paste(Y[1])), ylab=expression(paste(Y[2])), pch=20, cex=0.5, col="grey60")

tau= seq(0.01, 0.9, length.out = 6)
tau = c(0.25,0.5,0.75)
##### averaged contours
for (i in 1:length(tau)) {
  
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  
  for(k in 1:simulationsize){
    U_centered_corresponding = map_list[[k]][U_number]
    Y_tau = Y_cap_list[[k]][U_centered_corresponding,]
    if(k==1){Y_tau_avg=Y_tau} else 
    {Y_tau_avg = rbind(Y_tau_avg, Y_tau)}
  }
  Y_tau_avg = unique(Y_tau_avg)
  ahull.obj = ahull(Y_tau_avg[,1], Y_tau_avg[,2], 1000)
  plot(ahull.obj,lwd=c(1,1,2), add=T, wpoints=F)
  
  
}




############ Flexible quantile contours with non convex margins ########
t_length=10
t =seq(0.1,1,length.out = t_length)
library(splines)
## intercept is not included by default.
# Number of knots = df-degree(3)-1(if intercept)
basis <- bs(x = t, df=7, degree = 3,
            Boundary.knots = c(0,1), intercept = TRUE)


### regression coeffcients
## since we have generated b-spline matrix we provide gammas instead of betas
gamma1 = c(1, 0.7, 0.8, 0.9, 0.7, 0.8, 0.9)
gamma2 = c(1, 0.8, 0.7, 0.8, 0.9, 0.7, 0.8)

#gamma1 = rep(1,5)
#gamma2 = rep(1,5)
n_each = 300
Y1=matrix(0, nrow=length(t), ncol = n_each); mean1=matrix(0, nrow=length(t), ncol = n_each)
Y2=matrix(0, nrow=length(t), ncol = n_each); mean2=matrix(0, nrow=length(t), ncol = n_each)
for(i in 1:length(t)){
  #X = runif(n_each, -1, 1)
  #Z = runif(n_each, 0, 1)
  #phi = runif(n_each, 0, 2*pi)
  #R = 0.2*Z*(1+(1-abs(X))/2)
  
  X = runif(n_each, -0.5, 0.5)
  Z = runif(n_each, 0, 0.5)
  phi = runif(n_each, 0, 2*pi)
  R = 0.2*Z*(0.5+(0.5-abs(X))/2)
  
  epsilon1_t = X+ R*cos(phi)
  epsilon2_t = X^2 + R*sin(phi)
  
  mean1[i,] = t(basis[i,])%*%gamma1
  mean2[i,] = t(basis[i,])%*%gamma2
  
  Y1[i, ] = t(basis[i,])%*%gamma1 + epsilon1_t 
  Y2[i, ] = t(basis[i,])%*%gamma2 + epsilon2_t 
  
  #Y1[i, ] = 1 + epsilon1_t  ## model without slope (just intercept)
  #Y2[i, ] = 1 + epsilon2_t  
  
}

Y1_alltime = vec(t(Y1))
Y2_alltime = vec(t(Y2))
Y = cbind(Y1_alltime, Y2_alltime)

plot(Y)

plot(Y1[5,], Y2[5,])




library(scatterplot3d)
##plotting the simulated launches
s=scatterplot3d(Y1[,1],Y2[,1],t,color="grey60", #xlim = c(-0.5,1.5), ylim = c(0,1.2),
                xlab = "x", ylab = "y",zlab = "t",angle=30, cex.symbols = 0.2)

for(i in 2:n_each){
  s$points3d(Y1[,i],Y2[,i],t,col="grey60", cex=0.2)
}
s$points3d(mean1[,i],mean2[,i],t,col="red", cex=0.2, type = "l")

## Gurobi ##
library(gurobi)

##############################################
grid_size=21
m = grid_size^2

Y1_alltime = vec(t(Y1))
Y2_alltime = vec(t(Y2))
Y = cbind(Y1_alltime, Y2_alltime)

#plot(Y)



n = n_each*length(t)

df = 7
basis_model = bs(x = t, df= df, degree = 3,
                 Boundary.knots = c(0,1), intercept = TRUE)

vec_n= vector()
for(i in 1:t_length){
  vec_i = c(rep(basis_model[i,], times= n_each))
  vec_n = c(vec_n, vec_i)
}


X =  matrix(vec_n, nrow = n, ncol =df, byrow=T)

#X = rep(1,n)  ## without slope
nu = rep(1/n, n)
mu = rep(1/m, m)
# simulating U from U[0,1]^2
#U = cbind(runif(m, 0,1), runif(m,0,1))


## simulating U from a regular grid
x <- seq(0, 1, length.out = grid_size)
y <- seq(0, 1, length.out = grid_size)
U <- as.matrix(expand.grid(x = x, y = y))




######################## DUAL ############################
obj1 = vec(as.matrix(nu))
obj2 = vec(mu%*%t(nu)%*%X)

model_dual = list()

##objective function
model_dual$obj = c(obj1, obj2)

## A matrix of constraints (lhs)
A1_dual = Matrix(kronecker(X= diag(n), Y= rep(1,m)), sparse = T)
A2_dual = Matrix(kronecker(X= X , Y= diag(m)), sparse = T)



#object.size(A2_dual)
#class(A1_dual)


model_dual$A = cbind(A1_dual,A2_dual)

## b vector of constraints (rhs)
rhs_dual = vec(U%*%t(Y))
model_dual$rhs <- rhs_dual

model_dual$vtype <- 'C'  ## Continuous variable type (default)
model_dual$sense <- '>'
model_dual$modelsense <- 'min'

params <- list(OutputFlag=0)
result_dual <- gurobi(model_dual, params)


print(result_dual$objval)
#print(result_dual$x)
psi = result_dual$x[1:n]
vec_b = result_dual$x[-(1:n)]
b = matrix(vec_b, nrow = m, ncol = df) ## b constructed according to x

## obtaining the beta's using numererical differentiation
library(numDeriv)

#covariate = c(1, t[15], t[15]^2, t[15]^3)
## predicting at  t = t[15]
covariate = predict(basis_model, newx = t[5])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
plot(bTx)




finite.differences <- function(x, y) {
  
  stepsize = abs(x[1,1]-x[2,1])
  n <- dim(x)[1]
  
  # Initialize a vector of length n to enter the derivative approximations
  fdx_u <- vector(length = grid_size)
  fdx_f = rep(0,grid_size); fdx_b = rep(0,grid_size)
  # Iterate through the values using the forward differencing method
  fdx = matrix(0, nrow=grid_size, ncol=grid_size)
  
  for(u in 1:grid_size){
    for (i in 1:(grid_size-1)) {
      fdx_f[i] <- (y[i+1 +grid_size*(u-1)] - y[i +grid_size*(u-1)]) / stepsize
    }
    
    # Iterate through the values using the backward differencing method
    for (i in 2:grid_size) {
      fdx_b[i] <- (y[i +grid_size*(u-1)] - y[i-1 +grid_size*(u-1)]) / stepsize
    }
    
    fdx_u <- (fdx_f+fdx_b)/2
    fdx_u[1] <- fdx_f[1]
    fdx_u[grid_size] <- fdx_b[grid_size]
    
    fdx[,u] = fdx_u
  }
  
  return(fdx)
}




beta1Tx = finite.differences(U, bTx)


fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))

# Here U=(u2,u1) so we take the transpose
beta2Tx = t(finite.differences(U, b2Tx))

Y_cap =cbind(c(beta1Tx), c(beta2Tx))


norm_vec <- function(x) sqrt(sum(x^2))

U_centered = U-0.5


## transform uniform reference to spherical uniform
radius = seq(0.00001, 1, length.out = grid_size)

num <- grid_size # number of points you want on the unit circle

pts_all <- t(sapply(1:num,function(p) c(radius[1]*cos(2*p*pi/num),radius[1]*sin(2*p*pi/num))))
for(i in 2:grid_size){
  pts.circle <- t(sapply(1:num,function(p) c(radius[i]*cos(2*p*pi/num),radius[i]*sin(2*p*pi/num))))
  pts_all = rbind(pts_all, pts.circle)
}

U_n = pts_all

#plot(U_n)
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

library(pdist)
## U_n as rows, Y_n as columns
dist_mat = as.matrix(pdist(X=U_n, Y=U))
ranks = rank(c(dist_mat))

##defining cost matrix
cmat = matrix(ranks, nrow=grid_size^2, ncol=grid_size^2)

library(adagio)
matching = assignment(cmat)

perm = matching$perm


plot(Y1[5,],Y2[5,], xlab=expression(paste(Y[1])), ylab=expression(paste(Y[2])), pch=20, cex=0.5, col="grey60")

library(alphahull)
#plot(Y_cap, xlab="Y1_cap", ylab= "Y2_cap")
tau= seq(0.01,1,length.out=6)

for (i in 1:length(tau)) {
  
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 0.6)
  plot(ahull.obj,lwd=c(1,1,2), add=T, wpoints=F)
  
  
}



## connect contours over time i.e. plot estimated spline function for each u
Y1_cap = matrix(0, nrow = m, ncol = length(t)); Y2_cap = matrix(0, nrow = m, ncol = length(t))
for (i in 1:length(t)) {
  covariate = predict(basis_model, newx = t[i])[1,]
  bTx= apply(b, 1, function(x) t(x)%*%covariate)
  beta1Tx = finite.differences(U, bTx)
  
  fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
  b2Tx = vec(t(fmatrix))
  beta2Tx = t(finite.differences(U, b2Tx))
  
  Y1_cap[,i] = c(beta1Tx)
  Y2_cap[,i] = c(beta2Tx)
  
}

### plotting the estimated spline function
library(scatterplot3d)
s=scatterplot3d(Y1_cap[1,],Y2_cap[1,],t,type="l" , xlim=c(-0.2,1.6), ylim = c(0.5,2),
                cex.symbols =0.4, pch = 20,
                xlab = expression(paste(hat(Y)[1])),ylab = "", zlab = "t",
                color ="grey60")

for(i in 2:dim(U)[1]){
  s$points3d(Y1_cap[i,],Y2_cap[i,],t, type = "l", cex=0.4,pch=20,col="grey60")
}
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])-0.6
y <- dims[3]+ 0.08*diff(dims[3:4])-0.3
text(x,y,expression(paste(hat(Y)[2])),srt=10)

## observed spline function
s$points3d(mean1[,1],mean2[,1],t,col="red", cex=0.2, type = "l")

##median
#s$points3d(Y1_cap[145,],Y2_cap[145,],t, type = "l", cex=0.4,pch=20,col="blue")

#legend("topright", legend=c("Predicted", "True"),
#       col=c("black", "red"), lty=1, cex=0.8)






##### comparing flexible quantile contours with Kong's directional quantile envelopes #######
library(quantreg)
library(splines)
library(pracma)
library(Matrix)
library(mvtnorm)




####################### define necessary functions ##################

sav=function(v,a)
{
  a1=v-a
  d1=sign(a1)<=0
  a1[which(d1)]=0
  a2=-v-a
  d2=sign(a2)<=0
  a2[which(d2)]=0
  a1-a2
}



ADMM=function(X,Y,B,O,rho,tau0,lambda)
{
  
  pmd=dim(O)[1]
  
  eabs=10^-4
  erel=10^-2
  tau=rep(tau0,length(Y))
  one=rep(1,length(Y))
  #y=rep(y0,length(Y))
  sqrtly=sqrt(length(Y))
  sqrtY=sqrt(sum(Y^2))
  
  
  U=0
  resid=Y-X%*%B
  t=0
  tx=t(X)
  solve=solve(2*lambda/rho*O+tx%*%X)
  
  
  
  #A=Matrix(matrix(NA,nrow=dim(O),ncol=10),sparse=T)
  while (t<=10)
  {
    
    a=1/rho
    v=U+Y-X%*%B-(2*tau-one)/rho         
    rresid=sav(v,a)  
    
    BB=solve%*%tx%*%(Y+U-rresid) 
    
    s=rho*X%*%(BB-B)      
    r=-rresid-X%*%BB+Y
    #print(sum(abs(BB-B)))
    resid=rresid
    B=BB
    U=U-resid-X%*%B+Y
    yy=rho*U
    
    edual=sqrtly*eabs+erel*sqrt(sum(yy^2))
    epri=sqrtly*eabs+erel*max(sqrt(sum(resid^2)),sqrt(sum((X%*%B)^2)),sqrtY)
    #cat(norm(as.matrix(s),type="2"),edual,norm(as.matrix(r),type="2"),epri,"\n")
    if ((sqrt(sum(as.matrix(s)^2))<=edual)&(sqrt(sum(as.matrix(r)^2))<=epri)) break;
    
    t=t+1
    #A[,t-1]=B
  }
  
  #A
  B
  
}



qdir=function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2
  cbind(cos(ang),sin(ang))
}



kernel.g=function(t)
{
  dnorm(t,0,1)
}            #gaussian kernel




kernel.u=function(t)
{
  if(abs(t)<=1) 0.5
  else 0
}             #uniform kernel



kloc=function(r,r0)
{
  #h=h*2*pi/dirs
  s=sqrt(t(qd[r,]-qd[r0,])%*%(qd[r,]-qd[r0,]))
  kernel.g(1/hc*s)
}

kst=function(r,r0,s=sigma)
{
  bs=t(B[,r0]-B[,r])%*%solve(s)%*%(B[,r0]-B[,r])
  kernel.g(1/C*bs)
}


K=function(r,r0)
{
  #nk=1
  nor=1/((kernel.g(0))^2)
  nk=nor*as.numeric(kloc(r,r0)*kst(r,r0,s))
  nk
  #print(nk)
}    #obtain K



errorhat=function(X0,Y0,B)
{
  e=Y0-X0%*%B
}



Sigma=function(X,e=errorhat)
{
  
  ne=length(e)
  q=quantile(e,probs=tau[k])
  b=0.9*min(sd(e),(IQR(e))/1.34)*(ne^(-1/5))
  fq=1/(ne*b)*sum(kernel.g((q-e)/b))
  tau[k]*(1-tau[k])/(fq^2)*solve(t(X)%*%X)
}




Rr=function(r0,step)
{
  r00=r0+10
  rr=c(c(91:100),c(1:100),c(1:10))
  rr[(r00-floor(h^step)):(r00+floor(h^step))]
}






Xdesign=function(r0,step)
{
  x=matrix(NA,nrow=dim(simuX)[1],ncol=dim(simuX)[2])
  nop=floor(h^step)*2+1
  for (i in 1:nop)
  {
    if(i==1)
    {x=simuX*K(Rr0[1],r0)}
    else
      x=rbind(x,simuX*K(Rr0[i],r0))
    
  }    
  kron(x,b) 
  #print(K(Rr0[i],r0)) 
}


YwoK=function(step)
{
  y=vector()
  nop=floor(h^step)*2+1
  for (i in 1:nop)
  {
    y=rbind(y,simuY%*%qd[Rr0[i],])
  }                 
  y
  
}


Yresponse=function(r0,step)
{
  y=vector()
  nop=floor(h^step)*2+1
  for (i in 1:nop)
  {
    y=rbind(y,simuY%*%qd[Rr0[i],]*K(Rr0[i],r0))
  } 
  
  y
}	

Chisq=function(alpha,step,pm)
{
  x=alpha/(step)       #1-alpha/(c-3+1)
  qchisq(x,pm)
}




inln <- function(a,b) ## compute intersection of two lines (a point)
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}


actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}


abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}




eli=function(qln) 
{
  dirs=dim(qln)[1]
  clnid=1
  cln=1
  stop=FALSE
  filn=c(1:dirs)
  ss=0
  sln=0
  while (!stop)
  {
    abc=dertln(filn,clnid)
    a=qln[abc[1],]
    b=qln[abc[2],]
    c=qln[abc[3],]
    tt=actln(a,b,c)
    if (ss==0&tt==1) sln=filn[clnid]
    if (tt*ss==1&cln==sln) stop=TRUE
    #print('ss,tt,sln,cln')
    #print(c(ss,tt,sln,cln))
    ss=tt
    m=length(filn)
    if (tt==1)
    {
      if (clnid==m) {cln=filn[1] 
      clnid=1} else 
      {cln=filn[clnid+1] 
      clnid=clnid+1}
    } else { 
      clnidt=clnid
      if (clnid==1) {cln=filn[m] 
      clnid=m-1} else 
      {cln=filn[clnid-1] 
      clnid=clnid-1}
      filn=filn[-clnidt]
    }
  } 
  
  actlns=qln[filn,]
  
  
  
  m=dim(actlns)[1]
  
  
  
  intpts=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts[i,]=inln(actlns[i,],actlns[j,])
  }
  z=list(actlns=actlns,intpts=intpts)
  z 
  #print(z)                
}



PPlot=function(B,choice)         #plot the contour
{
  
  predmat = matrix(0,dirs)
  
  for (kk in 1:dirs)
  {
    bb=bs(xposition,df=NULL,intercept=T,knots=knots2,degree=3,Boundary.knots = range(simup))
    predmat[kk] = kron(xnew,bb)%*%B[,kk]
    
  }          # Fitting quantile regression
  
  
  
  qln=cbind(qd,-predmat)
  #print(qln)
  el=eli(qln)
  intpts=el$intpts 
  #el$actlns 
  #par(pty="s")
  if (k==1){
    plot(a[,2],a[,3], xlab=expression(paste(Y[1])), ylab=expression(paste(Y[2])), pch=20, cex=0.5, col="grey60")
    # points(amale[,3],amale[,4],col='grey60',pch=1,cex=0.5)
    #points(afemale[,3],afemale[,4],col='red',pch=1,cex=0.5)
    
  }
  if (choice==0)  
  {
    polygon(intpts,border="light blue",lwd=0.8,cex=.5)
  }  
  else 
  {
    
    polygon(intpts,border="black",lwd=2,cex=.5)
  } 
  
  
}  


Plot=function(B)         #plot the contour
{
  
  predmat = matrix(0,dirs)
  
  
  
  for (kk in 1:dirs)
  {
    bb=bs(xposition,df=NULL,knots=knots2,intercept=T,degree=3,Boundary.knots = range(simup))
    predmat[kk] = kron(xnew,bb)%*%B[,kk]
    
  }          # Fitting quantile regression
  
  
  
  qln=cbind(qd,-predmat)
  #print(qln)
  el=eli(qln)
  #intpts=el$intpts 
  el$actlns 
  
  
} 

tau = seq(0.006125,0.3, length.out = 6)

constant=seq(0.8,1.1, length.out = length(tau))

lamin = rep(0.001, length(tau)) ####### initial lambda

lam= rep(0.1, length(tau))     ####### ps lambda

knots2=seq(0.02,0.98,by=0.07)
knots2=seq(0.02,0.98,length.out = 6) ## for noncovex try




library(matrixcalc) ## for vec
library(Matrix) ## for kronecker product
library(gurobi)

t_length=10
t = simup= seq(0.01,1,length.out = t_length) ## match with simup
library(splines)
## intercept is not included by default.
# Number of knots = df-degree(3)-1(if intercept)
basis <- bs(x = t, df=7, degree = 3,
            Boundary.knots = c(0,1), intercept = TRUE)



n=300     ####### num of subjects
J=10    ####### num of positions
p=1     ####### num of covariates


tn=n*J

simuX=cbind(rep(1,n))    #design matrix X ##just intercept
simuXdata=kron(simuX,matrix(rep(1,J),ncol=1))

#simup=seq(0.01,1,length.out=J)     #position
simupos=rep(t,n)
b=bs(simup,df=NULL,intercept=T,knots=knots2, degree=3,Boundary.knots = range(simup))



## non convex data from in 1-t form
Y1_alltime = vec(Y1)
Y2_alltime = vec(Y2)
Y = cbind(Y1_alltime, Y2_alltime)
simuY = Y




dirs=100
qd = qdir(dirs)
b=bs(simup,df=NULL,intercept=T,knots=knots2,degree=3,Boundary.knots = range(simup))
p=dim(simuX)[2]
m=dim(b)[2]
J=length(simup)
x.ini=kron(simuX,b)

#obtain H'' and omega
L=dim(b)[1]
cc2=rep(c(1,-2,1,rep(0,(L-2))),L)
cc2=cc2[1:(L*L)]
cc2=matrix(cc2,nrow=L,byrow=T)
cc2[L,1]=0
ddH=cc2%*%b
d=diag(rep(1,p))
o=0

for (jj in 1:J)
{
  ott=kron(d,ddH[jj,]%*%t(ddH[jj,]))
  o=o+ott
} 


tau=seq(0.01, 0.28, length.out = 6)
#tau= c(0.01,0.29)
for (k in 1:length(tau)) {
  
  B=matrix(NA,ncol=dirs,nrow=p*m)
  #se=rep(NA,dim(X)[2])
  for (kk in 1:dirs)
  {
    fit=rq(simuY%*%qd[kk,]~bs(simupos,df=NULL,knots=knots2,degree=3,Boundary.knots = range(simupos)),tau=tau[k])
    B[,kk]=fit$coef
  }
  
  B0=B
  
  
  
  
  ########### Step 0: obtain initial beta in each diraction ###########
  
  B.initial=matrix(NA,ncol=dirs,nrow=p*m)
  
  
  for (i in 1:dirs)
  {
    y.ini=simuY%*%qd[i,]
    B.initial[,i]=ADMM(x.ini,y.ini,B0[,i],o,0.4,tau[k],lamin[k]) 
    
    #ADMM(x.ini,y.ini,B0[,i],o,0.4,tau[k],0.01)
    
    #QRRADMMCPP(x.ini,y.ini,o,B0[,i], 10 ,tau[k],0.4 ,0.01)    #x,y,omega,initial B, max iteration,tau,rho,lambda
  }            # obtain initial B, 100 directions
  
  B=B.initial
  
  
  
  
  tt = t[5]
  ########### Step 1-8 ###########
  xnew=matrix(c(1),nrow=1)
  xposition=tt
  yyydata=cbind(simupos,simuY)
  a=yyydata[(simupos>=xposition-0.0101)&(simupos<=xposition+0.0101),]
  
  
  alpha=0.8
  n=dim(simuX)[1]
  C=n^constant[k]*qchisq(alpha,1)
  h=1.15
  
  #PPlot(B.initial,0)

  
  B.update=matrix(NA,ncol=dirs,nrow=p*m)
  sig=matrix(NA,ncol=p*m,nrow=p*m*dirs)
  r0=seq(1,dirs,by=1)
  
  for (step in 1:8)
  {
    hc=h^step*2*pi/dirs
    nop=floor(h^step)*2+1
    X0=kron(matrix(rep(1,nop),ncol=1),x.ini)
    
    
    if(step>2)
    {
      if(length(list)!=0) 
      {
        r0=r0[-list]
      }
      
      if(length(r0)==0)  break
    }
    list=vector()
    print(step)
    print(r0)
    
    for (ii in 1:length(r0))
    {
      
      Rr0=Rr(r0[ii],step)
      
      if(step==1)
      {
        X00=x.ini
        Y00=simuY%*%qd[r0[ii],]
        e=errorhat(X00,Y00,B[,r0[ii]])
        s=Sigma(X00,e)
      }
      else
      {
        s=sig[(p*m*(r0[ii]-1)+1):(p*m*r0[ii]),]
      }
      
      X=Xdesign(r0[ii],step)
      Y=Yresponse(r0[ii],step)
      B.update[,r0[ii]]=ADMM(X,Y,B[,r0[ii]],o,0.4,tau[k],lam[k])	

      
      if(step>1)
      {
        dd=((B.update[,r0[ii]]-B1[,r0[ii]])^2)/diag(sig1[(p*m*(r0[ii]-1)+1):(p*m*r0[ii]),])
        cs=Chisq(alpha,step-1,1)	
        dcs=sum(as.numeric(dd>cs))
        
        if(dcs>1)       #Chisq(alpha,c,1)
        {
          list=c(list,ii)
          B.update[,r0[ii]]=B[,r0[ii]]	
        }
        
      }
      
      
      Y0=YwoK(step)
      e=errorhat(X0,Y0,B.update[,r0[ii]])
      
      sig[(p*m*(r0[ii]-1)+1):(p*m*r0[ii]),]=Sigma(X,e)
      
    }
  
    
    if(step==1)
    {
      B1=B.update
      sig1=sig
    }
    
    B=B.update
    
  }
  PPlot(B,1)
  
}







##################### Trivariate Functional Data ######################


### simulating from the data generating model
# Y_i(t) = X_i(t)^T \beta + epsilon_i(t) , i = 1,2,...,n   
# deterministic mean part plus the random convex error
set.seed(1512)
set.seed(1513)


t_length=6
t =seq(0.1,1,length.out = t_length)
library(splines)
## intercept is not included by default.
# Number of knots = df-degree(3)-1(if intercept)
basis <- bs(x = t, df=7, degree = 3,
            Boundary.knots = c(0,1), intercept = TRUE)


### regression coeffcients
## since we have generated b-spline matrix we provide gammas instead of betas
gamma1 = c(1, 0.7, 0.8, 0.9, 0.7, 0.8, 0.9)
gamma2 = c(1, 0.8, 0.7, 0.8, 0.9, 0.7, 0.6)
gamma3 = c(1, 0.7, 0.7, 0.8, 0.7, 0.9, 0.6)


n_each = 300
Y1=matrix(0, nrow=length(t), ncol = n_each); mean1=matrix(0, nrow=length(t), ncol = n_each)
Y2=matrix(0, nrow=length(t), ncol = n_each); mean2=matrix(0, nrow=length(t), ncol = n_each)
Y3=matrix(0, nrow=length(t), ncol = n_each); mean3=matrix(0, nrow=length(t), ncol = n_each)
for(i in 1:length(t)){
  
  X = runif(n_each, -1, 1)
  Z = runif(n_each, 0, 1)
  phi = runif(n_each, 0, 2*pi)
  R = 0.2*Z*(1+(1-abs(X))/2)
  
  epsilon1_t = X+ R*cos(phi)
  epsilon2_t = X^2 + R*sin(phi)
  epsilon3_t = X+ R 
  
  mean1[i,] = t(basis[i,])%*%gamma1
  mean2[i,] = t(basis[i,])%*%gamma2
  mean3[i,] = t(basis[i,])%*%gamma3
  
  Y1[i, ] = t(basis[i,])%*%gamma1 + epsilon1_t 
  Y2[i, ] = t(basis[i,])%*%gamma2 + epsilon2_t 
  Y3[i, ] = t(basis[i,])%*%gamma3 + epsilon3_t#+ rnorm(n_each, 0, 0.1) #+ epsilon3_t 
  
}

library(rgl)
bg3d("white")
rgl.points(Y1[4,],Y2[4,],Y3[4,], col="black")

library(scatterplot3d)
##plotting the simulated launches
s=scatterplot3d(Y1[3,],Y2[3,],Y3[3,],color="grey60", #xlim = c(-0.5,1.5), ylim = c(0,1.2),
                xlab = "x", ylab = "y",zlab = "z",angle=30, cex.symbols = 0.2)


library(scatterplot3d)
##plotting the simulated launches
s=scatterplot3d(Y1[,1],Y2[,1],t,color="grey60", #xlim = c(-0.5,1.5), ylim = c(0,1.2),
                xlab = "x", ylab = "y",zlab = "t",angle=30, cex.symbols = 0.2)

for(i in 2:n_each){
  s$points3d(Y1[,i],Y2[,i],t,col="grey60", cex=0.2)
}
s$points3d(mean1[,i],mean2[,i],t,col="red", cex=0.2, type = "l")


## Gurobi ##
library(gurobi)

##############################################
grid_size= 9
m = grid_size^3

Y1_alltime = vec(t(Y1))
Y2_alltime = vec(t(Y2))
Y3_alltime = vec(t(Y3))
Y= cbind(Y1_alltime, Y2_alltime, Y3_alltime)

n = n_each*length(t)

df = 7
basis_model = bs(x = t, df= df, degree = 3,
                 Boundary.knots = c(0,1), intercept = TRUE)

vec_n= vector()
for(i in 1:t_length){
  vec_i = c(rep(basis_model[i,], times= n_each))
  vec_n = c(vec_n, vec_i)
}


X =  matrix(vec_n, nrow = n, ncol =df, byrow=T)

#X = rep(1,n)  ## without slope
nu = rep(1/n, n)
mu = rep(1/m, m)
# simulating U from U[0,1]^2
#U = cbind(runif(m, 0,1), runif(m,0,1))


## simulating U from a 3d regular grid
x <- seq(0, 1, length.out = grid_size)
y <- seq(0, 1, length.out = grid_size)
z <- seq(0, 1, length.out = grid_size)

U <- as.matrix(expand.grid(x = x, y = y, z = z))


library(matrixcalc) ## for vec
library(Matrix) ## for kronecker product
library(Matrix) ## for constructing sparse matrix



######################## DUAL ############################
obj1 = vec(as.matrix(nu))
obj2 = vec(mu%*%t(nu)%*%X)

model_dual = list()

##objective function
model_dual$obj = c(obj1, obj2)

## A matrix of constraints (lhs)
A1_dual = Matrix(kronecker(X= diag(n), Y= rep(1,m)), sparse = T)
A2_dual = Matrix(kronecker(X= X , Y= diag(m)), sparse = T)

model_dual$A = cbind(A1_dual,A2_dual)

## b vector of constraints (rhs)
rhs_dual = vec(U%*%t(Y))
model_dual$rhs <- rhs_dual

model_dual$vtype <- 'C'  ## Continuous variable type (default)
model_dual$sense <- '>'
model_dual$modelsense <- 'min'

params <- list(OutputFlag=0)
result_dual <- gurobi(model_dual, params)


print(result_dual$objval)
#print(result_dual$x)
psi = result_dual$x[1:n]
vec_b = result_dual$x[-(1:n)]
b = matrix(vec_b, nrow = m, ncol = df) ## b constructed according to x

## obtaining the beta's using numererical differentiation
library(numDeriv)

#covariate = c(1, t[15], t[15]^2, t[15]^3)
## predicting at  t = t[15]
covariate = predict(basis_model, newx = t[3])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
#plot(bTx)


finite.differences.3d <- function(x, y) {
  
  stepsize = abs(x[1,1]-x[2,1])
  n <- dim(x)[1]
  fdx=vector()
  
  u=1
  while (u < n ) {
    
    fdx_f = rep(0,grid_size); fdx_b= rep(0,grid_size); rep(0,grid_size) 
    # Iterate through the values using the forward differencing method
    for (i in 1:(grid_size-1)) {
      fdx_f[i] <- (y[i+1 +(u-1)] - y[i + (u-1)]) / stepsize
    }
    
    # Iterate through the values using the backward differencing method
    for (i in 2:grid_size) {
      fdx_b[i] <- (y[i + (u-1)] - y[i-1 + (u-1)]) / stepsize
    }
    fdx_temp <- (fdx_f+fdx_b)/2
    fdx_temp[1] <- fdx_f[1]
    fdx_temp[grid_size] <- fdx_b[grid_size]
    
    fdx = c(fdx, fdx_temp)
    
    u = u + grid_size
  }
  return(fdx)
}




y1_cap = finite.differences.3d(U, bTx)



f= bTx

j=1
fu2 = vector()
for(i in 1:grid_size)
{
  for (u1 in 1:grid_size) {
    for (u2 in 1:grid_size) {
      fu2[j]  = f[grid_size^2*(i-1) + u1 + grid_size*(u2-1)]
      j=j+1
    }
  }
}

y2_cap_distortred = finite.differences.3d(U, fu2)

j=1
y2_cap = vector()
for(i in 1:grid_size)
{
  for (u1 in 1:grid_size) {
    for (u2 in 1:grid_size) {
      y2_cap[grid_size^2*(i-1) + u1 + grid_size*(u2-1)]= y2_cap_distortred[j]  
      j=j+1
    }
  }
}


j=1
fu3 = vector()
for (u1 in 1:grid_size^2) {
  for (u2 in 1:grid_size) {
    fu3[j]  = f[u1 + grid_size^2*(u2-1)]
    j=j+1
  }
}


y3_cap_distorted = finite.differences.3d(U, fu3)

j=1
y3_cap = vector()
for (u1 in 1:grid_size^2) {
  for (u2 in 1:grid_size) {
    y3_cap[u1 + grid_size^2*(u2-1)]=  y3_cap_distorted[j]  
    j=j+1
  }
}




library(scatterplot3d)
##plotting the simulated launches
s=scatterplot3d(Y1[3,],Y2[3,],Y3[3,],color="grey60", #xlim = c(-0.5,1.5), ylim = c(0,1.2),
                xlab = "x", ylab = "y",zlab = "z",angle=30, cex.symbols = 0.2)

s$points3d(y1_cap,y2_cap,y3_cap,col="red", cex=0.2)






### spherical uniform 3d distriution
radius = seq(0.00001, 1, length.out = grid_size*sqrt(grid_size))

num <- grid_size*sqrt(grid_size) # number of points you want on each sphere

theta <- runif(num,0,2*pi)
phi = runif(num,0,pi)
x=radius[1]*sin(phi)*cos(theta)
y=radius[1]*sin(phi)*sin(theta)
z=radius[1]*cos(phi)

pts_all <- cbind(x,y,z)
for(i in 2:length(radius)){
  theta <- runif(num,0,2*pi)
  phi = runif(num,0,pi)
  x=radius[i]*sin(phi)*cos(theta)
  y=radius[i]*sin(phi)*sin(theta)
  z=radius[i]*cos(phi)
  
  pts.sphere <- cbind(x,y,z)
  pts_all = rbind(pts_all, pts.sphere)
}

U_n = pts_all


euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

library(pdist)
## U_n as rows, Y_n as columns
dist_mat = as.matrix(pdist(X=U_n,Y=U))
ranks = rank(c(dist_mat))

##defining cost matrix
cmat = matrix(ranks, nrow=grid_size^3, ncol=grid_size^3)

library(adagio)
matching = assignment(cmat)

perm = matching$perm
ycap = cbind(y1_cap, y2_cap, y3_cap)


tau= seq(0.01,1,length.out=6)


library(alphashape3d)
#ashape.obj = ashape3d(ycap, alpha = 0.3)
#bg3d("white")
#rgl.points(Y1[3,],Y2[3,],Y3[3,], col="black")
#plot(ashape.obj,clear = F, transparency = 0.4, triangles = T, edges = F, vertices = F, col = c(4, 4, 4), indexAlpha = "all")

#rgl.points(y1_cap,y2_cap,y3_cap, col="red")
#plot(ashape.obj,clear = F, transparency = 0.4, triangles = T, edges = F, vertices = F, col = c(4, 4, 4), indexAlpha = "all")

bg3d("white")
rgl.points(Y1[3,],Y2[3,],Y3[3,], col="black")
tau= seq(0.1,0.9,length.out=6)
tau= 0.25
for (i in 1:length(tau)) {
  
  U_number = which(apply(U_n, 1, norm_vec) <= tau)
  U_centered_corresponding = perm[U_number]
  Y_tau = ycap[U_centered_corresponding,]
  
  bg3d("white")
  rgl.points(Y1[3,],Y2[3,],Y3[3,], col="black")
  ashape.obj = ashape3d(Y_tau, alpha = 0.25)
  plot(ashape.obj,clear = F, transparency = 0.4, triangles = T, edges = F, vertices = F, col = c(4, 4, 4), indexAlpha = "all")
  
  
}



######################### repeating simulation to compute mse #################

## for all time points
Y1_cap = matrix(0, nrow = m, ncol = length(t)); Y2_cap = matrix(0, nrow = m, ncol = length(t));  Y3_cap = matrix(0, nrow = m, ncol = length(t))
for (k in 1:length(t)) {
  covariate = predict(basis_model, newx = t[k])[1,]
  bTx= apply(b, 1, function(x) t(x)%*%covariate)
  
  
  
  y1_cap = finite.differences.3d(U, bTx)
  
  
  f= bTx
  
  j=1
  fu2 = vector()
  for(i in 1:grid_size)
  {
    for (u1 in 1:grid_size) {
      for (u2 in 1:grid_size) {
        fu2[j]  = f[grid_size^2*(i-1) + u1 + grid_size*(u2-1)]
        j=j+1
      }
    }
  }
  
  y2_cap_distortred = finite.differences.3d(U, fu2)
  
  j=1
  y2_cap = vector()
  for(i in 1:grid_size)
  {
    for (u1 in 1:grid_size) {
      for (u2 in 1:grid_size) {
        y2_cap[grid_size^2*(i-1) + u1 + grid_size*(u2-1)]= y2_cap_distortred[j]  
        j=j+1
      }
    }
  }
  
  
  j=1
  fu3 = vector()
  for (u1 in 1:grid_size^2) {
    for (u2 in 1:grid_size) {
      fu3[j]  = f[u1 + grid_size^2*(u2-1)]
      j=j+1
    }
  }
  
  
  y3_cap_distorted = finite.differences.3d(U, fu3)
  
  j=1
  y3_cap = vector()
  for (u1 in 1:grid_size^2) {
    for (u2 in 1:grid_size) {
      y3_cap[u1 + grid_size^2*(u2-1)]=  y3_cap_distorted[j]  
      j=j+1
    }
  }
  
  
  Y1_cap[,k] = y1_cap
  Y2_cap[,k] = y2_cap
  Y3_cap[,k] = y3_cap
  
}


euc_sim = vector()
simsize = 10
for (sim in 1:simsize) {
  
  t_length=6
  t =seq(0.1,1,length.out = t_length)
  library(splines)
  basis <- bs(x = t, df=7, degree = 3,
              Boundary.knots = c(0,1), intercept = TRUE)
  
  
  gamma1 = c(1, 0.7, 0.8, 0.9, 0.7, 0.8, 0.9)
  gamma2 = c(1, 0.8, 0.7, 0.8, 0.9, 0.7, 0.6)
  gamma3 = c(1, 0.7, 0.7, 0.8, 0.7, 0.9, 0.6)
  
  
  n_each = 200
  Y1=matrix(0, nrow=length(t), ncol = n_each); mean1=matrix(0, nrow=length(t), ncol = n_each)
  Y2=matrix(0, nrow=length(t), ncol = n_each); mean2=matrix(0, nrow=length(t), ncol = n_each)
  Y3=matrix(0, nrow=length(t), ncol = n_each); mean3=matrix(0, nrow=length(t), ncol = n_each)
  for(i in 1:length(t)){
    
    X = runif(n_each, -0.5, 0.5)
    Z = runif(n_each, 0, 0.5)
    phi = runif(n_each, 0, 2*pi)
    R = 0.2*Z*(0.5+(0.5-abs(X))/2)
    
    epsilon1_t = X+ R*cos(phi)
    epsilon2_t = X^2 + R*sin(phi)
    epsilon3_t =X+ R*sin(phi)
    
    mean1[i,] = t(basis[i,])%*%gamma1
    mean2[i,] = t(basis[i,])%*%gamma2
    mean3[i,] = t(basis[i,])%*%gamma3
    
    Y1[i, ] = t(basis[i,])%*%gamma1 + epsilon1_t 
    Y2[i, ] = t(basis[i,])%*%gamma2 + epsilon2_t 
    Y3[i, ] = t(basis[i,])%*%gamma3 + epsilon3_t#+ rnorm(n_each, 0, 0.1) #+ epsilon3_t 
    
  }
  
  
  ##############################################
  grid_size= 9
  m = grid_size^3
  
  Y1_alltime = vec(t(Y1))
  Y2_alltime = vec(t(Y2))
  Y3_alltime = vec(t(Y3))
  Y= cbind(Y1_alltime, Y2_alltime, Y3_alltime)
  
  n = n_each*length(t)
  
  df = 7
  basis_model = bs(x = t, df= df, degree = 3,
                   Boundary.knots = c(0,1), intercept = TRUE)
  
  vec_n= vector()
  for(i in 1:length(t)){
    vec_i = c(rep(basis_model[i,], times= n_each))
    vec_n = c(vec_n, vec_i)
  }
  
  
  X =  matrix(vec_n, nrow = n, ncol =df, byrow=T)
  
  nu = rep(1/n, n)
  mu = rep(1/m, m)
  
  ## simulating U from a 3d regular grid
  x <- seq(0, 1, length.out = grid_size)
  y <- seq(0, 1, length.out = grid_size)
  z <- seq(0, 1, length.out = grid_size)
  
  U <- as.matrix(expand.grid(x = x, y = y, z = z))
  
  
  ######################## DUAL ############################
  obj1 = vec(as.matrix(nu))
  obj2 = vec(mu%*%t(nu)%*%X)
  
  model_dual = list()
  
  ##objective function
  model_dual$obj = c(obj1, obj2)
  
  ## A matrix of constraints (lhs)
  A1_dual = Matrix(kronecker(X= diag(n), Y= rep(1,m)), sparse = T)
  A2_dual = Matrix(kronecker(X= X , Y= diag(m)), sparse = T)
  
  model_dual$A = cbind(A1_dual,A2_dual)
  
  ## b vector of constraints (rhs)
  rhs_dual = vec(U%*%t(Y))
  model_dual$rhs <- rhs_dual
  
  model_dual$vtype <- 'C'  ## Continuous variable type (default)
  model_dual$sense <- '>'
  model_dual$modelsense <- 'min'
  
  params <- list(OutputFlag=0)
  result_dual <- gurobi(model_dual, params)
  
  
  print(result_dual$objval)
  #print(result_dual$x)
  psi = result_dual$x[1:n]
  vec_b = result_dual$x[-(1:n)]
  b = matrix(vec_b, nrow = m, ncol = df) ## b constructed according to x
  
  
  
  ## for all time points
  Y1_cap = matrix(0, nrow = m, ncol = length(t)); Y2_cap = matrix(0, nrow = m, ncol = length(t));  Y3_cap = matrix(0, nrow = m, ncol = length(t))
  for (k in 1:length(t)) {
    covariate = predict(basis_model, newx = t[k])[1,]
    bTx= apply(b, 1, function(x) t(x)%*%covariate)
    
    y1_cap = finite.differences.3d(U, bTx)
    
    
    f= bTx
    
    j=1
    fu2 = vector()
    for(i in 1:grid_size)
    {
      for (u1 in 1:grid_size) {
        for (u2 in 1:grid_size) {
          fu2[j]  = f[grid_size^2*(i-1) + u1 + grid_size*(u2-1)]
          j=j+1
        }
      }
    }
    
    y2_cap_distortred = finite.differences.3d(U, fu2)
    
    j=1
    y2_cap = vector()
    for(i in 1:grid_size)
    {
      for (u1 in 1:grid_size) {
        for (u2 in 1:grid_size) {
          y2_cap[grid_size^2*(i-1) + u1 + grid_size*(u2-1)]= y2_cap_distortred[j]  
          j=j+1
        }
      }
    }
    
    
    j=1
    fu3 = vector()
    for (u1 in 1:grid_size^2) {
      for (u2 in 1:grid_size) {
        fu3[j]  = f[u1 + grid_size^2*(u2-1)]
        j=j+1
      }
    }
    
    
    y3_cap_distorted = finite.differences.3d(U, fu3)
    
    j=1
    y3_cap = vector()
    for (u1 in 1:grid_size^2) {
      for (u2 in 1:grid_size) {
        y3_cap[u1 + grid_size^2*(u2-1)]=  y3_cap_distorted[j]  
        j=j+1
      }
    }
    
    
    Y1_cap[,k] = y1_cap
    Y2_cap[,k] = y2_cap
    Y3_cap[,k] = y3_cap
    
  }
  
  
  
  ### mse
  
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  true_mean = cbind(mean1[,1],mean2[,1], mean3[,1])
  predicted = cbind(Y1_cap[365,],Y2_cap[365,], Y3_cap[365,])
  
  euc=0
  for (i in 1:length(t)) {
    euc= euc + euc.dist(true_mean[i,],predicted[i,])
  }
  
  euc_sim[sim] = euc/length(t)
  
  
}


mean(euc_sim)


######### Application to Air Pollution Data ####

#### Data preprocessing
library(rgdal)



# read raw data 2011 covariates####
grib_GDAL_201110 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20111001_0000_000.grb")
grib_GDAL_201111 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20111101_0000_000.grb")
grib_GDAL_201112 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20111201_0000_000.grb")

# read raw data 2012 covariates####
grib_GDAL_201201 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20120101_0000_000.grb")
grib_GDAL_201202 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20120201_0000_000.grb")
grib_GDAL_201203 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20120301_0000_000.grb")
grib_GDAL_201204 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20120401_0000_000.grb")
grib_GDAL_201205 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20120501_0000_000.grb")
grib_GDAL_201206 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20120601_0000_000.grb")
grib_GDAL_201207 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20120701_0000_000.grb")
grib_GDAL_201208 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20120801_0000_000.grb")
grib_GDAL_201209 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20120901_0000_000.grb")
grib_GDAL_201210 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20121001_0000_000.grb")
grib_GDAL_201211 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20121101_0000_000.grb")
grib_GDAL_201212 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20121201_0000_000.grb")

## 2013
grib_GDAL_201301 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20130101_0000_000.grb")
grib_GDAL_201302 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20130201_0000_000.grb")
grib_GDAL_201303 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20130301_0000_000.grb")
grib_GDAL_201304 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20130401_0000_000.grb")
grib_GDAL_201305 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20130501_0000_000.grb")
grib_GDAL_201306 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20130601_0000_000.grb")
grib_GDAL_201307 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20130701_0000_000.grb")
grib_GDAL_201308 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20130801_0000_000.grb")
grib_GDAL_201309 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20130901_0000_000.grb")
grib_GDAL_201310 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20131001_0000_000.grb")
grib_GDAL_201311 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20131101_0000_000.grb")
grib_GDAL_201312 <- readGDAL("/Users/agarwag/Desktop/KAUST/Project 3/datasets/monthly narr/narrmon-a_221_20131201_0000_000.grb")



# transform LCC coordinate to Long/Lat ####
longlat <- spTransform(SpatialPoints(coordinates(grib_GDAL_201201), proj4string=grib_GDAL_201201@proj4string), CRS("+proj=longlat +datum=WGS84"))
narr_201110 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201110@data[,91], TMP=grib_GDAL_201110@data[,263],
                          RH=grib_GDAL_201110@data[,156], TCDC=grib_GDAL_201110@data[,214],
                          APCP=grib_GDAL_201110@data[,4], HGT850=grib_GDAL_201110@data[,63],
                          HGT500=grib_GDAL_201110@data[,73],
                          UGRD=grib_GDAL_201110@data[,277], VGRD=grib_GDAL_201110@data[,323],
                          WS=sqrt(grib_GDAL_201110@data[,277]^2 + grib_GDAL_201110@data[,323]^2))


narr_201111 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201111@data[,91], TMP=grib_GDAL_201111@data[,263],
                          RH=grib_GDAL_201111@data[,156], TCDC=grib_GDAL_201111@data[,214],
                          APCP=grib_GDAL_201111@data[,4], HGT850=grib_GDAL_201111@data[,63],
                          HGT500=grib_GDAL_201111@data[,73],
                          UGRD=grib_GDAL_201111@data[,277], VGRD=grib_GDAL_201111@data[,323],
                          WS=sqrt(grib_GDAL_201111@data[,277]^2 + grib_GDAL_201111@data[,323]^2))


narr_201112 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201112@data[,91], TMP=grib_GDAL_201112@data[,263],
                          RH=grib_GDAL_201112@data[,156], TCDC=grib_GDAL_201112@data[,214],
                          APCP=grib_GDAL_201112@data[,4], HGT850=grib_GDAL_201112@data[,63],
                          HGT500=grib_GDAL_201112@data[,73],
                          UGRD=grib_GDAL_201112@data[,277], VGRD=grib_GDAL_201112@data[,323],
                          WS=sqrt(grib_GDAL_201112@data[,277]^2 + grib_GDAL_201112@data[,323]^2))



## 2012
narr_201201 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201201@data[,91], TMP=grib_GDAL_201201@data[,263],
                          RH=grib_GDAL_201201@data[,156], TCDC=grib_GDAL_201201@data[,214],
                          APCP=grib_GDAL_201201@data[,4], HGT850=grib_GDAL_201201@data[,63],
                          HGT500=grib_GDAL_201201@data[,73],
                          UGRD=grib_GDAL_201201@data[,277], VGRD=grib_GDAL_201201@data[,323],
                          WS=sqrt(grib_GDAL_201201@data[,277]^2 + grib_GDAL_201201@data[,323]^2))


narr_201202 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201202@data[,91], TMP=grib_GDAL_201202@data[,263],
                          RH=grib_GDAL_201202@data[,156], TCDC=grib_GDAL_201202@data[,214],
                          APCP=grib_GDAL_201202@data[,4], HGT850=grib_GDAL_201202@data[,63],
                          HGT500=grib_GDAL_201202@data[,73],
                          UGRD=grib_GDAL_201202@data[,277], VGRD=grib_GDAL_201202@data[,323],
                          WS=sqrt(grib_GDAL_201202@data[,277]^2 + grib_GDAL_201202@data[,323]^2))


narr_201203 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201203@data[,91], TMP=grib_GDAL_201203@data[,263],
                          RH=grib_GDAL_201203@data[,156], TCDC=grib_GDAL_201203@data[,214],
                          APCP=grib_GDAL_201203@data[,4], HGT850=grib_GDAL_201203@data[,63],
                          HGT500=grib_GDAL_201203@data[,73],
                          UGRD=grib_GDAL_201203@data[,277], VGRD=grib_GDAL_201203@data[,323],
                          WS=sqrt(grib_GDAL_201203@data[,277]^2 + grib_GDAL_201203@data[,323]^2))


narr_201204 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201204@data[,91], TMP=grib_GDAL_201204@data[,263],
                          RH=grib_GDAL_201204@data[,156], TCDC=grib_GDAL_201204@data[,214],
                          APCP=grib_GDAL_201204@data[,4], HGT850=grib_GDAL_201204@data[,63],
                          HGT500=grib_GDAL_201204@data[,73],
                          UGRD=grib_GDAL_201204@data[,277], VGRD=grib_GDAL_201204@data[,323],
                          WS=sqrt(grib_GDAL_201204@data[,277]^2 + grib_GDAL_201204@data[,323]^2))


narr_201205 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201205@data[,91], TMP=grib_GDAL_201205@data[,263],
                          RH=grib_GDAL_201205@data[,156], TCDC=grib_GDAL_201205@data[,214],
                          APCP=grib_GDAL_201205@data[,4], HGT850=grib_GDAL_201205@data[,63],
                          HGT500=grib_GDAL_201205@data[,73],
                          UGRD=grib_GDAL_201205@data[,277], VGRD=grib_GDAL_201205@data[,323],
                          WS=sqrt(grib_GDAL_201205@data[,277]^2 + grib_GDAL_201205@data[,323]^2))


narr_201206 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201206@data[,91], TMP=grib_GDAL_201206@data[,263],
                          RH=grib_GDAL_201206@data[,156], TCDC=grib_GDAL_201206@data[,214],
                          APCP=grib_GDAL_201206@data[,4], HGT850=grib_GDAL_201206@data[,63],
                          HGT500=grib_GDAL_201206@data[,73],
                          UGRD=grib_GDAL_201206@data[,277], VGRD=grib_GDAL_201206@data[,323],
                          WS=sqrt(grib_GDAL_201206@data[,277]^2 + grib_GDAL_201206@data[,323]^2))


narr_201207 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201207@data[,91], TMP=grib_GDAL_201207@data[,263],
                          RH=grib_GDAL_201207@data[,156], TCDC=grib_GDAL_201207@data[,214],
                          APCP=grib_GDAL_201207@data[,4], HGT850=grib_GDAL_201207@data[,63],
                          HGT500=grib_GDAL_201207@data[,73],
                          UGRD=grib_GDAL_201207@data[,277], VGRD=grib_GDAL_201207@data[,323],
                          WS=sqrt(grib_GDAL_201207@data[,277]^2 + grib_GDAL_201207@data[,323]^2))


narr_201208 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201208@data[,91], TMP=grib_GDAL_201208@data[,263],
                          RH=grib_GDAL_201208@data[,156], TCDC=grib_GDAL_201208@data[,214],
                          APCP=grib_GDAL_201208@data[,4], HGT850=grib_GDAL_201208@data[,63],
                          HGT500=grib_GDAL_201208@data[,73],
                          UGRD=grib_GDAL_201208@data[,277], VGRD=grib_GDAL_201208@data[,323],
                          WS=sqrt(grib_GDAL_201208@data[,277]^2 + grib_GDAL_201208@data[,323]^2))


narr_201209 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201209@data[,91], TMP=grib_GDAL_201209@data[,263],
                          RH=grib_GDAL_201209@data[,156], TCDC=grib_GDAL_201209@data[,214],
                          APCP=grib_GDAL_201209@data[,4], HGT850=grib_GDAL_201209@data[,63],
                          HGT500=grib_GDAL_201209@data[,73],
                          UGRD=grib_GDAL_201209@data[,277], VGRD=grib_GDAL_201209@data[,323],
                          WS=sqrt(grib_GDAL_201209@data[,277]^2 + grib_GDAL_201209@data[,323]^2))


narr_201210 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201210@data[,91], TMP=grib_GDAL_201210@data[,263],
                          RH=grib_GDAL_201210@data[,156], TCDC=grib_GDAL_201210@data[,214],
                          APCP=grib_GDAL_201210@data[,4], HGT850=grib_GDAL_201210@data[,63],
                          HGT500=grib_GDAL_201210@data[,73],
                          UGRD=grib_GDAL_201210@data[,277], VGRD=grib_GDAL_201210@data[,323],
                          WS=sqrt(grib_GDAL_201210@data[,277]^2 + grib_GDAL_201210@data[,323]^2))


narr_201211 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201211@data[,91], TMP=grib_GDAL_201211@data[,263],
                          RH=grib_GDAL_201211@data[,156], TCDC=grib_GDAL_201211@data[,214],
                          APCP=grib_GDAL_201211@data[,4], HGT850=grib_GDAL_201211@data[,63],
                          HGT500=grib_GDAL_201211@data[,73],
                          UGRD=grib_GDAL_201211@data[,277], VGRD=grib_GDAL_201211@data[,323],
                          WS=sqrt(grib_GDAL_201211@data[,277]^2 + grib_GDAL_201211@data[,323]^2))


narr_201212 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201212@data[,91], TMP=grib_GDAL_201212@data[,263],
                          RH=grib_GDAL_201212@data[,156], TCDC=grib_GDAL_201212@data[,214],
                          APCP=grib_GDAL_201212@data[,4], HGT850=grib_GDAL_201212@data[,63],
                          HGT500=grib_GDAL_201212@data[,73],
                          UGRD=grib_GDAL_201212@data[,277], VGRD=grib_GDAL_201212@data[,323],
                          WS=sqrt(grib_GDAL_201212@data[,277]^2 + grib_GDAL_201212@data[,323]^2))


### 2013
narr_201301 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201301@data[,91], TMP=grib_GDAL_201301@data[,263],
                          RH=grib_GDAL_201301@data[,156], TCDC=grib_GDAL_201301@data[,214],
                          APCP=grib_GDAL_201301@data[,4], HGT850=grib_GDAL_201301@data[,63],
                          HGT500=grib_GDAL_201301@data[,73],
                          UGRD=grib_GDAL_201301@data[,277], VGRD=grib_GDAL_201301@data[,323],
                          WS=sqrt(grib_GDAL_201301@data[,277]^2 + grib_GDAL_201301@data[,323]^2))


narr_201302 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201302@data[,91], TMP=grib_GDAL_201302@data[,263],
                          RH=grib_GDAL_201302@data[,156], TCDC=grib_GDAL_201302@data[,214],
                          APCP=grib_GDAL_201302@data[,4], HGT850=grib_GDAL_201302@data[,63],
                          HGT500=grib_GDAL_201302@data[,73],
                          UGRD=grib_GDAL_201302@data[,277], VGRD=grib_GDAL_201302@data[,323],
                          WS=sqrt(grib_GDAL_201302@data[,277]^2 + grib_GDAL_201302@data[,323]^2))


narr_201303 <- data.frame(Day=01, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          HPBL=grib_GDAL_201303@data[,91], TMP=grib_GDAL_201303@data[,263],
                          RH=grib_GDAL_201303@data[,156], TCDC=grib_GDAL_201303@data[,214],
                          APCP=grib_GDAL_201303@data[,4], HGT850=grib_GDAL_201303@data[,63],
                          HGT500=grib_GDAL_201303@data[,73],
                          UGRD=grib_GDAL_201303@data[,277], VGRD=grib_GDAL_201303@data[,323],
                          WS=sqrt(grib_GDAL_201303@data[,277]^2 + grib_GDAL_201303@data[,323]^2))







############################ read pm 2.5  ################################
pm25_day_2011 <- read.table("/Users/agarwag/Desktop/KAUST/Project 3/datasets/pm2.5/2011_pm25_daily_average.txt", header = TRUE, sep = ",")
pm25_day_2012 <- read.table("/Users/agarwag/Desktop/KAUST/Project 3/datasets/pm2.5/2012_pm25_daily_average.txt", header = TRUE, sep = ",")
pm25_day_2013 <- read.table("/Users/agarwag/Desktop/KAUST/Project 3/datasets/pm2.5/2013_pm25_daily_average.txt", header = TRUE, sep = ",")

names(pm25_day_2013)
pm25_day_2011$Month <- as.numeric(paste(substr(pm25_day_2011$Date,1,4),substr(pm25_day_2011$Date,6,7),sep=""))
pm25_mon_2011 <- aggregate(pm25_daily_average.ug.m3. ~ Month + FIPS + Longitude + Latitude, data=pm25_day_2011, mean)
names(pm25_mon_2011)[names(pm25_mon_2011)=="pm25_daily_average.ug.m3."] <- "pm25_monthly_mean"

pm25_day_2012$Month <- as.numeric(paste(substr(pm25_day_2012$Date,1,4),substr(pm25_day_2012$Date,6,7),sep=""))
pm25_mon_2012 <- aggregate(pm25_daily_average ~ Month + FIPS + Longitude + Latitude, data=pm25_day_2012, mean)
names(pm25_mon_2012)[names(pm25_mon_2012)=="pm25_daily_average"] <- "pm25_monthly_mean"

pm25_day_2013$Month <- as.numeric(paste(substr(pm25_day_2013$Date,1,4),substr(pm25_day_2013$Date,6,7),sep=""))
pm25_mon_2013 <- aggregate(pm25_daily_average ~ Month + FIPS + Longitude + Latitude, data=pm25_day_2013, mean)
names(pm25_mon_2013)[names(pm25_mon_2013)=="pm25_daily_average"] <- "pm25_monthly_mean"

#str(pm25_mon_2012); head(pm25_mon_2012) # 867396 obs
pm25_201110 <- subset(pm25_mon_2011, Month==201110); str(pm25_201110)
pm25_201111 <- subset(pm25_mon_2011, Month==201111); str(pm25_201111)
pm25_201112 <- subset(pm25_mon_2011, Month==201112); str(pm25_201112)

#2012
pm25_201201 <- subset(pm25_mon_2012, Month==201201); str(pm25_201201)
pm25_201202 <- subset(pm25_mon_2012, Month==201202); str(pm25_201202)
pm25_201203 <- subset(pm25_mon_2012, Month==201203); str(pm25_201203)
pm25_201204 <- subset(pm25_mon_2012, Month==201204); str(pm25_201204)
pm25_201205 <- subset(pm25_mon_2012, Month==201205); str(pm25_201205)
pm25_201206 <- subset(pm25_mon_2012, Month==201206); str(pm25_201206)
pm25_201207 <- subset(pm25_mon_2012, Month==201207); str(pm25_201207)
pm25_201208 <- subset(pm25_mon_2012, Month==201208); str(pm25_201208)
pm25_201209 <- subset(pm25_mon_2012, Month==201209); str(pm25_201209)
pm25_201210 <- subset(pm25_mon_2012, Month==201210); str(pm25_201210)
pm25_201211 <- subset(pm25_mon_2012, Month==201211); str(pm25_201211)
pm25_201212 <- subset(pm25_mon_2012, Month==201212); str(pm25_201212)


#2013
pm25_201301 <- subset(pm25_mon_2013, Month==201301); str(pm25_201301)
pm25_201302 <- subset(pm25_mon_2013, Month==201302); str(pm25_201302)
pm25_201303 <- subset(pm25_mon_2013, Month==201303); str(pm25_201303)



#### merging the cordinates
long_pm25<- unique(pm25_201201$Longitude)
lat_pm25 <- unique(pm25_201201$Latitude)

data_2011_10 <- subset(narr_201110, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2011_11 <- subset(narr_201111, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2011_12 <- subset(narr_201112, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))


data_2012_01 <- subset(narr_201201, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_02 <- subset(narr_201202, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_03 <- subset(narr_201203, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_04 <- subset(narr_201204, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_05 <- subset(narr_201205, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_06 <- subset(narr_201206, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_07 <- subset(narr_201207, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_08 <- subset(narr_201208, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_09 <- subset(narr_201209, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_10 <- subset(narr_201210, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_11 <- subset(narr_201211, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2012_12 <- subset(narr_201212, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))

data_2013_01 <- subset(narr_201301, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2013_02 <- subset(narr_201302, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))
data_2013_03 <- subset(narr_201303, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))


#==============================================================================================####
library(geosphere)

LL_pm25 <- cbind(pm25_201201$Longitude, pm25_201201$Latitude)
LL_narr <- cbind(data_2012_01$Longitude, data_2012_01$Latitude)
dim(LL_pm25); dim(LL_narr)

DMatrix <- distm(LL_pm25,LL_narr)    # 1398.215 sec
tmp <- apply(DMatrix, 1, min)                # 20.734 sec
min_ind <- which(DMatrix==tmp, arr.ind=T)         # 3.275 sec
id_data <- min_ind[order(min_ind[,1]),2]          # 0.002 sec

## pm2.5
pm25_201110$id_data <- id_data
data_pm25_201110 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201110, mean)

pm25_201111$id_data <- id_data
data_pm25_201111 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201111, mean)

pm25_201112$id_data <- id_data
data_pm25_201112 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201112, mean)



pm25_201201$id_data <- id_data
data_pm25_201201 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201201, mean)

pm25_201202$id_data <- id_data
data_pm25_201202 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201202, mean)

pm25_201203$id_data <- id_data
data_pm25_201203 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201203, mean)

pm25_201204$id_data <- id_data
data_pm25_201204 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201204, mean)

pm25_201205$id_data <- id_data
data_pm25_201205 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201205, mean)

pm25_201206$id_data <- id_data
data_pm25_201206 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201206, mean)

pm25_201207$id_data <- id_data
data_pm25_201207 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201207, mean)

pm25_201208$id_data <- id_data
data_pm25_201208 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201208, mean)

pm25_201209$id_data <- id_data
data_pm25_201209 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201209, mean)

pm25_201210$id_data <- id_data
data_pm25_201210 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201210, mean)

pm25_201211$id_data <- id_data
data_pm25_201211 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201211, mean)

pm25_201212$id_data <- id_data
data_pm25_201212 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201212, mean)

## 2013

pm25_201301$id_data <- id_data
data_pm25_201301 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201301, mean)

pm25_201302$id_data <- id_data
data_pm25_201302 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201302, mean)

pm25_201303$id_data <- id_data
data_pm25_201303 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201303, mean)

dim(data_pm25_201212)

## narr
n_narr <- dim(LL_narr)[1]
data_2011_10$id_data <- 1:n_narr
data_2011_11$id_data <- 1:n_narr
data_2011_12$id_data <- 1:n_narr

data_2012_01$id_data <- 1:n_narr
data_2012_02$id_data <- 1:n_narr
data_2012_03$id_data <- 1:n_narr
data_2012_04$id_data <- 1:n_narr
data_2012_05$id_data <- 1:n_narr
data_2012_06$id_data <- 1:n_narr
data_2012_07$id_data <- 1:n_narr
data_2012_08$id_data <- 1:n_narr
data_2012_09$id_data <- 1:n_narr
data_2012_10$id_data <- 1:n_narr
data_2012_11$id_data <- 1:n_narr
data_2012_12$id_data <- 1:n_narr

data_2013_01$id_data <- 1:n_narr
data_2013_02$id_data <- 1:n_narr
data_2013_03$id_data <- 1:n_narr

## merging
total_201110 <- merge(data_2011_10, data_pm25_201110, by="id_data"); str(total_201110)
total_201111 <- merge(data_2011_11, data_pm25_201111, by="id_data"); str(total_201111)
total_201112 <- merge(data_2011_12, data_pm25_201112, by="id_data"); str(total_201112)

total_201201 <- merge(data_2012_01, data_pm25_201201, by="id_data"); str(total_201201)
total_201202 <- merge(data_2012_02, data_pm25_201202, by="id_data"); str(total_201202)
total_201203 <- merge(data_2012_03, data_pm25_201203, by="id_data"); str(total_201203)
total_201204 <- merge(data_2012_04, data_pm25_201204, by="id_data"); str(total_201204)
total_201205 <- merge(data_2012_05, data_pm25_201205, by="id_data"); str(total_201205)
total_201206 <- merge(data_2012_06, data_pm25_201206, by="id_data"); str(total_201206)
total_201207 <- merge(data_2012_07, data_pm25_201207, by="id_data"); str(total_201207)
total_201208 <- merge(data_2012_08, data_pm25_201208, by="id_data"); str(total_201208)
total_201209 <- merge(data_2012_09, data_pm25_201209, by="id_data"); str(total_201209)
total_201210 <- merge(data_2012_10, data_pm25_201210, by="id_data"); str(total_201210)
total_201211 <- merge(data_2012_11, data_pm25_201211, by="id_data"); str(total_201211)
total_201212 <- merge(data_2012_12, data_pm25_201212, by="id_data"); str(total_201212)

total_201301 <- merge(data_2013_01, data_pm25_201301, by="id_data"); str(total_201301)
total_201302 <- merge(data_2013_02, data_pm25_201302, by="id_data"); str(total_201302)
total_201303 <- merge(data_2013_03, data_pm25_201303, by="id_data"); str(total_201303)

### Assign climatic region
library(sp)
library(maps)
library(maptools)
library(geosphere)

# The single argument to this function, pointsDF, is a data.frame in which:
#   - column 1 contains the longitude in degrees (negative in the US)
#   - column 2 contains the latitude in degrees

latlong2state <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('state', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}

# Assign States
coords <- as.data.frame(cbind(total_201201$Longitude,total_201201$Latitude))
State <- latlong2state(coords)

CR_NE <- c("maine", "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
           "pennsylvania", "new jersey", "delaware", "maryland")

CR_NE <- c( "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
            "pennsylvania","new jersey", "delaware")



CR <- rep(NA,length(State))
CR[State %in% CR_NE] <- "NE"

total_201110$State <- State; total_201110$CR <- CR; str(total_201110)
total_201111$State <- State; total_201111$CR <- CR; str(total_201111)
total_201112$State <- State; total_201112$CR <- CR; str(total_201112)

total_201201$State <- State; total_201201$CR <- CR; str(total_201201)
total_201202$State <- State; total_201202$CR <- CR; str(total_201202)
total_201203$State <- State; total_201203$CR <- CR; str(total_201203)
total_201204$State <- State; total_201204$CR <- CR; str(total_201204)
total_201205$State <- State; total_201205$CR <- CR; str(total_201205)
total_201206$State <- State; total_201206$CR <- CR; str(total_201206)
total_201207$State <- State; total_201207$CR <- CR; str(total_201207)
total_201208$State <- State; total_201208$CR <- CR; str(total_201208)
total_201209$State <- State; total_201209$CR <- CR; str(total_201209)
total_201210$State <- State; total_201210$CR <- CR; str(total_201210)
total_201211$State <- State; total_201211$CR <- CR; str(total_201211)
total_201212$State <- State; total_201212$CR <- CR; str(total_201212)

total_201301$State <- State; total_201301$CR <- CR; str(total_201301)
total_201302$State <- State; total_201302$CR <- CR; str(total_201302)
total_201303$State <- State; total_201303$CR <- CR; str(total_201303)


## subsetting climatic region
total_201110_ne = subset(total_201110, CR=="NE")
total_201111_ne = subset(total_201111, CR=="NE")
total_201112_ne = subset(total_201112, CR=="NE")

total_201201_ne = subset(total_201201, CR=="NE")
total_201202_ne = subset(total_201202, CR=="NE")
total_201203_ne = subset(total_201203, CR=="NE")
total_201204_ne = subset(total_201204, CR=="NE")
total_201205_ne = subset(total_201205, CR=="NE")
total_201206_ne = subset(total_201206, CR=="NE")
total_201207_ne = subset(total_201207, CR=="NE")
total_201208_ne = subset(total_201208, CR=="NE")
total_201209_ne = subset(total_201209, CR=="NE")
total_201210_ne = subset(total_201210, CR=="NE")
total_201211_ne = subset(total_201211, CR=="NE")
total_201212_ne = subset(total_201212, CR=="NE")

total_201301_ne = subset(total_201301, CR=="NE")
total_201302_ne = subset(total_201302, CR=="NE")
total_201303_ne = subset(total_201303, CR=="NE")


#################### Data Analysis ###########################
grid_size=25
m= grid_size^2
t_length = 6
t=seq(0,1, length.out = t_length)
n_each = dim(total_201210_ne)[1]


################################# 850 ##############################################
### 2011 - 2012 october-mar

Y1_alltime = c(total_201110_ne$HGT850,total_201111_ne$HGT850,total_201112_ne$HGT850,total_201201_ne$HGT850,total_201202_ne$HGT850,total_201203_ne$HGT850)
Y2_alltime = c(total_201110_ne$pm25_monthly_mean,total_201111_ne$pm25_monthly_mean,total_201112_ne$pm25_monthly_mean,total_201201_ne$pm25_monthly_mean,total_201202_ne$pm25_monthly_mean,total_201203_ne$pm25_monthly_mean)

### 2012 april-sept

Y1_alltime = c(total_201204_ne$HGT850,total_201205_ne$HGT850,total_201206_ne$HGT850,total_201207_ne$HGT850,total_201208_ne$HGT850,total_201209_ne$HGT850)
Y2_alltime = c(total_201204_ne$pm25_monthly_mean,total_201205_ne$pm25_monthly_mean,total_201206_ne$pm25_monthly_mean,total_201207_ne$pm25_monthly_mean,total_201208_ne$pm25_monthly_mean,total_201209_ne$pm25_monthly_mean)



Y= cbind(Y1_alltime,Y2_alltime)

df = 7
library(splines)
basis_model = bs(x = t, df= df, degree = 3,
                 Boundary.knots = c(0,1), intercept = TRUE)

vec_n= vector()
for(i in 1:t_length){
  vec_i = c(rep(basis_model[i,], times= n_each))
  vec_n = c(vec_n, vec_i)
}

n= n_each*length(t)

X =  matrix(vec_n, nrow = n, ncol =df, byrow=T)

nu = rep(1/n, n)
mu = rep(1/m, m)


## simulating U from a regular grid
x <- seq(0, 1, length.out = grid_size)
y <- seq(0, 1, length.out = grid_size)
U <- as.matrix(expand.grid(x = x, y = y))

library(matrixcalc) ## for vec
library(Matrix) ## for kronecker product ## for constructing sparse matrix
library(gurobi)
######################## DUAL ############################
obj1 = vec(as.matrix(nu))
obj2 = vec(mu%*%t(nu)%*%X)

model_dual = list()

##objective function
model_dual$obj = c(obj1, obj2)

## A matrix of constraints (lhs)
A1_dual = Matrix(kronecker(X= diag(n), Y= rep(1,m)), sparse = T)
A2_dual = Matrix(kronecker(X= X , Y= diag(m)), sparse = T)


model_dual$A = cbind(A1_dual,A2_dual)

## b vector of constraints (rhs)
rhs_dual = vec(U%*%t(Y))
model_dual$rhs <- rhs_dual

model_dual$vtype <- 'C'  ## Continuous variable type (default)
model_dual$sense <- '>'
model_dual$modelsense <- 'min'

params <- list(OutputFlag=0)
result_dual <- gurobi(model_dual, params)


print(result_dual$objval)
#print(result_dual$x)
psi = result_dual$x[1:n]
vec_b = result_dual$x[-(1:n)]
b = matrix(vec_b, nrow = m, ncol = df) ## b constructed according to x

## obtaining the beta's using numererical differentiation
library(numDeriv)

covariate = predict(basis_model, newx = t[5])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)

finite.differences <- function(x, y) {
  
  n <- dim(x)[1]
  
  # Initialize a vector of length n to enter the derivative approximations
  fdx_u <- vector(length = grid_size)
  fdx_f = rep(0,grid_size); fdx_b = rep(0,grid_size)
  # Iterate through the values using the forward differencing method
  fdx = matrix(0, nrow=grid_size, ncol=grid_size)
  
  for(u in 1:grid_size){
    for (i in 2:grid_size) {
      fdx_f[i-1] <- (y[i-1 +grid_size*(u-1)] - y[i +grid_size*(u-1)]) / (x[i-1+grid_size*(u-1),1] - x[i+grid_size*(u-1),1])
    }
    
    # Iterate through the values using the backward differencing method
    for (i in 1:(grid_size-1)) {
      fdx_b[i+1] <- (y[i +grid_size*(u-1)] - y[i+1 +grid_size*(u-1)]) / (x[i+grid_size*(u-1),1] - x[i+1+grid_size*(u-1),1])
    }
    
    fdx_u <- (fdx_f+fdx_b)/2
    fdx_u[1] <- fdx_f[1]
    fdx_u[grid_size] <- fdx_b[grid_size]
    
    fdx[,u] = fdx_u
  }
  
  return(fdx)
}

beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))

# Here U=(u2,u1) so we take the transpose
beta2Tx = t(finite.differences(U, b2Tx))


## transform uniform reference to spherical uniform
radius = seq(0.00001, 1, length.out = grid_size)

num <- grid_size # number of points you want on the unit circle

pts_all <- t(sapply(1:num,function(p) c(radius[1]*cos(2*p*pi/num),radius[1]*sin(2*p*pi/num))))
for(i in 2:grid_size){
  pts.circle <- t(sapply(1:num,function(p) c(radius[i]*cos(2*p*pi/num),radius[i]*sin(2*p*pi/num))))
  pts_all = rbind(pts_all, pts.circle)
}

U_n = pts_all

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

library(pdist)
## U_n as rows, Y_n as columns
dist_mat = as.matrix(pdist(X=U_n, Y=U))
ranks = rank(c(dist_mat))

##defining cost matrix
cmat = matrix(ranks, nrow=grid_size^2, ncol=grid_size^2)

library(adagio)
matching = assignment(cmat)

perm = matching$perm


library(alphahull)
tau= seq(0.1,0.9,length.out=6)

par(mfrow=c(2,3))
### oct -mar 2011

#oct
covariate = predict(basis_model, newx = t[1])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201110_ne$HGT850,total_201110_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main="Oct 2011", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5 )
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 30)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}

#nov
covariate = predict(basis_model, newx = t[2])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201111_ne$HGT850,total_201111_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main="Nov 2011", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 100)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}

#dec
covariate = predict(basis_model, newx = t[3])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201112_ne$HGT850,total_201112_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main = "Dec 2011", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 150)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}


#jan
covariate = predict(basis_model, newx = t[4])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201201_ne$HGT850,total_201201_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main ="Jan 2012", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 150)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}

#feb
covariate = predict(basis_model, newx = t[5])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201202_ne$HGT850,total_201202_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main ="Feb 2012", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 120)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}

#mar
covariate = predict(basis_model, newx = t[6])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201203_ne$HGT850,total_201203_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main="Mar 2012", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)

for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 40)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}


## connect contours over time i.e. plot estimated spline function for each u
Y1_cap_11_12 = matrix(0, nrow = m, ncol = length(t)); Y2_cap_11_12 = matrix(0, nrow = m, ncol = length(t))
for (i in 1:length(t)) {
  covariate = predict(basis_model, newx = t[i])[1,]
  bTx= apply(b, 1, function(x) t(x)%*%covariate)
  beta1Tx = finite.differences(U, bTx)
  
  fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
  b2Tx = vec(t(fmatrix))
  beta2Tx = t(finite.differences(U, b2Tx))
  
  Y1_cap_11_12[,i] = c(beta1Tx)
  Y2_cap_11_12[,i] = c(beta2Tx)
  
}
par(mfrow=c(1,1))
summary(c(Y1_cap_11_12))
summary(c(Y2_cap_11_12))

### plotting the estimated spline function
library(scatterplot3d)


s=scatterplot3d(Y1_cap_11_12[1,],Y2_cap_11_12[1,],t,type="l",  xlim=c(1350,1600), ylim = c(0,20),
                cex.symbols =0.4, pch = 20, z.ticklabs = NA,
                xlab = "Geopotential height", ylab = "",zlab = "", main = "",
                color ="grey60", cex.lab=1.5)

axis(2, at=seq(0,5,length.out = 6),las=1,tick=F,hadj = 0.5,labels=c("Oct 2011", "Nov 2011", "Dec 2011", "Jan 2012", "Feb 2012", "Mar 2012"), cex.axis=1)



for(i in 2:dim(U)[1]){
  s$points3d(Y1_cap_11_12[i,],Y2_cap_11_12[i,],t, type = "l", cex=0.4,pch=20,col="grey60")
}

s$points3d(Y1_cap_11_12[313,],Y2_cap_11_12[313,],t, type = "l", cex=0.4,pch=20,col="blue")

dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])-1
y <- dims[3]+ 0.08*diff(dims[3:4])-0.3
text(x,y,"PM2.5",srt=10, cex=1.5)




par(mfrow=c(2,3))
### april - sept 2012
#april
covariate = predict(basis_model, newx = t[1])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))
plot(total_201204_ne$HGT850,total_201204_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main="April 2012", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5 )
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 140)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}

#may
covariate = predict(basis_model, newx = t[2])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201205_ne$HGT850,total_201205_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main="May 2012", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 18)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}

#June
covariate = predict(basis_model, newx = t[3])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201206_ne$HGT850,total_201206_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main = "June 2012", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 30)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}


#July
covariate = predict(basis_model, newx = t[4])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201207_ne$HGT850,total_201207_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main ="July 2012", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 50)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}

#Aug
covariate = predict(basis_model, newx = t[5])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201208_ne$HGT850,total_201208_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main ="August 2012", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)
for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 20)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}

#sept
covariate = predict(basis_model, newx = t[6])[1,]
bTx= apply(b, 1, function(x) t(x)%*%covariate)
beta1Tx = finite.differences(U, bTx)

fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
b2Tx = vec(t(fmatrix))
beta2Tx = t(finite.differences(U, b2Tx))
Y_cap =cbind(c(beta1Tx), c(beta2Tx))

plot(total_201209_ne$HGT850,total_201209_ne$pm25_monthly_mean, xlab= "Geopotential Height", ylab = "PM2.5", main="September 2012", pch=20, col="grey60", cex=0.5, cex.lab=1.5, cex.main= 1.5)

for (i in 1:length(tau)) {
  U_number = which(apply(U_n, 1, norm_vec) < tau[i])
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1], Y_tau[,2], 40)
  plot(ahull.obj,lwd=c(1,1,1.5), add=T, wpoints=F)
  
}


## connect contours over time i.e. plot estimated spline function for each u
Y1_cap_11_12 = matrix(0, nrow = m, ncol = length(t)); Y2_cap_11_12 = matrix(0, nrow = m, ncol = length(t))
for (i in 1:length(t)) {
  covariate = predict(basis_model, newx = t[i])[1,]
  bTx= apply(b, 1, function(x) t(x)%*%covariate)
  beta1Tx = finite.differences(U, bTx)
  
  fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
  b2Tx = vec(t(fmatrix))
  beta2Tx = t(finite.differences(U, b2Tx))
  
  Y1_cap_11_12[,i] = c(beta1Tx)
  Y2_cap_11_12[,i] = c(beta2Tx)
  
}
par(mfrow=c(1,2))
### plotting the estimated spline function
library(scatterplot3d)
s=scatterplot3d(Y1_cap_11_12[1,],Y2_cap_11_12[1,],t,type="l",  xlim=c(1350,1600), ylim = c(0,20),
                cex.symbols =0.4, pch = 20, z.ticklabs = NA,
                xlab = "Geopotential height", ylab = "",zlab = "", main = "",
                color ="grey60", cex.lab = 1.5)

axis(2, at=seq(0,5,length.out = 6),las=1,tick=F,hadj = 0.5,labels=c("Apr 2012", "May 2012", "Jun 2012", "Jul 2012", "Aug 2012", "Sep 2012"), cex.axis=1)

for(i in 2:dim(U)[1]){
  s$points3d(Y1_cap_11_12[i,],Y2_cap_11_12[i,],t, type = "l", cex=0.4,pch=20,col="grey60")
}


s$points3d(Y1_cap_11_12[313,],Y2_cap_11_12[313,],t, type = "l", cex=0.4,pch=20,col="blue")

dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])-1
y <- dims[3]+ 0.08*diff(dims[3:4])-0.3
text(x,y,"PM2.5",srt=10, cex=1.5)









