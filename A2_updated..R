########################### 
# HW2: OLS and Probit
# Runling Wu 
# Date: Jan.26/2022
###########################

setwd("~/Dropbox/Duke/Econ 613/Assignment/A2/Data/")
getwd()

############################
# Exercise 1: OLS estimate
############################
rm(list=ls())
ID09 <- read.table("datind2009.csv", sep = ",", header = TRUE)

# Y: wage, X: age of individuals 

# 1.1 correlation between Y and X. 

cor(ID09$age, ID09$wage, use = "complete.obs")

# -0.1788512

# 1.2 Calculate the coefficients on this regression. 

# use the formula of beta 

# drop entire row if wage has missing data
ID09_m <- na.omit(ID09)

x <- as.matrix(ID09_m[, c("age")])
y <- as.matrix(ID09_m[, "wage"])

intercept <- rep(1, nrow(x))

# Coerce to matrices 
Y <- as.matrix(y)
X <- as.matrix(cbind(intercept, x))

beta_hat <- solve(t(X)%*%X)%*%t(X)%*%Y
beta_hat

# intercept 12568.7505 # slope 276.5593

#  checking
reg_ols = lm(wage ~ age, data = ID09_m )
reg_ols

# 1.3a Calculate the standard errors (standard formula of OLS)

# x matrix
nvar = ncol(X)
ndraw = nrow(ID09_m)

# standard errors
eps_hat     = Y - X %*% beta_hat
sigma_hat   = t(eps_hat) %*% eps_hat / (ndraw-nvar-1)

vcov = sigma_hat[1,1] * solve(t(X)%*%X)

# Estimate of standard errors
sd_beta =sqrt(diag(vcov))
sd_beta

# se for intercept 717.21727, se for slope 16.79163

# 1.3b(1) Calculate the standard errors (Using bootstrap)

library(AER)

R    = 49;                      # number of bootstraps
nind = nrow(ID09_m);            # number of individuals
nvar = 2  # number of variables

outs1 = mat.or.vec(R,nvar)
set.seed(123)

for (i in 1:R)
{
  samp     = sample(1:nind,nind,rep=TRUE)
  dat_samp = ID09_m[samp,]
  x0 = as.matrix(cbind(1, dat_samp$age))
  y0 = as.matrix(dat_samp$wage)
  b_beta0 = solve(t(x0)%*%x0)%*%t(x0)%*%y0
  outs1[i,] = b_beta0
}

sd_est1 = apply(outs1,2,sd)
sd_est1

# se for intercept is 566.19408 and se for slope is 15.75011 

# 1.3b(2) Calculate the standard errors (Using bootstrap)

R    = 499;                      # number of bootstraps
nind = nrow(ID09_m);            # number of individuals
nvar = 2  # number of variables

outs2 = mat.or.vec(R,nvar)
set.seed(123)

for (i in 1:R)
{
  samp     = sample(1:nind,nind,rep=TRUE)
  dat_samp = ID09_m[samp,]
  x0 = as.matrix(cbind(1, dat_samp$age))
  y0 = as.matrix(dat_samp$wage)
  b_beta0 = solve(t(x0)%*%x0)%*%t(x0)%*%y0
  outs2[i,] = b_beta0
}

sd_est1 = apply(outs2,2,sd)
sd_est1

# se for intercept is 583.55418   and se for slope is 15.73966. 

############################
# Exercise 2: Detrend Data
############################

# Read all individual datasets from 2005 to 2018. Append all these datasets.

data_path1 <- "~/Dropbox/Duke/Econ 613/Assignment/A2/Data/"
files <- dir(data_path1, pattern = "datind*")

ID0518 <- files %>%
  map(~ paste(data_path1, .x, sep = "")) %>%
  map(read_csv, col_types = cols(.default = "c"), col_select = !c(1)) %>%
  bind_rows

# data tidying work 
# transfer from strings to integers 

cols.num <- c("wage", "age", "year")
ID0518[cols.num] <- sapply(ID0518[cols.num],as.numeric)

sapply(ID0518, class) 

# 2.1 create a categorical variable which bins the age variables in 9 groups 
ID0518$age_group <- as.factor(ifelse(ID0518$age >= 17 & ID0518$age <= 25, '18-25',  
                  ifelse(ID0518$age >= 26 & ID0518$age <= 30, '26-30',
                  ifelse(ID0518$age >= 31 & ID0518$age <= 35, '31-35',
                  ifelse(ID0518$age >= 36 & ID0518$age <= 40, '36-40',
                  ifelse(ID0518$age >= 41 & ID0518$age <= 45, '41-45', 
                  ifelse(ID0518$age >= 46 & ID0518$age <= 50, '46-50', 
                  ifelse(ID0518$age >= 51 & ID0518$age <= 55, '51-55',  
                  ifelse(ID0518$age >= 56 & ID0518$age <= 60, '56-60',  
                  ifelse(ID0518$age >= 60,  '60+',0))))))))))


# 2.2 Plot the wage of each age group across years.

ID0518_u <-ID0518 %>% filter(wage != 0 & !is.na(wage))

# Method 1: use boxplot to plot the full distribution 

ggplot(ID0518_u  %>% filter(wage != 0 & !is.na(wage), wage < 200000), mapping = aes(x= age_group, y= wage, color = factor(year))) + 
  geom_boxplot()

# Method 2: plot the mean wage of different groups as a representative 

# find mean wage for each age for  each year 
meandf <- aggregate(ID0518_m$wage, list(ID0518_m$age_group, ID0518_m$year), mean)

ggplot(data=meandf, mapping = aes(x = Group.1, y = x, group = 1, color = factor(Group.2))) + 
  geom_point() + 
  xlab("Age") + ylab("Mean Wage")

# 2.3 including a time fixed effect, how do the estimated coefficients change?

###############################################################
#   The fixed effect model is as follows:                     #                       
#   Y(i,t) = bX(i,t) + gamma(t) +  e(it)                      #
# identification: using demean to identify beta, assume that  #
# attrition is exogenous, so that we wont have inconsistency  #
# we will drop those observation to make it a balanced panel  #
###############################################################

# Method 1: try to include the time in a matrix form unsucessful... 

nind = nrow(ID0518_u)
nper = 14 

IDage <- ID0518_u %>% select(age,year,idind)
IDwage <- ID0518_u %>% select(wage,year,idind)

# keep bl
# reshape long to wide 
IDage_w <- reshape(IDage, idvar = "idind", timevar = "year", 
                   v.names = c("age"), direction = "wide")

x = cbind(ID0518$age)

xMat = matrix(runif(nind*nper, 0,10),nind,nper)
eps = matrix(rnorm(nind*nper, 0, 1),nind,nper)
beta = runif(1)
fixE = runif(nind,0,0.25)
yMat = mat.or.vec(nind,nper)

for (i in 1:nper)
{
  yMat[,i] = fixE + beta*xMat[,i] + eps[,i]
}

reg1 = lm(c(yMat)~c(xMat))
summary(reg1)

# Method 2: Use year dummies 

library(plm)

ols_fetime <- plm(wage ~ age + 1, data=ID0519_m, index = "year", 
                  effect = "time", model = "within")
summary(ols_fetime) 

ols_fetime <- plm(wage ~ age + 1 + factor(year), data=ID0519_m)

ols_fetime <- plm(wage ~ age + factor(year), data=ID0518_m)
summary(ols_fetime)
# age 326.3130

texreg(list(ols_fetime, reg_ols)) 

########################################
# Exercise 3: Numerical Optimaization 
########################################

# data from 2007 
ID07 <- read.table("datind2007.csv", sep = ",", header = TRUE)

# data tidying work 
# transfer from strings to integers 

cols.num <- c("wage", "age")
ID07[cols.num] <- sapply(ID07[cols.num],as.numeric)

sapply(ID07, class) 

# 3.1 Exclude all individuals who are inactive.
em_ac <- ID07 %>% filter (empstat != "Inactive") 

# 3.2 likelihood of the probit of being employed.

em_ac$employed = ifelse(em_ac$empstat=="Employed", 1, 0)

X <- em_ac$age
Y <- em_ac$employed
K <- ncol(X)

flike = function(par,x1,yvar,model)
{
  xbeta           = par[1] + par[2]*x1  
  if (model=="plm")
  {
    xbeta           = par[1] + par[2]*x1 
    like            = log(dnorm(yvar-xbeta))
  }
  if (model=="probit")
  {
    pr              = pnorm(xbeta) # probit
    pr[pr>0.999999] = 0.999999
    pr[pr<0.000001] = 0.000001
    like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  }else 
  {
    pr              = exp(xbeta)/(1+exp(xbeta)) #logit
    pr[pr>0.999999] = 0.999999
    pr[pr<0.000001] = 0.000001
    like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  }
  return(-sum(like))
}

# test if the function is correct
probit = glm(employed ~ age, family = binomial(link = "probit"), data = em_ac)
summary(probit)

test_par = probit$coefficients
flike(test_par,X,Y, model = "probit") # 6582.155 
logLik(probit) # 'log Lik.' -6582.155 (df=2) 

summary(em_ac$age)

# 3.3 Optimize the model and interpret the coefficients.
# need to do many times to make sure it convergenes 
start = runif(2)
res  = optim(start,fn=flike,method="BFGS",control=list(trace=6,REPORT=1,maxit=1000),
             x1=em_ac$age,yvar=em_ac$employed,model="probit",hessian=TRUE)

# final  value 6582.155339  matched! 

fisher_info = solve(res$hessian)       # standard formula is -res$hessian but flike is return -like
prop_sigma  = sqrt(diag(fisher_info))
prop_sigma

est = cbind(probit$coefficients,res$par)
colnames(est) = c("R: GLM : est","R: own : est")
est
xtable(est)
# R: GLM : est R: own : est
# [1,]    3.8291873   3.82928037
# [2,]   -0.0678642  -0.06786624

# 3.4 we cannot include wage as a determinant of labor see in pdf

# a probit model is utilized to identify the factors affecting the labor market 
# participation. 

###############################
# Exercise 4: Discrete Choice
###############################

# Read all individual datasets from 2005 to 2015. Append all these datasets.

data_path2 <- "~/Dropbox/Duke/Econ 613/Assignment/A2/Data/0515/"
files <- dir(data_path2, pattern = "datind*")

ID0515 <- files %>%
  map(~ paste(data_path2, .x, sep = "")) %>%
  map(read_csv, col_types = cols(.default = "c"), col_select = !c(1)) %>%
  bind_rows

# data tidying work 
# transfer from strings to integers 

cols.num <- c("wage", "age")
ID0515[cols.num] <- sapply(ID0515[cols.num],as.numeric)

sapply(ID0515, class) 

# 4.1 Exclude all individuals who are inactive.

dis_ch <- ID0515 %>% filter (empstat != "Inactive") 

# 4.2 Write and optimize the probit, logit, and the linear probability models.

dis_ch$employed = ifelse(dis_ch$empstat=="Employed", 1, 0)

################################################
# 4.2 probit, logit, lpm  with time FE effects 
################################################

# recall 3.2 we have written probit, logit and lpm in one function.
# we will revise by adding year dummies on top of that. 
library(fastDummies)

dis_ch_d <- dummy_cols(dis_ch, select_columns = "year")

flike_FE = function(par,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,yvar,model)
{
  xbeta           = par[1]*x1+par[2]*x2+par[3]*x3+par[4]*x4+par[5]*x5+par[6]*x6+
    par[7]*x7+par[8]*x8+par[9]*x9+par[10]*x10+par[11]*x11+par[12]*x12
  if (model=="plm")
  {
    xbeta           = par[1]*x1+par[2]*x2+par[3]*x3+par[4]*x4+par[5]*x5+par[6]*x6+
      par[7]*x7+par[8]*x8+par[9]*x9+par[10]*x10+par[11]*x11+par[12]*x12
    like            = log(dnorm(yvar-xbeta))
  }
  if (model=="probit")
  {
    pr              = pnorm(xbeta) # probit
    pr[pr>0.999999] = 0.999999
    pr[pr<0.000001] = 0.000001
    like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  }else 
  {
    pr              = exp(xbeta)/(1+exp(xbeta)) #logit
    pr[pr>0.999999] = 0.999999
    pr[pr<0.000001] = 0.000001
    like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  }
  return(-sum(like))
}

##################################
# 4.2 (1) probit with time FE effects 
###################################

start2=runif(12)
res_probit  = optim(start2,fn=flike_FE,method="BFGS",control=list(trace=6,maxit=1000)
                    ,x1=dis_ch_d$age,yvar=dis_ch_d$employed,x2=dis_ch_d$year_2005,
                    x3=dis_ch_d$year_2006, x4=dis_ch_d$year_2007, 
                    x5=dis_ch_d$year_2008, x6=dis_ch_d$year_2009,
                    x7=dis_ch_d$year_2010, x8=dis_ch_d$year_2011,
                    x9=dis_ch_d$year_2012, x10=dis_ch_d$year_2013,
                    x11=dis_ch_d$year_2014, x12=dis_ch_d$year_2015,
                    model="probit",hessian=TRUE) 

# final  value 79892.070811 converged

# test if the function is correct 
reg1_probit = glm(employed~age + year_2005+year_2006+year_2007+year_2008+year_2009+
                    year_2010+year_2011+year_2012+year_2013+year_2014+year_2015-1, 
                  family=binomial(link = "probit"),
                  data = dis_ch_d)
test_reg1pr = reg1_probit$coefficients
logLik(reg1_probit) 
# log Lik.' -79892.06 (df=12) matched!! 

##################################
# 4.2 logit with time FE effects 
###################################

res_logit   = optim(start2,fn=flike_FE,method="BFGS",control=list(trace=6,maxit=1000)
                    ,x1=dis_ch_d$age,yvar=dis_ch_d$employed,x2=dis_ch_d$year_2005,
                    x3=dis_ch_d$year_2006, x4=dis_ch_d$year_2007, 
                    x5=dis_ch_d$year_2008, x6=dis_ch_d$year_2009,
                    x7=dis_ch_d$year_2010, x8=dis_ch_d$year_2011,
                    x9=dis_ch_d$year_2012, x10=dis_ch_d$year_2013,
                    x11=dis_ch_d$year_2014, x12=dis_ch_d$year_2015,
                    model="logit",hessian=TRUE)

# final  value 76900.350116 converged

# test if the function is correct 
reg1_logit = glm(employed~age + year_2005+year_2006+year_2007+year_2008+year_2009+
                    year_2010+year_2011+year_2012+year_2013+year_2014+year_2015-1, 
                  family=binomial(link = "logit"),
                  data = dis_ch_d)
test_reg1log = reg1_logit$coefficients
logLik(reg1_logit) 
#'log Lik.' -76900.35 (df=12) matched!! 

##################################
# 4.2(3) lpm with time FE effects 
###################################

res_plm     = optim(start2,fn=flike_FE,method="BFGS",control=list(trace=6,maxit=1000)
                          ,x1=dis_ch_d$age,yvar=dis_ch_d$employed,x2=dis_ch_d$year_2005,
                          x3=dis_ch_d$year_2006, x4=dis_ch_d$year_2007, 
                          x5=dis_ch_d$year_2008, x6=dis_ch_d$year_2009,
                          x7=dis_ch_d$year_2010, x8=dis_ch_d$year_2011,
                          x9=dis_ch_d$year_2012, x10=dis_ch_d$year_2013,
                          x11=dis_ch_d$year_2014, x12=dis_ch_d$year_2015,
                          ,model="plm",hessian=TRUE)

# may not work welll 

ntry = 100
out_ex4plm = mat.or.vec(ntry,12)
for (i0 in 1:ntry)
{
  start_ex4plm = runif(12)
  res_ex4plm      = optim(start_ex4plm, fn=flike_FE,method="BFGS",
                          control=list(trace=6),x1=dis_ch_d$age,
                          yvar=dis_ch_d$employed,x2=dis_ch_d$year_2005,
                          x3=dis_ch_d$year_2006, x4=dis_ch_d$year_2007, 
                          x5=dis_ch_d$year_2008, x6=dis_ch_d$year_2009,
                          x7=dis_ch_d$year_2010, x8=dis_ch_d$year_2011,
                          x9=dis_ch_d$year_2012, x10=dis_ch_d$year_2013,
                          x11=dis_ch_d$year_2014, x12=dis_ch_d$year_2015,
                          model="plm",hessian=TRUE)
  out_ex4plm[i0,] = res_ex4plm$par 
}

# final  value  76900.351580 converged 

# test if the function is correct 
reg1_lpm = lm(employed~age + year_2005+year_2006+year_2007+year_2008+year_2009+
                   year_2010+year_2011+year_2012+year_2013+year_2014+year_2015-1, 
                 data = dis_ch_d)
summary(reg1_lpm)
test_reg1lpm = reg1_lpm$coefficients
logLik(reg1_lpm) 

#''log Lik.' -76900.35 (df=12) 


ntry = 100
out_ex4pr = mat.or.vec(ntry,12)
for (i0 in 1:ntry)
{
  start_ex4pr = runif(12)
  res_ex4pr      = optim(start_ex4pr, fn=flike_FE,method="BFGS",
                        control=list(trace=6),x1=dis_ch_d$age,
                        yvar=dis_ch_d$employed,x2=dis_ch_d$year_2005,
                        x3=dis_ch_d$year_2006, x4=dis_ch_d$year_2007, 
                        x5=dis_ch_d$year_2008, x6=dis_ch_d$year_2009,
                        x7=dis_ch_d$year_2010, x8=dis_ch_d$year_2011,
                        x9=dis_ch_d$year_2012, x10=dis_ch_d$year_2013,
                        x11=dis_ch_d$year_2014, x12=dis_ch_d$year_2015,
                        model="probit",hessian=TRUE)
  out_ex4pr[i0,] = res_ex4pr$par 
}


ntry = 100
out_ex4lg = mat.or.vec(ntry,12)
for (i0 in 1:ntry)
{
  start_ex4lg = runif(12)
  res_ex4lg      = optim(start_ex4lg, fn=flike_FE,method="BFGS",
                        control=list(trace=6),x1=dis_ch_d$age,
                        yvar=dis_ch_d$employed,x2=dis_ch_d$year_2005,
                        x3=dis_ch_d$year_2006, x4=dis_ch_d$year_2007, 
                        x5=dis_ch_d$year_2008, x6=dis_ch_d$year_2009,
                        x7=dis_ch_d$year_2010, x8=dis_ch_d$year_2011,
                        x9=dis_ch_d$year_2012, x10=dis_ch_d$year_2013,
                        x11=dis_ch_d$year_2014, x12=dis_ch_d$year_2015,
                        model="logit",hessian=TRUE)
  out_ex4lg[i0,] = res_ex4lg$par 
}

###########################################
# 4.3 compare the estimated coefficients 
############################################ 

sink("discrete_choice_model.tex")
texreg(list(reg1_logit, reg1_probit, reg1_lpm), dcolumn = FALSE, 
       booktabs = TRUE,use.packages = FALSE, caption = "Models",float.pos = "hb!")
sink()

###############################
# Exercise 5: Marginal Effects 
###############################

library(mfx)

# 5.1 Compute the marginal effect of the previous probit and logit models.
# Marginal effect at the mean

flike_FE5 = function(par,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,yvar,model)
{
  xbeta           = par[1]+par[2]*x2+par[3]*x3+par[4]*x4+par[5]*x5+par[6]*x6+
    par[7]*x7+par[8]*x8+par[9]*x9+par[10]*x10+par[11]*x11+par[12]*x12
  if (model=="plm")
  {
    xbeta           = par[1]+par[2]*x2+par[3]*x3+par[4]*x4+par[5]*x5+par[6]*x6+
      par[7]*x7+par[8]*x8+par[9]*x9+par[10]*x10+par[11]*x11+par[12]*x12
    like            = log(dnorm(yvar-xbeta))
  }
  if (model=="probit")
  {
    pr              = pnorm(xbeta) # probit
    pr[pr>0.999999] = 0.999999
    pr[pr<0.000001] = 0.000001
    like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  }else 
  {
    pr              = exp(xbeta)/(1+exp(xbeta)) #logit
    pr[pr>0.999999] = 0.999999
    pr[pr<0.000001] = 0.000001
    like           = yvar*log(pr) + (1-yvar)*log(1-pr)
  }
  return(-sum(like))
}

start5=runif(12)
reg5_probit  = optim(start5,fn=flike_FE5,method="BFGS",control=list(trace=6,maxit=1000)
                    ,x2=dis_ch_d$age,yvar=dis_ch_d$employed,x3=dis_ch_d$year_2005,
                    x4=dis_ch_d$year_2006, x5=dis_ch_d$year_2007, 
                    x6=dis_ch_d$year_2008, x7=dis_ch_d$year_2009,
                    x8=dis_ch_d$year_2010, x9=dis_ch_d$year_2011,
                    x10=dis_ch_d$year_2012, x11=dis_ch_d$year_2013,
                    x12=dis_ch_d$year_2014,
                    model="probit",hessian=TRUE) 
# final  value 79892.071797 converged

reg5_logit  = optim(start5,fn=flike_FE5,method="BFGS",control=list(trace=6,maxit=1000)
                     ,x2=dis_ch_d$age,yvar=dis_ch_d$employed,x3=dis_ch_d$year_2005,
                     x4=dis_ch_d$year_2006, x5=dis_ch_d$year_2007, 
                     x6=dis_ch_d$year_2008, x7=dis_ch_d$year_2009,
                     x8=dis_ch_d$year_2010, x9=dis_ch_d$year_2011,
                     x10=dis_ch_d$year_2012, x11=dis_ch_d$year_2013,
                     x12=dis_ch_d$year_2014,
                     model="logit",hessian=TRUE) 

# final value  76900.355435 converged

margeff_eval = function(par,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,yvar,model)
{
  xbeta   = par[1]+ par[2]*x2+par[3]*x3+par[4]*x4+par[5]*x5
  +par[6]*x6+par[7]*x7+par[8]*x8+par[9]*x9+par[10]*x10+par[11]*x11+par[12]*x12
  pr      = ifelse(model=="probit",
                   dnorm(mean(xbeta)),dlogis(mean(xbeta)))
  marg            = par[-1]*pr
  return(marg)
}

 margeff_eval(reg5_probit$par,x2=dis_ch_d$age,
             yvar=dis_ch_d$employed,x3=dis_ch_d$year_2005,
             x4=dis_ch_d$year_2006, x5=dis_ch_d$year_2007, 
             x6=dis_ch_d$year_2008, x7=dis_ch_d$year_2009,
             x8=dis_ch_d$year_2010, x9=dis_ch_d$year_2011,
             x10=dis_ch_d$year_2012, x11=dis_ch_d$year_2013,
             x12=dis_ch_d$year_2014
             , model = "probit")


 # -0.023976957 -0.000598462  0.002701356  0.013084789  0.015356954 -0.003724729
#  -0.001630270  0.008565833  0.003505702 -0.007080582  0.004300036
 
 # checking 
 probitmfx(employed~age+year_2005+year_2006+year_2007+year_2008+year_2009
           +year_2010+year_2011+year_2012+year_2013+year_2014, 
           data=dis_ch_d, atmean= TRUE)
 
 # very similar results 
 
margeff_eval(reg5_logit$par,x2=dis_ch_d$age,
             yvar=dis_ch_d$employed,x3=dis_ch_d$year_2005,
             x4=dis_ch_d$year_2006, x5=dis_ch_d$year_2007, 
             x6=dis_ch_d$year_2008, x7=dis_ch_d$year_2009,
             x8=dis_ch_d$year_2010, x9=dis_ch_d$year_2011,
             x10=dis_ch_d$year_2012, x11=dis_ch_d$year_2013,
             x12=dis_ch_d$year_2014,
             model = "logit")

# -0.0270037760 -0.0123718775 -0.0077177930  0.0013091766  0.0053754071
# -0.0112710272 -0.0093489434  0.0028369538 -0.0002095949 -0.0095949843
# 0.0024206837

# checking 
logitmfx(employed~age+year_2005+year_2006+year_2007+year_2008+year_2009
          +year_2010+year_2011+year_2012+year_2013+year_2014, 
          data=dis_ch_d, atmean= TRUE)

# very similar results 

#=====================================================================
# 5.2 Construct the standard errors of the marginal effects (boostrap)
#=====================================================================
strat=runif(13)

x2=dis_ch_d$age
ydum=dis_ch_d$employed
x3=dis_ch_d$year_2005
x4=dis_ch_d$year_2006
x5=dis_ch_d$year_2007
x6=dis_ch_d$year_2008
x7=dis_ch_d$year_2009
x8=dis_ch_d$year_2010
x9=dis_ch_d$year_2011
x10=dis_ch_d$year_2012
x11=dis_ch_d$year_2013
x12=dis_ch_d$year_2014

margeff_eval(reg5_probit$par,x2=dis_ch_d$age,
             yvar=dis_ch_d$employed,x3=dis_ch_d$year_2005,
             x4=dis_ch_d$year_2006, x5=dis_ch_d$year_2007, 
             x6=dis_ch_d$year_2008, x7=dis_ch_d$year_2009,
             x8=dis_ch_d$year_2010, x9=dis_ch_d$year_2011,
             x10=dis_ch_d$year_2012, x11=dis_ch_d$year_2013,
             x12=dis_ch_d$year_2014
             , model = "probit")

nboot    = 4999;                      # number of bootstraps
nind = nrow(dis_ch_d);            # number of individuals

outs3 = mat.or.vec(11,R)

for (i in 1:R)
{
  indices = sample(1:nind,nind,replace=T)
  outs3[,i] = margeff_eval(reg5_probit$par, x2[indices],x3[indices], 
                           x4[indices],x5[indices],x6[indices],
                           x7[indices],x8[indices],x9[indices],x10[indices],
                           x11[indices],x12[indices],ydum[indices],"probit")
}

sd_est1 = apply(outs3,1,sd)
sd_est1

# checking 
probitmfx(employed~age+year_2005+year_2006+year_2007+year_2008+year_2009
          +year_2010+year_2011+year_2012+year_2013+year_2014, 
         data=dis_ch_d,atmean=TRUE)


margeff_eval(reg5_logit$par,x2=dis_ch_d$age,
             yvar=dis_ch_d$employed,x3=dis_ch_d$year_2005,
             x4=dis_ch_d$year_2006, x5=dis_ch_d$year_2007, 
             x6=dis_ch_d$year_2008, x7=dis_ch_d$year_2009,
             x8=dis_ch_d$year_2010, x9=dis_ch_d$year_2011,
             x10=dis_ch_d$year_2012, x11=dis_ch_d$year_2013,
             x12=dis_ch_d$year_2014
             , model = "logit")

nboot    = 4999;                      # number of bootstraps
nind = nrow(dis_ch_d);            # number of individuals

outs4 = mat.or.vec(11,R)

for (i in 1:R)
{
  indices = sample(1:nind,nind,replace=T)
  outs4[,i] = margeff_eval(reg5_probit$par, x2[indices],x3[indices], 
                           x4[indices],x5[indices],x6[indices],
                           x7[indices],x8[indices],x9[indices],x10[indices],
                           x11[indices],x12[indices],ydum[indices],"logit")
}

sd_est2 = apply(outs3,1,sd)
sd_est2

# checking 
logitmfx(employed~age+year_2005+year_2006+year_2007+year_2008+year_2009
         +year_2010+year_2011+year_2012+year_2013+year_2014, 
         data=dis_ch_d, atmean=TRUE)

#======================================end=======================================
