rm(list = ls())
food = read.csv("common_household_food.txt",header = T, row.names = NULL)
food = food[-962,]
food1 = food
rownames(food1) = food$Food
food1 = food1[,-1]
###
X = subset(food1,select = -KCal)

qqnorm(food1$KCal)
qqline(food1$KCal)

sum(food1$KCal == 0,na.rm = TRUE)

sum(is.na(food1$KCal))

###--- multicollinearity
panel.hist_line = function(x,...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2],0,1.5))
  h = hist(x,plot = FALSE)
  breaks = h$breaks;nB = length(breaks)
  y = h$counts;y = y/max(y)
  rect(breaks[-nB],0,breaks[-1],y,col = "cyan")
}
panel_cor = function(x,y,digits = 2,prefix = "",cex.cor,...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0,1,0,1))
  r = abs(cor(x,y))
  txt = format(c(r,0.123456789),digits = digits)[1]
  txt = paste0(prefix,txt)
  if(missing(cex.cor)) cex.cor = 0.8/strwidth(txt)
  text(0.5,0.5,txt,cex = 5)
}

pairs(food1[,c("Fat","SatFat","MonoUnSatFat","PolyUnSatFat")],upper.panel = panel.smooth,lower.panel = panel_cor,
      cex = 1.5,pch = 20,col = "dodgerblue3", bg = "navy blue",
      diag.panel = panel.hist_line, cex.labels = 2,font.labels = 2)

pairs(food1[,c("Ca","P","Fe","K","Na")],upper.panel = panel.smooth,lower.panel = panel_cor,
      cex = 1.5,pch = 20,col = "dodgerblue3", bg = "navy blue",
      diag.panel = panel.hist_line, cex.labels = 2,font.labels = 2)
# pairs(food1[,c("VitaA.IU.","VitaA.RE.","Thiamin","Riboflavin","Niacin","VitaC")],upper.panel = panel.smooth,lower.panel = panel_cor,
#       cex = 1.5,pch = 20,col = "dodgerblue3", bg = "navy blue",
#       diag.panel = panel.hist_line, cex.labels = 2,font.labels = 2)

mod = lm(KCal ~ ., data = food1)
library(car)
vif = vif(mod)
mean(vif)
max(vif)

summary(mod)
anova(mod,test = "Chisq")
summary(aov(mod))
Anova(mod,type = "II")

anova(mod)

###########  b
layout(1)
model1 = lm(KCal~Weight + Protein + Fat + Carb + Ca + P + Fe + K + Thiamin,data = food1)
par(cex = 1)
plot(model1$fitted.values,model1$residuals,cex.lab = 1.5)
plot(model1)

qqnorm(food1$KCal);qqline(food1$KCal)
## fit vs res ##
library(MASS)
plot(model1$fitted.values ,stdres(model1),xlab="Fitted values",
     ylab="Standardized residuals",pch = 16,
     main="Residual vs Fitted",cex.lab = 1.5) 
abline(h=0);abline(h=3,lty=2);abline(h=-3,lty=2)
## res QQ ##
qqnorm(stdres(model1),ylab="Standardized residuals",
       ylim=c(-3,3),cex.lab = 1.5);qqline(studres(model1))

plot(resid(model1),cex.lab = 1.5)

## outlier
plot(model1$fitted.values ,studres(model1),xlab="Fitted values",
     ylab="Studentized residuals",pch = 16,
     main="Residual vs Fitted") 
abline(h=0);abline(h=4,lty=2);abline(h=-4,lty=2)
text(model1$fitted.values[studres(model1)>4],
     studres(model1)[studres(model1)>4],
     labels = which(studres(model1)>4),pos = 2)
text(model1$fitted.values[studres(model1)<=-4],
     studres(model1)[studres(model1)<=-4],
     labels = which(studres(model1)<=-4),pos = 2)

which(abs(studres(model1))>=4)
outliers = seq(961)[abs(studres(model1))>=4]   ###  potential outliers
## DFFITS

lm.reg.dffits = dffits(model1)
plot(lm.reg.dffits, type = "h", ylab = "DFFITS", ylim = c(-3,3)) 
text(seq(961)[lm.reg.dffits> 2*sqrt(9/961)],lm.reg.dffits[lm.reg.dffits> 2*sqrt(9/961)], 
     labels = seq(961)[lm.reg.dffits> 2*sqrt(9/961)], cex = 0.8, pos = 2)
abline(h = c(-1,-2*sqrt(9/961), 0, 2*sqrt(9/961), 1), lty = 2) # specify your own h 

seq(961)[abs(lm.reg.dffits)> 2*sqrt(9/961)]

## cook's distance
lm.reg.cooksD = cooks.distance(model1)
plot(lm.reg.cooksD, type = "h", ylab="Cook's Distance",ylim=c(0,2)) 
text(lm.reg.cooksD, labels = index, cex = 1)
abline(h=qf(0.50,4,15), lty=2) #check whether D_i > f_0.5,p,n-p

seq(961)[lm.reg.cooksD>qf(0.5,9,952)]


influ = intersect(seq(961)[abs(studres(model1))>=4],union(seq(961)[abs(lm.reg.dffits)> 2*sqrt(9/961)],seq(961)[lm.reg.cooksD>qf(0.5,9,952)]))
setdiff(outliers,influ)

### transformation
library(MASS)
food2 = food1[-8,]
food2$KCal = food2$KCal+1
bc = boxcox(model1,data = food2[food2$KCal!=0,],lambda = seq(-2,2,0.01))
bc$x[which.max(bc$y)]

###-------------c   variables selection
rm(list = ls())
food = read.csv("common_household_food.txt",header = T, row.names = NULL)
food = food[-962,]
food1 = food
rownames(food1) = food$Food
food1 = food1[,-1]
food2 = food1[-8,]
X = subset(food2,select = -KCal)

library(leaps)
source("myregsub.R")
N = 1
my = my.regsub(X,y=food2$KCal,nbest=7,nvmax = 8,method="exhaustive")
my1 = my[43:49,]
library(stargazer)

stargazer(my1)


step(mod,direction = "backward",k = log(961))
step(mod,direction = "both",k = log(961))


modelf = lm(KCal ~ Fat + Protein + Carb + Chol + K + Thiamin + Weight,data = food2)
summary(modelf)
anova(modelf)

292

library(caret)
folds = createFolds(seq(960), k = 10, list = TRUE)
MSPE = seq(10)
for(i in 1:10)
{
  model =  lm(KCal ~ Fat + Protein + Carb + Chol + K + Thiamin + Weight,data = food2[-folds[[i]],])
  y_hat = predict(model, newdata = as.data.frame(X[folds[[i]],]))
  y_true = food2$KCal[folds[[i]]]
  MSPE[i] = sum((y_hat - y_true)^2)/96
}
MSPE
sum(MSPE>292)
mean(MSPE)

### d prediction
new_X = data.frame(Fat=1.5,
                   Protein = 3,
                   Carb = 26,
                   Chol = 0,
                   Thiamin = 0,
                   K = 95,
                   Weight = 33)
y_pre = predict(modelf,newdata = new_X,interval = "prediction",level = 0.95)
y_pre

####--------------------------------------------------------------------------------------
rm(list = ls())
aud = read.csv("audibility.csv")
aud1 = aud
aud1$Rayleigh[aud$Rayleigh < 0.05] = "P"
aud1$Rayleigh[aud$Rayleigh >= 0.05] = "N"
aud1$Ftest[aud$Ftest < 0.05] = "P"
aud1$Ftest[aud$Ftest >= 0.05] = "N"
aud1$SL[aud$SL < 0] = "N"
aud1$SL[aud$SL >= 0] = "P"

table(aud1$SL,aud1$Rayleigh)
table(aud1$SL,aud1$Ftest)
(107+355)/672
107/112
355/560

####  accuracy
aud2 = aud1[aud1$SPL == 20,]
#table(aud2$SL,aud2$Rayleigh)
table(aud2$SL,aud2$Ftest)




##  random effect
library(lme4)

## Default REML estimation
aud1$acc = rep(0,672)
aud1$acc[aud1$SL == aud1$Ftest & aud1$SL == aud1$Rayleigh] = 1

fit1 = glmer(acc ~  as.factor(Carrier) + (1|SPL) + (1|Participant),family = binomial,data = aud1)
summary(fit1)

aud2 = aud1
aud2$fre = rep(0,672)
aud2$fre[aud1$Carrier == 'a_F1' | aud1$Carrier == 'i_F1' | aud1$Carrier == 'u_F1'] = 'low'
aud2$fre[aud1$Carrier == 'a_F2' | aud1$Carrier == 'i_F2' | aud1$Carrier == 'u_F2'] = 'mid'
aud2$fre[aud1$Carrier == 's' | aud1$Carrier == 'sh'] = 'high'

fit2 = glmer(acc ~  as.factor(fre) + (1|SPL) + (1|Participant),family = binomial,data = aud2)
summary(fit2)

aud2 = aud1[aud1$fre == 'low',]
mean(aud2$Rayleigh == aud2$SL)
aud2 = aud1[aud1$fre == 'mid',]
mean(aud2$Rayleigh == aud2$SL)
aud2 = aud1[aud1$fre == 'high',]
mean(aud2$Rayleigh == aud2$SL)

aud2 = aud1[aud1$fre == 'low',]
mean(aud2$Ftest == aud2$SL)
aud2 = aud1[aud1$fre == 'mid',]
mean(aud2$Ftest == aud2$SL)
aud2 = aud1[aud1$fre == 'high',]
mean(aud2$Ftest == aud2$SL)


accu = data.frame(acc = c( 0.6388889,0.6785714,0.8928571,0.5873016,0.6746032, 0.8571429),
                  test = c(rep("Rayleigh",3),rep('Ftest',3)),
                  fre = rep(c('low','mid','high'),2))
accu
library(ggplot2)
ggplot(accu,aes(x = fre,y = acc,color = test,pch = test)) + geom_point(cex = 5) +
          theme(axis.text=element_text(size=20),
          axis.title=element_text(size=26,face="bold"))

###--------------  c
aud2 = aud
aud2$Rayleigh[aud$Rayleigh < 0.05] = "P"
aud2$Rayleigh[aud$Rayleigh >= 0.05] = "N"
aud2$Ftest[aud$Ftest < 0.05] = "P"
aud2$Ftest[aud$Ftest >= 0.05] = "N"

fit3 = glmer(as.factor(Rayleigh) ~  SL + (1|Participant),data = aud2,family = binomial('logit'))
summary(fit3)
coe = coef(fit3)
apply(coe$Participant,2,mean)

ran = ranef(fit3)
ran = ran$Participant
ran = ran$`(Intercept)`

new = data.frame(SL = seq(5,7,by = 0.1),Participant = user[1])
new = data.frame(SL = seq(5,7,by = 0.1),Participant = 'ANH4206')

predict(fit3,new,type = "response")
li = list()

for(i in 1:21)
{
  user = unique(aud2$Participant)
  new = data.frame(SL = seq(-10,50,by = 0.01),Participant = user[i])
  pre = predict(fit3,new,type = "response")
  li[[i]] = new$SL[pre>0.4998 & pre < 0.5002]
  SL = max(unlist(li))
}


## differ in carrier

fit4 = glmer(as.factor(Rayleigh) ~  SL + (1|Carrier) + (1|Participant),data = aud2,family = binomial('logit'))

max(unlist(li))

SL = rep(0,8)
for(j in 1:8){
  cari = unique(aud1$Carrier)
  
  for(i in 1:21)
  {
    user = unique(aud2$Participant)
    new = data.frame(SL = seq(-10,50,by = 0.01),Carrier = cari[j],Participant = user[i])
    pre = predict(fit4,new,type = "response")
    li[[i]] = new$SL[pre>0.4998 & pre < 0.5002]
    SL[j] = max(unlist(li))
  }
  
}




