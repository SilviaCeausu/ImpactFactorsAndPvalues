rm(list=ls())


#Necessary R packages
library("nlme")
library("arm")
library("RCurl")

#Download of p-value data from the figshare website. 
#There are one file for precise data and one file for inequality data (e.g., p<0.05) for each journal and for each year for which that journal was examined. 
url<-"https://ndownloader.figshare.com/articles/5786838/versions/1"

path1 <- tempfile(fileext = ".zip")
if (file.exists(path1))  'file alredy exists' else download.file(url, path1, mode="wb")
data<-unzip(zipfile = path1,exdir = tempdir())

all.files<-sapply(X=data, FUN=function(X) read.csv(X, header = FALSE, stringsAsFactors=FALSE))


file.names<-make.names(sapply(X=data, FUN=function(X) gsub(pattern = "*.csv$", 
                                                           replacement = "", 
                                                           x=tail(unlist(strsplit(X, "/")), 1))))
names(all.files)<-file.names

#Download of summary data for all journals examined. 
urlIF<-"https://ndownloader.figshare.com/files/10934903"

path2 <- tempfile(fileext = ".zip")
if (file.exists(path2))  'file alredy exists' else asd<-download.file(urlIF, path2, mode="wb")
IFvalues<-read.csv2(path2)

#################################################################################################
###PRECISE P-VALUES##############################################################################
#################################################################################################


#Selection of precise values.

precise.data<-unlist(lapply(all.files, is.numeric))

myfiles.precise<-all.files[precise.data]

#Organize data for modeling
organise.data<-function(files.list=myfiles.precise, if.data=IFvalues){
  data<-data.frame()
  for(i in 1:nrow(if.data)){
    no<-grep(if.data[i,3], as.character(names(files.list)))
    data<-rbind(data,
                cbind(Pvalues=files.list[[no]], if.data[i,], row.names = NULL))
  }
  return(data)
}

mydata<-organise.data()

#add all the relevant information
data.to.model.func<-function(data, Pcolumn=1, Journalcolumn=2, 
                            AcronymColumn=4, IFcolumn=5) {
  acronyms<-as.character(unique(data[,AcronymColumn]))
  per05<-c()    #percentage of precise p-values <0.05
  per05no<-c()  #number of precise instances of p-values<0.05
  TotCounts<-c()  #all precise p-values reported in one journal, in one year
  IF<-c()         #the impact factor of the journal that year
  journal<-c()    #the accronym of the journal
  year<-c()       #the year
  for (i in 1:length(acronyms)){
    rows.i<-data[,AcronymColumn]==acronyms[i]
    pvalues.i<-data[rows.i, Pcolumn]
    TotCounts[i]<-length(pvalues.i)
    per05no[i]<-length(pvalues.i[pvalues.i<0.05])
    per05[i]<-per05no[i]/TotCounts[i]#*100
    year[i]<-unique(data$Year[rows.i])
    IF[i]<-as.numeric(as.character(unique(data[rows.i,IFcolumn])))
    journal[i]<-unlist(strsplit(acronyms[i], as.character(year[i])))
  }
  data.to.model<-data.frame(cbind(per05, per05no, TotCounts, IF, journal, year))
  data.to.model[,1:4]<-apply(X = data.to.model[,1:4], 
                            MARGIN = 2, 
                            FUN = function(X) as.numeric(as.character(X)))
  return(data.to.model)
}

data.to.model2<-data.to.model.func(mydata)
data.to.model2<-data.to.model2[order(data.to.model2$IF),] #ordering the dataset according to increasing IF 


#Building of statistical models statistical models

create.model<-function(fixed, 
                       random, 
                       response=per05,
                       data = data.to.model2, 
                       family=c("binomial"), 
                       weights = data.to.model2$TotCounts) {
    model<-glmer(as.formula(paste(deparse(substitute(response)),
                                  "~", 
                                  fixed, 
                                  "+", 
                                  "(",
                                  random, 
                                  ")",
                                  sep = " ")), 
                 data=data, 
                 family=family, 
                 weights=weights)
}



glmer.journ0<-create.model(fixed=c(1), random=c("1 | journal"))
glmer.IFjourn1<-create.model(fixed = "IF", random = c("1 | journal"))
glmer.IFjourn2<-create.model(fixed = "IF", random = c("IF-1 | journal"))

#Comparing the models
all.models<-list(glmer.journ0,
                  glmer.IFjourn1, 
                  glmer.IFjourn2)

model.comparison<-data.frame( models=as.character(unlist(lapply(X = all.models, 
                                               FUN = function(X) X@call$formula))), 
                              AIC = unlist(lapply(X = all.models, 
                                                FUN = function(X) AIC(X))),
                              BIC = unlist(lapply(X = all.models, 
                                                FUN = function(X) BIC(X)))
                             )

#model.comparison<-model.comparison[order(model.comparison$AIC),] #order the models, lowest AIC first.

model.comparison

#all.models<-all.models[order(model.comparison$AIC)] #order the models, lowest AIC first.

#Coefficients and confidence intervals for all models
all.coef<-lapply(X=all.models, FUN = function(X) summary(X)$coefficients)



confidence<-lapply(X=all.models, FUN = function(X) confint(X))


################
#Selecting the best model 
best.model<-unlist(all.models[which(model.comparison$AIC==
                                               min(model.comparison$AIC))])

#predicting p-value proportions based only on the fixed effects of the best model
IFandPredicted<-cbind(data.to.model2$IF, 
                      invlogit( cbind(1, data.to.model2$IF) %*% 
                                  summary(best.model[[1]])$coefficients[,1]))

#Predicting p-value proportions based only on the fixed effects of the best model 
#lower limit of the 95% confidence interval (2.5%)
ci.best.model<-confint(best.model[[1]])
IFandPredicted025<-cbind(data.to.model2$IF, 
                         invlogit( cbind(1, data.to.model2$IF) %*% ci.best.model[2:3,1]))


#Predicting p-value proportions based only on the fixed effects of the best model
#the higher limit of the 95% confidence interval (97.5%)
IFandPredicted975<-cbind(data.to.model2$IF, 
                         invlogit( cbind(1, data.to.model2$IF) %*% ci.best.model[2:3,2]))





#calculating the average slopes for the fixed effect and the slopes at the edges of the confidence interval
#I am predicting the slopes by calculating the ratio between the predicted y values and the x values (IF)

SlopesAndCI<-function(Xmodel, data){
  ci<-confint(Xmodel[[1]])       #calculating the confidence intervals for the model coefficients
  coefAndCI<-rbind(summary(Xmodel[[1]])$coefficients[,1],
            t(ci)[,-1])
  slopes<-c(predicted = 0, ci025 = 0, ci975 = 0)
  for (i in 1:length(slopes)){
    IFandPred<-cbind(data$IF, 
                     invlogit( cbind(1, data$IF) %*% coefAndCI[i,]))  #calculating predicted values based on the coef of the best model
    temp.slopes<-c()
    for (j in 2:nrow(IFandPred)){
      temp.slopes[j-1]<-(IFandPred[j,2]-IFandPred[j-1,2])/ #calculating the slopes as a ration between the difference between predicted y values and
        (IFandPred[j,1]-IFandPred[j-1,1])   #the difference beween the respective x values (IF)
      finalSlope<-sum(temp.slopes)/j     #averaging the slope values obtained between different points. 
    }
    slopes[i]<-finalSlope
  }
  return(slopes)
}

reportedSlopes<-SlopesAndCI(best.model, data.to.model2)

#################################################################################################
###INEQUALITY P-VALUES###########################################################################
#################################################################################################


#Selection of inequality values.

inexact<-unlist(lapply(all.files, is.character))

myfiles.inexact<-all.files[inexact]

#organize data for modeling
organise.data.inexact<-function(files.list=myfiles.inexact, if.data=IFvalues){
  data<-data.frame()
  for(i in 1:nrow(if.data)){
    no<-grep(if.data[i,3], as.character(names(files.list)))
    data<-rbind(data,
                cbind(Pvalues=files.list[[no]], if.data[i,], row.names = NULL))
  }
  return(data)
}

inexact.data<-organise.data.inexact()

#add all relevant information

data.to.model.func.inexact<-function(data, Pcolumn=1, Journalcolumn=2, 
                                   AcronymColumn=4, IFcolumn=5) {
  data<-data[order(as.numeric(data$I.F)),]
  acronyms<-as.character(unique(data[,AcronymColumn]))
  per05<-c()          #percentage of inexact p<0.05
  TotCounts<-c()      #total number of inexct p-values reported
  IF<-c()             #impact factor
  journal<-c()        #journal acronym
  year<-c()           #year of publicatio
  per05no<-c()        #number of inexact p<0.05
  for (i in 1:length(acronyms)){
    rows.i<-data[,AcronymColumn]==acronyms[i]
    pvalues.i<-data[rows.i, Pcolumn]
    TotCounts[i]<-length(pvalues.i)
    per05no[i]<-length(pvalues.i[pvalues.i=="<0.005" | 
                                   pvalues.i=="<0.01" | 
                                   pvalues.i=="<0.05"])
    per05[i]<-per05no[i]/TotCounts[i]
    IF[i]<-as.numeric(as.character(unique(data[rows.i,IFcolumn])))
    year[i]<-unique(data$Year[rows.i])
    journal[i]<-unlist(strsplit(acronyms[i], as.character(year[i])))
  }
  data.to.model<-data.frame(cbind(per05, per05no, TotCounts, IF, journal, year))
  data.to.model[,1:4]<-apply(X = data.to.model[,1:4], 
                             MARGIN = 2, 
                             FUN = function(X) as.numeric(as.character(X)))
  return(data.to.model)
} 
data.to.model.inexact<-data.to.model.func.inexact(inexact.data)
data.to.model.inexact<-data.to.model.inexact[order(data.to.model.inexact$IF),]



#calculating: 1) the total number of p-values reported in each journal (precise + inexact)
data.to.model.inexact$AllPvalues<-data.to.model.inexact$TotCounts+data.to.model2$TotCounts
#2)the proportion on inexact p-values out of the total reported p-values
data.to.model.inexact$RatioinexactPrecise<-data.to.model.inexact$TotCounts/data.to.model.inexact$AllPvalues
#3) proportion of p-values below 0.05 overall for both precise and inexact p-values. 
data.to.model.inexact$per05all<-(data.to.model.inexact$per05no+data.to.model2$per05no)/data.to.model.inexact$AllPvalues


#statistical models for inexct p-values

create.model.inexact<-function(fixed, 
                             random, 
                             response=per05all,
                             data = data.to.model.inexact, 
                             family=c("binomial"), 
                             weights = data.to.model.inexact$AllPvalues) {
  model<-glmer(as.formula(paste(deparse(substitute(response)),
                                "~", 
                                fixed, 
                                "+", 
                                "(",
                                random, 
                                ")",
                                sep = " ")), 
               data=data, 
               family=family, 
               weights=weights)
}


glmer.journ0.inexact<-create.model.inexact(fixed=c(1), random=c("1 | journal"))
glmer.IFjourn1.inexact<-create.model.inexact(fixed = "IF", random = c("1 | journal"))
glmer.IFjourn2.inexact<-create.model.inexact(fixed = "IF", random = c("IF-1 | journal"))


#comparing the models for inexact p-values
all.models.inexact<-list(glmer.journ0.inexact, 
                       glmer.IFjourn1.inexact, 
                       glmer.IFjourn2.inexact)


model.comparison.inexact<-
  data.frame(models=as.character(unlist(lapply(X = all.models.inexact, 
                                               FUN = function(X) X@call$formula))), 
             AIC = unlist(lapply(X = all.models.inexact, 
                                 FUN = function(X) AIC(X))),
             BIC = unlist(lapply(X = all.models.inexact, 
                                 FUN = function(X) BIC(X)))
  )
all.models.inexact<-all.models.inexact[order(model.comparison.inexact$AIC)]
model.comparison.inexact<-model.comparison.inexact[order(model.comparison.inexact$AIC),]
model.comparison.inexact<-model.comparison.inexact[,-4]


#Coefficients and confidence intervals for all inexactp-values models
all.coef.inexact<-lapply(X=all.models.inexact, FUN = function(X) summary(X)$coefficients)

names(all.coef.inexact)<-unlist(lapply(X = all.models.inexact, FUN = function(X) X@call$formula))

confidence.inexact<-lapply(X=all.models.inexact, FUN = function(X) confint(X))
names(confidence.inexact)<-unlist(lapply(X = all.models.inexact, FUN = function(X) X@call$formula))


################
#CONFIDENCE INTERVALS FOR THE BEST MODEL 
best.model.inexact<-all.models.inexact[which(model.comparison.inexact$AIC==
                                           min(model.comparison.inexact$AIC))]


#Predicting p-value proportions based only on the fixed effects of the best model 
#lower limit of the 95% confidence interval (2.5%)
ci.inexact<-confint(best.model.inexact[[1]])
IFandPredicted025inexact<-cbind(data.to.model.inexact$IF,
                              invlogit( cbind(1, data.to.model.inexact$IF) %*% ci.inexact[2:3,1]))


#Predicting p-value proportions based only on the fixed effects of the best model
#the higher limit of the 95% confidence interval (97.5%)
IFandPredicted975inexact<-cbind(data.to.model.inexact$IF,
                              invlogit( cbind(1, data.to.model.inexact$IF) %*% ci.inexact[2:3,2]))


#calculating the average slopes for the fixed effect and the slopes at the edges of the confidence interval
#I am predicting the slopes by calculating the ratio between the predicted y values and the x values (IF)

SlopesAndCI<-function(Xmodel, data){
  ci<-confint(Xmodel[[1]])
  coefAndCI<-rbind(summary(Xmodel[[1]])$coefficients[,1],
                   t(ci)[,-1])
  slopes<-c(predicted = 0, ci025 = 0, ci975 = 0)
  for (i in 1:length(slopes)){
    IFandPred<-cbind(data$IF, 
                     invlogit( cbind(1, data$IF) %*% coefAndCI[i,]))
    temp.slopes<-c()
    for (j in 2:nrow(IFandPred)){
      temp.slopes[j-1]<-(IFandPred[j,2]-IFandPred[j-1,2])/      #calculating the slopes as a ration between the difference between predicted y values and
        (IFandPred[j,1]-IFandPred[j-1,1])                       #the difference beween the respective x values (IF)
      finalSlope<-sum(temp.slopes)/j                            ##averaging the slope values obtained between different points.
    }
    slopes[i]<-finalSlope
  }
  return(slopes)
}

SlopesAndCI(best.model.inexact, data.to.model.inexact)

#################################################################################################
#####FIGURES#####################################################################################
#################################################################################################
#MAIN PLOT - Figure 1a. Proportion of exact p-values reported across impact factors (IF). 

plot05<-function(data.to.model, all){
  data.to.model<-data.to.model[order(data.to.model$IF),]
  plot(data.to.model$IF, data.to.model$per05, pch=19, 
       xlab="Impact factor", ylab="Proportion of pvalues<0.05", 
       xaxt="n", bty="n", ylim = c(0, 1), type="n")
  par(cex=1)
  if(all) 
  {axis(1, cex.axis=0.8, lwd.ticks = 0, xlim = c(0,1), labels = seq(0.5, 18, 2.5), at=seq(0.5, 18, 2.5))
    polygon(x=c(IFandPredicted025[1,1], IFandPredicted975[1, 1],               #the shading showing the confidence interval of the slope
                IFandPredicted975[17, 1], IFandPredicted025[17,1]),
            y=c(IFandPredicted025[1,2], IFandPredicted975[1, 2], 
                IFandPredicted975[17,2], IFandPredicted025[17,2]),
            col="lightgrey", border = "lightgrey")
    points(data.to.model$IF, data.to.model$per05, pch=19)
  }
  else {points(data.to.model$IF, data.to.model$per05, pch=19)
    axis(1, labels=data.to.model$IF, at=data.to.model$IF, cex.axis=0.73, lwd.ticks = 0)
  }
}

dev.off()
par(mfrow=c(2,1))
plot05(data.to.model2, TRUE)
mtext(at=c(0.35), text = "A", cex=1.5)

#MAIN PLOT - Figure 1b. Proportion of all p-values (exact and inexact) reported across impact factors (IF).

plot05inexact<-function(data.to.model, all, resp.col){
  data.to.model<-data.to.model[order(data.to.model$IF),]
  plot(data.to.model$IF, data.to.model[,resp.col], pch=19, 
       xlab="Impact factor", ylab="Proportion of pvalues<0.05", 
       xaxt="n", bty="n", ylim = c(0, 1), type="n")
  par(cex=1)
  if(all) 
  {axis(1, cex.axis=0.8, lwd.ticks = 0, labels = seq(0.5, 18, 2.5), at=seq(0.5, 18, 2.5))
    polygon(x=c(IFandPredicted025inexact[1,1], IFandPredicted975inexact[1, 1],        #the shading showing the confidence interval of the slope
                IFandPredicted975inexact[17, 1], IFandPredicted025inexact[17,1]),
            y=c(IFandPredicted025inexact[1,2], IFandPredicted975inexact[1, 2],
                IFandPredicted975inexact[17,2], IFandPredicted025inexact[17,2]),
            col="lightgrey", border = "lightgrey")
    points(data.to.model$IF, data.to.model[,resp.col], pch=19)
  }
  else {points(data.to.model$IF, data.to.model[,resp.col], pch=19)
    axis(1, labels = seq(0.5, 18, 2.5), at=seq(0.5, 18, 2.5), cex.axis=0.73, lwd.ticks = 0)
  }
}


plot05inexact(data.to.model.inexact, TRUE, 9)
mtext(at=c(0.35), text = "B", cex=1.5)


#Supplimentary figure S1 a and b. Proportion of exact p-values below 0.05 reported in 
#a) 2012 and b) 2014 across the range of impact factors (IF).
par(mfrow=c(2,1))
plot05(data.to.model2[data.to.model2$year==2012,], FALSE)
mtext(at=c(0.35), text = "A", cex=1.5)
plot05(data.to.model2[data.to.model2$year==2014,], FALSE)
mtext(at=c(0.95), text = "B", cex=1.5)


#Figure S2 Proportion of inexact p-values out of the total reported values across the range of impact factors (IF).
dev.off()
#par(mfrow=c(1,1))
plot(data.to.model.inexact$IF, data.to.model.inexact$RatioinexactPrecise, 
     pch=19, xlab="Impact factor", 
     ylab="Proportion of pvalues reported as inexact values", xaxt="n", bty="n", ylim=c(0, 1))
axis(1, cex.axis=0.8, lwd.ticks = 0, 
     xlim = c(0,max(data.to.model.inexact$IF)), 
     labels = seq(0.5, 18, 2.5), at=seq(0.5, 18, 2.5))


