### R code from vignette source 'SimpleGuide.Rnw'

###################################################
### code chunk number 1: SimpleGuide.Rnw:85-86
###################################################
options(width=80,useFancyQuotes="UTF-8")


###################################################
### code chunk number 2: SimpleGuide.Rnw:92-97
###################################################
library(survival)
library(rstpm2)
brcancer <- transform(brcancer, recyear=rectime / 365.24)
fit.cox <- coxph(Surv(recyear,censrec==1)~hormon, data=brcancer)
summary(fit.cox)


###################################################
### code chunk number 3: SimpleGuide.Rnw:105-108
###################################################
fit <- stpm2(Surv(recyear,censrec==1)~hormon, data=brcancer, df=4)
summary(fit)
eform(fit)[2,]


###################################################
### code chunk number 4: SimpleGuide.Rnw:115-120
###################################################
plot(fit, newdata=data.frame(hormon=0), xlab="Time since diagnosis (years)")
lines(fit, newdata=data.frame(hormon=1), lty=2)
lines(survfit(Surv(recyear,censrec==1)~hormon, data=brcancer), col="blue", lty=1:2)
legend("topright", c("PH hormon=0","PH hormon=1","KM hormon=0","KM hormon=1"), 
       lty=1:2, col=c("black","black","blue","blue"))


###################################################
### code chunk number 5: SimpleGuide.Rnw:125-136
###################################################
library(ggplot2)
predHormon <- predict(fit, newdata=data.frame(hormon=0:1), 
                      type="hazard", grid=TRUE, full=TRUE, se.fit=TRUE)
predHormon <- transform(predHormon,Hormone=factor(hormon,labels=c("No","Yes")))
ggplot(predHormon, 
       aes(x=recyear,y=Estimate,ymin=lower,ymax=upper,fill=Hormone)) +
    facet_grid(~Hormone) + 
    xlab("Time since diagnosis (years)") + 
    ylab("Hazard") + 
    geom_ribbon() + 
    geom_line()


###################################################
### code chunk number 6: SimpleGuide.Rnw:139-145
###################################################
ggplot(predHormon, 
       aes(x=recyear,y=Estimate,ymin=lower,ymax=upper,fill=Hormone)) +
    xlab("Time since diagnosis (years)") + 
    ylab("Hazard") + 
    geom_ribbon(alpha=0.6) +
    geom_line()


###################################################
### code chunk number 7: SimpleGuide.Rnw:153-161
###################################################
par(mfrow=1:2)
plot(fit,newdata=data.frame(hormon=0), type="hdiff",
     exposed=function(data) transform(data, hormon=1),
     xlab="Time since diagnosis (years)")
plot(fit,newdata=data.frame(hormon=0), type="sdiff",
     var="hormon",
     xlab="Time since diagnosis (years)")
mtext("Effect of hormonal treatment", outer = TRUE, line=-3, cex=1.5, font=2)


