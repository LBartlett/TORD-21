# Script by Lewis J Bartlett - lewis.bartlett@uga.edu
# Maintained at https://github.com/LBartlett/TORD-21.git

# Load packages and data

library(afex)
library(frailtypack)

TO.Data <- read.csv(file = 'TOData.csv',
                    header = TRUE, stringsAsFactors = FALSE)

#####################

# Begin data manipulation for plotting and analysis

TO.PC.SSM <- data.frame('Colony' = sort(rep(unique(TO.Data$Colony), times = NROW(unique(TO.Data$Wave[which(!is.na(TO.Data$SSMites))])))), 
                        'Wave' = rep(unique(TO.Data$Wave[which(!is.na(TO.Data$SSMites))]), times = NROW(unique(TO.Data$Colony)))
)

TO.PC.SSM$Treatment <- NA
TO.PC.SSM$SSM <- NA
TO.PC.SSM$Yard <- NA
TO.PC.SSM$Inoc <- NA

for(A in 1:NROW(TO.PC.SSM)){
  
  TO.PC.SSM$Treatment[A] <- unique(TO.Data$Treatment[which(TO.Data$Colony == TO.PC.SSM$Colony[A])])
  
  TO.PC.SSM$Inoc[A] <- (unique(TO.Data$Inoculated[which(TO.Data$Colony == TO.PC.SSM$Colony[A])]))=='Yes'
  
  TO.PC.SSM$Yard[A] <- unique(TO.Data$Yard[which(TO.Data$Colony == TO.PC.SSM$Colony[A])])
  
  #TO.PC.SSM$SSM[A] <- sum(na.omit(TO.Data$SSMites[which(TO.Data$Colony == TO.PC.SSM$Colony[A] & TO.Data$Wave == TO.PC.SSM$Wave[A])]))
  
  TO.PC.SSM$SSM[A] <- TO.Data$SSMites[which(TO.Data$Colony == TO.PC.SSM$Colony[A] & TO.Data$Wave == TO.PC.SSM$Wave[A])]
  
  
}

TO.PC.SSM <- na.exclude(TO.PC.SSM)

TO.PC.SSM <- TO.PC.SSM[which(TO.PC.SSM$Wave <7),]

# Quick colour transparency function
Transpa <- function(color, percent) {
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  
  return(t.col)
}

#plot StickyScreen counts per colony through Yards

Shading <- c('lightblue2','orange1','purple1','red3','darkblue','gold2','green3','pink3','black','grey3')

par(mfrow = c(5,4),
    mar = c(4,6,2,2))

for(Y in unique(TO.PC.SSM$Yard)){
  
  for(X in 0:1){
    
    P1 <- 1:max(TO.PC.SSM$Wave)
    P2 <- log10(((P1-1)*((max(TO.PC.SSM$SSM))/(max(P1)-1)))+1)
    
    if(as.logical(X)){
      
      plot(P2 ~ P1, 
           main = paste0(Y, '  Inoculated'), ylab = expression(paste('(log'[10],'+1) Mite Counts'),sep='') ,xlab = 'Month',
           cex.axis = 1.2, cex.lab = 1.2, lwd = 1.8, pch = NA)
      
    }else{
      
      plot(P2 ~ P1, 
           main = paste0(Y, '  Exposed'), ylab = expression(paste('(log'[10],'+1) Mite Counts'),sep='') ,xlab = 'Month',
           cex.axis = 1.2, cex.lab = 1.2, lwd = 1.8, pch = NA)
    }
    
    CTrack <- 0
    
    for(C in unique(as.character(TO.PC.SSM$Colony[which(TO.PC.SSM$Yard == Y & TO.PC.SSM$Inoc == X)]))){
      
      CTrack <- CTrack + 1
      
      points(log10(TO.PC.SSM$SSM[which(TO.PC.SSM$Colony == C)]+1) ~ TO.PC.SSM$Wave[which(TO.PC.SSM$Colony == C)], 
             col = Transpa(Shading[CTrack],35),
             type = "b",
             lwd=2,
             cex = 2,
             pch = 20)
      
      
    }
  }
}

#####################

# Filter down data


PC.SSM <- TO.PC.SSM

PC.SSM <- PC.SSM[which(PC.SSM$Wave < 7),]


PS.TM <- data.frame('Yard' = sort(rep(unique(PC.SSM$Yard), times = NROW(unique(PC.SSM$Wave[which(!is.na(PC.SSM$SSM))])))), 
                    'Wave' = rep(unique(PC.SSM$Wave[which(!is.na(PC.SSM$SSM))]), times = NROW(unique(PC.SSM$Yard)))
)

PS.TM$Mites <- NA

for(A in 1:NROW(PS.TM)){
  
  PS.TM$Mites[A] <- sum(na.omit(PC.SSM$SSM[which(PC.SSM$Yard == PS.TM$Yard[A] & PC.SSM$Wave == PS.TM$Wave[A])]))
  
}

PS.TM$Mites[which(PS.TM$Mites == 0)] <- NA

PS.TM <- na.exclude(PS.TM)

TO.SV <- TO.Data[,c('Colony','Wave','Study','Yard','Inoculated','Alive')]
PC.SV <- TO.SV
PC.SV <- PC.SV[which(PC.SV$Wave < 7),]

PS.SV <- data.frame(Yard = sort(rep(unique(PC.SV$Yard), times = NROW(unique(PC.SV$Wave)))),
                    Wave = rep(unique(PC.SV$Wave), times = NROW(unique(PC.SV$Yard))))

PS.SV$TA <- NA

for(S in 1:NROW(PS.SV)){
  
  PS.SV$TA[S]  <- sum(PC.SV$Alive[which(PC.SV$Yard == PS.SV$Yard[S] & PC.SV$Wave == PS.SV$Wave[S] & PC.SV$Inoculated == 'No')])
  
}

GrowthCor <- data.frame(Yard = unique(PC.SSM$Yard))
GrowthCor$SG <- NA
GrowthCor$EG <- NA
GrowthCor$Total <- NA

# Models to estimate growth rates

for(C in 1:NROW(GrowthCor)){
  
  Yard <- GrowthCor$Yard[C]
  
  GrowthCor$SG[C] <- coef(glmer(SSM ~ Wave - 1 + (1|Colony ),
                                family = 'poisson',
                                data = PC.SSM[which(PC.SSM$Inoc == T & PC.SSM$Yard == Yard),]))[[1]][1,2]
  
  
  GrowthCor$EG[C] <- coef(glmer(SSM ~ Wave - 1 + (1|Colony ),
                                family = 'poisson',
                                data = PC.SSM[which(PC.SSM$Inoc == F & PC.SSM$Yard == Yard),]))[[1]][1,2]
  
  GrowthCor$Total[C] <- coef(glm(Mites ~ Wave - 1,
                                 family = 'poisson',
                                 data = PS.TM[which(PS.TM$Yard == Yard),]))[[1]]
  
}

# Plots and analyses

par(mfrow = c(1,1))

plot(GrowthCor$EG ~ GrowthCor$SG)
abline(lm(GrowthCor$EG ~ GrowthCor$SG))
anova(lm(GrowthCor$EG ~ GrowthCor$SG))
summary(lm(GrowthCor$EG ~ GrowthCor$SG))

plot(GrowthCor$Total ~ GrowthCor$SG)
abline(lm(GrowthCor$Total ~ GrowthCor$SG))
anova(lm(GrowthCor$Total ~ GrowthCor$SG))
summary(lm(GrowthCor$Total ~ GrowthCor$SG))


plot(GrowthCor$EG ~ GrowthCor$SG,
     pch = 20, cex = 1.8,
     ylab = 'Mite Growth Coefficient (Exposed)',
     xlab = 'Mite Growth Coefficient (Inoculated)',
     cex.lab = 1.25,
     cex.axis = 1.25)

lines(x= c(0.55, 1.4), 
      y = predict(lm(EG ~ SG, data = GrowthCor), newdata=data.frame(SG=c(0.55, 1.4))),
      lwd = 3,
      lty = 2)

plot(GrowthCor$Total ~ GrowthCor$SG,
     pch = 20, cex = 1.8,
     ylab = 'Mite Growth Coefficient (Apiary)',
     xlab = 'Mite Growth Coefficient (Inoculated)',
     cex.lab = 1.25,
     cex.axis = 1.25)

lines(x= c(0.55, 1.4), 
      y = predict(lm(Total ~ SG, data = GrowthCor), newdata=data.frame(SG=c(0.55, 1.4))),
      lwd = 3,
      lty = 2)



# Frailty

FrailDat <- data.frame('Colony' = unique(TO.Data$Colony[which(TO.Data$Inoculated == 'No')]))
FrailDat$Yard <- NA
FrailDat$Treatment <- NA
FrailDat$Died <- NA
FrailDat$TTE <- NA

for(N in 1:NROW(FrailDat)){
  
  FrailDat$Yard[N] <- as.character(unique(TO.Data$Yard[which(TO.Data$Colony == FrailDat$Colony[N])]))
  FrailDat$Treatment[N] <- GrowthCor$SG[which(GrowthCor$Yard == FrailDat$Yard[N])]
  FrailDat$Died[N] <- !(sum(TO.Data$Alive[which(TO.Data$Colony == FrailDat$Colony[N])]) == 24)
  FrailDat$TTE[N] <- sum(TO.Data$Alive[which(TO.Data$Colony == FrailDat$Colony[N])])
  
}

FrailDat$Died[which(FrailDat$TTE > 12)] <- FALSE
FrailDat$TTE[which(FrailDat$TTE > 12)] <- 12

library(frailtypack)

frailtyPenal(formula = Surv(TTE,Died, type = 'right') ~ Treatment + cluster(Yard),
             n.knots = 6, kappa = 10000,
             data=FrailDat)



# Survival plot (requires some rearranging)

TO.Surv.Frame <- data.frame('Yard' = rep(unique(TO.Data$Yard), times =23), 
                            'Wave' = rep(2:24, times = NROW(unique(TO.Data$Yard)))
)

TO.Surv.Frame$Treatment <- NA

TO.Surv.Frame$SurvWave <- NA
TO.Surv.Frame$DeadWave <- NA

for(M in 1:NROW(TO.Surv.Frame)){
  
  TO.Surv.Frame$SurvWave[M] <- sum(TO.Data$Alive[which(TO.Data$Wave == TO.Surv.Frame$Wave[M] & TO.Data$Yard == TO.Surv.Frame$Yard[M] & TO.Data$Inoculated == 'No')])
  
  TO.Surv.Frame$DeadWave[M] <- (sum(TO.Data$Alive[which(TO.Data$Wave == (TO.Surv.Frame$Wave[M] - 1) & TO.Data$Yard == TO.Surv.Frame$Yard[M] & TO.Data$Inoculated == 'No')])
                                -
                                  sum(TO.Data$Alive[which(TO.Data$Wave == TO.Surv.Frame$Wave[M] & TO.Data$Yard == TO.Surv.Frame$Yard[M] & TO.Data$Inoculated == 'No')])
  )
  
  TO.Surv.Frame$Treatment[M] <- unique(as.character(TO.Data$Treatment[which(TO.Data$Yard ==  TO.Surv.Frame$Yard[M])]))
}

TO.Surv.Frame$WaveF <- as.factor(TO.Surv.Frame$Wave)

TO.Surv.Frame$RV <- cbind(TO.Surv.Frame$SurvWave,TO.Surv.Frame$DeadWave)


# SURVIVAL STYLE PLOT

par(mfrow = c(1,1))

# Quick colour transparency function
Transpa <- function(color, percent) {
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  
  return(t.col)
  
}


TO.Surv.Frame$PropSurv <- TO.Surv.Frame$SurvWave/12

P1 <- (0:10)/10
P2 <- 2:12

plot(P1 ~ P2, 
     main = NA, ylab = 'Proportion Surviving' ,xlab = 'Month',
     cex.axis = 1.3, cex.lab = 1.5, lwd = 1.8, pch = NA)


for(Y in unique(TO.Surv.Frame$Yard)){
  
  
  YTCol <- rgb(red = (GrowthCor$SG[which(GrowthCor$Yard == Y)] - min(GrowthCor$SG))/(max(GrowthCor$SG) - min(GrowthCor$SG)),
               green = 0,
               blue = (1 - (GrowthCor$SG[which(GrowthCor$Yard == Y)] - min(GrowthCor$SG))/(max(GrowthCor$SG) - min(GrowthCor$SG))),
               alpha = 0.5)
  
  lines(y = sort(TO.Surv.Frame$PropSurv[which(TO.Surv.Frame$Yard == Y)], decreasing = T), 
        x = sort(TO.Surv.Frame$Wave[which(TO.Surv.Frame$Yard == Y)]), 
        col = YTCol, 
        type = 's', pch = 19, lwd = 6, cex = 1.4)
  
}


#### Additional alcohol wash bit

#
TO.AW <- TO.Data[which(!is.na(TO.Data$AWMites)),c('Colony','Yard','Inoculated','AWMites')]

TO.AW$PSS <- NA
TO.AW$PAB <- NA
TO.AW$PBC <- NA

for(N in 1:NROW(TO.AW)){
  
  TO.AW$PSS[N] <- TO.Data$SSMites[which(TO.Data$Wave == 6 & TO.Data$Colony == TO.AW$Colony[N])]
  
  TO.AW$PAB[N] <- TO.Data$ABees[which(TO.Data$Wave == 6 & TO.Data$Colony == TO.AW$Colony[N])]
  
  TO.AW$PBC[N] <- TO.Data$CBrood[which(TO.Data$Wave == 6 & TO.Data$Colony == TO.AW$Colony[N])]
  
}

TO.AW <- na.exclude(TO.AW)

TO.AW <- TO.AW[which(TO.AW$Inoculated == 'No'),]

par(mfrow = c(1,1), mar = c(6,6,2,2))

plot(log(TO.AW$AWMites+1) ~ log(TO.AW$PSS+1))

cor.test(x = log(TO.AW$AWMites+1), y = log(TO.AW$PSS+1))

cor.test(x = log(TO.AW$AWMites+1), y = log(TO.AW$PAB+1))
cor.test(x = log(TO.AW$AWMites+1), y = log(TO.AW$PBC+1))

boxplot(TO.AW$AWMites ~ TO.AW$Yard)

anova(glm(TO.AW$AWMites ~ TO.AW$Yard, family = 'poisson'), test = 'Chisq')

plot(log(TO.AW$AWMites+1) ~ log(TO.AW$PSS+1),
     ylab = 'Per-Capita Parasitism (log+1 transformed)',
     xlab = 'Mite Population Size (log+1 transformed)',
     pch = 20, cex = 1.5, cex.lab = 1.8, cex.axis = 1.5
)

