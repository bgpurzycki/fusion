###############################################################################
##### Identity Fusion, Outgroups, and Sacrifice:  A Cross-Cultural Test #######
###############################################################################

### Script written by Benjamin Grant Purzycki and Martin Lang
### contact email:  benjamin_purzycki@eva.mpg.de and martinlang@fas.harvard.edu

### Last updated December 21, 2018

setwd("")
set.seed(7)

## Load packages
library(psych)
library(brms)
library(lme4)

## Read data
d <- read.csv("d.csv")

## Transform data
d$t <- NA #add vector of total possible coins
d$t <- 30

o <- as.integer(d$EMOO) # outgroup relations
c <- as.integer(d$EMOC) # distant ingroup fusion
g <- as.integer(d$EMOI) # local ingroup fusion
y <- as.integer(d$COREL.S) # coins in distant ingroup cup
yx <- as.integer(d$SELF) # coins in players' cups
s <- as.integer(d$CORELSIM) # religious similarity of distant ingroup
site <- as.numeric(factor(d$SITE)) # field site
t <- as.integer(d$t) # total # of coins
chk <- as.integer(d$INGFIRST) # played other game first?
pr <- as.integer(d$TREATMENT) # treatment condition?

# Create subsets
newdata1 <- cbind(g, o, y, c, site, t)
newdata2 <- cbind(g, o, y, s, c, site, t)
newdata3 <- cbind(g, o, y, c, site, t, chk)
newdata4 <- cbind(g, o, y, c, site, t, pr)
n1 <- data.frame(na.omit(newdata1))
n2 <- data.frame(na.omit(newdata2))
n3 <- data.frame(na.omit(newdata3))
n4 <- data.frame(na.omit(newdata4))

## Statistics for Table 1
aggregate(d$EMOI, list(d$SITE), mean, na.rm = T)
aggregate(d$EMOI, list(d$SITE), sd, na.rm = T)

aggregate(d$EMOO, list(d$SITE), mean, na.rm = T)
aggregate(d$EMOO, list(d$SITE), sd, na.rm = T)

aggregate(d$SELF, list(d$SITE), mean, na.rm = T)
aggregate(d$SELF, list(d$SITE), sd, na.rm = T)

#####################
### Main Analyses ###
#####################

### Set priors
priorm1 <- c(set_prior("normal(0,1)", class = "b"),
             set_prior("cauchy(0,2)", class = "sd"),
             set_prior("lkj(4)", class = "cor"))

# Panels a and b in main plot
m1a <- brm(formula = y | trials(t) ~ mo(g)*mo(o) + (mo(g)*mo(o) | site), 
           data = n1, family = binomial("logit"),
           prior = priorm1,
           warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

m1b <- brm(formula = y | trials(t) ~ mo(s) + mo(g)*mo(o) + (mo(g)*mo(o) | site), 
           data = n2, family = binomial("logit"),
           prior = priorm1,
           warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

# Panels c and d in main plot
mg <- brm(formula = y | trials(t) ~ mo(g) + (mo(g) | site), 
          data = n1, family = binomial("logit"),
          prior = priorm1,
          warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

mo <- brm(formula = y | trials(t) ~ mo(o) + (mo(o) | site), 
          data = n1, family = binomial("logit"),
          prior = priorm1,
          warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

# Basic summaries, diagnostics, counterfactuals, and visualizatons
summary(m1a)
plot(m1a)
marginal_effects(m1a, surface = T)
ranef(m1a)
launch_shinystan(m1a)
stanplot(m1a, pars = "^b", type = "intervals")

logistic(fixef(m1a))
logistic(fixef(m1a)[1,1])
logistic(fixef(m1a)[1,1] + fixef(m1a)[2,1]*5 + fixef(m1a)[3,1]*1 + fixef(m1a)[4,1]*5*1)

logistic(fixef(m1b))
logistic(fixef(m1b)[1,1])
logistic(fixef(m1b)[1,1] + fixef(m1b)[2,1]) # add similarity
logistic(fixef(m1b)[1,1] + fixef(m1b)[2,1]*2) # add max similarity


logistic(fixef(mg))
logistic(fixef(mg)[1,1])
logistic(fixef(mg)[1,1] + fixef(mg)[2,1]*5)

logistic(fixef(mo))
logistic(fixef(mo)[1,1])
logistic(fixef(mo)[1,1] + fixef(mo)[2,1]*5)

# Posterior predictions
pp <- predict(m1a)
head(pp)

yob <- data.frame(g = 5, o = 1, t = 30, site = factor(c("1", "2", "3", "4", "5", "6", "7", "8")))
predict(m1a, newdata = yob)

meshuggah <- data.frame(g = 5, o = 5, t = 30, site = factor(c("1", "2", "3", "4", "5", "6", "7", "8")))
predict(m1a, newdata = meshuggah)

####################################
###### Supplementary Analyses ######
####################################

### Plots
mycol2 <- rgb(0,0,139, max=255, alpha = 100, names = "blue")
mycol3 <-  rgb(255, 140, 0, max = 255, alpha = 100, names = "orange")

# Group-level
denss <- density(s, na.rm = T)
densg <- density(g, na.rm = T)
denso <- density(o, na.rm = T)
plot(NA, xlab = "Ingroup and Outgroup Scores", ylab = "Density", xlim=c(0,6), ylim=c(0,.8), cex.lab = 1.3)
polygon(densg, col = mycol2)
polygon(denso, col = mycol3)
plot(NA, xlab = "Religious Similarity", ylab = "Density", xlim=c(-3,3), ylim=c(0,.5), cex.lab = 1.3)
polygon(denss, col = mycol2)

# Ingroup vs. Outgroup Densities by site
subvars <- c("SITE", "EMOI", "EMOO")
fac <- d[subvars]
validrows <- fac[complete.cases(fac),]
levels(validrows$SITE)[match("Pesqueiro", levels(validrows$SITE))] <- "Marajo"

par(mfrow=c(4,2), mar = c(4,4,3,3)) 
for (site in unique(validrows$SITE)){  
  sub <- subset(validrows, site == validrows$SITE) 
  densL <- density(sub$EMOI)
  densM <- density(sub$EMOO)
  demo <- data.frame(Ingroup=sub$EMOI, Outgroup=sub$EMOO)
  dens <- apply(demo, 2, density, na.rm=T)
  plot(NA, xlab = site, ylab = "Density", xlim=c(-1,6), ylim=c(0,1), cex.lab = 1.3)
  polygon(densM, col = mycol2)
  polygon(densL, col = mycol3)
  #legend("top", legend=names(dens), fill=c(mycol2, mycol3), bty = "n")
}

######################
## Bayesian extensions

ms <- brm(formula = y | trials(t) ~ mo(o)*mo(g) + (mo(o) | site), 
          data = n1, family = binomial("logit"),
          prior = priorm1,
          warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

mp <- brm(formula = y | trials(t) ~ mo(g)*mo(o) + pr + (mo(g)*mo(o) | site), 
          data = n4, family = binomial("logit"),
          prior = priorm1,
          warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

mt <- brm(formula = y | trials(t) ~ mo(g)*mo(o) + chk + (mo(g)*mo(o) | site), 
          data = n3, family = binomial("logit"),
          prior = priorm1,
          warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

mc1 <- brm(formula = y | trials(t) ~ mo(c)*mo(o) + (mo(c)*mo(o) | site), 
           data = n1, family = binomial("logit"),
           prior = priorm1,
           warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

mc2 <- brm(formula = y | trials(t) ~ mo(c) + (mo(c) | site), 
           data = n1, family = binomial("logit"),
           prior = priorm1,
           warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

### Tyvan Individualist Game

set.seed(7)

tyva <- read.csv("CERC_TYVA.csv")

# Transform data
tyva$t <- NA
tyva$t <- 30

o <- as.integer(tyva$OUTGREMO)
g <- as.integer(tyva$EMOI) # includes the data point with 4.5 converted to 5
y <- as.integer(tyva$INGROUP.1)
t <- as.integer(tyva$t)
chk <- as.integer(tyva$INGFIRST)
pr <- as.integer(tyva$TREATMENT)

newtyva <- cbind(o, g, y, t, chk, pr)
ntyva1 <- data.frame(na.omit(newtyva))

priortyva <- c(set_prior("normal(0,1)", class = "b"))

m1tyva <- brm(formula = y | trials(t) ~ mo(o)*mo(g) + mo(o) + mo(g), 
          data = ntyva1, family = binomial("logit"),
          prior = priortyva,
          warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

m2tyva <- brm(formula = y | trials(t) ~ mo(o) + mo(g), 
              data = ntyva1, family = binomial("logit"),
              prior = priortyva,
              warmup = 1000, iter = 5000, chains = 2, control = list(adapt_delta = .99), seed = 7)

####################################################################################################################
### IDENTITY FUSION - FREQUENTIST ANALYSIS #########################################################################
####################################################################################################################

# load libraries
libs <- c("car" ,"dplyr",  "DHARMa", "optimx", "psychometric", "Hmisc")
lapply(libs, require, character.only = TRUE)

# load data set
dat <- read.csv("d.csv", stringsAsFactors=FALSE)

# label variables
dat <- dat %>% 
  mutate(
    INGREMO = EMOI, # Ingroup fusion
    CORELEMO = EMOC, # Extended fusion
    OUTGREMO = EMOO) # Outgroup relations

# create data set w/o NAs
dat1<- dat %>%  
  filter(!is.na(OUTGREMO))


##############################
#### Fusion analysis #########

#First, check whether INGREMO and OUTGREMO negatively correlate
rcorr(dat$INGREMO, dat$OUTGREMO)
CIr(r=.28, n = 571, level = .95) # 95% CI for correlation

# what is the number of PPTs per cell of OUTGRMO and INGREMO interaction? Information for Table S1
table(dat$OUTGREMO, dat$INGREMO)



### Ingroup fusion models ################################

### Basic - check whether ingroup fusion predicts sacrifice ###
summary(fg <- glmer(cbind(COREL.S, 30-COREL.S) ~ INGREMO + (INGREMO|SITE), data = dat, family = 'binomial'))

# plot residuals
par(mfrow=c(3,1), mar=c(2,2,2,2))
plot(exp(predict(fg)),residuals(fg),  main = "Fitted vs Residuals")
qqnorm(residuals(fg), main = "QQ plot")
hist(residuals(fg), main = "Histogram")

# Assess homoscedasticity with DHARMa package
simulationOutput <- simulateResiduals(fittedModel = fg, n = nobs(fg))
plotSimulatedResiduals(simulationOutput = simulationOutput)


### Extended - check whether outgroup relations predicts sacrifice ###
summary(fo <- glmer(cbind(COREL.S, 30-COREL.S) ~ OUTGREMO + (OUTGREMO|SITE), data = dat, family = 'binomial'))

# plot residuals
par(mfrow=c(3,1), mar=c(2,2,2,2))
plot(exp(predict(fo)),residuals(fo),  main = "Fitted vs Residuals")
qqnorm(residuals(fo), main = "QQ plot")
hist(residuals(fo), main = "Histogram")

# Assess homoscedasticity with DHARMa package
simulationOutput <- simulateResiduals(fittedModel = fo, n = nobs(fo))
plotSimulatedResiduals(simulationOutput = simulationOutput)


### Interaction model - INGREMO and OUTGREMO should negatively interact ###
summary(fgo <- glmer(cbind(COREL.S, 30-COREL.S) ~ INGREMO * OUTGREMO + (INGREMO * OUTGREMO|SITE),
                     data = dat,
                     family = 'binomial',
                     control=glmerControl(optimizer="optimx",optCtrl  = list(method="bobyqa"))))
ss <- getME(fgo, c("theta","fixef"))
summary(update(fgo, start=ss,control=glmerControl(optCtrl=list(maxfun=2e4))))
# problems with convergence - we need to factor OUTGREMO

#trichotimize outgremo for plots
dat$OUTGREMO.RD <- dat$OUTGREMO; dat$OUTGREMO.RD[dat$OUTGREMO < 3] <- 0; dat$OUTGREMO.RD[dat$OUTGREMO == 3] <- 1;
dat$OUTGREMO.RD[dat$OUTGREMO == 4] <- 2; dat$OUTGREMO.RD[dat$OUTGREMO == 5] <- 2
dat$OUTGREMO.RD <- as.factor(dat$OUTGREMO.RD)

### Interaction model II ###
summary(fgo <- glmer(cbind(COREL.S, 30-COREL.S) ~ INGREMO * OUTGREMO.RD + (INGREMO * OUTGREMO.RD|SITE),
                     data = dat,
                     family = 'binomial',
                     control=glmerControl(optimizer="optimx",optCtrl  = list(method="bobyqa"))))
ss <- getME(fgo, c("theta","fixef"))
summary(update(fgo, start=ss,control=glmerControl(optCtrl=list(maxfun=2e4))))
# convergence almost OK

# plot residuals
par(mfrow=c(3,1), mar=c(2,2,2,2))
plot(exp(predict(fgo)),residuals(fgo),  main = "Fitted vs Residuals")
qqnorm(residuals(fgo), main = "QQ plot")
hist(residuals(fgo), main = "Histogram")

# Assess homoscedasticity with DHARMa package
simulationOutput <- simulateResiduals(fittedModel = fgo, n = nobs(fgo))
plotSimulatedResiduals(simulationOutput = simulationOutput)


# random-effects coefficients for sites
round(coef(fg)$SITE[1], digits = 3)
round(coef(fg)$SITE[2], digits = 3)
round(coef(fo)$SITE[1], digits = 3)
round(coef(fo)$SITE[2], digits = 3)
round(coef(fgo)$SITE[1], digits = 3)
round(coef(fgo)$SITE[3], digits = 3)
round(coef(fgo)$SITE[4], digits = 3)


### Distant fusion models ################################

### Basic - check whether ingroup fusion predicts sacrifice ###
summary(fc <- glmer(cbind(COREL.S, 30-COREL.S) ~ CORELEMO + (1|SITE), data = dat, family = 'binomial'))
# convergence problems, hence no CORELEMO random slopes

# plot residuals
par(mfrow=c(3,1), mar=c(2,2,2,2))
plot(exp(predict(fc)),residuals(fc),  main = "Fitted vs Residuals")
qqnorm(residuals(fc), main = "QQ plot")
hist(residuals(fc), main = "Histogram")

# Assess homoscedasticity with DHARMa package
simulationOutput <- simulateResiduals(fittedModel = fc, n = nobs(fc))
plotSimulatedResiduals(simulationOutput = simulationOutput)

### Interaction model II ###
summary(fco <- glmer(cbind(COREL.S, 30-COREL.S) ~ CORELEMO * OUTGREMO.RD + (CORELEMO * OUTGREMO.RD|SITE),
                     data = dat,
                     family = 'binomial',
                     control=glmerControl(optimizer="optimx",optCtrl  = list(method="bobyqa"))))
ss <- getME(fco, c("theta","fixef"))
summary(update(fco, start=ss,control=glmerControl(optCtrl=list(maxfun=2e4))))

# plot residuals
par(mfrow=c(3,1), mar=c(2,2,2,2))
plot(exp(predict(fco)),residuals(fco),  main = "Fitted vs Residuals")
qqnorm(residuals(fco), main = "QQ plot")
hist(residuals(fco), main = "Histogram")

# Assess homoscedasticity with DHARMa package
simulationOutput <- simulateResiduals(fittedModel = fco, n = nobs(fco))
plotSimulatedResiduals(simulationOutput = simulationOutput)


#  random-effects coefficients for sites
round(coef(fc)$SITE[1], digits = 3)
round(coef(fco)$SITE[1], digits = 3)
round(coef(fco)$SITE[3], digits = 3)
round(coef(fco)$SITE[4], digits = 3)


####################################################################################################################
### SUPPLEMENTARY PLOTS ############################################################################################
####################################################################################################################

#load libraries
libs <- c("ggplot2", "viridis", "akima", "rgl")
lapply(libs, require, character.only = TRUE)


### Figure S3 ############################
# Surface plot

# prepare dataset
dat.x <- na.omit(subset(dat, select = c(OUTGREMO, OUTGREMO, INGREMO, SITE, COREL.S)))

# Interpolate data fro smoother-looking plots
y=dat.x$INGREMO
x=dat.x$OUTGREMO
y = abs(y-6) # reverse INGREMO to correspond to Table S1
z=dat.x$COREL.S
s=akima::interp(x,y,z, duplicate = "median") # interpolation
dim(s$z) # extract dimensions from interpolated data

# a lame transformation of data for surface plots
{dat.c <- 0
  dat.c <- as.data.frame(s$x)
  colnames(dat.c) <- "x"
  dat.c2 <- matrix(nrow=40*40,ncol=3)
  for(i in 1:40){
    dat.c2[((40*i)-39):((40*i)),1] <- dat.c$x
  }
  dat.c$y <- s$y
  for(i in 1:40){
    dat.c2[((40*i)-39):((40*i)),2] <- dat.c$y[i]
  }
  
  for(i in 1:40){
    dat.c2[((40*i)-39):((40*i)),3] <- s$z[1:40,i]
  }
  dat.c2 <- as.data.frame(dat.c2)
  colnames(dat.c2) <- c("OUTGREMO","INGREMO","DISTANT_Cup")
}


# plot
ggplot(dat.c2, aes(rev(OUTGREMO), INGREMO, z = DISTANT_Cup)) +
  
  geom_contour(aes(colour = factor(round(DISTANT_Cup))), binwidth = 0.05) +
  
  scale_color_viridis(discrete=T) + 
  
  scale_x_continuous(limits=c(1,5), breaks = c(seq(1,5,1))) +
  
  ylab("Outgroup Relations") + xlab("Ingroup Fusion") +
  
  theme_bw() + 
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(2)),        
    axis.line = element_line(colour = "black"),
    legend.position = "", # you can plot legend here
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.title=element_text(size=rel(1.5)), 
    legend.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(2)),
    axis.text.y= element_text(size = rel(2)),
    axis.text.x= element_text(size = rel(2)),
    plot.margin=unit(c(1,1,1,1),"cm"),
    strip.text.x = element_text(size = rel(2))) 

# ggsave("Figure S3.tiff", width = 5.5 , height = 5.5, dpi = 300)

### Figure S4 ############################
# Surface plots

# Extract marginal effects from model m1a
dat.m <- marginal_effects(m1a, surface = T)
dat.m <- dat.m$`g:o`

# interpolate the marginal effects for smoother plot
y=dat.m$o
x=dat.m$g
x = abs(x-6) # reverse INGREMO values to correspond to Table S1 and Figure S3
z=dat.m$estimate__
s=akima::interp(x,y,z, duplicate = "median")

# use again the lame way to create a new data frame for interpolated data
{dat.c <- 0
  dat.c <- as.data.frame(s$x)
  colnames(dat.c) <- "x"
  dat.c2 <- matrix(nrow=40*40,ncol=3)
  for(i in 1:40){
    dat.c2[((40*i)-39):((40*i)),1] <- dat.c$x
  }
  dat.c$y <- s$y
  for(i in 1:40){
    dat.c2[((40*i)-39):((40*i)),2] <- dat.c$y[i]
  }
  
  for(i in 1:40){
    dat.c2[((40*i)-39):((40*i)),3] <- s$z[1:40,i]
  }
  dat.c2 <- as.data.frame(dat.c2)
  colnames(dat.c2) <- c("OUTGREMO","INGREMO","DISTANT_Cup")
}

# this step just increases the density of interpolated data for easier reading, can be set to 1
dat.c2$DISTANT_Cup <- dat.c2$DISTANT_Cup*3.5

ggplot(dat.c2, aes(rev(OUTGREMO), INGREMO, z = DISTANT_Cup)) + 
  geom_contour(aes(colour = factor(round(DISTANT_Cup))), binwidth = 0.05) +
  scale_color_viridis(discrete=T) + 
  theme_bw() + 
  scale_x_continuous(limits=c(1,5), breaks = c(seq(1,5,1))) +
  ylab("Outgroup Relations") + xlab("Ingroup Fusion") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(2)),        
    axis.line = element_line(colour = "black"),
    legend.position = "",
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.title=element_text(size=rel(1.5)), 
    legend.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(2)),
    axis.text.y= element_text(size = rel(2)),
    axis.text.x= element_text(size = rel(2)),
    plot.margin=unit(c(1,1,1,1),"cm"),
    strip.text.x = element_text(size = rel(2))) 

# ggsave("Figure S4a.tiff", width = 5.5 , height = 5.5, dpi = 300)


# Extract marginal effects from model m1b
dat.m <- marginal_effects(m1b, surface = T)
dat.m <- dat.m$`g:o`

# interpolate the marginal effects for smoother plot
y=dat.m$o
x=dat.m$g
x = abs(x-6) # reverse INGREMO values to correspond to Table S1 and Figure S3
z=dat.m$estimate__
s=akima::interp(x,y,z, duplicate = "median")

# create a new data frame for interpolated data
{dat.c <- 0
  dat.c <- as.data.frame(s$x)
  colnames(dat.c) <- "x"
  dat.c2 <- matrix(nrow=40*40,ncol=3)
  for(i in 1:40){
    dat.c2[((40*i)-39):((40*i)),1] <- dat.c$x
  }
  dat.c$y <- s$y
  for(i in 1:40){
    dat.c2[((40*i)-39):((40*i)),2] <- dat.c$y[i]
  }
  
  for(i in 1:40){
    dat.c2[((40*i)-39):((40*i)),3] <- s$z[1:40,i]
  }
  dat.c2 <- as.data.frame(dat.c2)
  colnames(dat.c2) <- c("OUTGREMO","INGREMO","DISTANT_Cup")
}

# this step just increases the density of interpolated data for easier reading
dat.c2$DISTANT_Cup <- dat.c2$DISTANT_Cup*3.5

ggplot(dat.c2, aes(rev(OUTGREMO), INGREMO, z = DISTANT_Cup)) + 
  geom_contour(aes(colour = factor(round(DISTANT_Cup))), binwidth = 0.05) +
  scale_color_viridis(discrete=T) + 
  theme_bw() + 
  scale_x_continuous(limits=c(1,5), breaks = c(seq(1,5,1))) +
  ylab("Outgroup Relations") + xlab("Ingroup Fusion") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(2)),        
    axis.line = element_line(colour = "black"),
    legend.position = "",
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.title=element_text(size=rel(1.5)), 
    legend.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(2)),
    axis.text.y= element_text(size = rel(2)),
    axis.text.x= element_text(size = rel(2)),
    plot.margin=unit(c(1,1,1,1),"cm"),
    strip.text.x = element_text(size = rel(2))) 

# ggsave("Figure S4b.tiff", width = 5.5 , height = 5.5, dpi = 300)

### Figure S5 ############################
# Points plot

# extract estimates with credibility intervals from model m1a
dat.m <- marginal_effects(m1a, surface = T)
dat.m <- dat.m$`g:o`
# order based on OUTGREMO values
dat.m2 <- dat.m[order(dat.m$o),]
# select only necessary values
dat.m2 <- dat.m2[,c(1,2,7,9,10)]
# select three levels of OUTGREMO
dat.m2 <- dat.m2[dat.m2$o==1 | dat.m2$o==3 |dat.m2$o == 5,]
# factor OUTGREMO
dat.m2$o <- factor(dat.m2$o)

dat.m2$variable <- 'labelx'


# plot
ggplot() +
  
  geom_hline(yintercept = 15, lty=2, lwd=1, colour="#736A62", alpha = 0.9) +
  
  geom_point(data = dat.m2,aes(x=g, y=estimate__, color=o, shape = o),
             position=position_dodge(0.25), size = 3) +
  
  geom_errorbar(data = dat.m2, aes(x=g,
                                   ymin=lower__,ymax=upper__, color =
                                     o), width=.2, size =1,
                position=position_dodge(0.25)) +
  
  scale_color_manual(name = "Outgroup relations",breaks = c("1","3","5"),
                     values = c("#999999", "#56B4E9", "#E69F00"),
                     labels = c("1", "3", "5")) +
  
  scale_shape_manual(name = "Outgroup relations",breaks = c("1","3","5"),
                     values = c(16, 17, 18),
                     labels = c("1", "3","5")) +
  
  scale_y_continuous(limits = c(0,30),
                     expand = c(0,0)) +
  
  ylab("Allocations to distant ingroup") + xlab("Ingroup Fusion") +
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(2)),       
    axis.line = element_line(colour = "black"),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.title=element_text(size=rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(2)),
    axis.text.y= element_text(size = rel(2)),
    axis.text.x= element_text(size = rel(2)),
    plot.margin=unit(c(1,1,1,1),"cm"),
    strip.text.x = element_text(size = rel(2)))

# ggsave("FigureS5.tiff", width = 7.5 ,height = 5.5, dpi = 300)



### Figure S6 ############################
# Line plot

## Re-create the fgo model first
# subset   
dat.x <- na.omit(subset(dat, select = c(OUTGREMO.RD, INGREMO, SITE, COREL.S)))
# model   
summary(fgo <- glmer(cbind(COREL.S, 30-COREL.S) ~ INGREMO * OUTGREMO.RD + (INGREMO * OUTGREMO.RD|SITE),
                     data = dat.x, family = 'binomial',
                     control=glmerControl(optimizer="optimx",optCtrl  = list(method="bobyqa"))))
ss <- getME(fgo,c("theta","fixef"))
summary(update(fgo,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4))))


# extract coefficients from model and logit transform
int1 = exp(fixef(fgo)[1])/(1+exp(fixef(fgo)[1])) # intercept for OUTGREMO.RD == 0
slp1 = exp(fixef(fgo)[2])/(1+exp(fixef(fgo)[2])) - 0.5 # slope for OUTGREMO.RD == 0
int2 = exp(fixef(fgo)[1] + fixef(fgo)[3])/(1 + exp(fixef(fgo)[1] + fixef(fgo)[3]))
slp2 = exp(fixef(fgo)[2] + fixef(fgo)[5])/(1 + exp(fixef(fgo)[2] + fixef(fgo)[5])) - 0.5
int3 = exp(fixef(fgo)[1] + fixef(fgo)[4])/(1 + exp(fixef(fgo)[1] + fixef(fgo)[4]))
slp3 = exp(fixef(fgo)[2] + fixef(fgo)[6])/(1 + exp(fixef(fgo)[2] + fixef(fgo)[6])) - 0.5

# INGREMO range for levels of OUTGREMO.RD
rmn1 = range(dat$INGREMO[dat$OUTGREMO.RD==0], na.rm = T)[1]
rmx1 = range(dat$INGREMO[dat$OUTGREMO.RD==0], na.rm = T)[2]
rmn2 = range(dat$INGREMO[dat$OUTGREMO.RD==1], na.rm = T)[1]
rmx2 = range(dat$INGREMO[dat$OUTGREMO.RD==1], na.rm = T)[2]
rmn3 = range(dat$INGREMO[dat$OUTGREMO.RD==2], na.rm = T)[1]
rmx3 = range(dat$INGREMO[dat$OUTGREMO.RD==2], na.rm = T)[2]


# plot figure
ggplot(data = dat.x, aes(x = INGREMO, y = COREL.S/30, color = OUTGREMO.RD)) +
  geom_line(alpha = 0.001, size = 1.7) + # this is just to set up legend
  
  geom_hline(yintercept = 0.5, lty=2, lwd=1, colour="#736A62", alpha = 0.9) +
  
  geom_segment(x = rmn1,
               xend = rmx1,
               y = int1 + slp1*rmn1,
               yend = int1 + slp1*rmx2,
               color = "#999999",size=1.7) +
  
  geom_segment(x = rmn2,
               xend = rmx2,
               y = int2 + slp2*rmn2,
               yend = int2 + slp2*rmx2,
               color = "#56B4E9",size=1.7) +
  
  geom_segment(x = rmn3,
               xend = rmx3,
               y = int3 + slp3*rmn3,
               yend = int3 + slp3*rmx3,
               color = "#E69F00",size=1.7) +
  
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  
  scale_color_manual(name = "Outgroup relations",
                     values = c("#999999", "#56B4E9","#E69F00"),
                     labels = c("1 & 2", "3 & 4", "5")) +
  
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(1,5), breaks = c(seq(1,5,1))) +
  
  theme_bw() +
  ylab("Probability of Allocation") + xlab("Ingroup Fusion") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(2)),       
    axis.line = element_line(colour = "black"),
    legend.position = c(0.99,0.97),
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.title=element_text(size=rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(2)),
    axis.text.y= element_text(size = rel(2)),
    axis.text.x= element_text(size = rel(2)),
    plot.margin=unit(c(1,1,1,1),"cm"),
    strip.text.x = element_text(size = rel(2)))

# ggsave("FigureS6.tiff", width = 7.5 ,height = 5.5, dpi = 300)
