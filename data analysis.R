setwd("C:/Users/Lenovo/Downloads/Thesis")
dataset <- read.csv("Audiomoth acoustics and water parameter_Master table.csv")
bat_activity <- dataset$total.bat.activity
species_richness <- dataset$species.richness
insect_abundance <- dataset$insect_abundance
mean_radiance <- dataset$mean_radiance
dist_citycenter <-dataset$dist_citycenter
size_category<-dataset$size_category
temp <- dataset$temp
rain <- dataset$rain
wind <-dataset$wind
cover_trees<-dataset$cover_trees
cover_bushes<-dataset$cover_bushes
cover_float<-dataset$cover_float
cover_upright<-dataset$cover_upright
cover_upright<-as.numeric(cover_upright)

total_phosphate<-dataset$total_phosphate
nitrate<-dataset$nitrate
carbon<-dataset$carbon
oxygen_content<-dataset$oxygen_content
salinity<-dataset$salinity
impervious_200m<-dataset$impervious_200m
impervious_200m<-as.numeric(impervious_200m)

plot(dataset)
str(dataset)
cor(dataset, method="spearman")
corrplot(dataset)
mydata.cor=cor(dataset)
corrplot(mydata.cor)
cor.mat <- dataset
cor.mat %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  corrplot(label = TRUE)


p_pip<-dataset$Ppip
p_pyg<-dataset$Ppyg
p_nat<-dataset$Pnat
Nyctaloid_grp<-dataset$Nyctaloid_grp
Myotis_grp<-dataset$Myotis_grp
Plecotus<-dataset$Plecotus
grassland_30m<-dataset$grassland_30m
cropland_30m<-dataset$crop_30m
Nnoc<-dataset$Nnoc
Eser<-dataset$Eser




pH<-dataset$pH
Tricoptera_abundance <-dataset$Trichoptera
Ephemeroptera_abundance<-dataset$Ephemeroptera
Ranatra_abundance<-dataset$Renatra
Corixidae_abundance<-dataset$Corixidae
Notonectidae_abundance<-dataset$Notonectidae
Ilyocoris_abundance<-dataset$Ilyocoris
Pleamin_abundance<-dataset$Plea_min
Coleoptera_abundance<-dataset$Coleoptera
Gerromorpha_abundance<-dataset$Gerromorpha
Chironomidae_abundance<-dataset$Chironomidae
Zygoptera_abundance<-dataset$Zygoptera
Plecoptera_abundance<-dataset$Plecoptera
insectorder_richness<-dataset$insectorder_richness
hist(bat_activity)
hist(species_richness)


kruskal.test(total_activity ~ water_bodytype, data = dataset)
boxplot(total_activity~water_bodytype)

Ppip.pca <- prcomp(dataset[,c(5,15:19,21:24, 27:29, 31, 32:35,37, 38, 39)], center = TRUE,scale. = TRUE)

#for model comparison
mean.var.plot = function(M12,M1){
  xb = predict(M1)
  g = cut(xb, breaks=unique(quantile(xb,seq(0,1,0.1))))
  m = tapply(M12$y, g, mean)
  v = tapply(M12$y, g, var)
  pr <- residuals(M12,"pearson")
  phi <- sum(pr^2)/df.residual(M12)
  x = seq(min(m),max(m),length.out = 100)
  line.data = data.frame(x=rep(x,2),y=c(x*phi,x*(1+x/M1$theta)),
                         model=c(rep("Q. Poisson",length(x)),rep("Neg. Binomial",length(x))))
  library(ggplot2)
  ggplot() + geom_point(aes(x=m,y=v)) + 
    geom_line(aes(x=x,y=y,linetype=model),data=line.data) + 
    theme_bw() + theme(panel.background = element_rect(rgb(.95,.95,.95))) +
    ylab("variance") + xlab("mean") +
    scale_linetype_manual(values = c("solid","dashed")) +
    ggtitle("Mean-Variance Relationship") 
}
#glmulti
res <- glmulti(total_activity ~ insect_abundance+insectorder_richness+cover_float+cover_upright+salinity+oxygen_content+phosphate+nitrate+carbon, data=dataset,
               level=1, fitfunction=rma.glmulti, crit="aicc", confsetsize=128)

#water parameters(total activity, drop1 method)

M1.nb<-glm.nb(species_richness~insect_abundance+size_category+cover_float+cover_upright+oxygen_content+total_phosphate+nitrate+carbon+temp+wind+rain, data=dataset)
warning(M1.nb)
res.nb <- residuals(M1.nb)
hist(res.nb)
summary(M1.nb)
plot(M1.nb)


M1<-glm(bat_activity~insect_abundance+size_category+cover_float+cover_upright+oxygen_content+total_phosphate+nitrate+carbon+temp+rain+wind, family="quasipoisson")
res1 <- residuals(M1)
hist(res1)
summary(M1)
plot(M1)
drop1(M1, test="F")
vif(M1)

drop1(M1, test="F")
M2<-glm(bat_activity~insect_abundance+size_category+cover_float+cover_upright+oxygen_content+total_phosphate+nitrate+carbon+temp+rain, family="quasipoisson")
res2<-residuals(M2)
hist(res2)
summary(M2)
plot(M2)
drop1(M2, test="F")

M3<-glm(bat_activity~insect_abundance+cover_float+cover_upright+oxygen_content+total_phosphate+nitrate+carbon+temp+rain, family="quasipoisson")
res3<-residuals(M3)
hist(res3)
summary(M3)
plot(M3)
drop1(M3, test="F")

M4<-glm(bat_activity~insect_abundance+cover_float+cover_upright+oxygen_content+total_phosphate+nitrate+temp+rain, family="quasipoisson")
res4<-residuals(M4)
hist(res4)
summary(M4)
plot(M4)
drop1(M4, test="F")

M5<-glm(bat_activity~insect_abundance+cover_float+cover_upright+total_phosphate+nitrate+temp+rain, family="quasipoisson")
res5<-residuals(M5)
hist(res5)
summary(M5)
plot(M5)
drop1(M5, test="F")

M6<-glm(bat_activity~insect_abundance+cover_float+total_phosphate+nitrate+temp+rain, family="quasipoisson")
res6<-residuals(M6)
hist(res6)
summary(M6)
plot(M6)
drop1(M6, test="F")


M7<-glm(bat_activity~insect_abundance+cover_float+total_phosphate+nitrate+temp, family="quasipoisson")
res7<-residuals(M7)
hist(res7)
summary(M7)
plot(M7)
drop1(M7, test="F")




anova (M1, M2, M3, M4, M5, M6,M7, test= "F")


ggplot(dataset,aes(y=bat_activity,x=insect_abundance))+geom_point()+geom_smooth(method="glm", grid="False")

ggplot(radial,aes(y=bat_activity,x=insect_abundance))+geom_point()+geom_smooth(method="lm")

ggPredict(M7,se=TRUE,interactive=TRUE)

mydf <- ggpredict(M7, term="insect_abundance")

plot(mydf)


mydf <-ggpredict(M7, c("nitrate", "bat_activity"))
plot(mydf)

ggplot(dataset, aes(bat_activity, insect_abundance)) + geom_point() +   geom_smooth(method="glm", method.args=list(family="quasipoisson"), fullrange=T, se=F) +   xlim(0, 700)

plot(ggpredict(M2, term="size_category"), rawdata = T)+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() 

                                                              predicted_glm <- predictvals(M7, xvar = "insect_abundace", yvar = "bat_activity")                                                         
                                                              
ggpredict(M7)                                                           ggpredict(M7, pred = insect_abundance, interval = TRUE)

plot<-ggplot(data=dataset, aes(x=x,y=y,col=species)) +
  geom_point(mapping=aes(x=x,y=y, shape=SEX), size=3)+
  scale_shape_manual(values=c(8, 20))+
  geom_smooth(mapping=aes(x=x,y=y),method=lm, se=FALSE)
plot+ theme_bw()

                                                                 
#water parameters species richness 
m.nb<-glm.nb(species_richness~insect_abundance+cover_float+cover_upright+salinity+oxygen_content+phosphate+nitrate+carbon+temp+wind+rain, link="log", data = dataset)
resm.nb <- residuals(m.nb)
hist(resm.nb)
summary(m.nb)
plot(m.nb)

m.p<-glm(species_richness~insect_abundance+cover_float+cover_upright+salinity+oxygen_content+phosphate+nitrate+carbon+temp+wind+rain, family="poisson", data = dataset)
resm.p <- residuals(m.p)
hist(resm.p)
summary(m.p)
plot(m.p)

m.z<-zeroinfl(species_richness~insect_abundance+cover_float+cover_upright+salinity+oxygen_content+phosphate+nitrate+carbon+temp+wind+rain, data = dataset)
resm.z <- residuals(m.z)
hist(resm.z)
summary(m.z)
plot(m.z)

mean.var.plot = function(m.z,m.nb){
  xb = predict(m.nb)
  g = cut(xb, breaks=unique(quantile(xb,seq(0,1,0.1))))
  m = tapply(m.z$y, g, mean)
  v = tapply(m.z$y, g, var)
  pr <- residuals(m.z,"pearson")
  phi <- sum(pr^2)/df.residual(m.z)
  x = seq(min(m),max(m),length.out = 100)
  line.data = data.frame(x=rep(x,2),y=c(x*phi,x*(1+x/m.nb$theta)),
                         model=c(rep("zero inflated",length(x)),rep("Neg. Binomial",length(x))))
  library(ggplot2)
  ggplot() + geom_point(aes(x=m,y=v)) + 
    geom_line(aes(x=x,y=y,linetype=model),data=line.data) + 
    theme_bw() + theme(panel.background = element_rect(rgb(.95,.95,.95))) +
    ylab("variance") + xlab("mean") +
    scale_linetype_manual(values = c("solid","dashed")) +
    ggtitle("Mean-Variance Relationship") 
}

mean.var.plot(m.z,m.nb)

ggPredict(M7, pred = insect_abundance, interval = TRUE, plot.points = TRUE)

ggeffect(M7, pred = insect_abundance, interval = TRUE, partial.residuals = TRUE,
            jitter = c(0.1, 0))



mydf <- ggpredict(M2, term="size_category", show.data = TRUE)
plot(mydf, rawdata = T)




ggpredict(M5, terms = "insect_abundance", xpos=0.75,vjust=1.5,show.error = TRUE)
plot(mydf, )
ggPredict(M5,xpos=0.75,vjust=1.5,show.error = TRUE)

mydf <- ggpredict(M1, terms = "insect_abundance")
ggplot(mydf, aes(x, predicted)) +
  geom_line() + geom_point()+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)




require(ggiraph)
require(ggiraphExtra)
require(plyr)
ggPredict(m1,se=TRUE)

require(moonBook)
require(ggplot2)
require(ggiraph)
require(plyr)


library (ggiraphExtra)
library(ggplot2)


install.packages("ellipsis")

require(devtools)


#species richness and water parameters

m1<-glm(species_richness~insect_abundance+size_category+cover_float+cover_upright+oxygen_content+total_phosphate+nitrate+carbon+temp+wind+rain, family="quasipoisson")
resm1 <- residuals(m1)
hist(resm1)
summary(m1)
plot(m1)


 
drop1(m1, test="F")
m2<-glm(species_richness~insect_abundance+size_category+cover_float+cover_upright+oxygen_content+total_phosphate+nitrate+temp+wind+rain, family="quasipoisson")
resm2 <- residuals(m2)
hist(resm2)
summary(m2)
plot(m2)

drop1(m2, test="F")
m3<-glm(species_richness~insect_abundance+size_category+cover_float+cover_upright+total_phosphate+nitrate+temp+rain+wind, family="quasipoisson")
resm3 <- residuals(m3)
hist(resm3)
summary(m3)
plot(m3)

drop1(m3, test="F")
m4<-glm(species_richness~insect_abundance+size_category+cover_float+cover_upright+total_phosphate+nitrate+temp+wind, family="quasipoisson")
resm4 <- residuals(m4)
hist(resm4)
summary(m4)
plot(m4)

drop1(m4, test="F")
m5<-glm(species_richness~insect_abundance+size_category+cover_float+total_phosphate+nitrate+temp+wind, family="quasipoisson")
resm5 <- residuals(m5)
hist(resm5)
summary(m5)
plot(m5)


drop1(m5, test="F")
m6<-glm(species_richness~size_category+cover_float+total_phosphate+nitrate+temp+wind, family="quasipoisson")
resm6 <- residuals(m6)
hist(resm6)
summary(m6)
plot(m6)


drop1(m6, test="F")

anova(m1, m2, m3, m4, m5, m6,m7, test="F")



m7<-glm(species_richness~size_category+cover_float+nitrate+temp+wind, family="quasipoisson")
resm7 <- residuals(m7)
hist(resm7)
summary(m7)
plot(m7)


#urban parameters total bat activity
l<-glm(bat_activity~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+grassland_30m+temp+wind+rain, family="quasipoisson", data=dataset)
resl <- residuals(l)
hist(resl)
summary(l)
plot(l)

drop1(l, test="F")
l1<-glm(bat_activity~dist_citycenter+mean_radiance+cover_trees+cover_bushes+grassland_30m+cropland_30m+temp+wind+rain, family="quasipoisson", data=dataset)
resl1 <- residuals(l1)
hist(resl1)
summary(l1)

drop1(l1, test="F")
l2<-glm(bat_activity~mean_radiance+cover_trees+cover_bushes+grassland_30m+cropland_30m+temp+rain+wind, family="quasipoisson", data=dataset)
resl2 <- residuals(l2)
hist(resl2)
summary(l2)

drop1(l2, test ="F")

l3<-glm(bat_activity~mean_radiance+cover_trees+cover_bushes+grassland_30m+cropland_30m+temp+rain, family="quasipoisson", data=dataset)
resl3 <- residuals(l3)
hist(resl3)
summary(l3)
plot(l3)

drop1(l3, test="F")
l4<-glm(bat_activity~mean_radiance+cover_trees+cover_bushes+grassland_30m+temp+rain, family="quasipoisson", data=dataset)
resl4 <- residuals(l4)
hist(resl4)
summary(l4)
plot(l4)

drop1(l4, test="F")
l5<-glm(bat_activity~mean_radiance+cover_trees+cover_bushes+grassland_30m+rain, family="quasipoisson", data=dataset)
resl5 <- residuals(l5)
hist(resl5)
summary(l5)
plot(l5)

drop1(l5, test="F")
l6<-glm(bat_activity~cover_trees+cover_bushes+grassland_30m+rain, family="quasipoisson", data=dataset)
resl6 <- residuals(l6)
hist(resl6)
summary(l6)
plot(l6)

l7<-glm(bat_activity~cover_trees+cover_bushes+grassland_30m, family="quasipoisson", data=dataset)
resl7 <- residuals(l7)
hist(resl7)
summary(l7)
plot(l7)


anova(l,l1,l2,l3,l4,l5,l6, test="F")


plot<- ggpredict(M1, terms = "temp")

#URBAN PARAMETERS AND SPECIES RICHNESS

plot(ggpredict(U, term="cover_trees"), rawdata = T)+ theme(panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank())


U.nb<-glm.nb(species_richness~mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+grassland_30m+temp+wind+rain, link = "log", data=dataset)
resU.nb <- residuals(U.nb)
hist(resU.nb)
summary(U.nb)
plot(U.nb)


U<-glm(species_richness~dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resU <- residuals(U)
hist(resU)
summary(U)
plot(U)
drop1(U, test="F")



backwards = step(U)
backwards = step(U,trace=0) 
formula(backwards)
summary(backwards)
forwards = step(U,
                scope=list(lower=formula(nothing),upper=formula(fullmod)), direction="forward")


drop1(U,test="F" )
U1<-glm(species_richness~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resU1 <- residuals(U1)
hist(resU1)
summary(U1)
plot(U1)

drop1(U1,test="F" )
U2<-glm(species_richness~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+temp+wind+rain, family = "quasipoisson", data=dataset)
resU2 <- residuals(U2)
hist(resU2)
summary(U2)
plot(U2)

drop1(U2,test="F" )
U3<-glm(species_richness~ dist_citycenter+impervious_200m+cover_trees+cover_bushes+temp+wind+rain, family = "quasipoisson", data=dataset)
resU3 <- residuals(U3)
hist(resU3)
summary(U3)
plot(U3)

drop1(U3,test="F" )
U4<-glm(species_richness~ dist_citycenter+impervious_200m+cover_trees+cover_bushes+temp+wind, family = "quasipoisson", data=dataset)
resU4 <- residuals(U4)
hist(resU4)
summary(U4)
plot(U4)

drop1(U4,test="F" )
U5<-glm(species_richness~dist_citycenter+cover_trees+cover_bushes+temp+wind, family = "quasipoisson", data=dataset)
resU5 <- residuals(U5)
hist(resU5)
summary(U5)
plot(U5)
anova(U, U1, U2, U3, U4, U5,U6, test="F")

drop1(U5,test="F" )
U6<-glm(species_richness~cover_trees+cover_bushes+temp+wind, family = "quasipoisson", data=dataset)
resU6 <- residuals(U6)
hist(resU6)
summary(U6)
plot(U6)



plot(ggpredict(U, term="cropland"), rawdata = T)+ theme(panel.grid.major = element_blank(), 
                                                             panel.grid.minor = element_blank())


#ppip
p<-glm(p_pip~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resp <- residuals(p)
hist(resp)
summary(p)
plot(p)

drop1(p, test="F")
p1<-glm(p_pip~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resp1 <- residuals(p1)
hist(resp1)
summary(p1)
plot(p1)

drop1(p1, test="F")
p2<-glm(p_pip~dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+temp+wind+rain, family = "quasipoisson", data=dataset)
resp2 <- residuals(p2)
hist(resp2)
summary(p2)
plot(p2)

drop1(p2, test="F")
p3<-glm(p_pip~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+temp+wind, family = "quasipoisson", data=dataset)
resp3 <- residuals(p3)
hist(resp3)
summary(p3)
plot(p3)

drop1(p3, test="F")
p4<-glm(p_pip~ mean_radiance+impervious_200m+cover_trees+cover_bushes+temp+wind, family = "quasipoisson", data=dataset)
resp4 <- residuals(p4)
hist(resp4)
summary(p4)
plot(p4)

drop1(p4, test="F")

anova(p,p1,p2,p3,p4,p5, test="F")

p5<-glm(p_pip~ mean_radiance+impervious_200m+cover_trees+cover_bushes+temp, family = "quasipoisson", data=dataset)
resp5 <- residuals(p5)
hist(resp5)
summary(p5)
plot(p5)


ggcoef(p4)

plot(ggpredict(n, term="mean_radiance"), rawdata = T)+ theme(panel.grid.major = element_blank(), 
                                                                    panel.grid.minor = element_blank())





m4.emm <- emmeans(p4, specs= "mean_radiance", "cover_bushes", "grassland" "temp" "wind")
m4.emm



gg.effects <- ggdotplot(data = p_pip,
                        x="contrast",
                        y="activity",
                        color = "steelblue",
                        fill = "steelblue",
                        size=0.5) +
  geom_errorbar(aes(x=contrast,
                    ymin=asymp.LCL,
                    ymax=asymp.UCL),
                width=0.15,
                color="steelblue") 


#ppyg
y<-glm(p_pyg~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resy <- residuals(y)
hist(resy)
summary(y)
plot(y)

drop1(y, test="F")
y1<-glm(p_pyg~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resy1 <- residuals(y1)
hist(resy1)
summary(y1)
plot(y1)


drop1(y1, test="F")
y2<-glm(p_pyg~ mean_radiance+impervious_200m+cover_trees+cover_bushes+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resy2 <- residuals(y2)
hist(resy2)
summary(y2)
plot(y2)

drop1(y2, test="F")
y3<-glm(p_pyg~ mean_radiance+cover_trees+cover_bushes+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resy3 <- residuals(y3)
hist(resy3)
summary(y3)
plot(y3)



drop1(y3,test="F")

y4<-glm(p_pyg~ mean_radiance+cover_trees+cover_bushes+grassland_30m+temp+wind, family = "quasipoisson", data=dataset)
resy4 <- residuals(y4)
hist(resy4)
summary(y4)
plot(y4)

drop1(y4, test="F")


anova(y, y1, y2, y3, y4,y5, test="F")

y5<-glm(p_pyg~mean_radiance+cover_trees+cover_bushes+temp+wind, family = "quasipoisson", data=dataset)
resy5 <- residuals(y5)
hist(resy5)
summary(y5)
plot(y5)
 
drop1(y5, test="F")


plot(ggpredict(G1, term="cover_trees"), rawdata = T)+ theme(panel.grid.major = element_blank(), 
                                                              panel.grid.minor = element_blank())




#PNAT
n<-glm(p_nat~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resn <- residuals(n)
hist(resn)
summary(n)

drop1(n, test="F")
n1<-glm(p_nat~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resn1 <- residuals(n1)
hist(resn1)
summary(n1)

drop1(n1, test="F")
n2<-glm(p_nat~ dist_citycenter+mean_radiance+cover_trees+cover_bushes+grassland_30m+temp+rain+wind, family = "quasipoisson", data=dataset)
resn2 <- residuals(n2)
hist(resn2)
summary(n2)
plot(n2)

drop1(n2, test="F")
n3<-glm(p_nat~ dist_citycenter+mean_radiance+cover_trees+cover_bushes+grassland_30m+temp+rain, family = "quasipoisson", data=dataset)
resn2 <- residuals(n2)
hist(resn2)
summary(n2)
plot(n2)

drop1(n3, test="F")
n4<-glm(p_nat~ dist_citycenter+mean_radiance+cover_trees+cover_bushes+temp, family = "quasipoisson")
resn4 <- residuals(n3)
hist(resn4)
summary(n4)
plot(n3)

anova(n,n1,n2,n3,n4, test="F")


plot(ggpredict(n, term="distance_citycenter"), rawdata = T)+ theme(panel.grid.major = element_blank(), 
                                                             panel.grid.minor = element_blank())

mydf <- ggpredict(n, term="distance_citycenter", show.data = TRUE)
plot(mydf, rawdata = T)



#Myotis 
G1<-glm(Myotis_grp~dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+grassland_30m+cropland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resG1 <- residuals(G1)
hist(resG1)
summary(G1)
plot(G1)


drop1(G1, test="F")
G2<-glm(Myotis_grp~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resG2 <- residuals(G2)
hist(resG2)
summary(G2)

drop1(G2, test="F")
G3<-glm(Myotis_grp~ dist_citycenter+mean_radiance+impervious_200m+cover_bushes+grassland_30m+temp+wind, family = "quasipoisson", data=dataset)
resG3 <- residuals(G3)
hist(resG3)
summary(G3)

drop1(G3, test="F")

anova(G1, G2, G3,G4, G5, test="F")

G4<-glm(Myotis_grp~ dist_citycenter+mean_radiance+impervious_200m+cover_bushes+grassland_30m+temp, family = "quasipoisson", data=dataset)
resG2 <- residuals(G2)
hist(resG2)
summary(G2)

drop1(G4, test="F")
G5<-glm(Myotis_grp~ dist_citycenter+mean_radiance+impervious_200m+cover_bushes+grassland_30m, family = "quasipoisson", data=dataset)
resG5 <- residuals(G5)
hist(resG5)
summary(G5)


plot(ggpredict(G, term="cover_bushes"), rawdata = T)+ theme(panel.grid.major = element_blank(), 
                                                             panel.grid.minor = element_blank())


#Nyctaloid

N<-glm(Nyctaloid_grp~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+grassland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resN <- residuals(N)
hist(resN)
summary(N)

drop1(N, test="F")
N1<-glm(Nyctaloid_grp~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+temp+wind+rain, family = "quasipoisson", data=dataset)
resN1 <- residuals(N1)
hist(resN1)
summary(N1)

drop1(N1, test="F")
N2<-glm(Nyctaloid_grp~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+temp+rain+wind, family = "quasipoisson", data=dataset)
resN2 <- residuals(N2)
hist(resN2)
summary(N2)

drop1(N2, test="F")
N3<-glm(Nyctaloid_grp~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+temp+rain, family = "quasipoisson", data=dataset)
resN3 <- residuals(N3)
hist(resN3)
summary(N3)

drop1(N3, test="F")
N4<-glm(Nyctaloid_grp~ dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+rain, family = "quasipoisson", data=dataset)
resN4 <- residuals(N4)
hist(resN4)
summary(N4)
drop1(N4, test="F")

N5<-glm(Nyctaloid_grp~ dist_citycenter+impervious_200m+cover_trees+cover_bushes+rain, family = "quasipoisson", data=dataset)
resN5 <- residuals(N5)
hist(resN5)
summary(N5)

drop1(N5, test="F")

N6<-glm(Nyctaloid_grp~ dist_citycenter+impervious_200m+cover_trees+cover_bushes, family = "quasipoisson", data=dataset)
resN6 <- residuals(N6)
hist(resN6)
summary(N6)


anova(N, N1, N2, N3, N4, N5, N6, test="F")

plot(ggpredict(G1, term="mean_radiance"), rawdata = F)+ theme(panel.grid.major = element_blank(), 
                                                             panel.grid.minor = element_blank())


#Plecotus
cmp.out = glm.cmp(Plecotus~distance_citycenter+mean_radiance+cover_trees+cover_bushes+cropland+grassland+temp+wind+rain, method = optim.method, control = optim.control) 
print(cmp.out)

P<-glm(Plecotus~dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+grassland_30m+temp+wind+rain ,family="quasipoisson", data=dataset)
resP <- residuals(P)
hist(resP)
summary(P)
plot(P)


drop1(P, test="F")
P1<-glm(Plecotus~dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+cropland_30m+temp+wind+rain, family= "quasipoisson", data=dataset)
resP1 <- residuals(P1)
hist(resP1)
summary(P1)
plot(P1)

drop1(P1, test="F")
P2<-glm(Plecotus~dist_citycenter+mean_radiance+cover_trees+cover_bushes+cropland_30m+temp+wind+rain, family= "quasipoisson", data=dataset)
resP2 <- residuals(P2)
hist(resP2)
summary(P2)

drop1(P2, test="F")
P3<-glm(Plecotus~dist_citycenter+cover_trees+cover_bushes+cropland_30m+temp+wind+rain, family= "quasipoisson", data=dataset)
resP3 <- residuals(P3)
hist(resP3)
summary(P3)



anova(P,P1, P2, P3, test="F")

drop1(P3, test="F")
P4<-glm(Plecotus~dist_citycenter+cover_trees+cropland_30m+temp+rain, family= "quasipoisson", data=dataset)
resP4 <- residuals(P4)
hist(resP4)
summary(P4)

drop1(P4, test="F")
P5<-glm(Plecotus~cover_trees+cropland_30m+temp+rain, family= "quasipoisson", data=dataset)
resP5 <- residuals(P5)
hist(resP5)
summary(P5)

drop1(P5, test="F")
P6<-glm(Plecotus~cover_trees+temp, family= "quasipoisson", data=dataset)
resP6 <- residuals(P6)
hist(resP6)
summary(P6)

anova(P, P1, P2, P3, P4, P5, P6, test = "F")

#insect abundance and other parameters

i<-glm( insect_abundance ~ + salinity + oxygen_content + nitrate + cover_float, 
    family = "quasipoisson")
res_i<-residuals(i)
hist(res_i)
summary(i)
vif(i)

i<-glm(insect_abundance~size_category+oxygen_content+carbon+mean_radiance+cover_bushes+cover_float+cover_trees+co+grassland, family= "quasipoisson")
res_i<-residuals(i)
hist(res_i)
summary(i)

drop1(i, test= "F")


p<-glm(species_richness~ insect_order+Tricoptera_abundance +Ephemeroptera_abundance+Ranatra_abundance+Corixidae_abundance+Notonectidae_abundance+Ilyocoris_abundance+Pleamin_abundance+Coleoptera_abundance+Gerromorpha_abundance+Chironomidae_abundance+Zygoptera_abundance+Plecoptera_abundance, family = "quasipoisson")
summary(p)


#negative binomial
M1nb<-glm.nb(total_activity~waterbody_length+insect_abundance+cover_float+cover_upright+salinity+conductivity+oxygen_content+phosphate+nitrate+carbon, link="log", data=dataset)
summary(M1nb)
vif(M1nb)

M.nb<-glm.nb(total_activity~waterbody_length+insect_abundance+cover_float+cover_upright+phosphate+nitrate, link="log")


M.nb<-glm.nb(total_activity~waterbody_length+insect_abundance+cover_float+cover_upright+phosphate+nitrate, link="log", control = glm.control(maxit = 2500), init.theta = 1.0)
summary(M.nb)

#species richness (negative binomial)
s.nb<-glm.nb(species_richness~waterbody_length+cover_float+cover_upright+phosphate+nitrate, link="log")
resnb <- residuals(s.nb)
hist(resnb)
summary(s.nb)
plot(s.nb)
vif(s.nb)

drop1(s.nb, test="Chi")



#PLOTS

ggplot(dataset, aes(x=mean_radiance, y=total_activity)) +
  geom_point(shape=1) +    
  geom_smooth(method=lm)   

ggplot(dataset, aes(x=mean_radiance, y=species_richness)) +
  geom_point(shape=1) +    
  geom_smooth(method=lm) 

ggplot(dataset, aes(x=distance_citycenter, y=Plecotus)) +
  geom_point(shape=1) +    
  geom_smooth(method=lm) 
ggplot(dataset, aes(x=waterbody_length, y=total_activity)) + geom_boxplot()
boxplot(total_activity~waterbody_length)

plot(dataset$Total.bat.activity,dataset$Length.m.,pch=16,col=dataset$Type.of.waterbody)


myplot<-ggplot(dataset,aes(x=size_category,y=bat_activity,col=water_bodytype, size=6))+geom_point()

myplot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(dataset, aes(x=size_category, y=bat_activity,)) + geom_point() 
   
boxplot(bat_activity~size_category)

library(ggplot2)
mydf <- ggpredict(N2, terms = "distance_citycenter")
ggplot(mydf, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = p_nat, ymax = distance_citycenter), alpha = .1)

install.packages("devtools")


a<-glm(insect_abundance~size_category+cover_upright+nitrate+dist_citycenter+mean_radiance+impervious_200m+cover_trees+cover_bushes+grassland_30m, family="quasipoisson")
summary(a)

drop1 (a, test="F")
install.packages("BiodiversityR")
library(BiodiversityR)
