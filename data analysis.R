--Setting up the working directory--

setwd("C:/Users/Lenovo/Downloads/Thesis")

--loading the dataset--
dataset <- read.csv("Audiomoth acoustics and water parameter_Master table.csv")
bat_activity <- dataset$total.bat.activity
species_richness <- dataset$species.richness
insect_abundance <- dataset$insect_abundance

mean_radiance <- dataset$mean_radiance
dist_citycenter <-dataset$dist_citycenter
size_category<-dataset$size_category
cover_trees<-dataset$cover_trees
cover_bushes<-dataset$cover_bushes
cover_float<-dataset$cover_float
cover_upright<-dataset$cover_upright
cover_upright<-as.numeric(cover_upright) -- changing the data into numeric--
grassland_30m<-dataset$grassland_30m
cropland_30m<-dataset$crop_30m

temp <- dataset$temp
rain <- dataset$rain
wind <-dataset$wind

total_phosphate<-dataset$total_phosphate
nitrate<-dataset$nitrate
carbon<-dataset$carbon
oxygen_content<-dataset$oxygen_content
salinity<-dataset$salinity
impervious_200m<-dataset$impervious_200m
impervious_200m<-as.numeric(impervious_200m)

--checking the summary of the dataset--
output <- summary(dataset)
print(output)

--Creating new dataframe from existing column--
df<-dataset[,c('mean_radiance', 'dist_citycenter','impervious_200m', 'insect_abundance', 'insect_orderrichness', 'cover_tress','cover_bushes', 'cover_float', 'cover_upright', 'size.category', 'salinity, 'oxygen_content', 'total_phosphate', 'nitrate', 'carbon')]

-- checking the colineratity of the data--
cor(df, method="spearman")
corrplot(df)

p_pip<-dataset$Ppip
p_pyg<-dataset$Ppyg
p_nat<-dataset$Pnat
Nyctaloid_grp<-dataset$Nyctaloid_grp
Myotis_grp<-dataset$Myotis_grp
Plecotus<-dataset$Plecotus
Nnoc<-dataset$Nnoc
Eser<-dataset$Eser

hist(bat_activity)
hist(species_richness)

library(ggplot2) --(getting the ggplot2 pacakage for plotting)--
require(ggiraph)
require(ggiraphExtra)
require(plyr)
#water parameters(total activity, drop1 method)

M1<-glm(bat_activity~insect_abundance+size_category+cover_float+cover_upright+oxygen_content+total_phosphate+nitrate+carbon+temp+rain+wind, family="quasipoisson")
res1 <- residuals(M1)
hist(res1)
summary(M1)
plot(M1)
drop1(M1, test="F")

--kept on droping until all the parametrs were significant)
--Then used the ANOVA command to compare the models with an analysis of deviance--

anova (M1, M2, M3, M4, M5, M6,M7, test= "F")

--Trying out different plotting method--
ggplot(dataset,aes(y=bat_activity,x=insect_abundance))+geom_point()+geom_smooth(method="glm", grid="False")

mydf <- ggpredict(M2, term="size_category", show.data = TRUE)
plot(mydf, rawdata = T)

ggPredict(M7,se=TRUE,interactive=TRUE)
mydf <- ggpredict(M7, term="insect_abundance")
plot(mydf)

---CONTNIUED THE SAME MODELING METHOD FOR SPECIES RICHNESS AND BAT ACTIVITY AMONG OTHER PARAMETERS--
                                                                 
--TRYING OTHER MODELING METHOD--
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


--for model comparison--
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
