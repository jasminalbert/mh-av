#working directing is /mhav/code
datLoc <- "../data/"
seedsLoc <- paste0(datLoc, "competition_all_seeds_2021.csv")
seeds <- read.csv(seedsLoc, row.names=1, stringsAsFactors=T)
head(seeds,10)

## notes from ashely ##
#output = seeds/stems
#mean_density_halfm2 = observed density of background competitor in the field (stem counts). some missing
#max_density_halfm2 = maximum possible density calculated based on seed mass (seeding at 8g/m2) of background
str(seeds)
unique(seeds$block) #1,2,3,4
#nut_trt(3): C-compost, F-fertlized, N-none?
#ppt_trt(3): D-drought, W-wet, XC-control
unique(seeds$plot) #1-36
unique(seeds[,c("block","plot")]) #9 plots per block, in order
unique(seeds$subplot) #1-7
dim(unique(seeds[,c("block","plot","subplot")])) #252 rows = 4*9*7
levels(seeds$background) #(7)AVBA, Control, ERBO, HOMU, LOMU, TACA, TRHI
levels(seeds$phyto) #(5)AVBA, ERBO, LOMU, TACA, TRHI,..no HOMU
dim(unique(seeds[,c("phytonum","phyto")])) #seems random, 30 
#role (3): competitor-intra, control-background, invader
seeds[seeds$role=="control",] #one subplot per plot, bg
seeds$id <- 1:nrow(seeds)

#i am interested in avba-taca dynamics in control plots only?
sp <- c("AVBA", "TACA")
spb <- c(sp,"Control")
mas <- seeds[seeds$phyto%in%sp&seeds$background%in%c(sp,"Control"),]
dim(mas) #216,15
3*9*4*2 #two phytos, 3bg==3sub, 9plots, 4 block
sum(is.na(mas$output))
mas[is.na(mas$output),]
#NAs in output and in some rows NAs in bg mean density

#treat NAs in phyto data as zero
startpdat <- which(names(mas)=="tot_stems")
endpdat <- which(names(mas)=="output")
length(mas[,startpdat:endpdat][is.na(mas[,startpdat:endpdat])])==sum(is.na(mas[,startpdat:endpdat]))
mas[,startpdat:endpdat][is.na(mas[,startpdat:endpdat])] <- 0
sum(is.na(mas[,startpdat:endpdat]))==0

#todo: number of NAs per role and phyto treatment
mas[colSums(apply(mas,1,is.na))>0,]
chck <- expand.grid(unique(mas$role),unique(mas$phyto))
for (i in rownames(chck)){
	tmp <-mas[mas$role==chck[i,1]& mas$phyto==chck[i,2],]
	print(chck[i,]);#print(tmp)
	print(nrow(tmp))
	nas <- tmp[rowSums(is.na(tmp))>0,]
	print(nrow(nas))
	print(nas)	
#only NAs in phyto dat and not bg dat when role=control..make sense bc zero	
}
#12 out of 36 NAs for all of four invader+competitor scenarios..12*4=48, 
#all in control nutrient plots?
dim(mas[mas$nut_trt=="N"&mas$role!="control",])
#ALL NA
3*3*4 #reps per invader+competitor
#nutrient,water,block
#reps per nutrient trt
4*1*3*3*2
#block,nutrient,water,bg,phyto

#ALL NO NUTRIENT PLOTS HAVE NA's for mean bg density
#otherwise, no other NA's..

#take out for now?
end<-which(names(mas)=="background")
rvars <- c(var="mean_density_halfm2",trt="ppt_trt",names(mas)[1:end])
pptden<-reshape(mas,timevar=rvars['trt'],idvar=names(mas)[1:end], direction="wide",drop=names(mas)[!names(mas)%in%c(rvars)])
head(pptden)
str(pptden)
colSums(is.na(pptden))
s<- "Control"
boxplot(pptden[pptden$background!=s,-(1:end-1)], names=unique(mas$ppt_trt), main="allBg!Control")
boxplot(pptden[,-(1:end-1)], names=unique(mas$ppt_trt), main="allBg");s="TACA"
boxplot(pptden[pptden$background==s,-(1:end-1)], names=unique(mas$ppt_trt), main=s);s<-"AVBA"
boxplot(pptden[pptden$background==s,-(1:end-1)], names=unique(mas$ppt_trt), main=s); s<- "Control"
boxplot(pptden[pptden$background==s,-(1:end-1)], names=unique(mas$ppt_trt), main=s);


rvars <- c(var="mean_density_halfm2",trt="nut_trt",names(mas)[1:end])
pptden<-reshape(mas,timevar=rvars['trt'],idvar=names(mas)[1:end], direction="wide",drop=names(mas)[!names(mas)%in%c(rvars)])
head(pptden)
str(pptden)
colSums(is.na(pptden))
s<- "Control"
boxplot(pptden[pptden$background!=s,-(1:end-1)], names=unique(mas$nut_trt), main="allBg!Control")
boxplot(pptden[,-(1:end-1)], names=unique(mas$nut_trt), main="allBg");s="TACA"
boxplot(pptden[pptden$background==s,-(1:end-1)], names=unique(mas$nut_trt), main=s);s<-"AVBA"
boxplot(pptden[pptden$background==s,-(1:end-1)], names=unique(mas$nut_trt), main=s); s<- "Control"
boxplot(pptden[pptden$background==s,-(1:end-1)], names=unique(mas$nut_trt), main=s);

# s6 <- seeds[seeds$background==c(sp,"Control"),]
# s5 <- seeds[seeds$phyto==sp & seeds$background==c(sp,"Control"),]
# s2 <-seeds[seeds$phyto==sp,]; dim(s2) #247 
# s3 <- (s2[s2$phyto==sp & s2$background==spb,])
# s4 <- (s2[s2$background==spb,])
# s7 <- seeds[seeds$phyto%in%sp,] ; dim(s7); #507
# unique(s2$background)
# seeds <- 
# str(seeds)
# int <- intersect(s2$id, s7$id)
# length(int)

# head(s2)
# head(s7[s7$id%in%int,])
# head(s7[!s7$id%in%int,])
# unique(s7$phyto);unique(s2$phyto)
# str(s2$phyto)
# str(s7$phyto)



# 1:10 %in% 3:7
# 1:10==3:7
# s <- seeds$phyto[1:20]
# cbind(as.character(s),sp,s==sp)

