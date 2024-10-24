#setwd("/Users/jasminalbert/Documents/HALLETTLAB/mhav/code")
library(minpack.lm)
datLoc <- "../data/"
filename <- "MHAVseeds.RDS"
mas <- readRDS(paste0(datLoc,filename))
#mas <- mas[mas$nut_trt!="N",] 	
	#if using max stems out, getting Ntrt is fine
str(mas)
max(mas$tot_stems)
#seedsIn =3, but max stems out is 4
#stems in stems out model or seedsind seedsout
#seeds in: have for phytos but background is in grams seeds
#seeds out:have for phytos but background is stems out
#stems in: no
#stems out: yes for both
#can i do seeds in seeds out but alpha are not symmetric units?
#so alpha_i,j is effects of stems of j - no seed in timet
#seeds in for phytos and bg to seeds out for phyto 
#then i wont even use the background data?

#seeds in for bg = 8g/m2
#subplots are 0.5x0.5 m2 = 0.25m2
#=8*.25 = 2g/m2
#MH one seed = 0.0058g (from JS)
2/0.0058 #=344.8276 seeds in
#AV: .015g
2/0.015 #=133.3333 seedsin
#AS expected 530.504 and 116.959 max stems out from seeding rate
#so that means she has one seed = 
2/530.504 #=0.00377 MH g/seed
2/116.959 #=0.01710001 AV g/seed  
#this is for the whole subplot but the plot was divided into six subsubplots?
#divide maxstemsout by 6 to get bg seeds in?




