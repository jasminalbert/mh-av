fitModels <- function(models,seed){
	MHoutput <- as.data.frame(matrix(nrow = 0, ncol = 6))
	names(MHoutput) = c("estimate", "se", "t", "p", "params", "species")

	m1out <- nlsLM(models[["mh"]], start=list(lambda=1, aiM = .01, aiA=.01), upper = c(200, 1, 1), control=nls.lm.control(maxiter=500), trace=T, data = seed[!is.na(seed$TACAOut),])
  
	outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  	names(outreport) = c("estimate", "se", "t", "p")
  	outreport$params <- row.names(outreport)
  	outreport$species <- "MH"
  	MHoutput <- rbind(MHoutput, outreport)
#took out constrant :    lower = c(0, 0, 0),
	AVoutput <- as.data.frame(matrix(nrow = 0, ncol = 6))	
	names(AVoutput) = c("estimate", "se", "t", "p", "params", "species")
	m1out <- nlsLM(models[["av"]], start=list(lambda=1, aiM = .01, aiA=.01), upper = c(200, 1, 1),control=nls.lm.control(maxiter=500), trace=T, data = seed[!is.na(seed$AVBAOut),])
  
 	outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  	names(outreport) = c("estimate", "se", "t", "p")
  	outreport$params <- row.names(outreport)
  	outreport$species <- "AV"
  	AVoutput <- rbind(AVoutput, outreport)
#MHoutput
	modelOut <- rbind(AVoutput,MHoutput)
	return(modelOut)
}