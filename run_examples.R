#### Lagged Exact Bayesian Online Changepoint Detection with Parameter Estimation
#### Code submission

# The file contains an example of running the algorithm and visualizing the results

source("LEXO_functions.R")

# Data Generation
x = data_gen_poisson(1,5,1000,4)
totalLength=length(x$Samples[[1]])
# Running algorithm
res = LEXO_poisson(x$Samples[[1]],H=1/250, ell=30)
R= res$RunLength

#pdf("Examples_plot.pdf")
m <- matrix(c(1,2,3,4,5,5),nrow=3,ncol=1,byrow=TRUE)
layout(mat = m,heights = c(0.6,0.6,0.6))

### finding MAP of run length at each time point
MAPmat= NULL
for (i in 1:(ell+1)){
  MAPmat = cbind(MAPmat, apply(R[[i]],2,which.max))
}
lag_to_plot = c(0,10,30) #(lag=0 -> EXO)
plot_colors <- c("black","blue","orange")
# MAP plot
plot(1:length(x$Samples[[1]]), MAPmat[,1], type="n",xlab="time",ylab="Run Length",main="MAP run-length")
for (i in 1:length(lag_to_plot)){
  points(1:(totalLength-lag_to_plot[i]), MAPmat[1:(totalLength-lag_to_plot[i]),lag_to_plot[i]+1],type="l",lwd=2, col=plot_colors[i])
}
# Posterior Mean and plot
MeanPosterior=res$PosMean$Mean
plot(1:length(x$Samples[[1]]), MeanPosterior[,1], type="n",xlab="time",ylab="Posterior Mean",main="Posterior Mean")
for (i in 1:length(lag_to_plot)){
  points(1:(totalLength-lag_to_plot[i]), MeanPosterior[1:(totalLength-lag_to_plot[i]),lag_to_plot[i]+1],type="l",lwd=2, col=plot_colors[i])
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("EXO", "LEXO-10", "LEXO-30"), 
       col=plot_colors, lwd=2, horiz = TRUE)
#dev.off()

