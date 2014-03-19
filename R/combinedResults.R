#
#
#

clusterResults = read.csv(dataFile("clusterRandomizedResults.csv"))
clusterResults$deviation = clusterResults$approxPower - clusterResults$empiricalPower
clusterResultsCount = length(clusterResults$deviation)

longitudinalResults = read.csv(dataFile("longitudinalResults.csv"))
longitudinalResults$deviation = longitudinalResults$approxPower - longitudinalResults$empiricalPower
longitudinalResultsCount = length(longitudinalResults$deviation)

combinedResults = data.frame(deviation=c(clusterResults$deviation,longitudinalResults$deviation))
combinedResultsCount = length(combinedResults$deviation)

powerResults = data.frame(deviation=c(clusterResults$deviation,longitudinalResults$deviation, 
                       c(clusterResults$deviation,longitudinalResults$deviation)),
           group=c(rep("Cluster Randomized", clusterResultsCount),
                   rep("Longitudinal", longitudinalResultsCount),
                   rep("All Designs", combinedResultsCount)))
mean(combinedResults$deviation)
range(combinedResults$deviation)
fivenum(combinedResults$deviation)

pdf(file="../inst/figures/PowerBoxPlot_overall.pdf")
boxplot(powerResults$deviation ~ powerResults$group, ylim=c(-0.1,0.1))
abline(h=0,lty=3)
dev.off()



