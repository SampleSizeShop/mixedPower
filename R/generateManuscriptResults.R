# 
#  Package mixedPower calculates power for the linear mixed model
#  Copyright (C) 2013 University of Colorado Denver.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Streetf, Fifth Floor, Boston, MA  02110-1301, USA.
#


#' generateManuscriptResults
#' 
#' This function reproduces the simulation study results for the manuscript
#' in the README file
#'
#' @param study.seed the random number seed (defaults to 7634)
#' @param study.data.dir the directory into which data files are written (defaults to
#' current working directory)
#' @param study.figures.dir the directory into which pdf figures are written (defaults
#' to the current working directory)
#' @param study.runEmpirical if true, empirical power values will be recalculated. If false,
#' existing empirical power values will be loaded from the R package.
#' 
#' @note 
#' The empirical power calculations may take several hours to run
#'
runSimulationStudy <- function(study.data.dir=getwd(), 
                               study.figures.dir=getwd(),
                               study.runEmpirical=FALSE) {
  
  # run longitudinal designs
  cat("### Calculating power for longitudinal designs\n")
  calculatePower.longitudinal(data.dir=study.data.dir, 
                              figures.dir=study.figures.dir,
                              runEmpirical=study.runEmpirical)
  cat("### Summarizing longitudinal results\n")
  summarizeResults.longitudinal(data.dir=study.data.dir, 
                                figures.dir=study.figures.dir)
  
  # run cluster randomized designs
  cat("### Calculating power for cluster randomized designs\n")
  calculatePower.clusterRandomized(data.dir=study.data.dir, 
                              figures.dir=study.figures.dir,
                              runEmpirical=study.runEmpirical)
  cat("### Summarizing cluster randomized results\n")
  summarizeResults.clusterRandomized(data.dir=study.data.dir, 
                                figures.dir=study.figures.dir)

  
  cat("### Generating combined results")
  powerResults.clusterRandomized = 
    read.csv(file.path(study.data.dir, 
                       "clusterRandomizedResults.csv"))[,c("empiricalPower",
                                                           "approxPower")]
  powerResults.longitudinal = 
    read.csv(file.path(study.data.dir, 
                       "longitudinalResults.csv"))[,c("empiricalPower",
                                                           "approxPower")]
  powerResults.combined = rbind(
    powerResults.clusterRandomized,
    powerResults.longitudinal,
    powerResults.clusterRandomized,
    powerResults.longitudinal
  )
  powerResults.combined$method = c(
    rep("All Designs", nrow(powerResults.clusterRandomized) +
          nrow(powerResults.longitudinal)),
    rep("Cluster Randomized", nrow(powerResults.clusterRandomized)),
    rep("Longitudinal", nrow(powerResults.longitudinal))
  )
  powerResults.combined$deviation = 
    powerResults.combined$approxPower - powerResults.combined$empiricalPower
  
  # Plot deviation from empirical across all designs
  pdf(file=file.path(study.figures.dir, "PowerBoxPlot_Overall.pdf"), family="Times")
  par(mfrow=c(1,1), lab=c(3,3,7))
  boxplot(deviation ~ method, data=powerResults.combined, range=0,
          las=1, ylim=c(-0.1,0.1),
          ylab="Deviation from Empirical Power")
  abline(h=0,lty=3)
  dev.off()
  
  tiff(file.path(study.figures.dir, "PowerBoxPlot_Overall.tiff"), units="in", 
       width=7, height=7, res=300, compression = 'lzw')
  par(mfrow=c(1,1), lab=c(3,3,7))
  boxplot(deviation ~ method, data=powerResults.combined, range=0,
          las=1, ylim=c(-0.1,0.1),
          ylab="Deviation from Empirical Power")
  abline(h=0,lty=3)
  dev.off()
  
  summary_df = data.frame(
    designType=c(rep("All Designs",5), rep("Cluster Randomized", 5),
                 rep("Longitudinal", 5)),
    stat=rep(c("Median", "Min", "1st quartile", "3rd quartile", "Max"), 3),
    value=c(
      round(quantile(powerResults.combined[powerResults.combined$method=="All Designs",]$deviation, probs=c(0.5,0,0.25,0.75,1)), 3),
      round(quantile(powerResults.combined[powerResults.combined$method=="Cluster Randomized",]$deviation, probs=c(0.5,0,0.25,0.75,1)), 3),
      round(quantile(powerResults.combined[powerResults.combined$method=="Longitudinal",]$deviation, probs=c(0.5,0,0.25,0.75,1)), 3)
    )
  )
  write.csv(summary_df, file.path(study.data.dir, "powerDeviationSummary.csv"), row.names=F)
}

#' runAppliedExample
#' 
#' This function runs the applied example in the manuscript seen in the
#' README. It produces a power value as output
#'
runAppliedExample <- function() {
  cat("### Running applied example")
  hennrikusExample()
}

