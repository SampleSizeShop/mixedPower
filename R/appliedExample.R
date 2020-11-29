# 
#  Package mixedPower calculates power for the linear mixed model
#  Copyright (C) 2013 Sarah Kreidler.
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
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#
# Applied example for the manuscript
#
# A proposed cluster-randomized trial in which
# worksites are randomized to 2 smoking cessation programs
#
hennrikusExample = function() {
  
  completeSize = 30
  clusterComplete = matrix(rep(1,completeSize))
  incompleteSize = 20
  clusterIncomplete = matrix(rep(1,incompleteSize))
  
  group1X = matrix(c(1,0), nrow=1)
  group2X = matrix(c(0,1), nrow=1)
  
  design = new("design.mixed", name = "Hennrikus Examle", 
               description = "Cluster randomized trial of smoking cessation intervention",
               xPatternList = c(
                 new("missingDataPattern", group=1, observations=1:completeSize, size=25,
                     designMatrix=clusterComplete %x% group1X),
                 new("missingDataPattern", group=1, observations=1:incompleteSize, size=15,
                     designMatrix=clusterIncomplete %x% group1X),
                 new("missingDataPattern", group=1, observations=1:completeSize, size=25,
                     designMatrix=clusterComplete %x% group2X),
                 new("missingDataPattern", group=1, observations=1:incompleteSize, size=15,
                     designMatrix=clusterIncomplete %x% group2X)
               ),
               beta = matrix(c(25,0)),
               Sigma = (125^2) * (0.04 * (matrix(rep(1,30)) %*% t(matrix(rep(1,30)))) + 
                                    diag(30) * (1 - 0.04))
  )
  # get the appropriate hypothesis
  glh = new("glh.mixed",
            alpha = 0.05,
            fixedContrast = matrix(c(1,-1), nrow=1),
            thetaNull = matrix(0),
            test = "Wald, KR ddf")
  
  return(mixedPower(design, glh))
  
}

