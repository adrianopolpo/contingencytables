# code for R (http://www.R-project.org)
# Copyright (C) 2018 Natalia Lombardi de Oliveira,
#                    Marcio Alves Diniz,
#                    Carlos Alberto de Bragança Pereira,
#                    Adriano Polpo (polpo@ufscar.br).
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    For the GNU General Public License see <http://www.gnu.org/licenses/>.
#
#    Please see the paper:
#    DOI: 10.1371/journal.pone.0199102 (To cite this code, please cite the paper).

source("progs.r")

if (FALSE) {
  ## This may take 2 days running...
  set.seed(1234)
  homog.2x2.10.power  <- power.homog.2x2(n=c(10,10))
  homog.2x2.10        <- homog.2x2.10.power$p
  homog.2x2.30        <- homog.2x2(30,30)
  homog.2x2.100.power <- power.homog.2x2(n=c(100,100))
  homog.2x2.100       <- homog.2x2.100.power$p
  homog.2x3.30        <- homog.2x3(30,30)
  homog.3x3.15        <- homog.3x3(15,15,15)

  indep.2x2.30 <- indep.2x2(30)
  indep.2x3.30 <- indep.2x3(30)
  indep.3x3.15 <- indep.3x3(15)
  indep.3x3.25 <- indep.3x3(25)

  hw10.power  <- power.hw(n=10)
  hw10        <- hw10.power$p
  hw30        <- HW(30)
  hw100.power <- power.hw(n=100)
  hw100       <- hw100.power$p

  save.image()
}

if (TRUE) {
set.height <- 10
set.width  <- 10
set.cex    <- 0.6
set.cex.labels <- 0.6

# Homogeneity 2x2
homog.labels <- rbind(
  c(         " ","asymptotic",  "exact",  "exact","asymptotic",  "exact",  "exact"),
  c("Chi-Square",       "LRT",    "LRT",   "FBST",      "FBST", "Fisher","Barnard"),
  c(   "p-value",   "p-value","P-value","e-value",   "e-value","p-value","p-value"))

name <- "fig_homo2x2_30"
pairs2(homog.2x2.30[,c(12,10,9,15,13,17,16)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

name <- "fig_homo2x2_100"
pairs2(homog.2x2.100[,c(12,10,9,15,13,17,16)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

####
homog.labels <- rbind(
  c(         " ","asymptotic",  "exact",  "exact","asymptotic"),
  c("Chi-Square",       "LRT",    "LRT",   "FBST",      "FBST"),
  c(   "p-value",   "p-value","P-value","e-value",   "e-value"))

# Homogeneity 2x3
name <- "fig_homo2x3_30"
pairs2(homog.2x3.30[,c(14,12,11,17,15)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

# Homogeneity 3x3
name <- "fig_homo3x3_15"
pairs2(homog.3x3.15[,c(17,15,14,20,18)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

# Power - Homogeneity 2x2
homog.labels <- rbind(
  c("asymptotic","exact",         " ",  "exact", "exact"),
  c(       "LRT",  "LRT","Chi-Square","Barnard","Fisher"),
  c(         " ",    " ",         " ",      " ",     " "))

name <- "fig_homo2x2_10-power"
pairs2(homog.2x2.10.power$power[,c(4,3,6,5,7)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

name <- "fig_homo2x2_100-power"
pairs2(homog.2x2.100.power$power[,c(4,3,6,5,7)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

plot.power3d(homog.2x2.10.power$power,x="LRT",y="P",a1=140,name="homog10-power_LRT-P")
plot.power3d(homog.2x2.10.power$power,x="P",y="Chi",a1=140,name="homog10-power_P-Chi")
plot.power3d(homog.2x2.10.power$power,x="Chi",y="Bar",a1=140,name="homog10-power_Chi-Bar")
plot.power3d(homog.2x2.10.power$power,x="Bar",y="Fish",a1=140,name="homog10-power_Bar-Fish")


###
indep.labels <- rbind(
  c(         " ","asymptotic",  "exact",  "exact","asymptotic"),
  c("Chi-Square",       "LRT",    "LRT",   "FBST",      "FBST"),
  c(   "p-value",   "p-value","P-value","e-value",   "e-value"))

# Independence 2x2
name <- "fig_indep2x2_30"
pairs2(indep.2x2.30[,c(11,9,8,14,12)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

# Independence 2x3
name <- "fig_indep2x3_30"
pairs2(indep.2x3.30[,c(13,11,10,16,14)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

# Independence 3x3
name <- "fig_indep3x3_15"
pairs2(indep.3x3.15[,c(16,14,13,19,17)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

name <- "fig_indep3x3_25"
pairs2(indep.3x3.25[,c(16,14,13,19,17)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()


# Hardy-Weiberg
hw.labels <- rbind(
  c(         " ","asymptotic",  "exact",  "exact","asymptotic",  "exact"),
  c("Chi-Square",       "LRT",    "LRT",   "FBST",      "FBST","Barnard"),
  c(   "p-value",   "p-value","P-value","e-value",   "e-value","p-value"))


name <- "fig_hw_30"
pairs2(hw30[,c(10,8,7,13,11,14)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

name <- "fig_hw_100"
pairs2(hw100[,c(10,8,7,13,11,14)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

# Power - Hardy-Weiberg
hw.labels <- rbind(
  c("asymptotic","exact",         " ",  "exact"),
  c(       "LRT",  "LRT","Chi-Square","Barnard"),
  c(         " ",    " ",         " ",      " "))


name <- "fig_hw_10-power"
pairs2(hw10.power$power[,c(4,3,6,5)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

name <- "fig_hw_100-power"
pairs2(hw100.power$power[,c(4,3,6,5)],upper.panel=NULL,
       labels=homog.labels,cex.labels=set.cex.labels,pch=18,cex=set.cex,family="Arial")
dev2bitmap(paste(name,".tiff",sep=""),height=set.height,width=set.width,units ='cm',
           type="tiff24nc",res=300,method="pdf")
dev.off()

plot.power3d(hw10.power$power,x="LRT",y="P",name="hw10-power_LRT-P")
plot.power3d(hw10.power$power,x="P",y="Chi",name="hw10-power_P-Chi")
plot.power3d(hw10.power$power,x="Chi",y="Bar",name="hw10-power_Chi-Bar")

}


