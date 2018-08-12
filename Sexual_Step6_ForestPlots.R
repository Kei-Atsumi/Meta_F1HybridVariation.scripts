############################################
# Forest plots Type models
############################################

rm(list=ls())

library(tidyverse)
library(gridExtra)

#################################################################################
# Reciplocal
#################################################################################
# Plot Type model

setwd("/media/sf_Dropbox/Meta_F1HybridVariation/Analysis")

recip.lnRR <- read.csv("recip.lnRR.csv")
plot.recip <- ggplot(recip.lnRR, aes(y = estimate, x = factors))+
  geom_rect(aes(ymin= 0.5, ymax= 1.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  geom_rect(aes(ymin=-0.5, ymax=-1.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                colour="black", width=.0001) +
  geom_point(size=8, colour="black", aes(shape=type))+
  scale_shape_manual(values=c(17, 16))+
  geom_hline(yintercept = 0, linetype="dashed", size = 0.3)+ 
  theme (panel.grid.major = element_blank())+
  theme_classic((base_size = 20))+ 
  theme(plot.title=element_text(size=17, face="bold", vjust=2, lineheight=1))+ 
  xlab(" ") + ylab("Effect size (lnRR)")+     # Title of x-axis and y-axis
  scale_x_discrete(limits = rev(recip.lnRR$factors))+ # reverse x-axis
  theme(axis.ticks.y=element_blank()) +
  theme(axis.line.y = element_blank())+
  theme(axis.title.x=element_text(face="plain", size=35, vjust=-0.35), axis.text.x= element_text(face="plain", size=25))+
  theme(axis.title.y=element_text(face="plain", size=35, vjust=1), axis.text.y= element_text(face="plain", size=25))+
  theme(legend.position="none") +
  coord_flip()
plot.recip
ggsave(file = "plot.recip.lnRR.eps", plot = plot.recip, height = 13, width = 9)
ggsave(file = "plot.recip.lnRR.png", plot = plot.recip, height = 13, width = 9)



recip.SMD <- read.csv("recip.SMD.csv")
plot.recip <- ggplot(recip.SMD, aes(y = estimate, x = factors))+
  geom_rect(aes(ymin= 0.5, ymax= 1.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  geom_rect(aes(ymin= 1.5, ymax= 1.0,xmin=-Inf,xmax=Inf),fill="Grey 75", inherit.aes = FALSE)+
  geom_rect(aes(ymin=-0.5, ymax=-1.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                colour="black", width=.0001) +
  geom_point(size=8, colour="black", aes(shape=type))+
  scale_shape_manual(values=c(17, 16))+
  geom_hline(yintercept = 0, linetype="dashed", size = 0.3)+ 
  theme (panel.grid.major = element_blank())+
  theme_classic((base_size = 20))+ 
  theme(plot.title=element_text(size=17, face="bold", vjust=2, lineheight=1))+ 
  xlab(" ") + ylab("Effect size (SMD)")+     # Title of x-axis and y-axis
  scale_x_discrete(limits = rev(recip.SMD$factors))+ # reverse x-axis
  theme(axis.ticks.y=element_blank()) +
  theme(axis.line.y = element_blank())+
  theme(axis.title.x=element_text(face="plain", size=35, vjust=-0.35), axis.text.x= element_text(face="plain", size=25))+
  theme(axis.title.y=element_text(face="plain", size=35, vjust=1), axis.text.y= element_text(face="plain", size=25))+
  theme(legend.position="none") +
  coord_flip()
plot.recip
ggsave(file = "plot.recip.SMD.eps", plot = plot.recip, height = 13, width = 9)
ggsave(file = "plot.recip.SMD.png", plot = plot.recip, height = 13, width = 9)



#################################################################################
# Epistasis
#################################################################################
# Plot Type model

setwd("/media/sf_Dropbox/Meta_F1HybridVariation/Analysis")

epi <- read.csv("epistasis.csv")
plot.epi <- ggplot(epi, aes(y = estimate, x = factors))+
  geom_rect(aes(ymin= 0.5, ymax= 1.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  geom_rect(aes(ymin=-0.5, ymax=-1.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  geom_rect(aes(ymin= 1.5, ymax= 1.0,xmin=-Inf,xmax=Inf),fill="Grey 75", inherit.aes = FALSE)+
  geom_rect(aes(ymin=-1.5, ymax=-1.0,xmin=-Inf,xmax=Inf),fill="Grey 75", inherit.aes = FALSE)+
  geom_rect(aes(ymin= 1.5, ymax= 2.0,xmin=-Inf,xmax=Inf),fill="Grey 65", inherit.aes = FALSE)+
  geom_rect(aes(ymin=-1.5, ymax=-2.0,xmin=-Inf,xmax=Inf),fill="Grey 65", inherit.aes = FALSE)+
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                colour="black", width=.0001) +
  geom_point(size=8, colour="black", aes(shape=type)) + 
  scale_shape_manual(values=c(17, 16))+
  geom_hline(yintercept = 0, linetype="dashed", size = 0.3)+ 
  theme (panel.grid.major = element_blank())+
  theme_classic((base_size = 20))+ 
  theme(plot.title=element_text(size=17, face="bold", vjust=2, lineheight=1))+ 
  xlab(" ") + ylab("Effect size (lnCVR)")+     # Title of x-axis and y-axis
  scale_x_discrete(limits = rev(epi$factors))+ # reverse x-axis
  theme(axis.ticks.y=element_blank()) +
  theme(axis.line.y = element_blank())+
  theme(axis.title.x=element_text(face="plain", size=35, vjust=-0.35), axis.text.x= element_text(face="plain", size=25))+
  theme(axis.title.y=element_text(face="plain", size=35, vjust=1), axis.text.y= element_text(face="plain", size=25))+
  theme(legend.position="none") +
  coord_flip()
plot.epi
ggsave(file = "plot.epistasis.eps", plot = plot.epi, height = 13, width = 10)
ggsave(file = "plot.epistasis.png", plot = plot.epi, height = 13, width = 10)




############################################################################
# recip old
############################################################################

meta.recip <- read.csv("random_recip.csv")
plot.meta.recip <- ggplot(meta.recip, aes(y = estimate, x = "intrcpt"))+
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                colour="black", width=.0001) +
  geom_point(size=5, colour="black", shape=19) + 
  geom_hline(yintercept = 0, linetype="dashed", size = 0.3)+ 
  theme (panel.grid.major = element_blank())+
  theme_classic((base_size = 20))+ 
  theme(plot.title=element_text(size=17, face="bold", vjust=2, lineheight=1))+ 
  xlab(" ") + ylab("Effect size (lnRR)")+     # Title of x-axis and y-axis
  scale_y_continuous(limits = c(-0.1, 0.2))+
  theme(axis.ticks.y=element_blank()) +
  theme(axis.line.y = element_blank())+
  theme(axis.title.x=element_text(face="plain", size=35, vjust=-0.35), axis.text.x= element_text(face="plain", size=25))+
  theme(axis.title.y=element_text(face="plain", size=35, vjust=1), axis.text.y= element_text(face="plain", size=25))+
  coord_flip()
plot.meta.recip
ggsave(file = "plot.meta.recip.eps", plot = plot.meta.recip, height = 3, width = 9)

reg.recip <- read.csv("reg.recip.csv")
plot.reg.recip <- ggplot(reg.recip, aes(y = estimate, x = factors))+
  #  geom_rect(aes(ymin= 0.5, ymax= 1.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  #  geom_rect(aes(ymin=-0.5, ymax=-1.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  #  geom_rect(aes(ymin= 1.5, ymax= 1.0,xmin=-Inf,xmax=Inf),fill="Grey 75", inherit.aes = FALSE)+
  #  geom_rect(aes(ymin=-1.5, ymax=-1.0,xmin=-Inf,xmax=Inf),fill="Grey 75", inherit.aes = FALSE)+
  #  geom_rect(aes(ymin= 1.5, ymax= 2.0,xmin=-Inf,xmax=Inf),fill="Grey 65", inherit.aes = FALSE)+
  #  geom_rect(aes(ymin=-1.5, ymax=-2.0,xmin=-Inf,xmax=Inf),fill="Grey 65", inherit.aes = FALSE)+
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                colour="black", width=.0001) +
  geom_point(size=5, colour="black", shape=19) + 
  geom_hline(yintercept = 0, linetype="dashed", size = 0.3)+ 
  theme (panel.grid.major = element_blank())+
  theme_classic((base_size = 20))+ 
  theme(plot.title=element_text(size=17, face="bold", vjust=2, lineheight=1))+ 
  xlab(" ") + ylab("Effect size (lnRR)")+     # Title of x-axis and y-axis
  scale_x_discrete(limits = dat.recip$factors)+
  theme(axis.ticks.y=element_blank()) +
  theme(axis.line.y = element_blank())+
  theme(axis.title.x=element_text(face="plain", size=35, vjust=-0.35), axis.text.x= element_text(face="plain", size=25))+
  theme(axis.title.y=element_text(face="plain", size=35, vjust=1), axis.text.y= element_text(face="plain", size=25))+
  coord_flip()
plot.recip
ggsave(file = "plot.reg.recip.eps", plot = plot.reg.recip, height = 5, width = 9)

############################################################################
# epi old
############################################################################


reg.epi <- read.csv("reg.epistasis.csv")

plot.epi <- ggplot(reg.epi, aes(y = estimate, x = factors))+
  #  geom_rect(aes(ymin= 1.0, ymax= 5.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  #  geom_rect(aes(ymin=-1.0, ymax=-5.0,xmin=-Inf,xmax=Inf),fill="Grey 90", inherit.aes = FALSE)+
  #  geom_rect(aes(ymin=-5.0, ymax=-10.0,xmin=-Inf,xmax=Inf),fill="Grey 75", inherit.aes = FALSE)+
  #  geom_rect(aes(ymin=-10.0, ymax=-15.0,xmin=-Inf,xmax=Inf),fill="Grey 65", inherit.aes = FALSE)+
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                colour="black", width=.0001) +
  geom_point(size=5, colour="black", shape=19) + 
  geom_hline(yintercept = 0, linetype="dashed", size = 0.3)+ 
  theme (panel.grid.major = element_blank())+
  theme_classic((base_size = 20))+ 
  theme(plot.title=element_text(size=17, face="bold", vjust=2, lineheight=1))+ 
  xlab(" ") + ylab("Effect size (lnCVR)")+     # Title of x-axis and y-axis
  scale_x_discrete(limits = reg.epi$factors, trans = "reverse")+
  theme(axis.ticks.y=element_blank()) +
  theme(axis.line.y = element_blank())+
  theme(axis.title.x=element_text(face="plain", size=35, vjust=-0.30), axis.text.x= element_text(face="plain", size=25))+
  theme(axis.title.y=element_text(face="plain", size=30, vjust=1), axis.text.y= element_text(face="plain", size=25))+coord_flip()
plot.epi
ggsave(file = "plot.reg.epistasis.eps", plot = plot.epi, height = 6.5, width = 9)
