

#  This function was used and modified under the terms of the GNU General Public License as published by
#  the Free Software Foundation
#  This function relies on R/ggplot2 https://ggplot2.tidyverse.org/

forestplot=function(rrdata, orientation='horizontal'){
  library(ggplot2)

  if (orientation=="vertical") {
    p = ggplot(data=rrdata,
    aes(x = Group,y = OddsRatio, ymin = LowerLimit, ymax = UpperLimit ))+
    scale_x_discrete(limits=rrdata$Group) + 
    geom_pointrange(aes(col=Group))+
    geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
    #xlab('Group')+ 
    ylab("Odds Ratio (95% Confidence Interval)")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.5,cex=1)+
    facet_wrap(~Condition,strip.position="left",nrow=9,scales = "free_y") +
    guides(fill=FALSE) +
    theme(plot.title=element_text(size=16,face="bold"),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(face="bold"),
    axis.title=element_text(size=12,face="bold"),
    strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
    theme(plot.margin = unit(c(1,1,0.5,0.5), "cm"))+
    coord_flip()
  }
 if (orientation=="horizontal") {
    p = ggplot(data=rrdata,
    aes(x = Group,y = OddsRatio, ymin = LowerLimit, ymax = UpperLimit ))+
    scale_x_discrete(limits=rrdata$Group) + 
    geom_pointrange(aes(col=Group))+
    geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
    #xlab('Group')+ 
    ylab("Odds Ratio (95% Confidence Interval)")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.2,cex=1)+
    facet_wrap(~Condition,strip.position="top",nrow=1,scales = "free_x") +
    guides(fill=FALSE) +
    theme(plot.title=element_text(size=16,face="bold"),
    axis.text.x=element_text(face="bold"),
    axis.title=element_text(size=12,face="bold"))+
    theme(plot.margin = unit(c(0.5, 0.5, 1, 1), "cm"))+
    scale_y_log10(breaks=c(0.5,1,2))
 }   
 return(p)
}

 