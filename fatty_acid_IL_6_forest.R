rm(list = ls())
setwd ("C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/Year Three/fatty acids and inflammation") # Set the file path to your local directory


library(ggplot2)
library(ggpubr)

#read in results
MR_results <- read.csv("C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/Year Three/fatty acids and inflammation/fatty_acid_IL_6.csv")
colnames(MR_results)[1]<- "Exposure"
MR_results[c(1,5,9,13), "Method"]<- "Weighted Mode"
MR_results[c(2,6,10,14), "Method"]<- "Weighted Median"
MR_results[c(3,7,11,15), "Method"]<- "MR-Egger"
MR_results[c(4,8,12,16), "Method"]<- "Inverse Variance Weighted"
MR_results[c(1:16), "Outcome"] <- " IL_6"
#Generate relationship variable
MR_results$Relationships <- paste(MR_results$Exposure, "->",MR_results$Outcome)
#set the factor levels of the method variable - will make sure they are plotted in the right order
MR_results$Method= factor(MR_results$Method, levels = unique(MR_results$Method))


#Where DHA is the exposure
DHA_exp=MR_results[MR_results$Exposure=="DHA",]


#set your preferred colour or colours
colours=c(rep("#4f549c", times=9))

#build the plot
fp <- ggplot(data=DHA_exp, aes(y=Method, x=beta, xmin=CIL, xmax=CIU,colour=Method)) +
  geom_pointrange() + 
  geom_errorbarh(height=.1) +
  geom_vline(xintercept=0, lty=2) + 
  xlab("Beta (95% CI)") +
  theme_bw() + 
  scale_color_manual(values = c("Inverse Variance Weighted" = "#ffd966" ,
                                "MR-Egger" = "#cc0000",
                                "Weighted Median" = "#6fa8dc",
                                "Weighted Mode"= "#8395a3" ))

print(fp)

fp2 <- fp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(legend.position = "none")


print (fp2)
##You might need to change the hjust value here to position the name of each facet correctly
fp3 <- fp2 + ggtitle("Relationship") + theme(plot.title = element_text(size=11, hjust =2.6)) 

plot1=fp3+facet_grid(Relationships~., scales="free", labeller = label_wrap_gen(20))+
  theme(strip.text.y = element_text(size = 11,  angle = 0))+
  theme(strip.background =element_rect(fill="white"))

print(fp3)
print(fp)
print(plot1)



#Where LA is the outcome
LA_exp=MR_results[MR_results$Exposure=="LA",]
#set your preferred colour or colours
colours=c(rep("#2f9c2b", times=9))
#build the plot
fp <- ggplot(data=LA_exp, aes(y=Method, x=(beta), xmin=(CIL), xmax=(CIU),colour=Method)) +
  geom_pointrange() + 
  geom_errorbarh(height=.1) +
  geom_vline(xintercept=0, lty=2) + 
  xlab("beta (95% CI)") +
  theme_bw() + 
  scale_color_manual(values = c("Inverse Variance Weighted" = "#ffd966" ,
                                "MR-Egger" = "#cc0000",
                                "Weighted Median" = "#6fa8dc",
                                "Weighted Mode"= "#8395a3" )) 

print(fp)
fp2 <- fp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(legend.position = "none")

print(fp2)

##You might need to change the hjust value here to position the name of each facet correctly
fp3 <- fp2 + ggtitle("Relationship") + theme(plot.title = element_text(size=11, hjust =2.8)) 

print(fp3)

plot2=fp3+facet_grid(Relationships~., scales="free", labeller = label_wrap_gen(20))+
  theme(strip.text.y = element_text(size = 11,  angle = 0))+
  theme(strip.background =element_rect(fill="white"))
final_plot<-ggarrange(plot1, plot2, ncol = 2, nrow = 1)
print(final_plot)
print(plot2)



#Where Omega3 is the outcome
omega_3_exp=MR_results[MR_results$Exposure=="Omega-3",]
#set your preferred colour or colours
colours=c(rep("#2f9c2b", times=9))
#build the plot
fp <- ggplot(data=omega_3_exp, aes(y=Method, x=(beta), xmin=(CIL), xmax=(CIU),colour=Method)) +
  geom_pointrange() + 
  geom_errorbarh(height=.1) +
  geom_vline(xintercept=0, lty=2) + 
  xlab("beta (95% CI)") +
  theme_bw() + 
  scale_color_manual(values = c("Inverse Variance Weighted" = "#ffd966" ,
                                "MR-Egger" = "#cc0000",
                                "Weighted Median" = "#6fa8dc",
                                "Weighted Mode"= "#8395a3" )) 

print(fp)
fp2 <- fp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(legend.position = "none")

print(fp2)

##You might need to change the hjust value here to position the name of each facet correctly
fp3 <- fp2 + ggtitle("Relationship") + theme(plot.title = element_text(size=11, hjust =2.8)) 

print(fp3)

plot3=fp3+facet_grid(Relationships~., scales="free", labeller = label_wrap_gen(20))+
  theme(strip.text.y = element_text(size = 11,  angle = 0))+
  theme(strip.background =element_rect(fill="white"))


print(plot3)


#Where LA is the outcome
LA_exp=MR_results[MR_results$Exposure=="LA",]
#set your preferred colour or colours
colours=c(rep("#2f9c2b", times=9))
#build the plot
fp <- ggplot(data=LA_exp, aes(y=Method, x=(beta), xmin=(CIL), xmax=(CIU),colour=Method)) +
  geom_pointrange() + 
  geom_errorbarh(height=.1) +
  geom_vline(xintercept=0, lty=2) + 
  xlab("beta (95% CI)") +
  theme_bw() + 
  scale_color_manual(values = c("Inverse Variance Weighted" = "#ffd966" ,
                                "MR-Egger" = "#cc0000",
                                "Weighted Median" = "#6fa8dc",
                                "Weighted Mode"= "#8395a3" ))

print(fp)
fp2 <- fp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(legend.position = "none")

print(fp2)

##You might need to change the hjust value here to position the name of each facet correctly
fp3 <- fp2 + ggtitle("Relationship") + theme(plot.title = element_text(size=11, hjust =2.8)) 

print(fp3)

plot2=fp3+facet_grid(Relationships~., scales="free", labeller = label_wrap_gen(20))+
  theme(strip.text.y = element_text(size = 11,  angle = 0))+
  theme(strip.background =element_rect(fill="white"))
final_plot<-ggarrange(plot1, plot2, ncol = 2, nrow = 1)
print(final_plot)
print(plot2)



#Where Omega6 is the outcome
omega_6_exp=MR_results[MR_results$Exposure=="Omega-6",]
#set your preferred colour or colours
colours=c(rep("#2f9c2b", times=9))
#build the plot
fp <- ggplot(data=omega_6_exp, aes(y=Method, x=(beta), xmin=(CIL), xmax=(CIU),colour=Method)) +
  geom_pointrange() + 
  geom_errorbarh(height=.1) +
  geom_vline(xintercept=0, lty=2) + 
  xlab("beta (95% CI)") +
  theme_bw() + 
  scale_color_manual(values = c("Inverse Variance Weighted" = "#ffd966" ,
                                "MR-Egger" = "#cc0000",
                                "Weighted Median" = "#6fa8dc",
                                "Weighted Mode"= "#8395a3" ))

print(fp)
fp2 <- fp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(legend.position = "none")

print(fp2)

##You might need to change the hjust value here to position the name of each facet correctly
fp3 <- fp2 + ggtitle("Relationship") + theme(plot.title = element_text(size=11, hjust =2.8)) 

print(fp3)

plot4=fp3+facet_grid(Relationships~., scales="free", labeller = label_wrap_gen(20))+
  theme(strip.text.y = element_text(size = 11,  angle = 0))+
  theme(strip.background =element_rect(fill="white"))


print(plot4)


png('C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/Year Three/fatty acids and inflammation/fatty_acids_IL_6.png', res=400, height=2000, width=4000)
plot.new()
ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)
dev.off()




final_plot<-ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)

print(final_plot)
