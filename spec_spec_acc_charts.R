# Libraries to load
library(ggplot2)
library(ROCR)
library(readr)
library(ggpubr)
# Load in main file 
outter <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/ACGSMeeting2018/ACGS_outter_sr_sensandSpec.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)
acc <- ggplot(data=outter, aes(x= Tool, y=Accuracy, fill=Threshold)) + geom_bar(stat = "identity", position = 'dodge') + ylab("Accuracy (%)") +labs(title = "Accuracy") + scale_x_discrete(limits=c("Alamut3", "ADA", "RF","SSF","MES","NNS", "GS" ))
sen <- ggplot(data=outter, aes(x= Tool, y=Sensitivity, fill=Threshold)) + geom_bar(stat = "identity", position = 'dodge') + ylab("Sensitivity (%)")  +labs(title = "Sensitivity") + scale_x_discrete(limits=c("Alamut3", "ADA", "RF","SSF","MES","NNS", "GS" ))
spec <-ggplot(data=outter, aes(x= Tool, y=Specificity, fill=Threshold)) + geom_bar(stat = "identity", position = 'dodge') + ylab("Specificity (%)") +labs(title = "Specificity")+ scale_x_discrete(limits=c("Alamut3", "ADA", "RF","SSF","MES","NNS", "GS" ))

plot(acc)
plot(sen)
plot(spec)

# Plot all 5' and 3' locations on in one page  
ggarrange(sen,spec, acc,
          ncol = 1, nrow = 3)
