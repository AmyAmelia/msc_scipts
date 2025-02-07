# Libraries to load
library(ggplot2)
library(ROCR)
library(readr)
library(ggpubr)
# Load in main file 
Lit_splice_all <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/SQL/Lit_splice_reduc_all.txt", 
                                                         "\t", escape_double = FALSE, na = "NA", 
                                                         trim_ws = TRUE)
# Varaint distribution in relation to VEP location charaterisation - count
# want to add but need to format  geom_text(stat='count', aes(label=..count..), vjust=-1) + # Varaint distribution in relation to VEP location charaterisation - %
# way 5 not working https://sebastiansauer.github.io/percentage_plot_ggplot2_V2/
region <- ggplot(Lit_splice_all, aes(Location, fill = Lit_splice_all$`Variants Effect on splicing`)) + geom_bar(position = position_dodge(width = 1))+ labs( title="Variant location in relation to transcripts (as predicted by VEP)") +
   xlab("")+  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + labs(fill= "Effect of varaint\n on splicing")
#+ scale_fill_discrete("Effect of variant \non splicing", breaks=c("Effects Splicing", ))


# Distribution of varaints -ALL
# dist to splice site ± 100bp note: excludes 172 variants outside of range
loc <- ggplot(Lit_splice_all, aes(distNearestSS, fill = Lit_splice_all$`Variants Effect on splicing`)) + geom_histogram(alpha = 0.7, bins =200, ) + xlim(range(-100:100)) + ylim(range(0:500)) + labs( title="Distribution of splice site variants in relation to the nearest splice site") + xlab("Distance to nearest splice site (bp)") + theme_bw()  + labs(fill= "Effect of varaint\n on splicing")
plot(loc)
# plot location figures on same page
ggarrange(region,loc, 
          labels = c("A", ""),ncol = 1, nrow = 2) 



################# Individual tools ######################
# Load truncated files for per tool location display then generate histograms per tool.
# SSF
Lit_splice_SSF <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/SQL/Lit_splice_reduc_SSF.txt", 
                                   "\t", escape_double = FALSE, na = "NA", 
                                   trim_ws = TRUE)
# dist to splice site ± 100bp
ssf_loc <-ggplot(Lit_splice_SSF, aes(distNearestSS, fill = Effect)) + geom_histogram(alpha = 0.7, bins =
                                                                             200,show.legend = FALSE ) + xlim(range(-100:100)) + ylim(range(0:500)) + labs( title="Splice Site Finder") + xlab("Nearest splice site (bp)") + theme_bw()

# MES
Lit_splice_MES <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/SQL/Lit_splice_reduc_MES.txt", 
                                   "\t", escape_double = FALSE, na = "NA", 
                                   trim_ws = TRUE)
# dist to splice site ± 100bp: 168 rows excluded
mes_loc <-ggplot(Lit_splice_MES, aes(distNearestSS, fill = Effect)) + geom_histogram(alpha = 0.7, bins =
  200, show.legend = FALSE)  + xlab("Nearest splice site (bp)") + xlim(range(-100:100)) + ylim(range(0:500)) + labs( 
    title="MaxEntScan") + theme_bw()

# NNS
Lit_splice_NNS <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/SQL/Lit_splice_reduc_NNS.txt", 
                                   "\t", escape_double = FALSE, na = "NA", 
                                   trim_ws = TRUE)
# dist to splice site ± 100bp
nns_loc<-ggplot(Lit_splice_NNS, aes(distNearestSS, fill = Effect)) + geom_histogram(alpha = 0.7, bins =
                                                                             200,show.legend = FALSE) + xlim(range(-100:100)) + ylim(range(0:500)) + labs( title="NNsplice") + xlab("Nearest splice site (bp)") + theme_bw()

#GS
Lit_splice_GS <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/SQL/Lit_splice_reduc_GS.txt", 
                                   "\t", escape_double = FALSE, na = "NA", 
                                   trim_ws = TRUE)
# dist to splice site ± 100bp

gs_loc <-ggplot(Lit_splice_GS, aes(distNearestSS, fill = Effect)) + geom_histogram(alpha = 0.7, bins =
                                                                             200, show.legend = FALSE) + xlim(range(-100:100)) + ylim(range(0:500)) + labs( title="Gene Splicer") + xlab("Nearest splice site (bp)") + theme_bw()

#ada
Lit_splice_ada <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/SQL/Lit_splice_reduc_ada.txt", 
                                   "\t", escape_double = FALSE, na = "NA", 
                                   trim_ws = TRUE)
# dist to splice site ± 100bp
ada_loc <- ggplot(Lit_splice_ada, aes(distNearestSS, fill = Effect)) + geom_histogram(alpha = 0.7, bins =
                                                                             200,show.legend = FALSE ) + xlim(range(-100:100)) + ylim(range(0:500)) + labs( title="Ada boost") + xlab("Nearest splice site (bp)") + theme_bw()

#rf
Lit_splice_rf <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/SQL/Lit_splice_reduc_rf.txt", 
                                   "\t", escape_double = FALSE, na = "NA", 
                                   trim_ws = TRUE)
# dist to splice site ± 100bp
rf_loc <- ggplot(Lit_splice_rf, aes(distNearestSS, fill = Effect)) + geom_histogram(alpha = 0.7, bins =
                                                                             200, show.legend = FALSE) + xlim(range(-100:100)) + ylim(range(0:500)) + labs( title="Random forest") + xlab("Nearest splice site (bp)") + theme_bw()
#spidex
Lit_splice_spidex <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/SQL/Lit_splice_reduc_spidex.txt", 
                                   "\t", escape_double = FALSE, na = "NA", 
                                   trim_ws = TRUE)
# dist to splice site ± 100bp
spi_loc <- ggplot(Lit_splice_spidex, aes(distNearestSS, fill = Effect)) + geom_histogram(alpha = 0.7, bins =
                                                                             200,show.legend = FALSE ) + xlim(range(-100:100)) + ylim(range(0:500)) + labs( title="Spidex") + xlab("Nearest splice site (bp)") + theme_bw()

# Plot all location graphs in one page  
ggarrange(ssf_loc, mes_loc, nns_loc, gs_loc , ada_loc, rf_loc, spi_loc, 
          labels = c("A", "B", "C","D","E", "F","G"),
          ncol = 3, nrow = 4) 
