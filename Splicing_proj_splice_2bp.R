# Libraries to load
library(ggplot2)
library(ROCR)
library(readr)
library(ggpubr)
theme_set(theme_classic())

# Load in file for splice region vartaints. 
Lit_splice_2bp <- read_delim("~/Google Drive/Splicing_Proj/02_Data_out/ACGSMeeting2018/Truthset_splicesite_var.txt", 
                                         "\t", escape_double = FALSE, na = "NA", 
                                         trim_ws = TRUE)
# Distribution of varaints -Splice site
# dist to splice site ± 20bp note: excludes 172 variants outside of range
ggplot(Lit_splice_2bp, aes(distNearestSS, fill = Effect)) + geom_histogram(alpha = 0.4, bins =25, palette = "Spectral") + xlim(range(-3:3)) + ylim(range(0:500)) + labs( title="Distribution of splice site variants in relation to the nearest splice site") + xlab("Distance to nearest splice site (bp)") + theme_bw()

# SSF score distribution
ssf_dist <- ggplot(Lit_splice_2bp, aes(SSF_PercentChange, fill = Effect)) + geom_density(alpha = 0.4,show.legend = FALSE)+ xlim(range(-25:125))+ theme_bw()  + labs( title="Splice Site Finder")
# MES score distribution
MES_dist <- ggplot(Lit_splice_2bp, aes(MES_PercentChange, fill = Effect)) + geom_density(alpha = 0.4,show.legend = FALSE)+ xlim(range(-50:125))+ theme_bw() + labs( title="MaxEntScan")
# NNS score distribution
NNS_dist <- ggplot(Lit_splice_2bp, aes(NNS_PercentChange, fill = Effect)) + geom_density(alpha = 0.4,show.legend = FALSE)+ xlim(range(-50:125))+ theme_bw() + labs( title="NNSplice")
# GS score distribution
GS_dist <- ggplot(Lit_splice_2bp, aes(GS_PercentChange, fill = Effect)) + geom_density(alpha = 0.4,show.legend = FALSE)+ xlim(range(-50:125))+ theme_bw() + labs( title="GeneSplicer")
# ada score distribution
ada_dist <- ggplot(Lit_splice_2bp, aes(ada_score, fill = Effect)) + geom_density(alpha = 0.4,show.legend = FALSE,adjust=1000)+ xlim(range(0:1))+ theme_bw() + labs( title="Ada Boost")
plot(ada_dist)
# rf score distribution
rf_dist <- ggplot(Lit_splice_2bp, aes(rf_score, fill = Effect)) + geom_density(alpha = 0.4,show.legend = FALSE)+ xlim(range(0:1))+ theme_bw() + labs( title="Random Forest")
# spidex score distribution
spidex_dist <- ggplot(Lit_splice_2bp, aes(inv_spidex_zscore, fill = Effect)) + geom_density(alpha = 0.4,show.legend = FALSE)+ xlim(range(-6:6))+ theme_bw() + labs( title="Spidex")
 plot(GS_dist)

# Plot all density graphs in one page  
ggarrange(ssf_dist, MES_dist, NNS_dist, GS_dist, ada_dist, rf_dist, spidex_dist, 
          labels = c("A", "B", "C","D","E", "F","G"),
          ncol = 3, nrow = 4)


######################################## ROC curves ###############################################
############### spidex ############### <- grey46
pred_spi <- prediction(Lit_splice_2bp$inv_spidex_zscore, Lit_splice_2bp$int_splicegroup)
# true pos v false neg
perf_spi <- performance(pred_spi, measure = "tpr", x.measure = "fpr")
# plot SSF roc
plot(perf_spi, col='grey46',  main = "Performance of individual tools within \n the canonical splice site",xlab = 
       "1- Specificity", ylab = "Sensitivity")

############### GS ############### <- magenta
pred_GS <- prediction(Lit_splice_2bp$GS_PercentChange, Lit_splice_2bp$int_splicegroup)
# true pos v false neg
perf_GS <- performance(pred_GS, measure = "tpr", x.measure = "fpr")
# plot SSF roc
#plot(perf_GS, col='magenta',add = TRUE)
plot(perf_GS, col='grey46',  main = "Performance of individual tools within \n the canonical splice site",xlab = 
       "1- Specificity", ylab = "Sensitivity")

############### MES ############### <- red
pred_MES <- prediction(Lit_splice_2bp$MES_PercentChange, Lit_splice_2bp$int_splicegroup)
# true pos v false neg
perf_MES <- performance(pred_MES, measure = "tpr", x.measure = "fpr")
# MES precission recal
recMES <- performance(pred_MES, measure = "prec", x.measure = "rec")
# MES accuracy
acc_MES <- performance(pred_MES,measure = "acc")
# plot MES roc
plot(perf_MES, col='red',add = TRUE)
# plot MES accuracy plot
# plot(acc_MES)

############### NNS ############### <- green4
pred_NNS <- prediction(Lit_splice_2bp$NNS_PercentChange, Lit_splice_2bp$int_splicegroup)
# true pos v false neg
perf_NNS <- performance(pred_NNS, measure = "tpr", x.measure = "fpr")
# plot SSF roc
plot(perf_NNS, col='green3',add = TRUE)

############### SSF ############### <- blue
pred_SSF <- prediction(Lit_splice_2bp$SSF_PercentChange, Lit_splice_2bp$int_splicegroup)
# true pos v false neg
perf_SSF <- performance(pred_SSF, measure = "tpr", x.measure = "fpr")
# plot SSF roc
plot(perf_SSF, col='blue', add=TRUE)

############### ADA ############### <- cyan
pred_ada <-prediction(Lit_splice_2bp$ada_score, Lit_splice_2bp$int_splicegroup)
perf_ada <- performance(pred_ada, measure = "tpr", x.measure = "fpr")
plot(perf_ada, col= "cyan",add = TRUE) # add = TRUE
# ada accuracy
acc_ada = performance(pred_ada, measure = "acc")

############### rf ############### <- gold
pred_rf <- prediction(Lit_splice_2bp$rf_score, Lit_splice_2bp$int_splicegroup)
# true pos v false neg
perf_rf <- performance(pred_rf, measure = "tpr", x.measure = "fpr")
# plot SSF roc
plot(perf_rf, col='gold',add = TRUE)

# add benchmark line
abline(a=0,b=1)

# Calculate area under curve (AUC)
#MES 
mesauc <- performance(pred_MES, "auc")
mesauc <- unlist(slot(mesauc,"y.values"))
mesauc <- round(mesauc,4)
txt_mesauc <- paste(c("MES: "), mesauc,sep="") # make veriable to combine text and AUC score.
#GS
gsauc <- performance(pred_GS, "auc")
gsauc <- unlist(slot(gsauc,"y.values"))
gsauc <- round(gsauc,4)
txt_gsauc <- paste(c("GS: "), gsauc,sep="") # make veriable to combine text and AUC score.
#SSF
ssfauc <- performance(pred_SSF, "auc")
ssfauc <- unlist(slot(ssfauc,"y.values"))
ssfauc <- round(ssfauc,4)
txt_ssfauc <- paste(c("SSF: "), ssfauc,sep="") # make veriable to combine text and AUC score.
#Ada
adaauc <- performance(pred_ada, "auc")
adaauc <- unlist(slot(adaauc,"y.values"))
adaauc <- round(adaauc,4)
txt_adaauc <- paste(c("Ada: "), adaauc,sep="") # make veriable to combine text and AUC score.
#RF
rfauc <- performance(pred_rf, "auc")
rfauc <- unlist(slot(rfauc,"y.values"))
rfauc <- round(rfauc,4)
txt_rfauc <- paste(c("RF: "), rfauc,sep="") # make veriable to combine text and AUC score.
#NNS
nnsauc <- performance(pred_NNS, "auc")
nnsauc <- unlist(slot(nnsauc,"y.values"))
nnsauc <- round(nnsauc,4)
txt_nnsauc <- paste(c("NNS: "), nnsauc,sep="") # make veriable to combine text and AUC score.
#Spidex
spiauc <- performance(pred_spi, "auc")
spiauc <- unlist(slot(spiauc,"y.values"))
spiauc <- round(spiauc,4)
txt_spiauc <- paste(c("Spidex: "), spiauc,sep="") # make veriable to combine text and AUC score.

# MES <- red | GS <- magenta | SSF <- blue |ADA <- cyan |rf <- gold |NNS <- green4 |spidex <- grey46
# Add AUC as legend
legend(.5,.5, c(txt_mesauc, txt_nnsauc, txt_ssfauc,  txt_gsauc, txt_adaauc, txt_rfauc),col=c(
  'red','green3','blue','grey64','cyan', 'gold'),lty=1,bty="n", title = "AUC", cex=1.2, y.intersp=0.5 )



# PAPER PMID: 27313609 : HSF, 2%; MES, 10%; NNSplice, 5%; and ASSP, 10%. 
################ Optimum cutoff/ sensativity/ spec per tool ################ 
############### SSF ############### 
opt.cut = function(perf_SSF, pred_SSF) {
  cut.ind = mapply(FUN=function(x,y,p){
    d=(x-0)^2 + (y-1)^2
    ind = which(d ==min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf_SSF@x.values, perf_SSF@y.values, pred_SSF@cutoffs)
}
print('SSF')
print(opt.cut(perf_SSF, pred_SSF))

# calaculate sens and spec for a given cut off values SSF =>  5%  
# calculate the spec and sens
spec_SSF <- performance(pred_SSF, measure = "sens", x.measure = "spec")
cutoff <-5
ix <- which.min(abs(spec_MES@alpha.values[[1]] - cutoff)) #good enough in our case
ssfsensitivity <- spec_SSF@y.values[[1]][ix] #note the order of arguments to `perfomance` and of x and y in `perf`
ssfspecificity <- spec_SSF@x.values[[1]][ix]
print("SSF")
print(ssfsensitivity)
print(ssfspecificity)

############### MES ############### 
opt.cut = function(perf_MES, pred_MES) {
  cut.ind = mapply(FUN=function(x,y,p){
    d=(x-0)^2 + (y-1)^2
    ind = which(d ==min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf_MES@x.values, perf_MES@y.values, pred_MES@cutoffs)
}
print('MES')
print(opt.cut(perf_MES, pred_MES))

# calaculate sens and spec for a given cut off values MES =>  10%
# calculate the spec and sens
spec_MES <- performance(pred_MES, measure = "sens", x.measure = "spec")
cutoff <-10
ix <- which.min(abs(spec_MES@alpha.values[[1]] - cutoff)) #good enough in our case
sensitivity <- spec_MES@y.values[[1]][ix] #note the order of arguments to `perfomance` and of x and y in `perf`
specificity <- spec_MES@x.values[[1]][ix]
print(sensitivity)
print(specificity)

############### NNS ############### 
opt.cut = function(perf_NNS, pred_NNS) {
  cut.ind = mapply(FUN=function(x,y,p){
    d=(x-0)^2 + (y-1)^2
    ind = which(d ==min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf_NNS@x.values, perf_NNS@y.values, pred_NNS@cutoffs)
}
print('NNS')
print(opt.cut(perf_NNS, pred_NNS))

# calaculate sens and spec for a given cut off values NNS =>  5%  
# calculate the spec and sens
spec_NNS <- performance(pred_NNS, measure = "sens", x.measure = "spec")
cutoff <-5
ix <- which.min(abs(spec_NNS@alpha.values[[1]] - cutoff)) #good enough in our case
NNSsensitivity <- spec_NNS@y.values[[1]][ix] #note the order of arguments to `perfomance` and of x and y in `perf`
NNSspecificity <- spec_NNS@x.values[[1]][ix]
print("NNS")
print(NNSsensitivity)
print(NNSspecificity)

############### GS ############### 
opt.cut = function(perf_GS, pred_GS) {
  cut.ind = mapply(FUN=function(x,y,p){
    d=(x-0)^2 + (y-1)^2
    ind = which(d ==min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf_GS@x.values, perf_GS@y.values, pred_GS@cutoffs)
}
print('GS')
print(opt.cut(perf_GS, pred_GS))

# calaculate sens and spec for a given cut off values GS =>  5%  
# calculate the spec and sens
spec_GS <- performance(pred_GS, measure = "sens", x.measure = "spec")
cutoff <-10
ix <- which.min(abs(spec_GS@alpha.values[[1]] - cutoff)) #good enough in our case
GSsensitivity <- spec_GS@y.values[[1]][ix] #note the order of arguments to `perfomance` and of x and y in `perf`
GSspecificity <- spec_GS@x.values[[1]][ix]
print("GS")
print(GSsensitivity)
print(GSspecificity)

############### Ada ############### 
opt.cut = function(perf_ada, pred_ada) {
  cut.ind = mapply(FUN=function(x,y,p){
    d=(x-0)^2 + (y-1)^2
    ind = which(d ==min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf_ada@x.values, perf_ada@y.values, pred_ada@cutoffs)
}
print('ada')
print(opt.cut(perf_ada, pred_ada))

# calaculate sens and spec for a given cut off values ada =>  5%  
# calculate the spec and sens
spec_ada <- performance(pred_ada, measure = "sens", x.measure = "spec")
cutoff <-0.4
ix <- which.min(abs(spec_ada@alpha.values[[1]] - cutoff)) #good enough in our case
adasensitivity <- spec_ada@y.values[[1]][ix] #note the order of arguments to `perfomance` and of x and y in `perf`
adaspecificity <- spec_ada@x.values[[1]][ix]
print("ada")
print(adasensitivity)
print(adaspecificity)

############### rf ############### 
opt.cut = function(perf_rf, pred_rf) {
  cut.ind = mapply(FUN=function(x,y,p){
    d=(x-0)^2 + (y-1)^2
    ind = which(d ==min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf_rf@x.values, perf_rf@y.values, pred_rf@cutoffs)
}
print('rf')
print(opt.cut(perf_rf, pred_rf))

# calaculate sens and spec for a given cut off values rf =>  5%  
# calculate the spec and sens
spec_rf <- performance(pred_rf, measure = "sens", x.measure = "spec")
cutoff <-0.4
ix <- which.min(abs(spec_rf@alpha.values[[1]] - cutoff)) #good enough in our case
rfsensitivity <- spec_rf@y.values[[1]][ix] #note the order of arguments to `perfomance` and of x and y in `perf`
rfspecificity <- spec_rf@x.values[[1]][ix]
print("rf")
print(rfsensitivity)
print(rfspecificity)

############### spidex ############### 
opt.cut = function(perf_spi, pred_spi) {
  cut.ind = mapply(FUN=function(x,y,p){
    d=(x-0)^2 + (y-1)^2
    ind = which(d ==min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf_spi@x.values, perf_spi@y.values, pred_spi@cutoffs)
}
print('spi')
print(opt.cut(perf_spi, pred_spi))

# calaculate sens and spec for a given cut off values rf =>  5%  
# calculate the spec and sens
spec_spi <- performance(pred_spi, measure = "sens", x.measure = "spec")
cutoff <-2
ix <- which.min(abs(spec_spi@alpha.values[[1]] - cutoff)) #good enough in our case
spisensitivity <- spec_spi@y.values[[1]][ix] #note the order of arguments to `perfomance` and of x and y in `perf`
spispecificity <- spec_spi@x.values[[1]][ix]
print("spi")
print(spisensitivity)
print(spispecificity)