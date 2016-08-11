#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)


args = commandArgs(trailingOnly=TRUE)

inputFile = args[1]
# inputFile = "/Users/charles/programmation/perl/bench_leon/example/report.tab"


Unicorn = read.table(inputFile,header=TRUE)


########## TREAT for size compression
mdat2 = melt(Unicorn, id.vars=c("File_ID","Size_fastQ"),
             measure.vars=c("Size_comp_leon_lossy", "Size_comp_leon_lossless", "Size_comp_gzip1", "Size_comp_gzip6", "Size_comp_gzip9"))

mdat2$ID = paste(mdat2$File_ID," - ", round(mdat2$Size_fastQ/1024/1024/1024,2), "Go")

out <- split(mdat2,f = mdat2$ID)
for(i in 1:length(out)) { out[[i]]$PC = paste(round(100-out[[i]]$value*100/out[[i]]$Size_fastQ,2),"%")}

combined = vector('list')
for(i in 1:length(out)){ combined = rbind(combined,out[[i]]) }

# BOX PLOT

g <- ggplot(combined, aes(x=variable,y=100-(value*100/Size_fastQ), colour = variable))
boxplot <- g + geom_boxplot(notch = FALSE) + 
  scale_y_continuous(breaks = seq(0,100, by=5)) + 
  geom_jitter(width = 0.2) + 
  facet_grid(.~variable) + 
  facet_grid(.~variable, scales = "free") + 
  theme(legend.position="none") +
  labs(list(title = paste("Compression rate between LEON\nand GZIP with",length(out),"fastQ of different size"), x = "Compression mode", y = "Ratio of compression (%)"))

#boxplot
# linear graph (compression function of the size)
g <- ggplot(combined, aes(x=Size_fastQ,y=100-(value*100/Size_fastQ), colour = variable))
linear_compression <- g + stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange") + 
  coord_cartesian(ylim=c(0,100)) + 
  scale_y_continuous(breaks = seq(0,100, by=10)) +
  scale_x_continuous(breaks = c((1024^2)*100,(1024^3),(1024^3)*5,(1024^3)*10,(1024^3)*15,(1024^3)*20,(1024^3)*25), labels = c("100 Mo","1 Go","5 Go","10 Go","15 Go","20 Go","25 Go")) +
  stat_summary(fun.y = mean, geom = "line") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1.2, size =7)) +
  labs(list(title = "Compression rate between LEON\nand GZIP depending on the size of FASTQ", x = "Size of initial FastQ", y = "Ratio of compression (%)"))


#linear_compression
############### TIME
#Unicorn = read.table("/Users/charles/Documents/time_compress.tab",header=TRUE,sep = "\t")


mdat2 = melt(Unicorn, id.vars=c("File_ID","Size_fastQ"),
             measure.vars=c("Time_comp_leon_lossy", "Time_comp_leon_lossless",	"Time_comp_gzip1", "Time_comp_gzip6", "Time_comp_gzip9", "Time_uncomp_leon_lossy", "Time_uncomp_leon_lossless",	"Time_uncomp_gzip1", "Time_uncomp_gzip6", "Time_uncomp_gzip9"))


mdat2$ID = paste(mdat2$File_ID,round(mdat2$Size_fastQ/1024/1024/1024,2), "Go")
mdat2$type = gsub("_\\w*","",gsub("Time_uncomp", "uncompress", mdat2$variable))
mdat2$type = gsub("Time", "compress", mdat2$type)
mdat2$soft = gsub("\\w*_", "", mdat2$variable)

#########################

g <- ggplot(mdat2, aes(x=Size_fastQ,y=(value/60),colour=soft))

time2 = g + facet_wrap(~type, scales = "free") + 
  scale_x_continuous(breaks = c((1024^2)*100,(1024^3),(1024^3)*5,(1024^3)*10,(1024^3)*15,(1024^3)*20,(1024^3)*25), labels = c("100 Mo","1 Go","5 Go","10 Go","15 Go","20 Go","25 Go")) + 
  stat_summary(fun.y = mean, geom = "line") + 
  stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
  theme(axis.text.x = element_text(angle = 50, hjust = 1.2, size =7)) +
  labs(list(title = "Time for compress/uncompress with LEON\nand GZIP depending on the size of FASTQ", x = "Size of initial FastQ", y = "Time (in minutes)"))

#time
###########################
# Save Graph
png(filename="/Users/charles/programmation/perl/bench_leon/example/boxplot_compression.png")
plot(boxplot)
dev.off()

png(filename="/Users/charles/programmation/perl/bench_leon/example/point_compression.png")
plot(linear_compression)
dev.off()

png(filename="/Users/charles/programmation/perl/bench_leon/example/point_time.png")
plot(time)
dev.off()
