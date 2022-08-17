#!/usr/bin/Rscript
# usage of this script: ./02_Rscript_rzooroh.R --chrom chrom_name OR ./02_Rscript_rzooroh.R -c chrom_name

# specify library path in R and install/load packages
.libPaths("/mnt/scratch/nodelete/c1322153/Rpackages4")
#install.packages("RZooRoH", lib = "/mnt/scratch/nodelete/c1322153/Rpackages4", repos = "https://www.stats.bris.ac.uk/R/")
#install.packages("R.utils", lib = "/mnt/scratch/nodelete/c1322153/Rpackages4", repos = "https://www.stats.bris.ac.uk/R/")
#install.packages("optparse", lib = "/mnt/scratch/nodelete/c1322153/Rpackages4", repos = "https://www.stats.bris.ac.uk/R/")
library(RZooRoH)
library(R.utils)
library(dplyr)
library(optparse)

# reading in extra variable -c, and it becomes opt$c
option_list = list(
  make_option(c("-c", "--chrom"), action="store", default=NA, type='character',
              help="name of chromosome"))
opt = parse_args(OptionParser(option_list=option_list))

#### Reading in the data ####
east_in <- paste(c("east_mLutLut_dp_q.vcf.gz_",opt$c,".gen.gz"), collapse="")
east_samples <- paste(c("east_mLutLut_dp_q.vcf.gz_",opt$c,".samples"), collapse="")
east <- zoodata(genofile = east_in, zformat = "gp",
                 samplefile = east_samples)
rm(east_in, east_samples)

wales_in <- paste(c("wales_mLutLut_dp_q.vcf.gz_",opt$c,".gen.gz"), collapse="")
wales_samples <- paste(c("wales_mLutLut_dp_q.vcf.gz_",opt$c,".samples"), collapse="")
wales <- zoodata(genofile = wales_in, zformat = "gp",
                samplefile = wales_samples)
rm(wales_in, wales_samples)

sweng_in <- paste(c("sweng_mLutLut_dp_q.vcf.gz_",opt$c,".gen.gz"), collapse="")
sweng_samples <- paste(c("sweng_mLutLut_dp_q.vcf.gz_",opt$c,".samples"), collapse="")
sweng <- zoodata(genofile = sweng_in, zformat = "gp",
                samplefile = sweng_samples)
rm(sweng_in, sweng_samples)

north_in <- paste(c("north_mLutLut_dp_q.vcf.gz_",opt$c,".gen.gz"), collapse="")
north_samples <- paste(c("north_mLutLut_dp_q.vcf.gz_",opt$c,".samples"), collapse="")
north <- zoodata(genofile = north_in, zformat = "gp",
                samplefile = north_samples)
rm(north_in, north_samples)

#### Defining the model ####
model5 <- zoomodel(K=10, krates=c(8, 16, 32, 64, 128, 256, 512, 1024, 2048, 5120))

#### Running zoorun ####
wales_model5 <- zoorun(model5, wales)
east_model5 <- zoorun(model5, east)
sweng_model5 <- zoorun(model5, sweng)
north_model5 <- zoorun(model5, north)

#### Saving/writing results ####
# Save BICs
bic <- c(wales_model5@modbic, north_model5@modbic, east_model5@modbic, sweng_model5@modbic)
write.csv(bic, file=paste("bic_model5_",opt$c,".csv", collpase=""))

# Plots of results
pdf(file=paste("plot1_model5_",opt$c,".pdf", collapse=""),height=6,width=10)
zooplot_prophbd(list(Wales=wales_model5,North=north_model5,
                     East=east_model5, SWEng=sweng_model5),style='barplot')
dev.off()

pdf(file=paste("plot2_model5_",opt$c,".pdf", collapse=""),height=6,width=10)
zooplot_prophbd(list(Wales=wales_model5,North=north_model5,
                     East=east_model5, SWEng=sweng_model5),style='lines', cumulative = TRUE)
dev.off()

pdf(file=paste("plot3_model5_",opt$c,".pdf", collapse=""),height=12,width=12)
zooplot_individuals(list(Wales=wales_model5,North=north_model5,
                         East=east_model5, SWEng=sweng_model5), cumulative = T)
dev.off()

pdf(file=paste("plot4_model5_",opt$c,".pdf", collapse=""),height=6,width=10)
zooplot_partitioning(list(Wales=wales_model5,North=north_model5,
                          East=east_model5, SWEng=sweng_model5), nonhbd = FALSE,
                     col=c("yellow", "orange", "pink", "green", "purple", "blue", "red",
                           "brown", "cyan", "coral"))
dev.off()

pdf(file=paste("plot5_model5_",opt$c,".pdf", collapse=""),height=12,width=12)
zooplot_hbdseg(list(Wales=wales_model5,North=north_model5,
                    East=east_model5, SWEng=sweng_model5), chr = 1, coord=c(0,40000000))
dev.off()

pdf(file=paste("plot5a_model5_",opt$c,".pdf", collapse=""),height=12,width=12)
zooplot_hbdseg(list(Wales=wales_model5,North=north_model5,
                    East=east_model5, SWEng=sweng_model5), chr = 1, coord=c(0,10000000))
dev.off()

# Re-structuring and writing results as csv
df1_V2 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,1]))
df1_V3 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,2]))
df1_V4 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,3]))
df1_V5 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,4]))
df1_V6 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,5]))
df1_V7 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,6]))
df1_V8 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,7]))
df1_V9 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,8]))
df1_V10 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,9]))
df1_V11 <- as.data.frame(cbind(wales_model5@sampleids, wales_model5@realized[,10]))

df1_V2$kclass <- "8"
df1_V3$kclass <- "16"
df1_V4$kclass <- "32"
df1_V5$kclass <- "64"
df1_V6$kclass <- "128"
df1_V7$kclass <- "256"
df1_V8$kclass <- "512"
df1_V9$kclass <- "1024"
df1_V10$kclass <- "2048"
df1_V11$kclass <- "non-HBD"

df1_wales <- rbind(df1_V2, df1_V3, df1_V4, df1_V5, df1_V6, df1_V7, df1_V8, df1_V9, df1_V10, df1_V11)
df1_wales$pop <- "Wales"
df1_wales$v2 <- as.numeric(df1_wales$V2)
df1_wales$kclass <- as.factor(df1_wales$kclass)
summary(df1_wales)


df1_V2 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,1]))
df1_V3 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,2]))
df1_V4 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,3]))
df1_V5 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,4]))
df1_V6 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,5]))
df1_V7 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,6]))
df1_V8 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,7]))
df1_V9 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,8]))
df1_V10 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,9]))
df1_V11 <- as.data.frame(cbind(east_model5@sampleids, east_model5@realized[,10]))

df1_V2$kclass <- "8"
df1_V3$kclass <- "16"
df1_V4$kclass <- "32"
df1_V5$kclass <- "64"
df1_V6$kclass <- "128"
df1_V7$kclass <- "256"
df1_V8$kclass <- "512"
df1_V9$kclass <- "1024"
df1_V10$kclass <- "2048"
df1_V11$kclass <- "non-HBD"

df1_east <- rbind(df1_V2, df1_V3, df1_V4, df1_V5, df1_V6, df1_V7, df1_V8, df1_V9, df1_V10, df1_V11)
df1_east$pop <- "east"
df1_east$v2 <- as.numeric(df1_east$V2)
df1_east$kclass <- as.factor(df1_east$kclass)
summary(df1_east)


df1_V2 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,1]))
df1_V3 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,2]))
df1_V4 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,3]))
df1_V5 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,4]))
df1_V6 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,5]))
df1_V7 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,6]))
df1_V8 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,7]))
df1_V9 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,8]))
df1_V10 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,9]))
df1_V11 <- as.data.frame(cbind(north_model5@sampleids, north_model5@realized[,10]))

df1_V2$kclass <- "8"
df1_V3$kclass <- "16"
df1_V4$kclass <- "32"
df1_V5$kclass <- "64"
df1_V6$kclass <- "128"
df1_V7$kclass <- "256"
df1_V8$kclass <- "512"
df1_V9$kclass <- "1024"
df1_V10$kclass <- "2048"
df1_V11$kclass <- "non-HBD"

df1_north <- rbind(df1_V2, df1_V3, df1_V4, df1_V5, df1_V6, df1_V7, df1_V8, df1_V9, df1_V10, df1_V11)
df1_north$pop <- "north"
df1_north$v2 <- as.numeric(df1_north$V2)
df1_north$kclass <- as.factor(df1_north$kclass)
summary(df1_north)


df1_V2 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,1]))
df1_V3 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,2]))
df1_V4 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,3]))
df1_V5 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,4]))
df1_V6 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,5]))
df1_V7 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,6]))
df1_V8 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,7]))
df1_V9 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,8]))
df1_V10 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,9]))
df1_V11 <- as.data.frame(cbind(sweng_model5@sampleids, sweng_model5@realized[,10]))

df1_V2$kclass <- "8"
df1_V3$kclass <- "16"
df1_V4$kclass <- "32"
df1_V5$kclass <- "64"
df1_V6$kclass <- "128"
df1_V7$kclass <- "256"
df1_V8$kclass <- "512"
df1_V9$kclass <- "1024"
df1_V10$kclass <- "2048"
df1_V11$kclass <- "non-HBD"

df1_sweng <- rbind(df1_V2, df1_V3, df1_V4, df1_V5, df1_V6, df1_V7, df1_V8, df1_V9, df1_V10, df1_V11)
df1_sweng$pop <- "sweng"
df1_sweng$v2 <- as.numeric(df1_sweng$V2)
df1_sweng$kclass <- as.factor(df1_sweng$kclass)
summary(df1_sweng)


df1_all <- rbind(df1_wales, df1_east, df1_sweng, df1_north)
df1_all$id <- as.factor(df1_all$V1)
df1_all$pop <- as.factor(df1_all$pop)
df1_all$roh <- as.numeric(df1_all$V2)
towrite <- select(df1_all, id, pop, kclass, roh)
summary(towrite)
write.csv(towrite, file=paste("rzooroh_froh_model5_",opt$c,".csv",collapse=""))
