#! /bin/bash


module load samtools
module add parallel-20160622

# you need read depth at each position of reference, processing for each scaffold in parallel
parallel -j 20 "samtools depth *bam -r {} > ./_depthOfSamples/samtools.{}.depthOfSamples; echo {}.finished" :::: ../scafoldList



cd ./_depthOfSamples
# cd /auto/pruhonice1-ibot/home/mslenker/Projects/_DATA/RadSeq/RadCard3/_depthOfSamples


# Continue in R

module add r/r-4.1.1-intel-19.0.4-ilb46fy # R version 4.1.1


R

library(data.table,lib.loc="/auto/pruhonice1-ibot/home/mslenker/R/Rpacakges-R-4.1.1-intel-19.0.4-ilb46fy")

# args
minCoverage = 4  # reqquired minimal coverage of each site
percentMiss = 0.75 # acceptable anount of missing data for each site


list.of.files <- file.info(list.files(path=".", pattern=".depthOfSamples$"))
sizes <- file.info(list.files(path=".", pattern=".depthOfSamples$"))$size  # get the size for each file
inputFiles <- rownames(list.of.files)[which(sizes != 0)]   # subset the files that have non-zero size

# file="samtools.150.depthOfSamples"
for (file in inputFiles) {
    cat(file, "\n")
    dt <- fread(file,sep="\t", header=F)
    header=dt[,1:2]
    data=dt[,3:ncol(dt)]
    data$toTake = apply(data, 1, function(x) (sum(x >= minCoverage)/length(x) > percentMiss))  # coverage at least "minCoverage" in more than 75% of samples
    header = header[data$toTake,]
    #colnames(header) = c("CHROM", "POS")
    
    write.table(header, paste(file,".passedRegions", sep=""), row.names = F, quote = F, col.names = F, sep = "\t")
}

list.of.files <- file.info(list.files(path=".", pattern=".passedRegions$"))
sizes <- file.info(list.files(path=".", pattern=".passedRegions$"))$size  # get the size for each file
inputFiles <- rownames(list.of.files)[which(sizes != 0)]   # subset the files that have non-zero size
# f="samtools.101.depthOfSamples.passedRegions"

# histogram of distances among regions
distances = numeric()
for (f in inputFiles) {
    file = read.delim(f, header=F)
    dd = cbind(file$V2[-length(file$V2)], file$V2[-1]) # combine position with position+1
    dists = dd[,2] -  dd[,1]
    distances = c(distances, dists[dists > 1])  # remember distances, but omit continuous sites (distance == 1)
}

pdf()
hist(distances, breaks = seq(0, max(distances)+10, 1), ylim = c(0, 80), xlim = c(0, 3000))
dev.off()

betweenLociDist = 1000  # decision is made according the histogram. If two regions were within a distance of 1000 base pairs or less, they were deemed a single RAD locus. However, if the distance between them exceeded 1000 base pairs, the regions were treated as distinct and unlinked RAD loci.


for (f in inputFiles) {
    cat(f, "\n")
    file = read.delim(f, header=F)

    begin = NA
    regions=character()
    for (pos in file[,2]){
        if (is.na(begin)) {begin = pos; inside=pos}
        if ((pos - inside) > betweenLociDist) {
           end = inside 
           if (begin == end) {end= end + 1}
           regions = c(regions,paste(begin, "-", end, sep = ""))
           begin = pos
           inside = pos
        }
        if ((pos - inside) <= betweenLociDist) {
           inside = pos 
        }
    }
    regions = c(regions,paste(begin, "-", pos, sep = "")) # the last position

    write.table(cbind(f, regions), paste(f, ".regions.betweenLociDist", betweenLociDist, ".txt", sep = "" ),  row.names = F, quote = F, col.names = F, sep = "\t")
}


q()
# QUIT R

# Continue in BASH

cat *betweenLociDist1000.txt > ../regions1000apart.txt

sed -i 's/samtools.//' ../regions1000apart.txt
sed -i 's/.depthOfSamples.passedRegions//' ../regions1000apart.txt

















