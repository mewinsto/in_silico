##IN SILICO APEK1 DIGEST
library(ggplot2)


##ECITON
fragz_Ebur = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/frag_lengths_Ebur")
ra_Ebur = fragz_Ebur$V1

low_range = 300
high_range = 800
dig_length_Ebur = length(ra_Ebur)
count_Ebur = 0

for (i in 1:dig_length_Ebur){
  a = ifelse(low_range<ra_Ebur[i],1,0)
  b = ifelse(high_range>ra_Ebur[i],1,0)
  count_Ebur = count_Ebur + a*b
}

plot(density(log(ra_Ebur)),main="Fragment Distribution: In silico digest (ApeK1) of Eciton burchellii genome",xlab="Log Fragment Size", ylab="Frequency",sub="N = 552,698",xlim=c(0,10))
abline(v=log(low_range),lty=3)
abline(v=log(high_range),lty=3)

prop_Ebur = count_Ebur/dig_length_Ebur
count_Ebur

cut_freq_Ebur = dig_length_Ebur/sum(ra_Ebur)
G_size_Ebur = sum(ra_Ebur)

##Plot with ggplot (Eciton)
#Can adjust bandwidth with stat_density(adjust = XX), but current res is good.
density_plot_Ebur = ggplot(fragz_Ebur, aes(x = log(V1))) + geom_density(fill="gray") + xlab("Log Fragment Size") + ylab("Frequency") + ggtitle("Fragment Distribution: In silico digest (ApeK1) of Eciton burchellii genome") + geom_vline(aes(xintercept=log(low_range)),linetype="dashed", size = 1) + geom_vline(aes(xintercept=log(high_range)),linetype="dashed",size = 1)

##Atta Cephalotes
fragz_Acep = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/frag_lengths_Acep")
ra_Acep = fragz_Acep$V1

low_range = 300
high_range = 800
dig_length_Acep = length(ra_Acep)
count_Acep = 0

for (i in 1:dig_length_Acep){
  a = ifelse(low_range<ra_Acep[i],1,0)
  b = ifelse(high_range>ra_Acep[i],1,0)
  count_Acep = count_Acep + a*b
}

plot(density(log(ra_Acep)),main="Fragment Distribution: In silico digest (ApeK1) of Atta cephalotes genome",xlab="Log Fragment Size", ylab="Frequency",sub="N = 552,698",xlim=c(0,10))
abline(v=log(low_range),lty=3)
abline(v=log(high_range),lty=3)

prop_Acep = count_Acep/dig_length_Acep
count_Acep

cut_freq_Acep = dig_length_Acep/sum(ra_Acep)
G_size_Acep = sum(ra_Acep)

##Plot with ggplot (Atta Cephalotes)
density_plot_Acep = ggplot(fragz_Acep, aes(x = log(V1))) + geom_density(fill="gray") + xlab("Log Fragment Size") + ylab("Frequency") + ggtitle("Fragment Distribution: In silico digest (ApeK1) of Atta cephalotes genome") + geom_vline(aes(xintercept=log(low_range)),linetype="dashed", size = 1) + geom_vline(aes(xintercept=log(high_range)),linetype="dashed",size = 1)

##Atta echiniator
fragz_Aech = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/frag_lengths_Aech")
ra_Aech = fragz_Aech$V1

low_range = 300
high_range = 800
dig_length_Aech = length(ra_Aech)
count_Aech = 0

for (i in 1:dig_length_Aech){
  a = ifelse(low_range<ra_Aech[i],1,0)
  b = ifelse(high_range>ra_Aech[i],1,0)
  count_Aech = count_Aech + a*b
}

plot(density(log(ra_Aech)),main="Fragment Distribution: In silico digest (ApeK1) of Acromyrmex echiniator genome",xlab="Log Fragment Size", ylab="Frequency",sub="N = 552,698",xlim=c(0,10))
abline(v=log(low_range),lty=3)
abline(v=log(high_range),lty=3)

prop_Aech = count_Aech/dig_length_Aech
count_Aech

cut_freq_Aech = dig_length_Aech/sum(ra_Aech)
G_size_Aech = sum(ra_Aech)

##Plot with ggplot (Acromyrmex echiniator)
density_plot_Aech = ggplot(fragz_Aech, aes(x = log(V1))) + geom_density(fill="gray") + xlab("Log Fragment Size") + ylab("Frequency") + ggtitle("Fragment Distribution: In silico digest (ApeK1) of Acromyrmex echiniator genome") + geom_vline(aes(xintercept=log(low_range)),linetype="dashed", size = 1) + geom_vline(aes(xintercept=log(high_range)),linetype="dashed",size = 1)

##Camponotus floridianus
fragz_Cflo = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/frag_lengths_Cflo")
ra_Cflo = fragz_Cflo$V1

low_range = 300
high_range = 800
dig_length_Cflo = length(ra_Cflo)
count_Cflo = 0

for (i in 1:dig_length_Cflo){
  a = ifelse(low_range<ra_Cflo[i],1,0)
  b = ifelse(high_range>ra_Cflo[i],1,0)
  count_Cflo = count_Cflo + a*b
}

plot(density(log(ra_Cflo)),main="Fragment Distribution: In silico digest (ApeK1) of Camponotus floridianus genome",xlab="Log Fragment Size", ylab="Frequency",sub="N = 552,698",xlim=c(0,10))
abline(v=log(low_range),lty=3)
abline(v=log(high_range),lty=3)

prop_Cflo = count_Cflo/dig_length_Cflo
count_Cflo

cut_freq_Cflo = dig_length_Cflo/sum(ra_Cflo)
G_size_Cflo = sum(ra_Cflo)

##Plot with ggplot (Camponotus floridianus)
density_plot_Cflo = ggplot(fragz_Cflo, aes(x = log(V1))) + geom_density(fill="gray") + xlab("Log Fragment Size") + ylab("Frequency") + ggtitle("Fragment Distribution: In silico digest (ApeK1) of Camponotus floridianus genome") + geom_vline(aes(xintercept=log(low_range)),linetype="dashed", size = 1) + geom_vline(aes(xintercept=log(high_range)),linetype="dashed",size = 1)

##COBS
fragz_Cobs = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/frag_lengths_Cobs")
ra_Cobs = fragz_Cobs$V1

low_range = 300
high_range = 800
dig_length_Cobs = length(ra_Cobs)
count_Cobs = 0

for (i in 1:dig_length_Cobs){
  a = ifelse(low_range<ra_Cobs[i],1,0)
  b = ifelse(high_range>ra_Cobs[i],1,0)
  count_Cobs = count_Cobs + a*b
}

plot(density(log(ra_Cobs)),main="Fragment Distribution: In silico digest (ApeK1) of COBS genome",xlab="Log Fragment Size", ylab="Frequency",sub="N = 552,698",xlim=c(0,10))
abline(v=log(low_range),lty=3)
abline(v=log(high_range),lty=3)

prop_Cobs = count_Cobs/dig_length_Cobs
count_Cobs

cut_freq_Cobs = dig_length_Cobs/sum(ra_Cobs)
G_size_Cobs = sum(ra_Cobs)

##Plot with ggplot (COBS)
density_plot_Cobs = ggplot(fragz_Cobs, aes(x = log(V1))) + geom_density(fill="gray") + xlab("Log Fragment Size") + ylab("Frequency") + ggtitle("Fragment Distribution: In silico digest (ApeK1) of COBS genome") + geom_vline(aes(xintercept=log(low_range)),linetype="dashed", size = 1) + geom_vline(aes(xintercept=log(high_range)),linetype="dashed",size = 1)

##Harpegnathos saltator
fragz_Hsal = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/frag_lengths_Hsal")
ra_Hsal = fragz_Hsal$V1

low_range = 300
high_range = 800
dig_length_Hsal = length(ra_Hsal)
count_Hsal = 0

for (i in 1:dig_length_Hsal){
  a = ifelse(low_range<ra_Hsal[i],1,0)
  b = ifelse(high_range>ra_Hsal[i],1,0)
  count_Hsal = count_Hsal + a*b
}

plot(density(log(ra_Hsal)),main="Fragment Distribution: In silico digest (ApeK1) of Harpegnathos saltator genome",xlab="Log Fragment Size", ylab="Frequency",sub="N = 552,698",xlim=c(0,10))
abline(v=log(low_range),lty=3)
abline(v=log(high_range),lty=3)

prop_Hsal = count_Hsal/dig_length_Hsal
count_Hsal

cut_freq_Hsal = dig_length_Hsal/sum(ra_Hsal)
G_size_Hsal = sum(ra_Hsal)

##Plot with ggplot (Harpegnathos saltator)
density_plot_Hsal = ggplot(fragz_Hsal, aes(x = log(V1))) + geom_density(fill="gray") + xlab("Log Fragment Size") + ylab("Frequency") + ggtitle("Fragment Distribution: In silico digest (ApeK1) of Harpegnathos saltator genome") + geom_vline(aes(xintercept=log(low_range)),linetype="dashed", size = 1) + geom_vline(aes(xintercept=log(high_range)),linetype="dashed",size = 1)

##Linepithema humile
fragz_Lhum = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/frag_lengths_Lhum")
ra_Lhum = fragz_Lhum$V1

low_range = 300
high_range = 800
dig_length_Lhum = length(ra_Lhum)
count_Lhum = 0

for (i in 1:dig_length_Lhum){
  a = ifelse(low_range<ra_Lhum[i],1,0)
  b = ifelse(high_range>ra_Lhum[i],1,0)
  count_Lhum = count_Lhum + a*b
}

plot(density(log(ra_Lhum)),main="Fragment Distribution: In silico digest (ApeK1) of Linepithema humile genome",xlab="Log Fragment Size", ylab="Frequency",sub="N = 552,698",xlim=c(0,10))
abline(v=log(low_range),lty=3)
abline(v=log(high_range),lty=3)

prop_Lhum = count_Lhum/dig_length_Lhum
count_Lhum

cut_freq_Lhum = dig_length_Lhum/sum(ra_Lhum)
G_size_Lhum = sum(ra_Lhum)

##Plot with ggplot (Linepithema humile)
density_plot_Lhum = ggplot(fragz_Lhum, aes(x = log(V1))) + geom_density(fill="gray") + xlab("Log Fragment Size") + ylab("Frequency") + ggtitle("Fragment Distribution: In silico digest (ApeK1) of Linepithema humile genome") + geom_vline(aes(xintercept=log(low_range)),linetype="dashed", size = 1) + geom_vline(aes(xintercept=log(high_range)),linetype="dashed",size = 1)

##Pogonomyrmex barbatus
fragz_Pbar = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/frag_lengths_Pbar")
ra_Pbar = fragz_Pbar$V1

low_range = 300
high_range = 800
dig_length_Pbar = length(ra_Pbar)
count_Pbar = 0

for (i in 1:dig_length_Pbar){
  a = ifelse(low_range<ra_Pbar[i],1,0)
  b = ifelse(high_range>ra_Pbar[i],1,0)
  count_Pbar = count_Pbar + a*b
}

plot(density(log(ra_Pbar)),main="Fragment Distribution: In silico digest (ApeK1) of Pogonomyrmex barbatus genome",xlab="Log Fragment Size", ylab="Frequency",sub="N = 552,698",xlim=c(0,10))
abline(v=log(low_range),lty=3)
abline(v=log(high_range),lty=3)

prop_Pbar = count_Pbar/dig_length_Pbar
count_Pbar

cut_freq_Pbar = dig_length_Pbar/sum(ra_Pbar)
G_size_Pbar = sum(ra_Pbar)

##Plot with ggplot (Pogonomyrmex barbatus)
density_plot_Pbar = ggplot(fragz_Pbar, aes(x = log(V1))) + geom_density(fill="gray") + xlab("Log Fragment Size") + ylab("Frequency") + ggtitle("Fragment Distribution: In silico digest (ApeK1) of Pogonomyrmex barbatus genome") + geom_vline(aes(xintercept=log(low_range)),linetype="dashed", size = 1) + geom_vline(aes(xintercept=log(high_range)),linetype="dashed",size = 1)

##Solenopsis invicta
fragz_Sinv = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/frag_lengths_Sinv")
ra_Sinv = fragz_Sinv$V1

low_range = 300
high_range = 800
dig_length_Sinv = length(ra_Sinv)
count_Sinv = 0

for (i in 1:dig_length_Sinv){
  a = ifelse(low_range<ra_Sinv[i],1,0)
  b = ifelse(high_range>ra_Sinv[i],1,0)
  count_Sinv = count_Sinv + a*b
}

plot(density(log(ra_Sinv)),main="Fragment Distribution: In silico digest (ApeK1) of Solenopsis invicta genome",xlab="Log Fragment Size", ylab="Frequency",sub="N = 552,698",xlim=c(0,10))
abline(v=log(low_range),lty=3)
abline(v=log(high_range),lty=3)

prop_Sinv = count_Sinv/dig_length_Sinv
count_Sinv

cut_freq_Sinv = dig_length_Sinv/sum(ra_Sinv)
G_size_Sinv = sum(ra_Sinv)

##Plot with ggplot (Solenopsis invicta)
density_plot_Sinv = ggplot(fragz_Sinv, aes(x = log(V1))) + geom_density(fill="gray") + xlab("Log Fragment Size") + ylab("Frequency") + ggtitle("Fragment Distribution: In silico digest (ApeK1) of Solenopsis invicta genome") + geom_vline(aes(xintercept=log(low_range)),linetype="dashed", size = 1) + geom_vline(aes(xintercept=log(high_range)),linetype="dashed",size = 1)

######################
## COLLATE ALL DATA ##
######################

names = c("Acep","Aech","Cflo","Cobs","Ebur","Hsal","Lhum","Pbar","Sinv")
cut_freq_all = c(cut_freq_Acep,cut_freq_Aech,cut_freq_Cflo,cut_freq_Cobs,cut_freq_Ebur,cut_freq_Hsal,cut_freq_Lhum,cut_freq_Pbar,cut_freq_Sinv)
count_all = c(count_Acep,count_Aech,count_Cflo,count_Cobs,count_Ebur,count_Hsal,count_Lhum,count_Pbar,count_Sinv)
G_size_all = c(G_size_Acep,G_size_Aech,G_size_Cflo,G_size_Cobs,G_size_Ebur,G_size_Hsal,G_size_Lhum,G_size_Pbar,G_size_Sinv)
prop_all = c(prop_Acep,prop_Aech,prop_Cflo,prop_Cobs,prop_Ebur,prop_Hsal,prop_Lhum,prop_Pbar,prop_Sinv)

stats_all = cbind(names,cut_freq_all,count_all,G_size_all,prop_all)
stats_all = as.data.frame(stats_all)

#write.table(stats_all,file="Digest_stats")

digest_stats = read.table("~/Desktop/Current Projects/Eciton GBS Phylogeny/in_silico/Digest_stats",header=TRUE)
digest_stats_est = subset(digest_stats,names != "Ebur")

plot(density(digest_stats_est$count_all),main="Distribution of Maximum Expected Loci from ApeK1 Digest",xlab="")

estimate_names = c("Cut_Freq","Selected_Loci","Genome_Size","Loci_Prop")
mean_estimates = c(mean(digest_stats_est$cut_freq_all),mean(digest_stats_est$count_all),mean(digest_stats_est$G_size_all),mean(digest_stats_est$prop_all))
sd_estimates = c(sd(digest_stats_est$cut_freq_all),sd(digest_stats_est$count_all),sd(digest_stats_est$G_size_all),sd(digest_stats_est$prop_all))

estimates = rbind(mean_estimates,sd_estimates)

##PLOT ESTIMATES
plot(density(rnorm(10000,mean=estimates[1,1],sd=estimates[2,1])),main="Distribution of Cut Frequency Estimates from Ant Genomes",xlab="Cut Frequency")
abline(v=estimates[1,1],lty=2)
abline(v=digest_stats$cut_freq_all[5],col=2)
plot(density(rnorm(10000,mean=estimates[1,2],sd=estimates[2,2])),main="Distribution of Number of Size-Selected Loci Estimates from Ant Genomes",xlab="Number of Size-Selected Loci")
abline(v=estimates[1,2],lty=2)
abline(v=digest_stats$count_all[5],col=2)
plot(density(rnorm(10000,mean=estimates[1,3],sd=estimates[2,3])),main="Genome Sizes")
plot(density(rnorm(10000,mean=estimates[1,4],sd=estimates[2,4])),main="Distribution of Proportion of Loci Size_Selected Estimates from Ant Genomes")
abline(v=estimates[1,4],lty=2)
abline(v=digest_stats$prop_all[5],col=2)

##NUMBER OF READS
all_stats_s2 = read.csv("~/Desktop/all_stats_files/all_stats_s2.csv",header=TRUE)
run = c(rep(1,48),rep(2,48),rep(3,48),rep(4,48),rep(5,48),rep(6,48),rep(7,48),rep(8,48),rep(9,48))
plot(log(all_stats_s2$Nreads),log(all_stats_s2$passed.total),col=run,main="Log Reads per Sample by Passed Reads",xlab="Log Total Reads",ylab="Log Passed Reads")
abline(a=0,b=1)
plot(log(all_stats_s2$Nreads),all_stats_s2$passed.total/all_stats_s2$Nreads,col=run,main="Log Reads per Sample by Passed Reads (Normalized)",xlab="Log Total Reads",ylab="Log Passed Reads")
