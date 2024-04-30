#pipeline dada2 

#setwd("C:/Users/GIULIA/Desktop/phyto&zoo") #lavorare in questa cartella

#library("dada2") #caricare il pacchetto


# path <- "C:/Users/GIULIA/Desktop/phyto&zoo" # CHANGE ME to the directory containing the fastq files after unzipping, per dirgli che deve lavorare in quella cartella
# list.files(path) # digli di farti vedere tutti i file in quella cartella
# 
# # Forward and reverse fastq filenames have format: SAMPLENAME_22_001.fastq and SAMPLENAME_R2_001.fastq
# fnFs <- sort(list.files(path, pattern="1.fastq.gz", full.names = TRUE)) # gli dico di caricare tutti i file  con il nome per esempio 1.fastq.gz
# fnRs <- sort(list.files(path, pattern="2.fastq.gz", full.names = TRUE))
# # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# sample.names <- sapply(strsplit(basename(fnFs), "_1"), `[`,1)
# 
# sample.names
# 
# plotQualityProfile(fnFs[1:5])
# 
# plotQualityProfile(fnRs[1:5])
# 
# # Place filtered files in filtered/ subdirectory nuovo file Assegnare i nomi dei file per i file fastq.gz filtrati.
# filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
# filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# names(filtFs) <- sample.names
# names(filtRs) <- sample.names
# 
# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(50,50), trimRight = c(0,25),
#                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                      compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
# head(out) # ti fa vedere i primi elementi della tabella
# View(out) # per vedere tutto ciò che c'è nella tabella
# 
# errF <- learnErrors(filtFs, multithread=TRUE)
# 
# errR <- learnErrors(filtRs, multithread=TRUE)
# 
# save.image("C:/Users/GIULIA/Desktop/phyto&zoo/dada.RData")
# 
# plotErrors(errF, nominalQ=TRUE)
# 
# dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
# 
# dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# 
# mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# # Inspect the merger data.frame from the first sample
# head(mergers[[1]])
# 
# seqtab <- makeSequenceTable(mergers)
# dim(seqtab)
# 
# # Inspect distribution of sequence lengths
# table(nchar(getSequences(seqtab)))
# 
# seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seqtab.nochim)
# 
# #seqtab.nochim2 <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% 338:367]
# 
# seqtab.nochim2 <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% 338:367]
# 
# sum(seqtab.nochim2)/sum(seqtab)
# 
# getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim),rowSums(seqtab.nochim2))
# # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
# colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "lenghtreduced")
# rownames(track) <- sample.names
# head(track)
# 
# write.csv(seqtab.nochim, "sequence_tab.csv")
# 
# taxa <- assignTaxonomy(seqtab.nochim2, "C:/Users/GIULIA/Desktop/schiume/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
# 
# write.csv(taxa, "taxonomy.csv")
# taxa2 <- addSpecies(taxa, "C:/Users/GIULIA/Desktop/schiume/silva_species_assignment_v138.1.fa.gz")
# 
# taxa.print <- taxa # Removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print)
# 
load("C:/Users/Dell/Desktop/CNR/Daphnia/16Sdaphnia_gliulia/16Sdaphnia/dada.RData")
setwd("C:/Users/Dell/Desktop/CNR/Daphnia/16Sdaphnia_gliulia/16Sdaphnia")


rm(dadaFs,dadaRs,dna,errF,errR,mergers, asvs,seqtab,seqtab.nochim)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings") 
library("ggplot2")

theme_set(theme_bw())

#seqtab.nochim2<-read.csv2("sequence_tab.csv")

samples.out <- rownames(seqtab.nochim2)
vari<-as.data.frame(samples.out)
vari$Treatment<-c("0 daphnids", "0 daphnids", "0 daphnids","12 daphnids", "12 daphnids", "12 daphnids", "15 daphnids", "15 daphnids", "15 daphnids",  "3 daphnids", "3 daphnids", "3 daphnids", "6 daphnids", "6 daphnids", "6 daphnids", "9 daphnids", "9 daphnids", "9 daphnids", "WWTP effluent", "WWTP effluent", "WWTP effluent","Lake Water", "Lake Water", "Lake Water" )
vari$Replicate<-c("1", "2", "3")
vari$SampleType<-c("exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","T0","T0","T0","T0","T0","T0")

rownames(vari) <- samples.out
vari<-vari[,-1]
vari$Treatment<-factor(vari$Treatment, levels=c("0 daphnids", "3 daphnids", "6 daphnids", "9 daphnids", "12 daphnids", "15 daphnids", "Lake Water", "WWTP effluent"))
ps <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE), 
               sample_data(vari), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))

names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

table(phyloseq::tax_table(ps)[, "Phylum"], exclude = NULL) # check phyla
table(phyloseq::tax_table(ps)[, "Class"], exclude = NULL) # check classes
table(phyloseq::tax_table(ps)[, "Order"], exclude = NULL) # check orders, chloroplasts

ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) # 
table(phyloseq::tax_table(ps0)[, "Phylum"], exclude = NULL) # check phyla

ps1 <- subset_taxa(ps0, !Order %in% c("Chloroplast")) # 
table(phyloseq::tax_table(ps1)[, "Order"], exclude = NULL) # check orders  ##

plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"))

saveRDS(ps, file="ps.RDS")
library("ANCOMBC")
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
# 
colNew2<-c("turquoise4","lightblue3","lightblue4", "turquoise","gray65","gray5", "darkblue", "orange")   
ordi<-plot_ordination(ps, ord.nmds.bray, color="Treatment", title="A) nmds")+
      geom_point(size=5, alpha=0.3) +
       scale_color_manual(values=colNew2)



sample_data(ps)<-as.data.frame(c("exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","exp","T0","T0","T0","T0","T0","T0"))
psT0<-subset_samples(ps, SampleType=="T0")

#out = ancombc2(data = psT0,fix_formula="Treatment", tax_level = "Genus")

load("C:/Users/Dell/Desktop/CNR/Daphnia/16Sdaphnia_gliulia/16Sdaphnia/data_ancomBC_nclu.RData")
res = out$res

write.csv2(res, "ancombc_res.csv")
res_global = out$res_global
resDB<-as.data.frame(res)
TrueF<-subset(resDB, resDB$`diff_TreatmentWWTP effluent`=="TRUE")
tax_split<-stringr::str_split_fixed(TrueF$taxon, ":", 2)
TruN<-cbind(tax_split, TrueF)
TruG<-subset(TruN, TruN$`1`=="Genus")
TruWW<-subset(TruG, TruG$`lfc_TreatmentWWTP effluent`>0)
TruLW<-subset(TruG, TruG$`lfc_TreatmentWWTP effluent`<0)
library("microbiome")

genus_data = aggregate_taxa(ps, "Genus")
genus_table<-as.data.frame(genus_data@otu_table)
tgenus<-t(genus_table)

genusda<-tgenus[1:18,]
tgen<-t(genusda)
tgen2<-subset(tgen, rowSums(tgen)!=0)
genusda2<-t(tgen2)
vari2<-vari[1:18,]
vari2$dap<-c(0,0,0,12,12,12,15,15,15,3,3,3,6,6,6,9,9,9)
coris<-as.data.frame(cor(log(genusda2+1), vari2$dap))
sig<--0.7
neg<-subset(coris, coris$V1<sig)
library("reshape2")
tgen<-subset(genus_table, rownames(genus_table)%in%rownames(neg))
gen<-t(tgen[,1:18])
mgen<-melt(as.matrix(gen))
mgen$vari<-rep(vari2$dap, rep=nrow(mgen)/nrow(vari))
mgen$origin<-"neutral"
mgen$origin[mgen$Var2%in%TruWW$`2`]<-"WW"
mgen$origin[mgen$Var2%in%TruLW$`2`]<-"LW"

colNewBP<-c("darkblue", "gray5", "orange")
ggplot(data=mgen, aes(x=factor(vari), y=value, fill=origin, color=origin)) +
  geom_point(alpha=0.2)+
  scale_fill_manual(values=colNewBP)+
  scale_color_manual(values=colNewBP)+
  labs(title="", y="", x="Number of reads")+
  theme_minimal()+
  facet_wrap(~Var2, scales="free")

A<-ggplot(data=mgen, aes(x=factor(vari), y=value, fill=origin, color=origin)) +
  geom_point(alpha=0.2, size=3)+
  scale_fill_manual(values=colNewBP)+
  scale_color_manual(values=colNewBP)+
  labs(title="", y="genus abundance [read Nr.]", x="Number of daphnids")+
  theme_minimal()+
  facet_wrap(origin~Var2, scale="free")

A+ggh4x::facet_nested_wrap(~ origin + Var2, nest_line = element_line(linetype = 2), scale="free") +
  theme(strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "grey33"))
  
  #facet_wrap2(origin~Var2, scale="free_x", remove_labels ="y")
genus_names<-c("Beijerinckiaceae-alpha cluster", )


library(ggh4x)

qPCR<-read.csv2("elab_exp.csv")
qPCR$numero.daphnie<-factor(qPCR$numero.daphnie, levels=c("T0", "0", "3", "6", "9", "12", "15"))
qPCR<-qPCR[1:18,]


library("vegan")
library("usedist")
betabray<-vegdist(seqtab.nochim2,method="bray")
plot(hclust(betabray, method="average"),hang=-1, main='Bray-Curtis Bacteria', sub='', xlab='', cex=1) #plot cluster analysis of betapair
beta<-as.matrix(betabray)
betadf<-as.data.frame(beta)
groupd<-dist_groups(betabray, vari$Treatment)
WW<-subset(groupd, groupd$Group2=="WWTP effluent")
WW2<-subset(WW, WW$Group1!="WWTP effluent")
WW3<-subset(WW2, WW2$Group1!="Lake Water")


LW<-subset(groupd, groupd$Group2=="Lake Water")
LW2<-subset(LW, LW$Group1!="Lake Water")


LW3<-subset(LW2, LW2$Group1!="WWTP effluent")
nr<-as.data.frame(stringr::str_split_fixed(WW3$Group1, " ", 2))
lmWWdist<-lm(WW3$Distance~as.numeric(nr$V1))
summary(lmWWdist)
nr<-as.data.frame(stringr::str_split_fixed(LW3$Group1, " ", 2))
lmLWdist<-lm(LW3$Distance~as.numeric(nr$V1))
summary(lmLWdist)

library("stringr")

plot(WW$Group1, WW$Distance)
plot(LW$Group1, LW$Distance)
dist<-rbind(WW3,LW2)

ggplot(WW, aes(x=factor(Group1), y=Distance))  + geom_boxplot()
ggplot(LW, aes(x=factor(Group1), y=Distance))  + geom_boxplot()
paOTU<-as.data.frame(seqtab.nochim2)
aOTU<-as.data.frame(seqtab.nochim2)

paOTU[paOTU>1]=1
betabray2<-vegdist(paOTU,method="bray", binary=TRUE)
plot(hclust(betabray2, method="average"),hang=-1, main='Bray-Curtis Bacteria', sub='', xlab='', cex=1) #plot cluster analysis of betapair
beta2<-as.matrix(betabray2)
beta2df<-as.data.frame(beta2)
groupd<-dist_groups(betabray2, vari$Treatment)
WW<-subset(groupd, groupd$Group2=="WWTP effluent")
WW2<-subset(WW, WW$Group1!="WWTP effluent")
WW3<-subset(WW2, WW2$Group1!="Lake Water")
LW<-subset(groupd, groupd$Group2=="Lake Water")
LW2<-subset(LW, LW$Group1!="Lake Water")
LW3<-subset(LW2, LW2$Group1!="WWTP effluent")
nr<-as.data.frame(stringr::str_split_fixed(WW3$Group1, " ", 2))
lmWWdist<-lm(WW3$Distance~as.numeric(nr$V1))
summary(lmWWdist)
nr<-as.data.frame(stringr::str_split_fixed(LW3$Group1, " ", 2))
lmLWdist<-lm(LW3$Distance~as.numeric(nr$V1))
summary(lmLWdist)
colNewG<-c("darkblue", "orange")

boxplot(factor(WW$Group1), WW$distance)
plot(LW$Group1, LW$distance)
dist2<-rbind(WW3,LW2)

abudist<-ggplot(dist, aes(x=factor(Group1), y=Distance, fill=Group2,color=Group2))  + 
  geom_boxplot(alpha=0.2)+
  theme_minimal()+
  scale_fill_manual(values=colNewG)+
  scale_color_manual(values=colNewG)+
  guides(color="none")+
  theme(legend.position = c(0.9, 0.3),legend.background = element_rect(fill="grey97",linetype="solid", colour ="white") ,axis.title.y = ggtext::element_markdown())+
  labs(fill="distance to", title = "", x="")
abudist
dist$NR<-str_replace(dist$Group1, " daphnids", "")
distWW<-subset(dist, dist$Group2=="WWTP effluent" | dist$Group1!="WWTP effluent")
distLW<-subset(dist, dist$Group2=="Lake Water" | dist$Group1!="WWTP effluent")

lmWWab<-lm(distWW$Distance~as.numeric(distWW$NR))


PAdist<-ggplot(dist2, aes(x=factor(Group1), y=Distance, fill=Group2, color=Group2))  + 
  geom_boxplot(alpha=0.2)+
  theme_minimal()+
  scale_fill_manual(values=colNewG)+
  scale_color_manual(values=colNewG)+
  guides(color="none")+
  labs(fill="distance to", title = "B) presence absence based", x="")


dist2<-cowplot::plot_grid(abudist, PAdist,nrow = 2)
cowplot::plot_grid(ordi, PAdist, ncol = 2, rel_widths = c(2,3))
paOTU<-as.data.frame(seqtab.nochim2)
aOTU<-as.data.frame(seqtab.nochim2)
paOTU[paOTU>1]=1
numb<-taxa_names(ps)
colnames(paOTU)<-numb

# colnames(aOTU)<-numb
# subpa<-as.data.frame(paOTU[19:24,])
# tpa<-as.data.frame(t(subpa))
# tpa<-subset(tpa, rowSums(tpa)>0)
# pa<-as.data.frame(t(tpa))
# varis<-vari[19:24,]
# treat<-c("WWTP", "WWTP", "WWTP","Lake Water", "Lake Water", "Lake Water")
# sum<-aggregate(pa, by=list(treat),mean)
# rownames(sum)<-sum[,1]
# sum<-sum[,-1]
# tsum<-as.data.frame(t(sum))
# WWotus<-rownames(subset(tsum, tsum$`Lake Water`==0))
# LWotus<-rownames(subset(tsum, tsum$WWTP==0))
# taOTU<-t(aOTU)
# lake<-as.data.frame(t(subset(taOTU, rownames(taOTU)%in%LWotus)))
# WWTP<-as.data.frame(t(subset(taOTU, rownames(taOTU)%in%WWotus)))
# 
# meltWW<-melt(as.matrix(WWTP))
# meltWW$vari<-rep(vari$Treatment, 2374)
# 
# ggplot(data=meltWW, aes(x=value, y=factor(vari))) +
#   geom_boxplot()+
#   scale_x_log10()
# 
# WWTPda<-WWTP[1:18,]
# tww<-t(WWTPda)
# tww2<-subset(tww, rowSums(tww)!=0)
# WWTPda2<-t(tww2)
# vari2<-vari[1:18,]
# vari2$dap<-c(0,0,0,12,12,12,15,15,15,3,3,3,6,6,6,9,9,9)
# coris<-as.data.frame(cor(log(WWTPda2+1), vari2$dap))
# sig<--0.7
# neg<-filter(coris, coris<sig)
# tWW<-t(WWTP)
# tWW<-subset(tWW, rownames(tWW)%in%rownames(neg))
# WW<-t(tWW)
# mWW<-melt(as.matrix(WW))
# mWW$vari<-rep(vari$Treatment, rep=53)
# ggplot(data=mWW, aes(x=value, y=factor(vari))) +
#   geom_boxplot()+
#   facet_wrap(~Var2, scales="free")
# 
# 
# meltLake<-melt(as.matrix(lake))
# meltLake$vari<-rep(vari$Treatment, 514)
# 
# ggplot(data=meltLake, aes(x=value, y=factor(vari))) +
#   geom_boxplot()+
#   scale_x_log10()
# 
# LakeTPda<-lake[1:18,]
# tLake<-t(LakeTPda)
# tLake2<-subset(tLake, rowSums(tLake)!=0)
# LakeTPda2<-t(tLake2)
# vari2<-vari[1:18,]
# vari2$dap<-c(0,0,0,12,12,12,15,15,15,3,3,3,6,6,6,9,9,9)
# coris<-as.data.frame(cor(LakeTPda2, vari2$dap))
# sig2<--0.7
# neg<-filter(coris, coris<sig2)
# tLake<-t(lake)
# tLake<-subset(tLake, rownames(tLake)%in%rownames(neg))
# Lake<-t(tLake)
# mLake<-melt(as.matrix(Lake))
# mLake$vari<-rep(vari$Treatment, rep=17)
# ggplot(data=mLake, aes(x=value, y=factor(vari))) +
#   geom_boxplot()+
#   facet_wrap(~Var2, scales="free")
# 
# 
# 
# 
# 
# 
# aOTUda<-aOTU[1:18,]
# tOTU<-t(aOTUda)
# tOTU2<-subset(tOTU, rowSums(tOTU)!=0)
# aOTUda2<-t(tOTU2)
# vari2<-vari[1:18,]
# vari2$dap<-c(0,0,0,12,12,12,15,15,15,3,3,3,6,6,6,9,9,9)
# coris<-as.data.frame(cor(log(aOTUda2+1), vari2$dap))
# sig<--0.7
# neg<-filter(coris, coris<sig)
# tOTU<-t(aOTU)
# tOTU<-subset(tOTU, rownames(tOTU)%in%rownames(neg))
# OTU<-t(tOTU)
# 
# tOTU1<-subset(tOTU, rowSums(tOTU)>2000)
# OTU1<-t(tOTU1)
# mOTU<-melt(as.matrix(OTU1))
# 
# 
# 
# mOTU$vari<-rep(vari$Treatment, rep=46)
# mOTU$origin<-"neural"
# mOTU$origin[mOTU$Var2%in%colnames(WW)]="WW"
# mOTU$origin[mOTU$Var2%in%colnames(Lake)]="LW"
# mOTU$taxa<-mOTU$Var2
# mOTU$taxa[mOTU$taxa]
# mutate(mOTU)
# mOTU$taxa[mOTU$Var2]
# ggplot(data=mOTU, aes(x=value, y=factor(vari), color=origin)) +
#   geom_boxplot()+
#   facet_wrap(~Var2, scales="free")



####MICROBIOTA DAPHNIA ED INTI1


setwd("C:/Users/Dell/Desktop/CNR/Daphnia/16Sdaphnia_gliulia/16Sdaphnia")

library("reshape2")



library("ggplot2")
library("dplyr")
library("cowplot")
colNew2<-c("darkblue","gray25","deepskyblue1","gray65","lightblue1","gray85","orange")   

qPCR<-read.csv2("elab_exp.csv")
qPCR$numero.daphnie<-factor(qPCR$numero.daphnie, levels=c("T0", "0", "3", "6", "9", "12", "15"))
qPCR<-qPCR[1:18,]

LWWW<-ggplot(data=qPCR, aes(x = nr, y = int1.16s)) +
  geom_point(size=1, alpha=0.55) +
  xlab("# Daphnia")+
  ylab("intI1 / 16S rRNA gene copy")+
  theme_minimal()+
  theme(axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
  labs(y="*int*I1 per 16S rRNA gene copy", x="Number of Daphnia")+
  geom_smooth(method='lm', color="grey56")
lml<-lm(qPCR$int1.16s~qPCR$nr)
summary(lml)

qPCR_lake<-read.csv2("daphnia_lago.csv")
LW<-ggplot(data=qPCR_lake, aes(x = dap, y = intI1_16S/10)) +
    geom_point(size=1,alpha=0.55) +
  xlab("# Daphnia")+
  ylab("intI1 / 16S rRNA gene copy")+
  theme_minimal()+
    theme(axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10),axis.title.y = ggtext::element_markdown())+
  labs(y="*int*I1 per 16S rRNA gene copy", x="")+
  geom_smooth(method='lm', color="grey56")

lml<-lm(qPCR_lake$intI1_16S~qPCR_lake$dap)
summary(lml)

cowplot::plot_grid(LW, LWWW, ncol=1,labels = c("A","B"))



qPCRcip<-read.csv2("elab_cipais_new.csv")

y_title<-expression(paste(italic("int"), "I1 per 16S rRNA gene copy"))
names(labels1)<-c("fraction >0.2um","fraction >1.2um", "fraction >50um")
#season<-
colNewG<-c("black", "darkblue", "grey50")
labely<-as_labeller(c(">0.2um"=">0.2 \U00B5m",">1.2um"=">1.2 \U00B5m",">50um"=">50 \U00B5m"))
 dots<- ggplot(data=qPCRcip, aes(x = factor(month), y = int1_16s, color = as.factor(year))) +
  geom_point(size=3, alpha=0.4) +
  scale_color_manual(values=colNewG)+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = c(0.9, 0.3),legend.background = element_rect(fill="grey97",linetype="solid", colour ="white"),axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
    labs(y="*int*I1 per 16S rRNA gene copy", x="month", color="year")+
  facet_grid(.~fraction,scales="free", labeller=labely)

box<-ggplot(data=qPCRcip, aes(x=fraction, y=int1_16s, fill=place))+
  geom_boxplot(alpha=0.2)+
  theme_minimal()+
  scale_fill_manual(values=colNew2)+
    geom_jitter(alpha=0.5, width = 0.1)+
  theme(legend.position = c(0.8, 0.3),legend.background = element_rect(fill="grey97",linetype="solid", colour ="white"),axis.title.y = ggtext::element_markdown())+
  labs(y="*int*I1 per 16S rRNA gene copy", x=NULL, )+
  scale_x_discrete(position = "top", labels=c(">0.2 \U00B5m",">1.2 \U00B5m", ">50 \U00B5m" )) +
  scale_y_log10()
# #temp<-ggplot(data=qPCRcip, aes(x=factor(year), y=int1_16s, color=fraction))+
#   geom_boxplot()+
#   geom_jitter(alpha=0.5, position=position_dodge(0.7))+
#   scale_y_log10()
cowplot::plot_grid(box, dots,rel_widths = c(3,6), rel_heights = c(1), align = "h", axis = "bt", labels = c("A","B"))











chem<-read.csv2("chem.csv")
cor(chem)
dfTOC<-read.csv2("TOCTEP.csv")
m5<-subset(dfTOC, dfTOC$m==5)
dfTEP<-read.csv2("TEP.csv")
m5TEP<-subset(dfTEP, dfTEP$m==5)
chem$TOC<-m5$TOC[-c(1:12)]
lmc<-glm(cbind(int1,X16S)~m5$TOC+Temp.+N.NO3+RP+pH, family = binomial)
lmc<-glm(cbind(int1,X16S)~Temp. + O2 + Cond..20..C + TP +TOC,data=chem, family = binomial)
summary(lmc)

#nb<-glm(chem$X~m5$TOC+chem$Temp.+chem$N.NO3+chem$RP+chem$pH, family = binomial)
Lmc<-lm(asin(sqrt(int1_16S))~TN + O2. + Cond..20..C + TP +TOC, data=chem)
Lmc<-lm(log(int1_16S+1)~TN + O2. + Cond..20..C + TP +TOC, data=chem)

summary(Lmc)
performance::model_performance(Lmc)
performance::check_model(Lmc)
plot(chem$TOC,log(chem$int1_16S))
plot(chem$O2.,chem$int1_16S)
abline(lm(chem$int1_16S~chem$O2.))

Lmp<-lm(asin(sqrt(int_fito))~TN + O2. + Cond..20..C + TP +TOC, data=chem)
summary(Lmp)
performance::check_model(Lmp)
plot(chem$TOC,asin(sqrt(chem$int1_16S)))
abline(lm(chem$int1_16S~chem$O2.))
#colNewG<-c("lightblue3", "aquamarine3", "grey")
#O2<-
colNewC<-c("darkblue", "grey50")

O2<-ggplot(data=chem, aes(x = factor(month), y = O2.,group = as.factor(Year), color = as.factor(Year))) +
  geom_point(size=3, alpha=0.4)+
  scale_color_manual(values=colNewC)+
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9),legend.background = element_rect(fill="grey97",linetype="solid", colour ="white"), axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
  labs(y="O2 % saturation", x="month", color="year")
  #facet_grid(.~fraction,scales="free")
TOC<-ggplot(data=chem, aes(x = factor(month), y = TOC, group = as.factor(Year),color = as.factor(Year))) +
  geom_point(size=3, alpha=0.4) +
  scale_color_manual(values=colNewC)+
  theme_minimal()+
  theme(legend.position = "none", axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
  labs(y="TOC ug L-1", x="", color="year")

c1<-cowplot::plot_grid(TOC,O2, ncol=1,labels = c("A","B"))

fr0_2<-as.data.frame(cbind(chem$int1_16S, chem$TOC))
fr0_2$fraction<-"0.2"
fr1_2<-as.data.frame(cbind(chem$int_fito, chem$TOC))
fr1_2$fraction<-"1.2"
chemfr<-rbind(fr0_2,fr1_2)
cTOC<-ggplot(data=chemfr, aes(y = V1, x = V2, color=fraction)) +
  geom_point(size=3, alpha=0.4) +
  scale_color_manual(values=colNewG)+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = c(0.9, 0.8),legend.background = element_rect(fill="grey97",linetype="solid", colour ="white") ,axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
  labs(x="TOC ug L-1", y="intI1 16S", color=NULL)+
  geom_smooth(method = "lm", se = TRUE, alpha=0.15)

cO2<-ggplot(data=chem, aes(y = int1_16S, x = O2., color="aquamarine3")) +
  geom_point(size=3, alpha=0.4) +
  scale_color_manual(values=colNewG)+
  scale_y_log10()+
  theme_minimal()+
  theme(axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown(), legend.position = "none")+
  labs(x="O2 %", y="intI1 16S", color=NULL)+
  geom_smooth(method = "lm", se = TRUE, alpha=0.15)


c2<-cowplot::plot_grid(cTOC,cO2, nrow=1, labels=c("C","D") )
cowplot::plot_grid(c1,c2, ncol=1, rel_heights = c(3,2))


cond<-ggplot(data=chem, aes(x = factor(month), y = Cond..20..C,group = as.factor(Year), color = as.factor(Year))) +
  geom_point(size=3, alpha=0.4)+
  scale_color_manual(values=colNewC)+
  theme_minimal()+
  theme(legend.position = c(0.9, 0.9),legend.background = element_rect(fill="grey97",linetype="solid", colour ="white"), axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
  labs(y="uS", x="month", color="year")
#facet_grid(.~fraction,scales="free")
TP<-ggplot(data=chem, aes(x = factor(month), y = TP, group = as.factor(Year),color = as.factor(Year))) +
  geom_point(size=3, alpha=0.4) +
  scale_color_manual(values=colNewC)+
  theme_minimal()+
  theme(legend.position = "none", axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
  labs(y="TP ug L-1", x="", color="year")
#facet_grid(.~fraction,scales="free")
TN<-ggplot(data=chem, aes(x = factor(month), y = TN, group = as.factor(Year),color = as.factor(Year))) +
  geom_point(size=3, alpha=0.4) +
  scale_color_manual(values=colNewC)+
  theme_minimal()+
  theme(legend.position = "none", axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
  labs(y="TN ug L-1", x="", color="year")

c3<-cowplot::plot_grid(cond,TP,TN, ncol=1,labels = c("E","F","G"))
cowplot::plot_grid(c1,c2,ncol=1, rel_heights = c(2,1))

