setwd("C:/Users/Dell/Desktop/CNR/Daphnia/16Sdaphnia_gliulia/16Sdaphnia")



qPCRcip<-read.csv2("elab_cipais_new.csv")

y_title<-expression(paste(italic("int"), "I1 per 16S rRNA gene copy"))
names(labels1)<-c("fraction >0.2um","fraction >1.2um", "fraction >50um")
#season<-
colNewG<-c("black", "darkblue", "grey50")
labely<-as_labeller(c(">0.2um"="community",">1.2um"="particle-associated",">50um"="zooplankton associated"))
colNew2<-c("darkblue","gray5","deepskyblue1","gray65","lightblue1","gray85","orange")   
library("car")
library(emmeans)
lmC<-lm(int1_16s~fraction*place*factor(year)  ,data=qPCRcip)
Anova(lmC)
emmeans(lmC, pairwise ~ fraction)


lmC<-lm(int1_01~fraction*place*factor(year)  ,data=qPCRcip)
Anova(lmC)
emmeans(lmC, pairwise ~ fraction)



dots<- ggplot(data=qPCRcip, aes(x = factor(month), y = int1_16s, color = as.factor(year))) +
  geom_point(size=3, alpha=0.4) +
  scale_color_manual(values=colNewG)+
  scale_y_log10()+
  theme_minimal()+
  theme(axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
  labs(y="*int*I1 per 16S rRNA gene copy", x="month", color="year")+
  facet_grid(.~fraction,scales="free", labeller=labely)

box<-ggplot(data=qPCRcip, aes(x=fraction, y=int1_16s, fill=place))+
  geom_boxplot(alpha=0.2)+
  theme_minimal()+
  scale_fill_manual(values=colNew2)+
  geom_jitter(alpha=0.5, width = 0.1)+
  theme(axis.title.y = ggtext::element_markdown())+
  theme(legend.position = c(0.75, 0.2),legend.background = element_rect(fill="grey97",linetype="solid", colour ="white") ,axis.title.y = ggtext::element_markdown())+
   labs(y="*int*I1 per 16S rRNA gene copy", x=NULL, )+
  scale_x_discrete(position = "top", labels=c("community","particle-associated", "zooplankton-associated" )) +
  scale_y_log10()
box
cowplot::plot_grid(box, dots,rel_widths = c(3,5.5), rel_heights = c(1), align = "h", axis = "bt", labels = c("A","B"))

aggregate(qPCRcip$int1_01, by=list(qPCRcip$fraction), FUN=sum)
table(qPCRcip$fraction)
table(qPCRcip$int1_16s>0, by=list(qPCRcip$fraction))
library ("dplyr")
lapply(qPCRcip, count>0)
my_summary <- qPCRcip %>%
  count(fraction, int1_16s>0, sort = FALSE) 

setwd("C:/Users/Dell/Desktop/CNR/Daphnia/16Sdaphnia_gliulia/16Sdaphnia")
fw<-read.csv2("food_web.csv")
library("lattice")

cor(fw[,c(3,7,8)])
#plot(aov1)
library("lme4")
library("performance")
library("ggplot2")
lm1<-lm(log(copiRea+1)~as.factor(flg)+as.factor(dap)+as.factor(rot)+as.factor(flg):as.factor(dap)+as.factor(rot):as.factor(flg)+as.factor(rot):as.factor(dap), data=fw)
glm1<-glm(copiRea~flg+dap+rot+flg:dap+rot:flg+rot:dap, family=binomial, data=fw)
plot(lm1)
summary(lm1)
check_model(glm1)
check_model(lm1)

plot(as.factor(fw$FOODWEB), fw$copiRea)
colNew2<-c("darkblue","gray5","orange","gray65","lightblue1","gray85","orange")   

names2<-c("Bacteria", "Bacteria & Flagellate", "Bacteria, Flagellate & Rotifer", "Bacteria & Rotifer")
ggplot(data=fw, aes(x=factor(Daphnia), y = copiRea, ))+
  geom_boxplot(alpha=0.3)+
  theme_minimal()+
  scale_fill_manual(values="orange")+
  geom_jitter(aes(color=factor(Treatment)), size=4, alpha=0.5, width=0.15)+
  theme(legend.position = c(0.75, 0.8),legend.background = element_rect(fill="grey97",linetype="solid", colour ="white") ,axis.title.y = ggtext::element_markdown())+
  scale_color_manual(values=colNew2, labels=names2)+
  labs(y="*int*I1 per 16S rRNA gene copy", x=NULL, color="treatment" )+
  scale_x_discrete(position = "bottom", labels=c("no daphnids","+ daphnids" )) 
colNew2<-c("darkblue","gray5","orange","gray65","lightblue1","gray85","orange3", "lightblue4")   

ggplot(data=fw, aes(x=factor(Daphnia), y = cells_bact, ))+
  geom_boxplot(alpha=0.3)+
  theme_minimal()+
  scale_fill_manual(values="orange")+
  geom_jitter(aes(color=factor(Treatment)), size=4, alpha=0.5, width=0.15)+
  theme(legend.position = c(0.75, 0.8),legend.background = element_rect(fill="grey97",linetype="solid", colour ="white") ,axis.title.y = ggtext::element_markdown())+
  scale_color_manual(values=colNew2, labels=names2)+
  labs(y="*int*I1 per 16S rRNA gene copy", x=NULL, color="treatment" )+
  scale_x_discrete(position = "bottom", labels=c("no Daphnia","+ Daphnia" )) 

Bact<-ggplot(data=fw, aes(x=factor(FOODWEB), y = cells_bact))+
  geom_point(alpha=0.8, size=3)+
  theme_minimal()+
  labs(y="cells mL-1", x=NULL, color="treatment" )+
  guides(color=FALSE)+
  scale_x_discrete(position = "bottom") 

Flag<-ggplot(data=fw, aes(x=factor(FOODWEB), y = cells_Flag))+
  geom_point(alpha=0.8, size=3)+
  theme_minimal()+
  labs(y="cells mL-1", x=NULL, color="treatment" )+
  guides(color=FALSE)+
  scale_x_discrete(position = "bottom") 

cowplot::plot_grid(Bact, Flag, ncol=1,labels = c("A","B"))

qPCR<-read.csv2("elab_exp.csv")
qPCR$numero.daphnie<-factor(qPCR$numero.daphnie, levels=c("T0", "0", "3", "6", "9", "12", "15"))
qPCR<-qPCR[1:18,]

LWWW<-ggplot(data=qPCR, aes(x = nr, y = int1.16s)) +
  geom_point(size=2, alpha=0.55) +
  xlab("# Daphnia")+
  ylab("intI1 / 16S rRNA gene copy")+
  theme_minimal()+
  theme(axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10), axis.title.y = ggtext::element_markdown())+
  labs(y="*int*I1 per 16S rRNA gene copy", x="Number of daphnids")+
  geom_smooth(method='lm', color="grey5", alpha=0.2)
lml<-lm(qPCR$int1.16s~qPCR$nr)
summary(lml)

qPCR_lake<-read.csv2("daphnia_lago.csv")
LW<-ggplot(data=qPCR_lake, aes(x = dap, y = intI1_16S/10)) +
  geom_point(size=2,alpha=0.55) +
    ylab("intI1 / 16S rRNA gene copy")+
  theme_minimal()+
  theme(axis.text=element_text(size=10),plot.title = element_text(face="bold",hjust =0.015, vjust = -10),axis.title.y = ggtext::element_markdown())+
  labs(y="*int*I1 per 16S rRNA gene copy", x="")+
  geom_smooth(method='lm', color="grey5", alpha=0.2)

lml<-lm(qPCR_lake$intI1_16S/10~qPCR_lake$dap)
summary(lml)

cowplot::plot_grid(LW, LWWW, ncol=1,labels = c("A","B"))

