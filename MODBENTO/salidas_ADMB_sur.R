rm(list=ls(all=TRUE))  # Borra todas los objetos creados
setwd("C:/Users/mauricio.mardones/Documents/IFOP/Eval_Stock/ERIZO/2019/XSUR/MODBENTO")            # Prepare the environment
data.dir <- "C:/Users/mauricio.mardones/Documents/IFOP/Eval_Stock/ERIZO/2019/XSUR/MODBENTO"
source("C:/Users/mauricio.mardones/Documents/IFOP/Eval_Stock/ERIZO/2019/XSUR/MODBENTO/read.admb.R")
library("ggplot2")

#============================================================================================
system("tpl2cpp -bounds Modbento");system("adcomp -s Modbento");system("adlink -s Modbento")
system("Modbento")

#system("Modbento -mcmc 1000 -mcsave 100")
#system("Modbento -mceval")
#============================================================================================
dir()
datos.X <-read.rep('Modbento.rep')

# like<-matrix(ncol=6,nrow=5)
# for(i in 1:6){
#   datos.X1 <-read.rep(paste('s',i,'.rep',sep=""))
#   like[,i]<-datos.X1$CPUE_Desemb
# }


names(datos.X)
datos.X

year<-datos.X$Años
year

CPUE_obs<-datos.X$CPUE_obs[1,] 
CPUE_est<-datos.X$CPUE_obs[2,] 
length(CPUE_obs)
length(CPUE_est)


a1<-(CPUE_obs==0)
CPUE_obs[a1]<-NA

Desemb_obs<-datos.X$Desemb_obs[1,]
Desemb_est<-datos.X$Desemb_obs[2,]

Lmed_obs<-datos.X$Lmed_obs[1,]
Lmed_est<-datos.X$Lmed_obs[2,]
length(Lmed_obs)
a2<-(Lmed_obs==0)
Lmed_obs[a2]<-NA

Biomasa_total<-datos.X$Biomasa_total

Biomasa_explotable<-datos.X$Biomasa_explotable

Biomasa_desovante<-datos.X$Biomasa_desovante
BD<-Biomasa_desovante
BD

R_est<-datos.X$Reclutamiento[1,]
R_pred<-datos.X$Reclutamiento[2,]
R_est
R_pred

Mortalidad_pesca<-datos.X$F

Reducción<-datos.X$Reducción
Reducción

Selectividad<-datos.X$Selectividad
Selectividad
# Selectividad_1<-datos.X$Selectividad[1,]; Selectividad_1
# Selectividad_2<-datos.X$Selectividad[50,]; Selectividad_2
# Selectividad_3<-datos.X$Selectividad[54,]; Selectividad_3


#################### GRAFICOS ############################################
# CPUE
x11()
par(mfrow=c(2,2))
#qplot(year, datos.X$Desemb_obs[1,]) +
#  geom_line()

plot(year, CPUE_est,xaxp=c(year[1],year[length(year)]+1,length(year)), type="l", col="black", axes=T, ann=F, lwd=2, lty=1,las=1)
points(year, CPUE_obs, type ="p", lwd=1)
title(main="CPUE", col.main="black", font.main=2)
title(xlab="Años", col.lab="black")
title(ylab="Kg/hora", col.lab="black")
legend("topright", c("CPUE_obs", "CPUE_est"), lty=c(NA,1), pch=c(1,NA), bty="n", cex=0.8)


# DESEMBARQUES
#x11()
plot(year,Desemb_obs,xaxp=c(year[1],year[length(year)]+1,length(year)), type ="l", col = "black", axes=T, ann=F, lwd=2,las=1)
points(year,Desemb_est, type ="p") #Volvemos a dibujar camiones.
title(main="Desembarques(ton)", col.main="black", font.main=2) #título
title(xlab="Años", col.lab="black")
title(ylab="Toneladas", col.lab="black")
legend("topright", c("Desemb_obs", "Desemb_est"), lty=c(1,NA), pch=c(NA,1), bty="n", cex=0.8)

# BIOMASAS
#x11()
plot(year,Biomasa_explotable, xaxp=c(year[1],year[length(year)]+1,length(year)/2),type ="l", col="black", axes=T, ann=F, lwd=2, ylim=c(0,100000),las=1)
points(year, Biomasa_total, type ="l", lty=1, col="blue", lwd=2)
points(year, BD, type ="l", lty=1, col ="green", lwd=2)
title(main="Biomasas", col.main="black", font.main=2)
title(xlab="Años", col.lab="black")
title(ylab="Toneladas", col.lab="black")
legend("topright", c("Biom_expl", "Biom_tot", "Biom_Desov"), lty=c(1,1,1), col=c("black","blue", "green"), bty="n", cex=0.8)


dev.off()

############################ OTRA VENTANA DE GRAFICOS#####################
# RECLUTAMIENTOS
#x11()
par(mfrow=c(2,2))
plot(year,R_est,xaxp=c(year[1],year[length(year)]+1,length(year)), type ="l", col = "black", axes=T, ann=F, lwd=2, ylim=c(0,400),las=1)
points(year,R_pred, type ="l", lty=2, lwd=1, pch=1) #Volvemos a dibujar camiones.
title(main="Reclutamiento", col.main="black", font.main=2) #título
title(xlab="Años", col.lab="black")
title(ylab="Individuos", col.lab="black")
legend("bottomleft", c("R_est", "R_pred"), lty=c(1,3), bty="n", cex=0.8)


# TALLAS  MEDIAS
plot(year,Lmed_obs, xaxp=c(year[1],year[length(year)]+1,length(year)), type ="l", col = "black", axes=T, ann=F, lwd=2,las=1, ylim=c(40, 100))
points(year,Lmed_est, type ="p", lty=2) #Volvemos a dibujar camiones.
title(main="Lmed", col.main="black", font.main=2) #título
title(xlab="Años", col.lab="black")
title(ylab="cm", col.lab="black")
legend("bottomleft", c("Lmed_est", "Lmed_obs"), lty=c(1,NA), pch=c(NA,1), bty="n", cex=0.8)

# MORTALIDAD POR PESCA
#x11()
plot(year, Mortalidad_pesca, xaxp=c(year[1],year[length(year)]+1,length(year)*1),type ="l", col="black", axes=T, ann=F, lwd=2,las=1)
#title(main="Mortalidad por pesca", col.main="black", font.main=2) #título
title(xlab="Años", col.lab="black")
title(ylab="Mortalidad por pesca (F)", col.lab="black")
abline(h=0.19,col="red", lwd = "3")
legend("topleft", "F40 = 0.19", lty=1, bty="n", cex=1.2, col = "red", lwd = "2")

# REDUCCIOND DE BIOMASA
#x11()
plot(year, Reducción, xaxp=c(year[1],year[length(year)]+1,length(year)*1),type ="l", col="black", axes=T, ann=F, las=1,lwd=2)
#title(main="Reducción de BD", col.main="black", font.main=2) #título
title(xlab="Años", col.lab="black")
title(ylab="BD/BDo", col.lab="black")
abline(h = 0.4, col = "red", lwd = "3") #PBR??
legend("topright", "BD/BDo 40%", lty=1, bty="n", cex=1.2, col = "red", lwd = "2")

# otra ventana de graficos
# SELECTIVIDAD 

x11()
edad<-1:12
plot(edad, Selectividad_1, type ="l",col = "black", axes=T, ann=F,las=1)
points(edad,Selectividad_2, type ="l", lty=1, col="blue", lwd=2) #Volvemos a dibujar camiones.
points(edad, Selectividad_3, type = "l", lty=1, col="red")
#title(main="Selectividad", col.main="black", font.main=2) #título
title(xlab="Edad", col.lab="black")
title(ylab="S", col.lab="black")
legend("bottomright", c("Sel_60_01", "Sel_02_09", "Sel_10_13"), horiz = TRUE, lty=c(1,1,1), col=c("black","blue", "red"), bty="n", cex=1)



########################################################################
year<-datos.X$Años
nyears <- length(datos.X$Años)  

#Residuos                                                                       
Res_Cpue   <-log(CPUE_obs)-log(CPUE_est)                                             
Res_Des   <-log(Desemb_obs)-log(Desemb_est)     
Res_Lmed   <-log(Lmed_obs)-log(Lmed_est) 
Res_Recl   <- log(R_est)-log(R_pred)


cvcpue<-rep(0.2,nyears)
cvdes<-rep(0.05,nyears)
cvlmed <- rep(0.02,nyears)
cvrecl <- rep(0.05,nyears)
#cvevadir <- rep(0.20, nyears) #cuando tenga otro indice de abundancia

obsCpue95i <- CPUE_obs*exp(-1.96*cvcpue);obsCpue95s <-CPUE_obs*exp(1.96*cvcpue)
obsD95i <- Desemb_obs*exp(-1.96*cvdes);obsD95s <- Desemb_obs *exp(1.96*cvdes)
obsLmed95i <- Lmed_obs*exp(-1.96*cvlmed);obsLmed95s <- Lmed_obs*exp(1.96*cvlmed)
obsRecl95i <- R_est*exp(-1.96*cvrecl);obsRecl95s <- R_est*exp(1.96*cvrecl)


#==================================================================================
# I FIGURAS
#==================================================================================
# AJUSTE CPUE
x11()
par(mfrow=c(2,2))
par(mar=c(4,4,1,1)+0.5)
plot(year,CPUE_est,type="l",cex.axis=1,lwd=2,xaxp=c(year[1],year[length(year)]+1,length(year)*1),
     ylim=c(0,max(CPUE_est)+2), xlim = c(1960, max (year)+1),xaxs= "i",yaxs= "i",ylab="Indice Relativo CPUE",las=1,xlab="Años",cex.lab=1.2)
arrows(x0=year,y0=obsCpue95i,x1=year,y1=obsCpue95s,length=0.05,angle=90,lty=1,code=3)
points(year,CPUE_obs,cex=1.2,bg="black",pch=21)

#AJUSTE DESEMBARQUES

#AJUSTE LONGITUD MEDIA
x11()
plot(year,Desemb_est,type="l",cex.axis=1,lwd=2,xaxp=c(year[1],year[length(year)]+1,length(year)*1),
     ylim=c(0,max(Desemb_est)+1000),xlim = c(1960, max (year)+1),xaxs= "i",yaxs= "i",ylab="Desembarques (t.)",las=1,xlab="Años",cex.lab=1.2)
arrows(x0=year,y0=obsD95i,x1=year,y1=obsD95s,length=0.05,angle=90,lty=1,code=3)
points(year,Desemb_obs ,cex=1.2,bg="black",pch=21, lwd=2)
x11()
par(mar=c(4,4,1,1)+0.5)
plot(year,Lmed_est ,type="l",cex.axis=1,lwd=2,xaxp=c(year[1],year[length(year)]+5,(length(year)+4)*1),
     ylim=c(60,max(Lmed_est)+5),xlim = c(1960, max (year)+1), xaxs= "i",yaxs= "i",ylab="Longitud media (cm)",las=1,xlab="Años",cex.lab=1.2)
arrows(x0=year,y0=obsLmed95i,x1=year,y1=obsLmed95s ,length=0.05,angle=90,col=4,lty=1,code=3)
points(year,Lmed_obs,cex=1.2,bg="blue",pch=21,col=4)
legend("topright", c("L_Pred", "L_Obs"), horiz = TRUE, lty=c(1,1), lwd = 2, col = c("black","blue"), bty="n", cex=1.2)


#AJUSTE RECLUTAMIENTO
x11()
par(mar=c(4,4,1,1)+0.5)
plot(year,R_est,type="l",cex.axis=1,lwd=2,xaxp=c(year[1],year[length(year)]+5,(length(year)+4)*1),
     ylim=c(0,max(R_est)+350),xlim = c(1960, max (year)),xaxs= "i",yaxs= "i",ylab="Reclutamiento (t)",las=1,xlab="Años",cex.lab=1.2)
#arrows(x0=year,y0=obsRecl95i,x1=year,y1=obsRecl95s,length=0.05,angle=90,col=4,lty=1,code=3)
lines(year,R_pred,cex=1.2,bg="blue",pch=21,col=4, lwd=2)
legend("topright", c("R_est", "R_pred"), horiz = TRUE, lty=c(1,1), lwd = 2, col = c("black","blue"), bty="n", cex=1.2)



#==================================================================================
# RESIDUALES CPUE
x11()
par(mfrow=c(2,2),mar=c(2,4,2,1)+0.5)
plot(year,Res_Cpue,cex.axis=0.8,type="h",main="CPUE",ylab="Residuales (escala log)",xlab="")
abline(h=0,col="darkgray")
plot(log(CPUE_est),Res_Cpue, main="Residuales vs ajustado",ylab="Residuales",xlab="Valor ajustado")
abline(h=0,col="darkgray")
hist(Res_Cpue,xlab="Residuales",ylab="Frecuencia",main="Histograma de Residuos")
qqnorm(Res_Cpue); qqline(Res_Cpue, col = 2)

# RESIDUALES DESEMBARQUES
x11()
par(mfrow=c(2,2),mar=c(2,4,2,1)+0.5)
plot(year,Res_Des,cex.axis=0.8,type="h",main="Desembarques",ylab="Residuales (escala log)",xlab="")
abline(h=0,col="darkgray")
plot(log(Desemb_est),Res_Des, main="Residuales vs ajustado",ylab="Residuales",xlab="Valor ajustado")
abline(h=0,col="darkgray")
hist(Res_Des,xlab="Residuales",ylab="Frecuencia",main="Histograma de Residuos")
qqnorm(Res_Des); qqline(Res_Des, col = 2)

# RESIDUALES LONGITUD MEDIA
x11()
par(mfrow=c(2,2),mar=c(2,4,2,1)+0.5)
plot(year,Res_Lmed,cex.axis=0.8,type="h",main="Longitud media (cm)",ylab="Residuales (escala log)",xlab="")
abline(h=0,col="darkgray")
plot(log(Lmed_est),Res_Lmed, main="Residuales vs ajustado",ylab="Residuales",xlab="Valor ajustado")
abline(h=0,col="darkgray")
hist(Res_Lmed,xlab="Residuales",ylab="Frecuencia",main="Histograma de Residuos")
qqnorm(Res_Lmed); qqline(Res_Lmed, col = 2)

#RESIDUALES RECLUTAMIENTO

x11()
par(mfrow=c(2,2),mar=c(2,4,2,1)+0.5)
plot(year,Res_Recl ,cex.axis=0.8,type="h",main="Reclutamiento",ylab="Residuales (escala log)",xlab="")
abline(h=0,col="darkgray")
plot(log(Lmed_est),Res_Recl , main="Residuales vs ajustado",ylab="Residuales",xlab="Valor ajustado")
abline(h=0,col="darkgray")
hist(Res_Recl ,xlab="Residuales",ylab="Frecuencia",main="Histograma de Residuos")
qqnorm(Res_Recl ); qqline(Res_Recl , col = 2)


#============================================================#
# II. COMPOSICIÓN EDAD DE LAS CAPTURAS                       #
#============================================================#
names(datos.X)
year<-datos.X$Años
nyears <- length(datos.X$Años)    
tallas <- seq(40, 146,2)
ntalla <- length(tallas)
#Proporción observada                                        
pobsF<-datos.X$pobs                                         
#Proporción predicha                                         
ppredF<-datos.X$ppred

par(mfrow=c(6,4),mar=c(1,1,1,1)+1)
# par(mfrow = c(6, 4))
# par(cex = 0.6)
# par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
# par(tcl = -0.25)
# par(mgp = c(2, 0.6, 0))
for(i in 37:nyears){
  plot(tallas,pobsF[i,],ylab="",xlab="",cex.lab=2, cex.axis=2,ylim=c(0,0.18),type="h",lwd=3, cex.axis=2, axes = T)
  lines(tallas,ppredF[i,],col=2,lwd=2)
  text(55,0.15,year[i], cex = 2)
  box(col = "grey60")
}
mtext("x axis", side = 1, outer = TRUE, cex = 0.7, line = 2.2, col = "grey20")
mtext("y axis", side = 2, outer = TRUE, cex = 0.7, line = 2.2, col = "grey20")



#--------------------
# GRAFICO RESIDUALES TALLAS 
#---------------------

anos<-year[37:nyears]
obs <-pobsF[37:nyears,]
pre <-ppredF[37:nyears,]  
res <-obs-pre

rng <-range(res,na.rm=T)
dd  <-dim(res)
est <-matrix(NA,nrow=dd[1],ncol=dd[2])

for(j in 1:dd[1]){for(k in 1:dd[2]){val<-res[j,k]
if(val>0){est[j,k]<-val/rng[2]}
else{est[j,k]<-val/rng[1]*-1}}}
#id<-which(years%in%anos==FALSE)
#est <-est[-id,]
#png(paste(getwd(),"/resid_tallaF.png",sep=""),width=350,height=550)   
#x11(width=350,height=550)
par(mar=c(5.4,6.7,2,1),cex.axis=1,cex.lab=1.1)
image(tallas,anos,t(est),col=0,yaxt="n",xlab="",ylab="")
ee  <-dim(est)
for(n in 1:ee[1]){for(m in 1:ee[2]){vol<-est[n,m]
if(is.na(vol)==FALSE){
	if(vol>0){points(tallas[m],anos[n],pch=19,cex=6.82*sqrt(vol),col=1)}
	if(vol<0){points(tallas[m],anos[n],pch=1,cex=6.82*sqrt(vol*-1),col=1)}
}}}
mtext("Tallas",side=1,line=3.2,cex=1.1);posi<-seq(1,57,by=4)
axis(2,at=anos,labels=anos,las=2)
mtext("Años",side=2,line=4.7,cex=1.1)
box()
#dev.off()


#===================================================================================
#============================================================#
# III. Abundabcia a la edad por año                          #
#============================================================#
year<-datos.X$Años
nyears <- length(datos.X$Años)    
     
edades <- seq(1,12)
nedades<- length(edades)
#Abundancia a la edad                                       
Abuned<-datos.X$Abundancia
max(Abuned)


x11()
par(mfrow=c(8,8),mar=c(1,2,1,1)+0.5)
for(i in 1:nyears){
  plot(edades,Abuned[i,],ylab="",xlab="",ylim=c(0,150),type="h",lwd=2, col="18")
  text(8,40,year[i])
}

#===================================================================================
#============================================================#
# IV. Selectividad a la edad por año                          #
#============================================================#

year<-datos.X$Años
nyears <- length(datos.X$Años)    

edades <- seq(1,12)
nedades<- length(edades)
#ASelectividad a la edad                                       
Selectividad<-datos.X$Selectividad
max(Selectividad)


#x11()
par(mfrow=c(8,8),mar=c(1,2,1,1)+0.5)
for(i in 1:nyears){
  plot(edades,Selectividad[i,],ylab="",xlab="",ylim=c(0,1.2),type="l",lwd=2, col="18")
  text(8,0.2,year[i])
}
#==================================================================#
#probabilidad de captura a la talla
#==================================================================#
x11()
tallas <- seq(40, 146,2)

plot(tallas,datos.X$Prob_talla[1,], type ="n", ylab="Probabilidad", xlab="Tallas(mm)")
for(i in 1:12){
  lines(tallas,datos.X$Prob_talla[i,], col=i, lwd=3)
}

#====================================================================================
#falta desbloquearlo del codigo en YP
#============================================================#
# V. Biomasas Proyectadas                                    #
#============================================================#
names(datos.X)
BiomProyF0 <- datos.X$Biomasa_desovante_proyectada_para_cada_mF[1,]
BiomProyF0.5 <- datos.X$Biomasa_desovante_proyectada_para_cada_mF[2,]
BiomProyF1 <- datos.X$Biomasa_desovante_proyectada_para_cada_mF[3,]
BiomProyF1.5 <- datos.X$Biomasa_desovante_proyectada_para_cada_mF[4,]


#x11()
par(mfrow=c(2,1))
plot(c(year, seq(2019, 2028, 1)), c(datos.X$Biomasa_desovante,rep(NA,10)), type="l", lwd = "3", col="1",  
     ylab = "Biomasa Desovante (ton)", xlab = "años proyectados")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Biomasa_desovante,rep(BiomProyF0)), type="l", col="3")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Biomasa_desovante,rep(BiomProyF0.5)), type="l", col="4")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Biomasa_desovante,rep(BiomProyF1)), type="l", col="5")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Biomasa_desovante,rep(BiomProyF1.5)), type="l", col="6")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Biomasa_desovante,rep(NA,10)), lwd ="2", type="l", col="1")
legend(1975, 3000, c("F0", "F0.5", "F1", "F1.5"), lty=c(1,1,1,1,1), col=c("3","4", "5", "6","7"), 
       bty="n", cex=0.8, horiz =TRUE, yjust = 0.5, xjust = 0.5, lwd = "3")

#============================================================#
# V. Capturas Proyectadas                                    #
#============================================================#
##Revisar

names(datos.X)
CaptProyF0 <- datos.X$Capturas[1,]
CaptProyF0.5 <- datos.X$Capturas[2,]
CaptProyF1 <- datos.X$Capturas[3,]
CaptProyF1.5 <- datos.X$Capturas[4,]
#CaptProyF2 <- datos.X$YP_cada_mF[5,]


#x11()
plot(c(year, seq(2019, 2028, 1)), c(datos.X$Desemb_obs[1,],rep(NA,10)), type="l", lwd = "2", col="1", ylim = c(0, 20000), 
     ylab = "Capturas (ton)", xlab = "Años Proyectados")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Desemb_obs[1,],rep(CaptProyF0)), type="l", col="3")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Desemb_obs[1,],rep(CaptProyF0.5)), type="l", col="4")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Desemb_obs[1,],rep(CaptProyF1)), type="l", col="5")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Desemb_obs[1,],rep(CaptProyF1.5)), type="l", col="6")
#lines(c(year, seq(2016, 2025, 1)), c(datos.X$Desemb_obs[1,],rep(CaptProyF2)),type="l", col="7")
lines(c(year, seq(2019, 2028, 1)), c(datos.X$Desemb_obs[1,],rep(NA,10)), type="l", col="1", lwd = "3")
legend("topleft", c("F0", "F0.5", "F1", "F1.5"), lty=c(1,1,1,1,1), col=c("3","4", "5", "6","7"), 
       bty="n", cex=0.8, horiz =TRUE, yjust = 2, xjust = 2, lwd = "3")


#================================================================================#
# BIOMASA TOTAL Y DESOVANTE
#================================================================================#
std     <-read.table("Modbento.std",header=T,sep="",na="NA",fill=T) 
x  <-c(year,rev(year))
x1 <-c(year[1],year[nyears]+1,nyears+1/2) #xaxp
x2 <-c(year[1]-1,year[nyears]+1) #xlim

BT       <- subset(std,name=="BT")$value
BTstd    <- subset(std,name=="BT")$std
SSBt     <- subset(std,name=="BD")$value
SSBtstd  <- subset(std,name=="BD")$std

bt   <- c((BT-1.96*BTstd),rev((BT+1.96*BTstd)))
ssbt <- c((SSBt-1.96*SSBtstd),rev((SSBt+1.96*SSBtstd)))

#png(paste(dir.0,"/Fig14.png",sep=""))
#x11()
par(mfrow=c(2,1),mar=c(0,4,2,1)+1)

plot(x, bt,type="n",ylim=c(0,max(bt)),cex.axis=0.8,xaxs="i",yaxs="i",xlim=x2,xaxp=x1,
	ylab="Biomasa total (t)",las=1,xlab="Año",cex.lab=1.1, mgp=c(4,1,0))
polygon(x,bt, col="gray",border="gray")
lines(year,BT,lwd=2)

plot(x,ssbt,type="n",ylim=c(0,max(ssbt)),cex.axis=0.8,xaxs="i",yaxs="i",xlim=x2,xaxp=x1,                                                               
	ylab="Biomasa desovante (t)",las=1,xlab="Año",cex.lab=1.1, mgp=c(4,1,0))                
polygon(x,ssbt,col="gray", border="gray")
lines(year,SSBt,lwd=2)   
#lines(year,rep(rep$SSBpbr[3],nyears),lwd=2,col=2)
#text(2003,rep$SSBpbr[3]*10^-6+0.1,"BDRMS",cex=1.2)
#dev.off()

#================================================================================#
# RECLUTAMIENTOS Y DESVIOS
#================================================================================#
Reclutas       <- subset(std,name=="Rest")$value
Reclutasstd    <- subset(std,name=="Rest")$std
logdesvRt      <- subset(std,name=="dev_log_Ro")$value
logdesvRtstd   <- subset(std,name=="dev_log_Ro")$std

rt     <- c((Reclutas-1.96*Reclutasstd),rev(Reclutas+1.96*Reclutasstd))
logdrt <- c((logdesvRt-1.96*logdesvRtstd),rev(logdesvRt+1.96*logdesvRtstd))

#png(paste(dir.0,"/Fig15.png",sep=""))
#x11()
par(mfrow=c(2,1),mar=c(2,4,2,1)+0.5)

plot(x,rt , type="n", xaxp=x1,cex.axis=0.8,xaxs= "i",yaxs= "i",
	xlim=x2,ylab="Reclutamientos",las=1,xlab="Año",cex.lab=1.1)
polygon(x, rt , col="gray", border = "gray");lines(year,Reclutas,lwd=2)
abline(h=exp(std$log_Ro+0.5*0.6^2)*10^-6,col=2,lty=2)

plot(x, logdrt, type="n", xaxp=x1,cex.axis=0.8,xaxs= "i",yaxs= "i",
	xlim=x2,ylab="Desvios de los Reclutamientos",las=1,xlab="Año",cex.lab=1.1)
polygon(x, logdrt, col="gray", border = "gray");lines(year,logdesvRt,lwd=2)
abline(h=0,lty=2,col="darkgray")
#dev.off()


#================================================================================#
# MORTALIDAD POR PESCA
#================================================================================#
Ft       <- subset(std,name=="log_F")$value
Ftstd    <- subset(std,name=="log_F")$std

ft  <- c(exp((Ft)-1.96*(Ftstd)),rev(exp((Ft)+1.96*(Ftstd)))) 

#png(paste("/Fig16.png",sep=""),width=600,height=400)
#x11()
par(mfrow=c(1,1),mar=c(2.5,4,1,1)+0.5,oma=c(1,1,2,0))
plot(x, ft, xaxp=x1,cex.axis=0.8,xaxs= "i",yaxs= "i",
	xlim=x2,type="n", ylab="Mortalidad por pesca (F)",las=1,xlab="Año",cex.lab=1.1)
polygon(x, ft, col="gray", border = "gray");lines(year,exp(Ft),lwd=2)
#lines(year,rep(1.0,nyears),lty=2) #mortalidad natural
lines(year,rep(rep$Fs[2],nyears),lwd=2,col=2) #Frms
text(c(2012,2004),c(1.0+0.05,rep$Fs[2]+0.05),c("M","FRMS"),cex=1.2)
abline(h=0.19, col="red", lwd=2)
legend("topleft", "F40 = 0.19", lty=1, bty="n", cex=1.2, col = "red", lwd = "2")
dev.off()
#================================================================================#
# RPR
#================================================================================#
RPRt     <- subset(std,name=="RPR")$value
RPRtstd  <- subset(std,name=="RPR")$std
rprt   <- c((RPRt-1.96*RPRtstd),rev((RPRt+1.96*RPRtstd)))

x11()
par(mfrow=c(1,1),mar=c(2,4,2,1)+0.5)
plot(x,rprt,type="n",ylim=c(0,max(rprt)),cex.axis=0.8,xaxs="i",yaxs="i",xlim=x2,xaxp=x1,
	ylab="RPR",las=1,xlab="Año",cex.lab=1.1)
polygon(x,rprt, col="gray",border="gray")
abline(h = 0.4, col = "yellow", lwd = "3")#PBR??
abline(h= 0.2, col = "red", lwd =3)
lines(year,RPRt,lwd=2)




#==============================================================
# DIAGRAMA DE FASE PARA BENTONICOS
#==============================================================

datos.X <-read.rep('Modbento.rep')
names(datos.X)
datos.X

year<-datos.X$Años
length(year)


RPRlp<-datos.X$Reducción

Fb1<-datos.X$F[1:length(year)]


F40_b1<-0.19774


BD_BD40<-RPRlp/0.4
F_F40 <- c(Fb1/F40_b1)

#===============================================================
# Calculate confidence intervals
#===============================================================
std     <-read.table("Modbento.std",header=T,sep="",na="NA",fill=T) 

SSBt     <- subset(std,name=="BD")$value
SSBtstd  <- subset(std,name=="BD")$std
Ft       <- subset(std,name=="log_F")$value
Ftstd    <- subset(std,name=="log_F")$std

SSBo     <- datos.X$BD_virginal
BD40     <- SSBo*0.4
Fval     <- exp(tail(Ft,1))/F40_b1
lastB    <- tail(SSBt,1)
lastF    <- tail(Fval,1)
Qmult    <- -qnorm((1-(80/100))/2.0)

sbSE     <- tail(SSBtstd,1)
sb95     <- c(lastB-Qmult*sbSE,lastB+Qmult*sbSE)
B95      <- sb95/BD40
FvSE     <- tail(Ftstd,1)
F95      <- c(lastF*exp(-Qmult*FvSE),lastF*exp(Qmult*FvSE))
#===============================================================
# GRAFICO
#===============================================================
cx    <-c(0,0.5,0.5,0,0)        #colapso
cy    <-c(-0.3,-0.3,9,9,-0.3)   #colapso
sex   <-c(0.5,1,1,0.5,0.5)      #sobre-explotación
sey   <-c(-0.3,-0.3,9,9,-0.3)   #sobre-explotación
spx   <-c(1,3,3,1,1)            #sobrepesca
spy   <-c(1,1,9,9,1)            #sobrepesca
subex <-c(1,3,3,1,1)            #subexplotado
subey <-c(-0.3,-0.3,1,1,-0.3)   #subexplotado



#x11()
par(mar=c(4,4,1,1)+0.5)
plot(BD_BD40,F_F40,type="o",xlab="BD/BD40",ylab="F/F40",las=1,yaxp=c(0,9,9),ylim = c(0,8), xlim = c(0.3, 2.8), xaxp=c(-0.5,2.5,2))
polygon(cx,cy,col="tomato")
polygon(sex,sey,col="tan 1")
polygon(spx,spy,col="khaki1")
polygon(subex,subey,col="greenyellow")
arrows(x0=B95[1],y0=lastF,x1=B95[2],y1=lastF,length=0.05,angle=90,col=4,lwd=2,code=3)
arrows(x0=tail(BD_BD40,1),y0=F95[1],x1=tail(BD_BD40,1),y1=F95[2],length=0.05,angle=90,col=4,lwd=2,code=3)
lines(BD_BD40,F_F40,type="o",pch=21,bg="white",cex=1.1)
abline(h=1,v=c(0.5,1),lty=2)
points(c(BD_BD40[1],BD_BD40[length(year)]),c(F_F40[1],F_F40[length(year)]),col=c(1,1),bg=c(3,2),pch=21,cex=1.8)
text(BD_BD40[c(1,35:58)],F_F40[c(1,35:58)]+0.1,year[c(1,35:56)],cex=0.7)
text(1.8,0.5,"SUB-EXPLOTADO", cex = 1.2)
text(1.8,3,"SOBREPESCA", cex = 1.2)
text(0.7,3,"SOBRE-EXPLOTADO", cex = 1.2)
text(0.37,0.5,"COLAPSADO", cex = 1.2)
box()


#mas sencillo

plot(BD_BD40,F_F40,ylab="F/Frms",
     xlab="B/Brms",las=1,ylim=c(0,10),xlim=c(0,4),cex=1.2)
arrows(x0=B95[1],y0=lastF,x1=B95[2],y1=lastF,length=0.05,angle=90,col=4,lwd=2,code=3)
arrows(x0=tail(BD_BD40,1),y0=F95[1],x1=tail(BD_BD40,1),y1=F95[2],length=0.05,angle=90,col=4,lwd=2,code=3)
#lines(BD_BD40,F_F40,type="o",pch=21,bg="white",cex=1.1)
lines(BD_BD40,F_F40,lty=2,lwd=2)
points(c(BD_BD40[1],BD_BD40[length(year)]),c(F_F40[1],F_F40[length(year)]),col=c(1,1),bg=c(3,2),pch=21,cex=1.8)
#points(rpr[length(lc$ano)],frel[length(lc$ano)],pch=19) #aÃ±o mÃ¡s reciente
abline(v=0.5,col="red", lwd=2)
lines(c(1,1),c(0,1),col="green", lwd=2)
lines(c(1,4),c(1,1),col="green", lwd=2)
text(BD_BD40[c(1,35:58)],F_F40[c(1,35:58)]+0.2,year[c(1,35:58)],cex=0.7)
text(1.8,0.5,"SUB-EXPLOTADO", cex = 1.2)
text(1.8,5,"SOBREPESCA", cex = 1.2)
#text(0.7,5,"SOBRE-EXPLOTADO", cex = 1.2)
text(0.2,0.6,"COLAPSADO", cex = 1.2)


#=====================================================================
# CALCULO DE LA CBA CON FRMS=F40proxy
#=====================================================================

YTPp<-std[std$name=="YTPp",,]$value
YTPpstd<-std[std$name=="YTPp",,]$std

q     <- seq(0.1,0.5,0.1)  # niveles de riesgo (cuantiles)                                
nq    <- length(q)                                                                                   
CBA  <- matrix(ncol=nq,nrow=4)

for(i in 1:4){
  for(j in 1:nq){
    CBA[i,j]<-qnorm(q[j],YTPp[i],YTPpstd[i])}}
round(CBA,1)

########################################################################
# ESTIMACIÓN DE CURVA SPR
########################################################################

source("D:/Mauricio/2017_xnorte/X_NOR/MODBENTO/Fn_pbr.R")

#datos de entrada
Dat<-list()
Dat$M		    <- 0.25                  # M se asume constate para edades y años
Dat$Tspw	  <- 0.91             # corresponde a la época de desove en sardina es 2/12 (año biológico)
Dat$Mad	    <- datos.X$msex_edad      # Paso 1:  Fmediana serie histórica 1990-2008
Dat$Wmed	  <- datos.X$Wmed_edad     # vector del promedio del peso medio a la edad (suma las columnas) !!! revisar esto!!!
Dat$Sel     <- datos.X$Selectividad[nyears,] 

Fmort 	    <- seq(0,3.5,0.02)           # vector de mortalidad por pesca para la curva de YPR y SPR
nf          <- length(Fmort) 
R0 		      <- 1  
Amax        <- 12

SPRcurv 		   <- SPRFmort(R0,Fmort,Amax,Dat) 

#---------------------------------------------------------------------------------------------------------------------------
x11()
par(mfrow=c(2,1),mar=c(4,4.5,2,1)+0.5)

plot(seq(1,Amax,1),datos.X$msex_edad,type="l", ylim=c(0,1.05),ylab="Madurez y selectividad",xlab="Edad (años)",las=1,col=3,lwd=2)
lines(seq(1,Amax),datos.X$Selectividad[nyears,],col=4,lwd=2)
legend(4.5,max(datos.X$Selectividad[nyears,])-0.5,c("Madurez", "Selectividad"),
       col=c(3,4),bty="n", lwd=c(2,2),lty=c(1,1),cex=1)

plot(SPRcurv[,1],SPRcurv[,4],type="l", ylab="%BDPR",xlab="Mortalidad por pesca (F)",lwd=2,las=1,col=4)
abline(h=0.4,v=0.26,col=5,lty=1)
text(0.9,0.62,"FRMS = 0.26", cex=0.8,col=4)
#---------------------------------------------------------------------------------------------------------------------------





