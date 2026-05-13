library(missMDA)
library(corrplot)
library(dplyr)
library(beeswarm)
library(openxlsx)
library(MASS)

# Prepare data ---------------------------------------------------------------------------------
{
  df <- read.xlsx("./NA_full_pheno.xlsx")
  df <- df %>%
    filter(!SampleID %in% c("LarGlaMVZ172654", "LarOccWymMVZ172421","LarGlaMVZ172495")) # removed bad samples
  df2=as.data.frame(df[,c("BillHueAngle1","BillHueAngle2","MiddleToe_cm",
                          "OrbitalRingHueCat", "TailLength_cm","Tarsus",
                          "IrisHueAngle","IrisPigHueAngle","TarsusHueAngle",
                          "IrisPig","Weight_g","WingLength_cm","WingSpan_cm",
                          "P10_cm","P9_cm","P8_cm","P7_cm","P6_cm","P5_cm",
                          "HeadLength_cm","BillLength_cm","LPN_cm","BIW_cm","GAPE_cm",
                          "BID_cm","BIDA_cm","BIDP_cm","Mantle","Tip")])
  
  df3 <- df2 %>%
    mutate(
      MiddleToe_R=ifelse(is.na(MiddleToe_cm) | is.na(Tarsus), NA, MiddleToe_cm / Tarsus),
      TailLength_R=ifelse(is.na(TailLength_cm) | is.na(WingLength_cm), NA, TailLength_cm / WingLength_cm),
      Tarsus_R=ifelse(is.na(WingLength_cm) | is.na(Tarsus), NA, Tarsus / WingLength_cm),
      Weight_R=ifelse(is.na(Weight_g) | is.na(WingLength_cm), NA, Weight_g / WingLength_cm),
      WingLength_R=ifelse(is.na(WingSpan_cm) | is.na(WingLength_cm), NA, WingLength_cm / WingSpan_cm),
      P10_R=ifelse(is.na(P10_cm) | is.na(WingLength_cm), NA, P10_cm / WingLength_cm),
      P9_R=ifelse(is.na(P9_cm) | is.na(WingLength_cm), NA, P9_cm / WingLength_cm),
      P8_R=ifelse(is.na(P8_cm) | is.na(WingLength_cm), NA, P8_cm / WingLength_cm),
      P7_R=ifelse(is.na(P7_cm) | is.na(WingLength_cm), NA, P7_cm / WingLength_cm),
      P6_R=ifelse(is.na(P6_cm) | is.na(WingLength_cm), NA, P6_cm / WingLength_cm),
      P5_R=ifelse(is.na(P5_cm) | is.na(WingLength_cm), NA, P5_cm / WingLength_cm),
      BillLength_R=ifelse(is.na(HeadLength_cm) | is.na(BillLength_cm), NA, HeadLength_cm / BillLength_cm),
      BillShape_R=ifelse(is.na(BillLength_cm) | is.na(BID_cm), NA, BID_cm / BillLength_cm),
      Bill1_R=ifelse(is.na(BillLength_cm) | is.na(LPN_cm), NA, LPN_cm / BillLength_cm),
      Bill2_R=ifelse(is.na(BillLength_cm) | is.na(BIW_cm), NA, BIW_cm / BillLength_cm),
      Bill3_R=ifelse(is.na(BillLength_cm) | is.na(GAPE_cm), NA, GAPE_cm / BillLength_cm),
      Bill4_R=ifelse(is.na(BillLength_cm) | is.na(BIDA_cm), NA, BIDA_cm / BillLength_cm),
      Bill5_R=ifelse(is.na(BillLength_cm) | is.na(BIDP_cm), NA, BIDP_cm / BillLength_cm)
    )
  
  # Used to check with trait is use to standardised, whichever trait is more linear is chosen
  #plot(df2$WT,df2$WIS)
  #text(df2$WT,df2$WIL, labels = df$SampleID, pos = 3)
  
  pca=as.data.frame(df3[,c("BillHueAngle1","BillHueAngle2","MiddleToe_R",
                          "OrbitalRingHueCat", "IrisHueAngle","IrisPigHueAngle",
                          "TarsusHueAngle", "TailLength_R","Tarsus_R",
                          "IrisPig","Weight_R","WingLength_R",
                          "P10_R","P9_R","P8_R","P7_R","P6_R","P5_R",
                          "BillLength_R","BillShape_R","Bill1_R","Bill2_R","Bill4_R",
                          "Bill3_R","Bill5_R","Mantle","Tip")])
  
  # remove rows/columns that have high missingness
  colMeans(is.na(pca))
  df.clean <- df[rowMeans(is.na(pca)) <= 0.4, ]
  pca.clean <- pca[rowMeans(is.na(pca)) <= 0.4, ]
  
  # impute missing data
  nb <- estim_ncpPCA(pca.clean)
  pca.imputed <- imputePCA(pca.clean, ncp = nb$ncp)$completeObs
  pca <- as.data.frame(pca.imputed)
  
  df.pca <- prcomp(pca, center = TRUE,scale. = TRUE)
  summary(df.pca)
  #0.2343 0.1220
  
  df.clean$PC1 <- df.pca$x[, 1]
  df.clean$PC2 <- df.pca$x[, 2]
}  
  
# LDA ---------------------------------------------------------------------------------  
{  
  for(i in 1:dim(pca)[2]){
    pca[,i]=c(scale(as.numeric(pca[,i])))
  }
    # Use same scaled / centered data subset as PCA
  lddat=pca
  
  # Conservative: glaucescens = populations 1-3, occidentalis = 21-33
  train.pops <- c(1:3,21:33)
  train.data <- lddat[df.clean$locality_ID %in% train.pops, ]
  test.data <- lddat[!(df.clean$locality_ID %in% train.pops), ]
  train.data$species <- 1+1*(df.clean$locality_ID[df.clean$locality_ID %in% train.pops]<4)
  
  # Compute LDA
  ld1=lda(species~.,data=train.data)
  
  # Make predictions
  predictions <- ld1 %>% predict(test.data)
  
  # Add LD1 to data frame
  df.clean$LD1<- NA
  df.clean$LD1[df.clean$locality_ID %in% train.pops] <- predict(ld1)$x
  df.clean$LD1[!(df.clean$locality_ID %in% train.pops)] <- predictions$x
  }  

write.table(df.clean,"./NA.LDA.txt")

# plot ---------------------------------------------------------------------------------
cols <- ifelse(df.clean$locality_ID >= 1 & df.clean$locality_ID <= 3, "#009E73",
               ifelse(df.clean$locality_ID >= 21 & df.clean$locality_ID <= 33, "#CC79A7",
                      "gray"))
group_colors <- c("#009E73", "gray","#CC79A7")  

dfgroup=(1+1*(df.clean$locality_ID>3)+1*(df.clean$locality_ID>20))
df.sample <- read.xlsx("./Sampling.xlsx")
df.clean$Sampling <- df.sample$Sampling[
  match(df.clean$SampleID, df.sample$SampleID)
]
id_match <- df.clean$SampleID %in% df.sample$SampleID

#female_idx.NA <- df.clean$SEX == "0"
#male_idx.NA   <- df.clean$SEX == "1"
xlim.NA <- range(df.pca$x[,1])
ylim.NA <- range(df.pca$x[,2])

border_col_match <- ifelse(df.clean$Sampling[id_match] == "Random", "black",
                           "red")

border_full <- rep(NA, nrow(df.clean))  # 608 long
border_full[id_match] <- border_col_match


pdf("./New/NA.Pheno_20260414_LDA_scaled_2.pdf",width=5,height=5) #LDA plot
{
par(mar=c(4,3,2,1), oma=c(0,0,0,2), mgp=c(2,0.7,0))
layout(matrix(c(
  1, 2), nrow = 1, byrow = TRUE),
widths = c(3, 1))
set.seed(1)
y <- df.clean$LD1
bs <- beeswarm(y ~ dfgroup, method = "compactswarm", pwcol = group_colors[dfgroup],
               cex = 0.8,pch = 1)
x <- bs$x
y <- bs$y
ord <- order(dfgroup)

cols <- group_colors[dfgroup[ord]]
match_ord <- id_match[ord]

points(bs$x[match_ord],
       bs$y[match_ord],
       pch = 16,
       col = cols[match_ord])

b.vals=sort(setNames(as.numeric(ld1$scaling), rownames(ld1$scaling)),
            decreasing = TRUE)
b=barplot(b.vals,horiz = F)
text(b, b.vals, labels = names(b.vals), srt = 90, cex = 0.6, adj = 0)
}
dev.off()


pdf("NA.Pheno_202604.pdf") #Supplement; clinal plot + PCA
{
par(mfrow=c(2,1), mar=c(4,3,2,1), oma=c(0,0,0,4), mgp=c(2,0.7,0)) 
plot(df.clean$locality_ID, df.clean$LD1,type = "n",xlab="",ylab="LD1",xlim=c(1,max(df.clean$locality_ID,na.rm=T)))
points(df.clean$locality_ID, df.clean$LD1, col=cols, pch=16)
points(df.clean$locality_ID[id_match],
       df.clean$LD1[id_match],
       pch = 21,
       col = border_col_match,  # border color
       bg  = cols[id_match],        # fill color
       lwd = 1.5)
segments(3.5, par("usr")[3], 3.5, par("usr")[4], col="orange", lty=2)
segments(20.4, par("usr")[3], 20.4, par("usr")[4], col="orange", lty=2)
segments(11.5, par("usr")[3], 11.5, par("usr")[4], col="red", lty=2)
segments(20.5, par("usr")[3], 20.5, par("usr")[4], col="red", lty=2)
mtext("Introgression zone",1,at=7.5,col=c("orange"),line=-1,cex=0.8)
mtext("Hybrid zone",1,at=16,col=c("red"),line=-1,cex=0.8)

plot(df.pca$x[,1], df.pca$x[,2],
     col = cols,
     pch = 16, xlim = xlim.NA, ylim = ylim.NA,
     xlab = "PC1 (23.5%)",
     ylab = "PC2 (12.1%)")
arrows(0, 0, 
       df.pca$rotation[, 1] * diff(par("usr")[1:2]) * 0.8,
       df.pca$rotation[, 2] * diff(par("usr")[3:4]) * 0.8,
       col = "#00000099", length = 0.1)
text(df.pca$rotation[, 1] * diff(par("usr")[1:2]) * 0.8 * 1.05, 
     df.pca$rotation[, 2] * diff(par("usr")[3:4]) * 0.8 * 1.05,
     rownames(df.pca$rotation), col = "#00000099",cex=0.7)

dev.off()
}


