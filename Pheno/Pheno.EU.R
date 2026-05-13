library(missMDA)
library(corrplot)
library(dplyr)
library(beeswarm)
library(openxlsx)
library(MASS)

# Prepare data ---------------------------------------------------------------------------------
df <- read.xlsx("./EU_full_pheno.xlsx")
df <- df %>%
  filter(!LaromicsID %in% "YP9ZZ")# removed bad sample
df2=as.data.frame(df[,c("Tarsus","WingLength_mm","HeadLength_mm","BillLength_mm","BillGonys_mm","BillMin_mm","Weight_g",
                        "NB","NM","NW","T10","T9","T5","T4","OrbitalRingHueAngle","IrisPigment",
                        "W10_mm","W9_mm","B10_mm","B9_mm","B8_mm","B7_mm","B6_mm","B5_mm","B4_mm")])

# Used to check with trait is use to standardised, whichever trait is more linear is chosen
#plot(df$BillLength_mm,df$BillGonys_mm)
#text(df$BillLength_mm,df$BillGonys_mm, labels = df$LaromicsID, pos = 3)

df2 <- df2 %>%
  mutate(
    Weight_g = as.numeric(Weight_g)
  )

df3 <- df2 %>%
  mutate(
    Tar_R=ifelse(is.na(Tarsus) | is.na(WingLength_mm), NA, Tarsus / WingLength_mm),
    W10_R=ifelse(is.na(W10_mm) | is.na(WingLength_mm), NA, W10_mm / WingLength_mm),
    W9_R=ifelse(is.na(W9_mm) | is.na(WingLength_mm), NA, W9_mm / WingLength_mm),
    B10_R=ifelse(is.na(B10_mm) | is.na(WingLength_mm), NA, B10_mm / WingLength_mm),
    B9_R=ifelse(is.na(B9_mm) | is.na(WingLength_mm), NA, B9_mm / WingLength_mm),
    B8_R=ifelse(is.na(B8_mm) | is.na(WingLength_mm), NA, B8_mm / WingLength_mm),
    B7_R=ifelse(is.na(B7_mm) | is.na(WingLength_mm), NA, B7_mm / WingLength_mm),
    B6_R=ifelse(is.na(B6_mm) | is.na(WingLength_mm), NA, B6_mm / WingLength_mm),
    B5_R=ifelse(is.na(B5_mm) | is.na(WingLength_mm), NA, W10_mm / WingLength_mm),
    B4_R=ifelse(is.na(B4_mm) | is.na(WingLength_mm), NA, B4_mm / WingLength_mm),
    P10C_R=ifelse(is.na(W10_mm) | is.na(B10_mm), NA, W10_mm / B10_mm),
    P9C_R=ifelse(is.na(W9_mm) | is.na(B9_mm), NA, W9_mm / B9_mm),
    BillLength_R=ifelse(is.na(HeadLength_mm) | is.na(BillLength_mm), NA, HeadLength_mm / BillLength_mm),
    WT_R=ifelse(is.na(Weight_g) | is.na(WingLength_mm), NA, Weight_g / WingLength_mm),
    BillShape_R=ifelse(is.na(BillLength_mm) | is.na(BillGonys_mm), NA, BillGonys_mm / BillLength_mm),
    Bill6_R=ifelse(is.na(BillLength_mm) | is.na(BillMin_mm), NA, BillMin_mm / BillLength_mm)
  )

pca=as.data.frame(df3[,c("NB","NM","NW","T10","T9","T5","T4","OrbitalRingHueAngle","IrisPigment",
                         "Tar_R","W10_R","W9_R","B10_R","B9_R","B8_R","B7_R",
                         "B6_R","B5_R","B4_R","P10C_R","P9C_R","BillLength_R","WT_R","BillShape_R","Bill6_R")])

# remove rows/columns that have high missingness
colMeans(is.na(pca)) #remove any traits with more than 50% missingness
pca2=as.data.frame(df3[,c("NB","T10","T9","T5","T4","OrbitalRingHueAngle","IrisPigment",
                          "Tar_R","W10_R","W9_R","B10_R","B9_R","B8_R","B7_R","B4_mm",
                          "B6_R","P10C_R","P9C_R","BillLength_R","BillShape_R","Bill6_R")])

df.clean <- df[rowMeans(is.na(pca2)) <= 0.4, ]
pca.clean <- pca2[rowMeans(is.na(pca2)) <= 0.4, ]

# impute missing data
pca.clean <- as.data.frame(lapply(pca.clean, as.numeric))
pca.clean[!is.finite(as.matrix(pca.clean))] <- NA
nb <- estim_ncpPCA(pca.clean)
pca.imputed <- imputePCA(pca.clean, ncp = nb$ncp)$completeObs
pca <- as.data.frame(pca.imputed)

df.pca <- prcomp(pca, center = TRUE,scale. = TRUE)
summary(df.pca)
#0.2675 0.1966

df.clean$PC1 <- df.pca$x[, 1]
df.clean$PC2 <- df.pca$x[, 2]

# LDA ---------------------------------------------------------------------------------  
{
  for(i in 1:dim(pca)[2]){
    pca[,i]=c(scale(as.numeric(pca[,i])))
  }
  
# Use same scaled / centered data subset as PCA

lddat=pca

# Conservative: argentatus = populations 1-3, cachninnans = 14-15
train.pops <- c(1:3,14:16)
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

write.table(df.clean,"./EU.LDA.txt")

# plot ---------------------------------------------------------------------------------
cols <- ifelse(df.clean$locality_ID >= 1 & df.clean$locality_ID <= 3, "#56B4E9",
               ifelse(df.clean$locality_ID >= 14 & df.clean$locality_ID <= 16, "#E69F00",
                      "gray"))
group_colors <- c("#56B4E9", "gray","#E69F00") 

dfgroup=(1+1*(df.clean$locality_ID>3)+1*(df.clean$locality_ID>13))
df.sample <- read.xlsx("./Sampling.xlsx")
df.clean$Sampling <- df.sample$Sampling[
  match(df.clean$SampleID, df.sample$SampleID)
]
id_match <- df.clean$SampleID %in% df.sample$SampleID

#female_idx.EU <- df.clean$Sex == "F"
#male_idx.EU   <- df.clean$Sex == "M"
xlim.EU <- range(df.pca$x[,1])
ylim.EU <- range(df.pca$x[,2])


border_col_match <- ifelse(df.clean$Sampling[id_match] == "Random", "black",
                           "red")
border_full <- rep(NA, nrow(df.clean))  # 608 long
border_full[id_match] <- border_col_match


pdf("./New/EU.Pheno_20260414_LDA_scaled_2.pdf",width=5,height=5) #Figure 1 plot
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


pdf("./New/EU.Pheno_202604.pdf") #Supplement
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
  segments(3.5, par("usr")[3], 3.5, par("usr")[4], col="red", lty=2)
  segments(13.5, par("usr")[3], 13.5, par("usr")[4], col="red", lty=2)
  mtext("Hybrid zone",1,at=6,col=c("red"),line=-1,cex=0.8)
  
  plot(df.pca$x[,1], df.pca$x[,2],
       col = cols,
       pch = 16, xlim = xlim.EU, ylim = ylim.EU,
       xlab = "PC1 (26.7%)",
       ylab = "PC2 (19.7%)")
  arrows(0, 0, 
         df.pca$rotation[, 1] * diff(par("usr")[1:2]) * 0.8,
         df.pca$rotation[, 2] * diff(par("usr")[3:4]) * 0.8,
         col = "#00000099", length = 0.1)
  text(df.pca$rotation[, 1] * diff(par("usr")[1:2]) * 0.8 * 1.05, 
       df.pca$rotation[, 2] * diff(par("usr")[3:4]) * 0.8 * 1.05,
       rownames(df.pca$rotation), col = "#00000099",cex=0.7)
  
  dev.off()
}

