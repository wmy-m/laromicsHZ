library(dplyr)
## NA --------------------------------------------------------------------------------------------------------------
 df <- read.xlsx("./NA_full_pheno.xlsx")
  df <- df %>%
    filter(!SampleID %in% c("LarGlaMVZ172654", "LarOccWymMVZ172421","LarGlaMVZ172495")) # removed bad samples


## Making GEMMA input
gemma=as.data.frame(df[,c("SampleID","sex","locality_short","locality_ID",
                        "BillHueAngle1","BillHueAngle2","MiddleToe_cm",
                        "OrbitalRingHueAngle", "TailLength_cm","Tarsus",
                        "IrisHueAngle","IrisPigHueAngle","TarsusHueAngle",
                        "IrisPig","Weight_g","WingLength_cm","WingSpan_cm",
                        "P10_cm","P9_cm","P8_cm","P7_cm","P6_cm","P5_cm",
                        "HeadLength_cm","BillLength_cm","LPN_cm","BIW_cm","GAPE_cm",
                        "BID_cm","BIDA_cm","BIDP_cm","Mantle","Tip","OrbitalRingHueCat","ER1")])

gemma2 <- gemma %>%
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

gemma3=as.data.frame(gemma2[,c("SampleID","sex","locality_short","locality_ID",
                          "BillHueAngle1","BillHueAngle2","MiddleToe_R",
                         "OrbitalRingHueAngle", "IrisHueAngle","IrisPigHueAngle",
                         "TarsusHueAngle", "TailLength_R","Tarsus_R",
                         "IrisPig","Weight_R","WingLength_R",
                         "P10_R","P9_R","P8_R","P7_R","P6_R","P5_R",
                         "BillLength_R","BillShape_R","Bill1_R","Bill2_R","Bill4_R",
                         "Bill3_R","Bill5_R","Mantle","Tip","OrbitalRingHueCat","ER1")])
gemma4=gemma3[gemma3$locality_short %in% c("COO20", "DES13", "ESI17", "GRA14", "TAT11", "YAQ19"), ]
write.xlsx(gemma4, "GEMMA_NA_INPUT_2.xlsx")

## Making HZAR Pheno input
write.xlsx(gemma3, "NA_HZAR_Pheno.xlsx")

## EU --------------------------------------------------------------------------------------------------------------
rm(list = ls())

df <- read.xlsx("./EU_full_pheno.xlsx")
df <- df %>%
  filter(!LaromicsID %in% "YP9ZZ")# removed bad sample
df2=as.data.frame(df[,c("Tarsus","WingLength_mm","HeadLength_mm","BillLength_mm","BillGonys_mm","BillMin_mm","Weight_g",
                        "NB","NM","NW","T10","T9","T5","T4","OrbitalRingHueAngle","IrisPigment",
                        "W10_mm","W9_mm","B10_mm","B9_mm","B8_mm","B7_mm","B6_mm","B5_mm","B4_mm")])

## Making GEMMA input
gemma=as.data.frame(df[,c("SampleID","Sex","Locality_short","locality_ID","Tarsus","WingLength_mm","HeadLength_mm","BillLength_mm","BillGonys_mm","BillMin_mm","Weight_g",
                          "NB","NM","NW","T10","T9","T5","T4","OrbitalRingHueAngle","IrisPigment",
                          "W10_mm","W9_mm","B10_mm","B9_mm","B8_mm","B7_mm","B6_mm","B5_mm","B4_mm")])
gemma <- gemma %>%
  mutate(
    Weight_g = as.numeric(Weight_g)
  )

gemma2 <- gemma %>%
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

gemma3=as.data.frame(gemma2[,c("SampleID","Sex","Locality_short","locality_ID",
                               "NB","T10","T9","T5","T4","OrbitalRingHueAngle","IrisPigment",
                               "Tar_R","W10_R","W9_R","B10_R","B9_R","B8_R","B7_R","B4_mm",
                               "B6_R","P10C_R","P9C_R","BillLength_R","BillShape_R","Bill6_R")])
gemma4=gemma3[gemma3$locality_ID %in% c("14", "8", "11", "13", 
                                           "10", "5","4","12","6","9"), ]
write.xlsx(gemma4, "GEMMA_EU_INPUT.xlsx")

## Making HZAR Pheno input
write.xlsx(gemma3, "EU_HZAR_Pheno.xlsx")
