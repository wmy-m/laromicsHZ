library(bgchm)
library(tidyr)
library(dplyr)
library(data.table)

setwd("./8_Cline/TraitRegions/BGChm/manyloci")

load("./NA.combinedClines.TR.rda")

df <- data.frame(
  locus = name_est,
  chr_class = "autosome",
  center = c_est[,1],
  center_low = c_est[,2],
  center_high = c_est[,3],
  gradient = v_est[,1],
  gradient_low = v_est[,2],
  gradient_high = v_est[,3]
)

df2 <- df %>%
  separate(locus, into = c("chr", "pos"), sep = "_") %>%
  mutate(pos = sub("\\..*$", "", pos))

traits <- read.table("../../TraitRegions.txt", header = FALSE)
colnames(traits)=c("trait","chr","start","end")

setDT(df2)
setDT(traits)

for (t in traits$trait) {
  df2[, (t) := FALSE]
}

# Fill TRUE where SNP is in the trait region
for (i in seq_len(nrow(traits))) {
  df2[
    chr == traits$chr[i] &
      pos >= traits$start[i] &
      pos <= traits$end[i],
    (traits$trait[i]) := TRUE
  ]
}
df2[, bg := !Reduce(`|`, .SD), .SDcols = unique(traits$trait)]


#sz
sz_out <- sum2zero(df2$center, df2$gradient, transform = TRUE)
sz_out_low <- sum2zero(df2$center_low, df2$gradient_low, transform = TRUE)

df2$sz_center   <- sz_out$center      # centered centers
df2$sz_gradient <- sz_out$gradient    # centered gradients
df2$sz_center_low   <- sz_out_low$center      # centered centers
df2$sz_gradient_low <- sz_out_low$gradient    # centered gradients
df2$steep       <- sz_out_low$gradient > 1  # TRUE if lower CI > 1

write.csv(df2, "NA_bgc_TraitRegions.csv", row.names = FALSE)

# calculating proportion
trait_cols <- c(unique(traits$trait), "bg")

# Calculate proportions for each trait
props <- rbindlist(lapply(trait_cols, function(tr) {
  
  subset <- df2[get(tr) == TRUE]
  
  counts <- table(subset$steep)
  props <- prop.table(counts)
  
  data.table(
    trait = tr,
    steep_FALSE = props["FALSE"],
    steep_TRUE_count  = table(subset$steep)["TRUE"],
    steep_TRUE_props  = props["TRUE"]
  )
}), fill = TRUE)
write.csv(props, "NA_bgc_TraitRegions_props.csv", row.names = FALSE)

# plot
pdf("NA.bgchm.corr.pdf", width=5, height=5)
par(
  mfrow = c(3, 3),
  mar = c(3, 3, 1, 1),   # inner margins
  oma = c(1, 1, 1, 1),   # outer margins
  mgp = c(1.5, 0.4, 0),  # axis title, labels, line
  tcl = -0.2,            # shorter ticks
  xaxs = "i",
  yaxs = "i"
)

plot(NULL, ylab="Cline gradient", xlab="Cline center",
#     xlim=c(0, 1), ylim=c(0, 2.5))
     xlim=c(0, 1), ylim=c(0, 3.7))
points(df2[steep == "FALSE" & NAorb1 == "TRUE"]$sz_center, 
       df2[steep == "FALSE" & NAorb1 == "TRUE"]$sz_gradient, col="#00640080", pch=16, cex=0.4)
points(df2[steep == "TRUE" & NAorb1 == "TRUE"]$sz_center, 
       df2[steep == "TRUE" & NAorb1 == "TRUE"]$sz_gradient, col="#00640080", pch=16, cex=0.6)

plot(NULL, ylab="Cline gradient", xlab="Cline center",
     xlim=c(0, 1), ylim=c(0, 3.7))
points(df2[steep == "FALSE" & NAorb2 == "TRUE"]$sz_center, 
       df2[steep == "FALSE" & NAorb2 == "TRUE"]$sz_gradient, col="#00640080", pch=16, cex=0.4)
points(df2[steep == "TRUE" & NAorb2 == "TRUE"]$sz_center, 
       df2[steep == "TRUE" & NAorb2 == "TRUE"]$sz_gradient, col="#00640080", pch=16, cex=0.6)

plot(NULL, ylab="Cline gradient", xlab="Cline center",
     xlim=c(0, 1), ylim=c(0, 3.7))
points(df2[steep == "FALSE" & NAorb3 == "TRUE"]$sz_center, 
       df2[steep == "FALSE" & NAorb3 == "TRUE"]$sz_gradient, col="#00640080", pch=16, cex=0.4)
points(df2[steep == "TRUE" & NAorb3 == "TRUE"]$sz_center, 
       df2[steep == "TRUE" & NAorb3 == "TRUE"]$sz_gradient, col="#00640080", pch=16, cex=0.6)

plot(NULL, ylab="Cline gradient", xlab="Cline center",
     xlim=c(0, 1), ylim=c(0, 3.7))
points(df2[steep == "FALSE" & EUorb1 == "TRUE"]$sz_center, 
       df2[steep == "FALSE" & EUorb1 == "TRUE"]$sz_gradient, col="#FF8C0080", pch=16, cex=0.4)
points(df2[steep == "TRUE" & EUorb1 == "TRUE"]$sz_center, 
       df2[steep == "TRUE" & EUorb1 == "TRUE"]$sz_gradient, col="#FF8C0080", pch=16, cex=0.6)

plot(NULL, ylab="Cline gradient", xlab="Cline center",
     xlim=c(0, 1), ylim=c(0, 3.7))
points(df2[steep == "FALSE" & EUorb2 == "TRUE"]$sz_center, 
       df2[steep == "FALSE" & EUorb2 == "TRUE"]$sz_gradient, col="#FF8C0080", pch=16, cex=0.4)
points(df2[steep == "TRUE" & EUorb2 == "TRUE"]$sz_center, 
       df2[steep == "TRUE" & EUorb2 == "TRUE"]$sz_gradient, col="#FF8C0080", pch=16, cex=0.6)

plot(NULL, ylab="Cline gradient", xlab="Cline center",
     xlim=c(0, 1), ylim=c(0, 3.7))
points(df2[steep == "FALSE" & EUorb3 == "TRUE"]$sz_center, 
       df2[steep == "FALSE" & EUorb3 == "TRUE"]$sz_gradient, col="#FF8C0080", pch=16, cex=0.4)
points(df2[steep == "TRUE" & EUorb3 == "TRUE"]$sz_center, 
       df2[steep == "TRUE" & EUorb3 == "TRUE"]$sz_gradient, col="#FF8C0080", pch=16, cex=0.6)

plot(NULL, ylab="Cline gradient", xlab="Cline center",
     xlim=c(0, 1), ylim=c(0, 3.7))
points(df2[steep == "FALSE" & NAtip1 == "TRUE"]$sz_center, 
       df2[steep == "FALSE" & NAtip1 == "TRUE"]$sz_gradient, col="#EF65B180", pch=16, cex=0.4)
points(df2[steep == "TRUE" & NAtip1 == "TRUE"]$sz_center, 
       df2[steep == "TRUE" & NAtip1 == "TRUE"]$sz_gradient, col="#EF65B180", pch=16, cex=0.6)

plot(NULL, ylab="Cline gradient", xlab="Cline center",
     xlim=c(0, 1), ylim=c(0, 3.7))
points(df2[steep == "FALSE" & bg == "TRUE"]$sz_center, 
       df2[steep == "FALSE" & bg == "TRUE"]$sz_gradient, col="#DEDEDE", pch=16, cex=0.4)
points(df2[steep == "TRUE" & bg == "TRUE"]$sz_center, 
       df2[steep == "TRUE" & bg == "TRUE"]$sz_gradient, col="#DEDEDE", pch=16, cex=0.6)

dev.off()

pdf("NA.bgchm.cline.pdf", width=5, height=5)
par(
  mfrow = c(3, 3),
  mar = c(3, 3, 1, 1),   # inner margins
  oma = c(1, 1, 1, 1),   # outer margins
  mgp = c(1.5, 0.4, 0),  # axis title, labels, line
  tcl = -0.2,            # shorter ticks
  xaxs = "i",
  yaxs = "i"
)

gencline_plot(center=df2[steep == "FALSE" & NAorb1 == "TRUE"]$sz_center,v=df2[steep == "FALSE" & NAorb1 == "TRUE"]$sz_gradient,col="#DEDEDE",pdf=FALSE)
par(new = TRUE)
gencline_plot(center=df2[steep == "TRUE" & NAorb1 == "TRUE"]$sz_center,v=df2[steep == "TRUE" & NAorb1 == "TRUE"]$sz_gradient,col="darkgreen",pdf=FALSE)

gencline_plot(center=df2[steep == "FALSE" & NAorb2 == "TRUE"]$sz_center,v=df2[steep == "FALSE" & NAorb2 == "TRUE"]$sz_gradient,col="#DEDEDE",pdf=FALSE)
par(new = TRUE)
gencline_plot(center=df2[steep == "TRUE" & NAorb2 == "TRUE"]$sz_center,v=df2[steep == "TRUE" & NAorb2 == "TRUE"]$sz_gradient,col="darkgreen",pdf=FALSE)

gencline_plot(center=df2[steep == "FALSE" & NAorb3 == "TRUE"]$sz_center,v=df2[steep == "FALSE" & NAorb3 == "TRUE"]$sz_gradient,col="#DEDEDE",pdf=FALSE)
par(new = TRUE)
gencline_plot(center=df2[steep == "TRUE" & NAorb3 == "TRUE"]$sz_center,v=df2[steep == "TRUE" & NAorb3 == "TRUE"]$sz_gradient,col="darkgreen",pdf=FALSE)

gencline_plot(center=df2[steep == "FALSE" & EUorb1 == "TRUE"]$sz_center,v=df2[steep == "FALSE" & EUorb1 == "TRUE"]$sz_gradient,col="#DEDEDE",pdf=FALSE)
par(new = TRUE)
gencline_plot(center=df2[steep == "TRUE" & EUorb1 == "TRUE"]$sz_center,v=df2[steep == "TRUE" & EUorb1 == "TRUE"]$sz_gradient,col="darkorange",pdf=FALSE)

gencline_plot(center=df2[steep == "FALSE" & EUorb2 == "TRUE"]$sz_center,v=df2[steep == "FALSE" & EUorb2 == "TRUE"]$sz_gradient,col="#DEDEDE",pdf=FALSE)
par(new = TRUE)
gencline_plot(center=df2[steep == "TRUE" & EUorb2 == "TRUE"]$sz_center,v=df2[steep == "TRUE" & EUorb2 == "TRUE"]$sz_gradient,col="darkorange",pdf=FALSE)

gencline_plot(center=df2[steep == "FALSE" & EUorb3 == "TRUE"]$sz_center,v=df2[steep == "FALSE" & EUorb3 == "TRUE"]$sz_gradient,col="#DEDEDE",pdf=FALSE)
par(new = TRUE)
gencline_plot(center=df2[steep == "TRUE" & EUorb3 == "TRUE"]$sz_center,v=df2[steep == "TRUE" & EUorb3 == "TRUE"]$sz_gradient,col="darkorange",pdf=FALSE)

gencline_plot(center=df2[steep == "FALSE" & NAtip1 == "TRUE"]$sz_center,v=df2[steep == "FALSE" & NAtip1 == "TRUE"]$sz_gradient,col="#DEDEDE",pdf=FALSE)
par(new = TRUE)
gencline_plot(center=df2[steep == "TRUE" & NAtip1 == "TRUE"]$sz_center,v=df2[steep == "TRUE" & NAtip1 == "TRUE"]$sz_gradient,col="#EF65B1",pdf=FALSE)

gencline_plot(center=df2[steep == "FALSE" & bg == "TRUE"]$sz_center,v=df2[steep == "FALSE" & bg == "TRUE"]$sz_gradient,col="#DEDEDE",pdf=FALSE)
par(new = TRUE)
gencline_plot(center=df2[steep == "TRUE" & bg == "TRUE"]$sz_center,v=df2[steep == "TRUE" & bg == "TRUE"]$sz_gradient,col="darkgrey",pdf=FALSE)
dev.off()

