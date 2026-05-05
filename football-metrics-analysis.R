# 0. Libraries and Data --------------------------------------------------------
library(dplyr)
library(sbgcop)
library(xtable)
library(readr)
library(plyr)
library(ggplot2)
library(Rgbp)

#Data 
gameLogs <- readr::read_csv()

allSeasonsData <- readr::read_csv()

# 1. Preliminary calculations for discrimination and stability -----------------

## 1.1. Bootstrap --------------------------------------------------------------

### Aggregators

basic_aggr <- list(
  # totals
  MP     = function(df) sum(df$MP,     na.rm = TRUE),
  Gls    = function(df) sum(df$Gls,    na.rm = TRUE),
  Ast    = function(df) sum(df$Ast,    na.rm = TRUE),
  xG     = function(df) sum(df$xG,     na.rm = TRUE),
  npxG   = function(df) sum(df$npxG,   na.rm = TRUE),
  SCA    = function(df) sum(df$SCA,    na.rm = TRUE),
  GCA    = function(df) sum(df$GCA,    na.rm = TRUE),
  Prg_P  = function(df) sum(df$Prg_P,  na.rm = TRUE),
  Prg_C  = function(df) sum(df$Prg_C,  na.rm = TRUE),
  xA     = function(df) sum(df$xA,     na.rm = TRUE),
  xAG    = function(df) sum(df$xAG,    na.rm = TRUE),
  KP     = function(df) sum(df$KP,     na.rm = TRUE),
  PPA    = function(df) sum(df$PPA,    na.rm = TRUE),
  Prg_R  = function(df) sum(df$Prg_R,  na.rm = TRUE),
  DefSum = function(df) sum(df$DefSum, na.rm = TRUE),
  
  # ratios / percentages 
  PassesPerc        = function(df) mean(df$PassesPerc,        na.rm = TRUE),
  TakeOns_SuccPerc  = function(df) mean(df$TakeOns_SuccPerc,  na.rm = TRUE),
  AerDuel_Perc      = function(df) mean(df$AerDuel_Perc,      na.rm = TRUE),
  Gls_Sh            = function(df) mean(df$Gls_Sh,            na.rm = TRUE),
  SoT_Sh            = function(df) mean(df$SoT_Sh,            na.rm = TRUE)
)


## 1.2.  Definition of the function --------------------------------------------

bootstrap_foot_season <- function(gameLogs, identity = FALSE, aggr_funs = basic_aggr) {
  
  # Extract teams
  teams <- unique(gameLogs$Tm) 
  
  # Resample dates for each team
  resampledDates <- list()
  for (tm in teams) {
    dates <- unique(gameLogs$Date[gameLogs$Tm == tm])
    resampledDates[[tm]] <- if (identity) dates else sample(dates, replace = TRUE)
  }
  
  # create the bootstrapped season
  resampled_season <- do.call(rbind, lapply(teams, function(tm) {
    team_data <- gameLogs[gameLogs$Tm == tm, ]
    do.call(rbind, lapply(resampledDates[[tm]], function(d) team_data[team_data$Date == d, ]))
  }))
  
  # Aggregate stats for each player
  player_stats <- ddply(resampled_season, .(Player), function(df_player) {
    sapply(aggr_funs, function(f) f(df_player))
  })
  
  return(player_stats)
}

## 1.3 Verification ------------------------------------------------------------

gameLogs_sub <- gameLogs[1:200, ] # Subsample of 200 rows to keep things lightweight

# Run the function without bootstrap (identity = TRUE)
reconstructed <- bootstrap_foot_season(gameLogs_sub, identity = TRUE)

# Use the same results as the reference value
season_totals <- reconstructed

# Compare the stats
common_vars <- setdiff(colnames(reconstructed), "Player")

# Visualization loop (scatterplot for each metric)
for (metric in common_vars) {
  plot(reconstructed[[metric]],
       season_totals[[metric]],
       xlab = paste("Reconstructed", metric),
       ylab = paste("Original", metric),
       main = metric)
  abline(a = 0, b = 1, col = "red")
  browser()
}


## 1.4. Bootstrap with 20 iterations -------------------------------------------

# Select players who have played at least 50 minutes
minMP <- 50
minPlayers <- allSeasonsData$Player[allSeasonsData$MP > minMP]
totals <- allSeasonsData[allSeasonsData$Player %in% minPlayers, ]



# Loop
nboot <- 20
gdfSubset <- gameLogs[gameLogs$Player %in% minPlayers, ]

repsList <- lapply(1:nboot, function(x) {
  print(sprintf("boostrapping %i of %i", x, nboot))
  bootstrap_foot_season(gdfSubset, identity = FALSE)
})
reps <- do.call(rbind, repsList) 
  

# Save reps
save(reps, file=file.path(dat.dir, "footballReps.Rdata"))


# 2. Calculation of stability and discrimination -------------------------------

## 2.1. Calculation of variances -----------------------------------------------

getVariances <- function(totals, reps, seasons="23_24", exclude=c(), minSeason=4) {
  
  exclude <- c(exclude, "Player", "Tm")
  commonMetrics <- setdiff(intersect(colnames(totals), colnames(reps)),
                           exclude)
  
  ## Total variance
  TV <- apply(totals, 2, function(w) var(w, na.rm=TRUE))[commonMetrics]
  
  ## Single season variance for calculating discrimination
  SV <- apply(totals[totals$Season %in% seasons, ], 2,
              function(w) var(w, na.rm=TRUE))[commonMetrics]
  
  bootstrapVars <- ddply(reps, .(Player), function(x) {
    if(nrow(x) < 10) {
      vec <- c()
    } else{ 
      vec <- apply(x, 2, function(w) var(w, na.rm=TRUE))
      vec["Player"] <- x$Player[1]
    }
    vec
  })
  
  ## Average bootstrap vars within a season
  BV <- colMeans(data.matrix(bootstrapVars[, commonMetrics]), na.rm=TRUE)
  
  withinPlayerVar <- ddply(totals, .(Player), function(x) {
    if(nrow(x) >= minSeason) {
      vec <- apply(x, 2, function(w) var(w, na.rm=TRUE))
    } else{
      vec <- rep(NA, ncol(x))
      names(vec) <- colnames(x)
    }
    vec["Player"] <- x$Player[1]
    vec
  })
  WV <- colMeans(data.matrix(withinPlayerVar[, commonMetrics]), na.rm=TRUE)
  
  intersectPlayers <- intersect(withinPlayerVar$Player, bootstrapVars$Player)
  commonCols <- intersect(colnames(withinPlayerVar), colnames(bootstrapVars))
  commonCols <- setdiff(commonCols, c("Player", "Year", "Tm"))
  diffVars <-
    data.matrix(withinPlayerVar[withinPlayerVar$Player %in% intersectPlayers, commonCols]) -
    data.matrix(bootstrapVars[bootstrapVars$Player %in% intersectPlayers, commonCols])
  DV <- colMeans(diffVars, na.rm=TRUE)[commonMetrics]
  
  discriminationScores <- 1 - BV / SV
  
  stabilityScores <- 1 - (WV - BV) / (TV - BV)
  
  list(TV=TV, BV=BV, WV=BV, SV=SV,
       discrimination=discriminationScores,
       stability=stabilityScores)
  
}

varList <- getVariances(allSeasonsData, reps, minSeason = 1, seasons = "23_24")

disc <- varList$discrimination
stab <- varList$stability


## 2.2. Discrimination vs. Stability plot --------------------------------------

pdf(file.path(fig.dir, "reliability.pdf"))

plot(disc, stab, cex=0,
     xlab="Discrimination", ylab="Stability",
     xlim=c(0.3, 1.0), ylim=c(0.3, 1.0),
     main="Metric Reliabilities", cex.axis=1.2, cex.main=1.5, cex.lab=1.2)

text(disc, stab, labels=names(stab), cex=0)

dev.off()




# 3. Independance --------------------------------------------------------------

## 3.1. Data formatting --------------------------------------------------------
row_names <- paste(allSeasonsData$League, allSeasonsData$Player, seq_len(nrow(allSeasonsData)), sep = " - ")
allSeasonsData_clean <- allSeasonsData[, !(names(allSeasonsData) %in% c("League", "Player", "Season"))]
rownames(allSeasonsData_clean) <- row_names 
allSeasonsData <- allSeasonsData_clean
rm(allSeasonsData_clean)

## 3.2. Gaussian copulas -------------------------------------------------------
res <- sbgcop::sbgcop.mcmc(allSeasonsData) # Gaussian copula 

C <- apply(res$C.psamp, c(1,2), mean) # Mean matrix

imputedTab <- res$Y.pmean

save(res, C, imputedTab,
     file=file.path(dat.dir, "fball-sbg.RData"))


## 3.3. Dendrogram -------------------------------------------------------------
pdf(file.path(fig.dir, "fbballDendro.pdf"), width=14)
plot(hclust(dist(abs(C))), main="", xlab="", axes=F, sub="", ylab="Dependencies between fotball Metrics", cex=0.75, cex.lab=1.5)
dev.off()


## 3.4. PCA --------------------------------------------------------------------

### Eigen decomposition 
eig <- eigen(C)
vals <- eig$values
vecs <- eig$vectors
head(cumsum(vals)/sum(vals), n=20)


### Variance explained figure
pdf(file.path(fig.dir, "fbballVarExplained.pdf"), width=4, height=4)
plot(cumsum(vals)/sum(vals),ylim=c(0,1),pch=19,ylab="Variance Explained",xlab="Number of Components", main="All football Metrics", type="l", lwd=3, cex.lab=1.2)
grid()
dev.off()



### Return conditional variance
condVar <- function(cov, varCols, condCols=NULL) {
  if(is.null(condCols)) {
    condCols <- setdiff(colnames(cov), varCols)
  }
  allCols <- union(varCols, condCols)
  cov <- cov[allCols, allCols]
  cov[varCols, varCols] - cov[varCols, condCols] %*% solve(cov[condCols, condCols]) %*% cov[condCols, varCols]
}


### PCA by hand, without any package
normTab <- apply(imputedTab, 2, function(col) {
  n <- length(col)
  F <- ecdf(col)
  qnorm(n/(n+1)*F(col))
})
rownames(normTab) <- rownames(imputedTab)

for(i in 1:3) {
  vi <- vecs[, i]
  names(vi) <- colnames(normTab)
  superstat <- as.matrix(normTab) %*% vi
  names(superstat) <- rownames(normTab)
  print(tail(sort(superstat), n=15))
  print(head(sort(superstat), n=15))
  print(head(sort(abs(vi), decreasing=TRUE), n=10))
  
}

v1 <- vecs[, 1]
v2 <- vecs[, 2]
v3 <- vecs[, 3]
v4 <- vecs[, 4]
v5 <- vecs[, 5]
names(v1) <- colnames(normTab)
superstat1 <- as.matrix(normTab) %*% v1
names(superstat1) <- rownames(normTab)
names(v2) <- colnames(normTab)
superstat2 <- -1 *as.matrix(normTab) %*% v2
names(v3) <- colnames(normTab)
superstat3 <- -1 *as.matrix(normTab) %*% v3
names(v4) <- colnames(normTab)
superstat4 <- -1 *as.matrix(normTab) %*% v4
names(v5) <- colnames(normTab)
superstat5 <- -1 *as.matrix(normTab) %*% v5
rnms <- rownames(normTab)
rnms <- sapply(rnms, function(nm) paste(unlist(strsplit(nm, " - "))[1:2], collapse=" - "))
names(superstat1) <- names(superstat2) <- names(superstat3) <-
  names(superstat4) <- names(superstat5) <- rnms

RsquaredMatrix <- MostSimilar <- matrix(1, nrow=ncol(C), ncol=ncol(C))
rownames(RsquaredMatrix) <- rownames(MostSimilar) <- colnames(C)
for(metric in colnames(C)) {
  print(metric)
  nextOut <- metric
  count <- ncol(C)
  removedSet <- c(metric)
  rsq <- condVar(C, metric, setdiff(colnames(C), c(removedSet)))
  RsquaredMatrix[metric, count] <- rsq
  MostSimilar[metric, count] <- nextOut
  count <- count - 1
  
  while(count > 1) {
    
    remaining <- setdiff(colnames(C), removedSet)
    idx <- which.max(sapply(remaining, function(x)
      condVar(C, metric, setdiff(colnames(C), c(removedSet, x)))))
    nextOut <- setdiff(colnames(C), removedSet)[idx]
    removedSet <- c(removedSet, nextOut)
    rsq <- condVar(C, metric, setdiff(colnames(C), removedSet))
    RsquaredMatrix[metric, count] <- rsq
    MostSimilar[metric, count] <- nextOut
    count <- count - 1
  }
  MostSimilar[metric, 1] <- setdiff(remaining, nextOut)
}


# 4. Plot ----------------------------------------------------------------------

## 4.1. Download function ------------------------------------------------------
save_plot <- function(plot, filename) {
  png(filename = file.path(fig.dir, filename), width = 1000, height = 800)
  print(plot)
  dev.off()
}

## 4.2. Goal metrics -----------------------------------------------------------
goal_metrics <- c("Gls", "Ast", "Sh", "SoT", "xG", "npXG", "xAG", "SCA", "GCA")

valid_metrics_goals <- goal_metrics[goal_metrics %in% rownames(RsquaredMatrix)]
submat_goals <- RsquaredMatrix[valid_metrics_goals, , drop = FALSE]
df_goals <- melt(submat_goals)
colnames(df_goals) <- c("Metric", "Included_Metrics", "Independence")
df_goals$Included_Metrics <- as.numeric(as.character(df_goals$Included_Metrics))

colors_goals <- c("#332288", "#6699CC", "#88CCEE", "#44AA99",
                  "#117733", "#999933", "#DDCC77", "#661100",
                  "#CC6677")

p_goals <- ggplot(df_goals, aes(x = Included_Metrics, y = Independence, color = Metric)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = colors_goals) +
  labs(title = "Independence Curves: Goal-Related Metrics",
       x = "Number of Included Metrics",
       y = "Independence") +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right", 
    legend.title = element_blank(),
    legend.text = element_text(size = 20),   
    legend.key.size = unit(1.8, "cm"),         
    legend.spacing.y = unit(0.4, "cm")           
  )

print(p_goals)
save_plot(p_goals, "independence_goals_ggplot.png")



## 4.2. Pass metrics -----------------------------------------------------------
pass_metrics <- c("xAG", "Prg_P", "Prg_R", "Prg_C", "KP", "PPA", "AerDuel_Perc", "PassesPerc", "TakeOns_SuccPerc", "DefSum")

valid_metrics_passes <- pass_metrics[pass_metrics %in% rownames(RsquaredMatrix)]
submat_passes <- RsquaredMatrix[valid_metrics_passes, , drop = FALSE]
df_passes <- melt(submat_passes)
colnames(df_passes) <- c("Metric", "Included_Metrics", "Independence")
df_passes$Included_Metrics <- as.numeric(as.character(df_passes$Included_Metrics))

colors_passes <- c("#332288", "#6699CC", "#88CCEE", "#44AA99",
                   "#117733", "#999933", "#DDCC77", "#661100",
                   "#CC6677", "#AA4466")

p_passes <- ggplot(df_passes, aes(x = Included_Metrics, y = Independence, color = Metric)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = colors_passes) +
  labs(title = "Independence Curves: Passing-Related Metrics",
       x = "Number of Included Metrics",
       y = "Independence") +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right", 
    legend.title = element_blank(),
    legend.text = element_text(size = 20),   
    legend.key.size = unit(1.8, "cm"),         
    legend.spacing.y = unit(0.4, "cm")           
  )

print(p_passes)
save_plot(p_passes, "independence_passes_ggplot.png")
