## data file names
dataFiles <- c(
    AGP = 'AGP_SGV_samples_time.RData',
    IGP = 'IGP_samples_time.RData',
    NGP = 'NGP_SGV_samples_time.RData',
    PP  = 'PP_samples_time.RData'
)

library(coda)
df <- data.frame(Model=character(), param=character(), ESS=numeric(), time=numeric())
for(f in dataFiles) {
    load(f)
    essInd <- 30001:40000  ## MCMC samples to use for calculation
    ess <- effectiveSize(samples[essInd,])
    df <- rbind(df, data.frame(
        Model = gsub('_.*', '', f),
        param = names(ess),
        ESS = as.numeric(ess),
        time = as.numeric(mcmc_time_sec)
    ))
}

## reorder Model factor
df$Model <- factor(df$Model, levels = c('NGP','IGP','AGP','PP'))

## calculate efficiency
df$Efficiency <- df$ESS / df$time

## remove (Intercept) from PP model
df <- df[-80,]

## boxplots
library(ggplot2)
ggplot(df, aes(x = Model, y = Efficiency)) + geom_boxplot()
ggsave('boxplots.pdf', width = 3, height = 2)

## calculate minimum and mean efficiency for each model
library(dplyr)
library(tidyr)
df %>%
    group_by(Model) %>%
        summarize(Minimum = min(Efficiency), Mean = mean(Efficiency)) %>%
            gather(Statistic, Efficiency, c(Minimum, Mean)) ->
                dfStats

## barplots of minimum and mean efficiencies
ggplot(dfStats, aes(x = Statistic, y = Efficiency, fill = Model)) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.7)
ggsave('barplots.pdf', width = 4, height = 2)

## generate and open PDF file with figures
tools::texi2pdf('comparisons.tex')
system('open -a "Google Chrome" comparisons.pdf')




## load('NGP_SGV_samples_time.RData')
## a <- abs(cor(samples))
## thresh <- 0.7
## sum(a > thresh & a < 1)  ## 94
## ai <- which(a > thresh & a < 1, arr.ind = TRUE)
## names <- t(apply(ai, 1, function(x) unlist(dimnames(a[x[1], x[2], drop = FALSE]))))
## dimnames(names) <- NULL
## names
      

