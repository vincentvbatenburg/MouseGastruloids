require(tidyverse)
require(mFilter)
require(lomb)
require(data.table)

# Intensity profiles; takes ImageJ measurement file

int_raw       <- as.data.frame(fread("lfng_profiles.csv") %>% select(-contains("X")))

# detrending

n_timepoints  <- 69 # enter number of timepoints
dt            <- 15 # time interval between timepoints in minutes

de.trend      <- function(index){
  
  trace  <- as.data.frame(raw_int)[!is.na(raw_int[,index]), index]
  column <- as.data.frame(cbind((trace - min(trace) * 100) / (max(trace) - min(trace))
                                seq(0, n_timepoints * dt, length.out = length(trace)),
                                hpfilter((trace - min(trace) * 100) / (max(trace) - min(trace)) , freq = length(trace) * 7)$cycle,
                                index))
  
  colnames(column) <- c("Int", "Time", "Int_cycle", "replicate")
  
  return(column) }

int_detrended <- map_dfr(1:(dim(int_raw)[2]), de.trend)

# frequency decomposition

freq_bins    <- seq(0, 2, length.out = 100)

lomb.scanner <- function(replicatee){
  
  osc_temp <- int_detrended %>% filter(replicate == replicatee)
  
  freq_temp <- lsp(osc_temp$Int_norm, times = osc_temp$Time, type = "frequency", plot = T, to = 0.5/dt, ofac = 3)
  
  return(as.data.frame(cbind(freq_temp$scanned, freq_temp$power, replicatee))) }

freq_all     <- map_dfr(1:(dim(int_raw)[2]), lomb.scanner)

freq_binned  <- freq_all %>% 
     mutate(   power_norm     = ifelse(is.na(V2), 0, V2),
               frequency      = V1 * 60) %>%
   group_by(   frequency) %>%
     mutate(   frequency_bin  = freq_bins[which.min(abs(freq_bins - (frequency)))]) %>%
   group_by(   frequency_bin, replicatee) %>%
  summarise(   power_norm_bin = sum(power_norm))

