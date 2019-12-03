require(tidyverse)
require(data.table)

# Intensity profiles; takes ImageJ measurement file

length_max <- "colname" # distance colname with highest number of observations  

Uncx       <- fread(   "Uncx.csv") %>% 
             select(   length_max, contains("Y")) %>% 
             gather(   value = Intensity, key = replicate, -length_max) %>%
             mutate(   gene  = "Uncx")

Tbx18      <- fread(   "Tbx18.csv") %>% 
             select(   length_max, contains("Y")) %>% 
             gather(   value = Intensity, key = replicate, -length_max) %>%
             mutate(   gene  = "Tbx18")

# Smoothing and peak calling

Intensity_prof <- 
  bind_rows(   Uncx, Tbx18) %>% 
     filter(   !is.na(Intensity)) %>%
   group_by(   gene, replicate) %>% 
     mutate(   Intensity_sc  = (Intensity - min(Intensity, na.rm = T)) / max(Intensity - min(Intensity, na.rm = T), na.rm = T),
               AP_distance   = X8) %>%
       nest(   -replicate) %>%
     mutate(   smoothing1    = map(data, function(df){lowess(df$Intensity_sc, f = 15 / length(df$Intensity_sc))$y})) %>%
     unnest(   ) %>%
   group_by(   replicate) %>%     
     mutate(   peaks         = (c(0,diff(diff(smoothing1) < 0),0) > 0),
               cump1         = cumsum(peaks),
               cump2         = ifelse(peaks > 0, cump1, NA)) %>%
    ungroup(   ) %>% 
     filter(   gene == "Uncx") %>%
  
     ggplot(   aes(   x = AP_distance, y = Intensity_sc, col = gene)) +
         geom_line(   size   = 1) +
         geom_line(   aes(y = smoothing1), col = "black") +
         geom_text(   aes(y = smoothing1, label = cump2), col = "blue") +
scale_alpha_manual(   values = c(0,1), na.value = 0) +
        facet_grid(   rows   = vars(replicate),
                      cols   = vars(gene))

# Manuel peak selection

Intensity_peak <- Intensity_prof %>% 
  filter(   (replicate == "Y0" & (cump2 %in% c(1,2,3,5))) |
            (replicate == "Y1" & (cump2 %in% c(1,2,3))) |
            (replicate == "Y2" & (cump2 %in% c(3,4,6,7,9,11,12,14))) |
            (replicate == "Y3" & (cump2 %in% c(2,4,5,7))) | 
            (replicate == "Y4" & (cump2 %in% c(12:15))) |
            (replicate == "Y5" & (cump2 %in% c(3,5,7))) |
            (replicate == "Y6" & (cump2 %in% c(3,6,9,11))) |
            (replicate == "Y7" & (cump2 %in% c(1:7,9,11))) |
            (replicate == "Y8" & (cump2 %in% c(5:10, 12:14)))) %>%
group_by(   replicate) %>%
  mutate(   peak_to_peak = c(diff(AP_distance), 0),
            somite_index = as.numeric(as.factor(cump2))) %>% 
  filter(   peak_to_peak > 0) %>%
 ungroup(   ) %>% 
  mutate(   group = ifelse(replicate == "Y8", 5, 
                           ifelse(replicate == "Y7", 4, 
                           ifelse(replicate == "Y0" | replicate == "Y1", 1, 
                           ifelse(replicate == "Y5" | replicate == "Y6", 3,2))))) %>%
  
  ggplot(   aes(   x = somite_index, y = peak_to_peak, col = as.factor(replicate))) +
      geom_line(   ) +
scale_y_continuous(   limits = c(0, 150)) 







