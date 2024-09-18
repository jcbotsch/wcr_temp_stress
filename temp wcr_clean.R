#====load====
library(tidyverse)
library(lme4)
library(segmented)

addC = function(x,...){format(paste0(x, "°C"))} # add degrees C to axis

#====aesthetics====
theme_set(theme_bw()+
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.position = "top",
                  text = element_text(size = 12),
                  axis.title = element_text(size = 13),
                  legend.spacing = unit(0,units = 'points'),
                  legend.margin = margin(c(1,5,5,5)),
                  legend.box.spacing = unit(0, units = "points")))

#====read preliminary data====
tcrit <- read_csv("wcr_temp2024.csv")

#====prepare CTmax data=====
tcrit <- tcrit %>% 
  mutate(ctmax = as.numeric(ctmax),
         pre_temp = as.numeric(pre_temp)) 

#====Do control rootworm collected from different latitudes show differences in CTmax?====
# subset data
ctrls <- tcrit %>% filter(control == 1)

#average for each line
ctrls_summary <- ctrls %>% group_by(line, lat, BIO1, BIO5, BIO7, BIO10) %>% summarise(mean = mean(ctmax), se = sd(ctmax)/sqrt(n()), pre_temp = NA)

# fit model
ctrlm <- lm(mean~lat, data = ctrls_summary)  
summary(ctrlm)

#====How does pre-assay temperature exposure affect CTmax?====
tc_summary <- tcrit %>% 
  filter(control != 1,
         !is.na(ctmax)) %>% 
  group_by(line, pre_temp, lat, BIO1, BIO10, BIO5, BIO7) %>% 
  summarise(se = sd(ctmax)/sqrt(n()),
            lower = quantile(ctmax, 1-0.75),
            upper = quantile(ctmax, 0.75),
            median = median(ctmax),
            mean = mean(ctmax)) %>% 
  ungroup %>% 
  mutate(se_weight = (max(se)-se)/max(se))

# breakpoint
lmod <- lm(mean~pre_temp*line, weight = se_weight, data = tc_summary)
davies.test(lmod, ~pre_temp) # non significant davies test suggests no evidence for breakpoint

# get model results from linear model
summary(lmod)
car::Anova(lmod)

#####FIGURES#####
#====plot 1====

p1 <- ctrls_summary %>% 
  ungroup %>% 
  mutate(hjust = c(-0.1, 1, -0.1, -0.1, 1, -0.1)) %>% 
  ggplot()+
  geom_pointrange(aes(x = lat, y = mean, col = lat, ymin = mean - se, ymax = mean + se), alpha = 0.5)+
  geom_text(aes(x = lat, y = mean, label = line, hjust = hjust), vjust = 0)+
  labs(x = "Latitude",
       y = expression("CT"["max"]),
       color = "Source latitude",
       fill = "Source latitude")+
  scale_y_continuous(limits = c(42, 45), labels = addC)+
  scale_x_continuous(limits = c(NA, 46))+
  scale_color_viridis_c(option = "magma", end = 0.75, begin = 0.1, direction = -1)+
  scale_fill_viridis_c(option = "magma", end = 0.75, begin = 0.1, direction = -1) +
  theme(legend.position = "none")

# ggsave(plot = p1, width = 3, height = 3,dpi = 650, units = "in", filename = "latCTmax.pdf")

# fig.2
p2 <- tc_summary %>% 
  mutate(fit = predict(lmod, newdata = .),
         fit.se = predict(lmod, newdata = ., se = TRUE)$se) %>% 
  ggplot()+
  facet_wrap(~fct_rev(line), ncol = 1) +
  # data
  geom_pointrange(aes(x = pre_temp, y = mean, col = lat, ymin = mean - se, ymax = mean + se), alpha = 0.5) +
  # linear model fit
  geom_line(aes(x = pre_temp, y = fit, col = lat)) +
  geom_ribbon(aes(x = pre_temp, y = fit, ymin = fit - fit.se, ymax = fit + fit.se, fill = lat), alpha = 0.5) +
  # starting CTmax
  geom_rect(aes(fill = lat, xmin = -Inf, xmax = Inf, ymin = mean - se, ymax = mean + se, ), alpha = 0.2, data = ctrls_summary) +
  geom_hline(aes(col = lat,  yintercept = mean), linetype = "dashed", alpha = 0.5, linewidth = 0.3, data = ctrls_summary) +
  # aesthetics
  labs(x = "Exposure temperature",
       y = expression("CT"["max"]),
       color = "Source latitude",
       fill = "Source latitude")+
  scale_x_continuous(breaks = scales::pretty_breaks(), labels = addC)+
  scale_y_continuous(limits = c(39, 46), labels = addC)+
  geom_rug(aes(x = pre_temp), length = unit(0.3, "lines"), sides = "b", col = "lightgray")+
  scale_color_viridis_c(option = "magma", end = 0.75, begin = 0.1, direction = -1)+
  scale_fill_viridis_c(option = "magma", end = 0.75, begin = 0.1, direction = -1) +
  theme(legend.position = "none",
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(vjust = 1))


# ggsave(plot = p2, width = 3, height = 8,dpi = 650, units = "in", filename = "exposure_plot.pdf")
