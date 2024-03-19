# Script for sensitivity analysis
library(forecast) # for ARIMA function
library(lmtest) # for ARIMA eigen value
library(ppcor) # for variable correlation
library(tsvr) # for synchrony 
library(tidyverse) # for synchrony 
library(ggplot2) 
library(ggthemes)
library(ggpattern)
library(wesanderson)
library(colorspace)
library(broom)
library(truncnorm)
library(ggpubr)
###########################################################################################################
# 1. Scale parameters for sensitivity analysis
###########################################################################################################

df.sim.sens.std.unscale <- tibble(function.int = c(df.stability.summary.small$function.int,df.stability.summary.small$function.int),
           function.name =c(df.stability.summary.small$function.name,df.stability.summary.small$function.name),
           stability =c(df.stability.summary.small$ratio.i,df.stability.summary.small$ratio.j),
           alpha.hat.intra = c(df.stability.summary.small$a_slope.i.i,df.stability.summary.small$a_slope.j.j),
           alpha.hat.inter = c(df.stability.summary.small$a_slope.i.j,df.stability.summary.small$a_slope.j.i),
           alpha.0.intra = c(df.stability.summary.small$a_initial.i.i,df.stability.summary.small$a_initial.j.j),
           alpha.0.inter = c(df.stability.summary.small$a_initial.i.j,df.stability.summary.small$a_initial.j.i),
           C.intra = c(df.stability.summary.small$c.i.i,df.stability.summary.small$c.j.j), 
           C.inter = c(df.stability.summary.small$c.i.j,df.stability.summary.small$c.j.i),
           N.opt.intra = c(df.stability.summary.small$Nmax.i.i,df.stability.summary.small$Nmax.j.j),
           N.opt.inter = c(df.stability.summary.small$Nmax.i.j,df.stability.summary.small$Nmax.j.i),
           comp.com =c(df.stability.summary.small$comp.com,df.stability.summary.small$comp.com),
           min.GR =c(df.stability.summary.small$min.GR.Ni,df.stability.summary.small$min.GR.Nj),
           min.abundance =c(df.stability.summary.small$min.abundance.Ni,df.stability.summary.small$min.abundance.Nj),
           external_factor = c(df.stability.summary.small$external_factor,df.stability.summary.small$external_factor))

ggplot(df.sim.sens.std.unscale,aes(x=comp.com)) +
  geom_bar(stat="count")

# only select the community with noise impacting the abundance, to avoid small variance which skew stability 
df.sim.sens.std.unscale <- df.sim.sens.std.unscale %>%
  filter(external_factor =="Noisy change" )


str(df.sim.sens.std.unscale)
ggplot(df.sim.sens.std.unscale, aes(y= log(stability ), 
                                    x=as.factor(function.int),
                                    color=function.int ))+ 
  geom_boxplot()  


stand.variable <- c("alpha.hat.intra","alpha.hat.inter",
                    "alpha.0.intra","alpha.0.inter",
                    "C.intra","C.inter",
                    "N.opt.intra","N.opt.inter")
names(df.sim.sens.std.unscale)
df.sim.sens.std <- df.sim.sens.std.unscale

for( n in stand.variable){
  
  df.sim.sens.std[,n] <- scale(df.sim.sens.std[,n])[,1]
}
#df.sim.sens.std[,stand.variable] <- lapply(df.sim.sens.std.unscale[,stand.variable],
                                     # scale)
str(df.sim.sens.std)
###########################################################################################################
# 4. Sensitivity analysis
##########################################################################################################
#---- Sensitivity analysis -----
stability.distribution <- ggarrange(
  ggplot(df.sim.sens.std, aes(x=log(stability))) + 
    geom_density() + labs(title="Skewed stability"),
  ggplot(df.sim.sens.std, aes(x=log(stability))) + 
    geom_density() +xlim(-7,7) +
    labs(title="Stability to have a normal distribution")
  
)
ggsave(stability.distribution,
       "figures/stability.distribution.pdf", 
       width = 10, height = 8)



df.sim.sens.std <- df.sim.sens.std %>%
  dplyr::filter(comp.com == "one-species community" | comp.com =="two-species community") %>%
  mutate(stability = log(stability), 
         function.int = as.factor(function.int)) %>%
  dplyr::filter(!is.infinite(stability) | !is.na(stability)) %>%
  dplyr::filter(stability < 7 & stability > -7)


df.sim.sens.std <- df.sim.sens.std[!is.infinite(df.sim.sens.std$stability),]
#view(df.sim.sens.std) 
#---- run glms for each function with specific parameters----
model.0 <- lm(stability ~ as.factor(function.int),
              df.sim.sens.std)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.0)
par(mfrow=c(1,1))


model.1 <- glm(formula(paste0("stability ~ ",
                             paste(c("alpha.0.intra","alpha.0.inter","alpha.0.intra:alpha.0.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==1),])

pdf(file = "figures/glm.funct1.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.1)
mtext(paste0("function 1 \n stability ~ ",
                paste(c("alpha.0.intra","alpha.0.inter",
                        "alpha.0.intra:alpha.0.inter"),
                      collapse = "+")),
      side = 3, line = - 2.5, outer = TRUE)
par(mfrow=c(1,1))
dev.off()


model.2 <- glm(formula(paste0("stability~ ",
                             paste(c("alpha.0.intra*alpha.0.inter",
                                     "alpha.hat.intra*alpha.hat.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==2),])


pdf(file = "figures/glm.funct2.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.2)
mtext(paste0("function 2 \n stability ~ ",
             paste(c("alpha.0.intra*alpha.0.inter",
                     "alpha.hat.intra*alpha.hat.inter"),
                   collapse = "+")),
      side = 3, line = - 2.5, outer = TRUE)
par(mfrow=c(1,1))
dev.off()


model.3 <- glm(formula(paste0("stability ~ ",
                             paste(c("alpha.0.intra*alpha.0.inter",
                                     "alpha.hat.intra*alpha.hat.inter",
                                     "C.intra*C.inter",
                                     "N.opt.intra*N.opt.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==3),])
summary(model.3)
pdf(file = "figures/glm.funct3.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.3)
mtext(paste0("function 3 \n stability ~ ",
             paste(c("alpha.0.intra*alpha.0.inter",
                     "alpha.hat.intra*alpha.hat.inter",
                     "C.intra*C.inter",
                     "N.opt.intra*N.opt.inter"),
                   collapse = "+")),
      side = 3, line = - 2.5, outer = TRUE)
par(mfrow=c(1,1))
dev.off()

model.4 <- glm(formula(paste0("stability~ ",
                             paste(c("alpha.0.intra*alpha.0.inter",
                                     "alpha.hat.intra*alpha.hat.inter",
                                     "C.intra*C.inter",
                                     "N.opt.intra*N.opt.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==4),])
summary(model.4)
pdf(file = "figures/glm.funct4.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.4)
mtext(paste0("function 4 \n stability ~ ",
             paste(c("alpha.0.intra*alpha.0.inter",
                     "alpha.hat.intra*alpha.hat.inter",
                     "C.intra*C.inter",
                     "N.opt.intra*N.opt.inter"),
                   collapse = "+")),
      side = 3, line = - 2.5, outer = TRUE)
par(mfrow=c(1,1))
dev.off()

#---- predict glms trends for each parameters----

predict.1 <- bind_cols(df.sim.sens.std[which(df.sim.sens.std$function.int==1),],
            data.frame(predict.stability=predict(model.1,df.sim.sens.std[which(df.sim.sens.std$function.int==1),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter"),key="term",value="value")



predict.2 <- bind_cols(df.sim.sens.std[which(df.sim.sens.std$function.int==2),],
                       data.frame(predict.stability=predict(model.2,df.sim.sens.std[which(df.sim.sens.std$function.int==2),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter",
           "alpha.hat.intra","alpha.hat.inter"),key="term",value="value")

predict.3 <- bind_cols(df.sim.sens.std[which(df.sim.sens.std$function.int==3),],
                       data.frame(predict.stability=predict(model.3,df.sim.sens.std[which(df.sim.sens.std$function.int==3),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter",
           "alpha.hat.intra","alpha.hat.inter",
           "C.intra","C.inter",
           "N.opt.intra","N.opt.inter"),key="term",value="value")

predict.4 <- bind_cols(df.sim.sens.std[which(df.sim.sens.std$function.int==4),],
                       data.frame(predict.stability=predict(model.4,df.sim.sens.std[which(df.sim.sens.std$function.int==4),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter",
           "alpha.hat.intra","alpha.hat.inter",
           "C.intra","C.inter",
           "N.opt.intra","N.opt.inter"),key="term",value="value")

predict <- bind_rows(predict.1,predict.2,predict.3,predict.4)

my_cols <- c("#AA4499", "#DDCC77","#88CCEE", "#44AA99")


predict_sensitivity_plot <- predict %>%
  ggplot(aes(x=value,color=function.name,fill=function.name)) +
  #geom_smooth(aes(y=predict.stability), se=FALSE)+
  geom_point(aes(y=stability),alpha=0.1,size=0.5) +
  geom_smooth(aes(y=predict.stability))+
  theme_bw() +
  scale_color_manual(values = darken(my_cols, amount = .1)) + 
  scale_fill_manual(values = my_cols) + 
  facet_wrap(~term,nrow = 2,dir="v",scales="free") +
  theme(legend.background = element_rect(color = "white"),
        legend.position ="right",
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, units = "in"))+
  #facet_wrap(~comp, nrow = 2) +
  labs(title="Prediction of parameters' effect \nfor each function on the stability of one species in 2-species communities",
       y = "Stability ratio (log+1) \nlow values equal lower stability and < 1 is unstable",
       x = "Range of value specific to the parameters, from low to high",
       fill = "function", color = "function") 

predict_sensitivity_plot
ggsave("figures/predict_sensitivity_plot.pdf", 
       width = 10, height = 8)

#---- See effect size of each parameters----


sens_out <- tidy(model.4) %>% 
  mutate(function.name = "4.Sigmoid") %>% 
  bind_rows(
    #(tidy(model.0) %>%
     #  mutate(term = "function",
     #         function.int = c(1,2,3,4))),
    (tidy(model.1) %>%
       mutate(function.name = "1.Constant")),
    (tidy(model.3) %>%
       mutate(function.name = "3.Exp")),
    (tidy(model.2) %>%
       mutate(function.name = "2.Linear"))
  ) %>%
  filter(term != "(Intercept)")

sens_out$term <- factor(sens_out$term, 
                        levels = 
                          c("alpha.0.intra","alpha.0.inter",
                            "alpha.hat.intra","alpha.hat.inter",
                            "C.intra","C.inter",
                            "N.opt.intra","N.opt.inter")
)
        


int_sensitivity_plot <- sens_out %>% 
  filter(term!="function") %>%
  ggplot(aes(x = estimate, xmin = (estimate - std.error), xmax = (estimate + std.error),
             color = as.factor(function.name), 
             fill = as.factor(function.name), y = as.factor(term))) +
  #geom_vline(xintercept = c(1:9+0.5), alpha = 0.5, color = "gray50", size = 0.2) +
  geom_bar(stat = "identity",width =0.6,
           position = position_dodge(), alpha = 0.7) +
  geom_errorbar(position = position_dodge(), width = 0.6, show.legend = FALSE) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  #scale_x_continuous(labes=("lowe"))
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(.87,.90),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.key.size = unit(0.1, units = "in"),
        title = element_text(size = 16),
        axis.text.y = element_text(size = 16, angle = 0, 
                                   hjust = 1, margin = margin(t =10)),
        
        axis.text.x = element_text(size = 16, angle = 0, 
                                   hjust = 1, margin = margin(t =10)),
        axis.title.x = element_text(size = 18)) +
  #facet_wrap(~comp, nrow = 2) +
  labs(title="Effect size of parameters for each function on the stability \nof one species in 2-species community",
       y = "", x = "Effect Size on Stability of i \n ratio mean/var", fill = "", color = "") 

int_sensitivity_plot
ggsave("figures/sensitivity_analysis_stability_filter.pdf", 
       width = 10, height = 8)

# plot probability of each function to have coexistence

df <- df.sim.sens.std

fits <- predict(model.0, newdata = df, se.fit = TRUE)

# Get odds from modrl
df$prediction <- fits$fit
df$upper <- fits$fit + 1.96 * fits$se.fit
df$lower <- fits$fit - 1.96 * fits$se.fit

# Plot probabilities

prediction.coexistence <- ggplot(df, aes(x=as.factor(function.name), 
                                         y=prediction,
                                         color=as.factor(function.name))) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.25, size = 1, position = position_dodge(width = 0.4)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.4)) +
  geom_hline(yintercept=0,color="black",linetype="dashed")+
  scale_color_manual(values=my_cols) +
  guides(color="none")+
  theme_light(base_size = 16) +
  labs(color="",x="",title="Stability distribution for one species in a 2-species community,\nfrom unstable (neg) to stable (positive)",
       y="Stability prediction")
prediction.coexistence
ggsave(prediction.coexistence,
       "figures/prediction.stability.pdf")