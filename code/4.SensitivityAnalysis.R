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
###########################################################################################################
# 1. Scale parameters for sensitivity analysis
###########################################################################################################

df.sim.sens.std.unscale <- data.frame(function.int = c(df.stability.summary.small$function.int,df.stability.summary.small$function.int),
           function.name =c(df.stability.summary.small$function.name,df.stability.summary.small$function.name),
           stability =c(df.stability.summary.small$ratio.i,df.stability.summary.small$ratio.j),
           alpha.decay.intra = c(df.stability.summary.small$a_slope.i.i,df.stability.summary.small$a_slope.j.j),
           alpha.decay.inter = c(df.stability.summary.small$a_slope.i.j,df.stability.summary.small$a_slope.j.i),
           alpha.0.intra = c(df.stability.summary.small$a_initial.i.i,df.stability.summary.small$a_initial.j.j),
           alpha.0.inter = c(df.stability.summary.small$a_initial.i.j,df.stability.summary.small$a_initial.j.i),
           C.intra = c(df.stability.summary.small$c.i.i,df.stability.summary.small$c.j.j), 
           C.inter = c(df.stability.summary.small$c.i.j,df.stability.summary.small$c.j.i),
           N.opt.intra = c(df.stability.summary.small$Nmax.i.i,df.stability.summary.small$Nmax.j.j),
           N.opt.inter = c(df.stability.summary.small$Nmax.i.j,df.stability.summary.small$Nmax.j.i),
           comp.com =c(df.stability.summary.small$comp.com,df.stability.summary.small$comp.com),
           min.GR =c(df.stability.summary.small$min.GR.Ni,df.stability.summary.small$min.GR.Nj),
           min.abundance =c(df.stability.summary.small$min.abundance.Ni,df.stability.summary.small$min.abundance.Nj))

ggplot(df.sim.sens.std.unscale, aes(y= stability , 
                                    x=as.factor(function.int),
                                    color=function.int ))+ 
  geom_boxplot()  +
  ylim(c(0,100))

ggplot(df.sim.sens.std.unscale, aes(y= stability , 
                                    x=as.factor(function.int),
                                    color=function.int ))+ 
  geom_boxplot()  +
  ylim(c(0,100))

stand.variable <- c("alpha.decay.intra","alpha.decay.inter",
                    "alpha.0.intra","alpha.0.inter",
                    "C.intra","C.inter",
                    "N.opt.intra","N.opt.inter")
names(df.sim.sens.std.unscale)
df.sim.sens.std <- df.sim.sens.std.unscale
df.sim.sens.std[,stand.variable] <- lapply(df.sim.sens.std.unscale[,stand.variable],
                                      scale)

###########################################################################################################
# 4. Sensitivity analysis
##########################################################################################################
#---- Sensitivity analysis -----
df.sim.sens.std <- df.sim.sens.std %>%
  mutate(stability = case_when(stability > 100 ~ 100,
                               is.na(stability) ~ 0,
                               T~stability),
         function.int = as.factor(function.int))

model.0 <- lm(stability ~ as.factor(function.int),
              df.sim.sens.std)

model.1 <- glm(formula(paste0("stability ~ ",
                             paste(c("alpha.0.intra","alpha.0.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==1),])


model.2 <- glm(formula(paste0("stability~ ",
                             paste(c("alpha.0.intra","alpha.0.inter",
                                     "alpha.decay.intra","alpha.decay.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==2),])

model.3 <- glm(formula(paste0("stability ~ ",
                             paste(c("alpha.0.intra","alpha.0.inter",
                                     "alpha.decay.intra","alpha.decay.inter",
                                     "C.intra","C.inter",
                                     "N.opt.intra","N.opt.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==3),])

model.4 <- glm(formula(paste0("stability~ ",
                             paste(c("alpha.0.intra","alpha.0.inter",
                                     "alpha.decay.intra","alpha.decay.inter",
                                     "C.intra","C.inter",
                                     "N.opt.intra","N.opt.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==4),])

predict.1 <- bind_cols(df.sim.sens.std.unscale[which(df.sim.sens.std.unscale$function.int==1),],
            data.frame(predict.stability=predict(model.1,df.sim.sens.std[which(df.sim.sens.std$function.int==1),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter"),key="term",value="value")



predict.2 <- bind_cols(df.sim.sens.std.unscale[which(df.sim.sens.std.unscale$function.int==2),],
                       data.frame(predict.stability=predict(model.2,df.sim.sens.std[which(df.sim.sens.std$function.int==2),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter",
           "alpha.decay.intra","alpha.decay.inter"),key="term",value="value")

predict.3 <- bind_cols(df.sim.sens.std.unscale[which(df.sim.sens.std.unscale$function.int==3),],
                       data.frame(predict.stability=predict(model.3,df.sim.sens.std[which(df.sim.sens.std$function.int==3),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter",
           "alpha.decay.intra","alpha.decay.inter",
           "C.intra","C.inter",
           "N.opt.intra","N.opt.inter"),key="term",value="value")

predict.4 <- bind_cols(df.sim.sens.std.unscale[which(df.sim.sens.std.unscale$function.int==4),],
                       data.frame(predict.stability=predict(model.4,df.sim.sens.std[which(df.sim.sens.std$function.int==4),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter",
           "alpha.decay.intra","alpha.decay.inter",
           "C.intra","C.inter",
           "N.opt.intra","N.opt.inter"),key="term",value="value")

predict <- bind_rows(predict.1,predict.2,predict.3,predict.4)

my_cols <- c("#AA4499", "#DDCC77","#88CCEE", "#44AA99")


predict_sensitivity_plot <- predict %>%
  ggplot(aes(x=value,color=function.name,fill=function.name)) +
  #geom_smooth(aes(y=predict.stability), se=FALSE)+
  geom_point(aes(y=log(stability+1)),alpha=0.1,size=0.5) +
  geom_smooth(aes(y=log(predict.stability+1)))+
  theme_bw() +
  scale_color_manual(values = darken(my_cols, amount = .1)) + 
  scale_fill_manual(values = my_cols) + 
  facet_wrap(~term,nrow = 2,dir="v",scales="free") +
  scale_y_continuous(limits=c(0,5)) +
  theme(legend.background = element_rect(color = "white"),
        legend.position ="right",
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, units = "in"))+
  #facet_wrap(~comp, nrow = 2) +
  labs(title="Prediction of parameters' effect for each function on the stability of one species in 2-species community",
       y = "Stability ratio (log+1) \nlow values equal lower stability and < 1 is unstable",
       x = "Range of value specific to the parameters, from low to high",
       fill = "function", color = "function") 

predict_sensitivity_plot
ggsave("figures/predict_sensitivity_plot.pdf", 
       width = 10, height = 8)


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
                            "alpha.decay.intra","alpha.decay.inter",
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
        title = element_text(size = 18),
        axis.text.y = element_text(size = 16, angle = 0, 
                                   hjust = 1, margin = margin(t =10)),
        
        axis.text.x = element_text(size = 16, angle = 0, 
                                   hjust = 1, margin = margin(t =10)),
        axis.title.x = element_text(size = 18)) +
  #facet_wrap(~comp, nrow = 2) +
  labs(title="Effect size of parameters for each function on the stability of one species\nin 2-species community",
       y = "", x = "Effect Size on Stability of i \n ratio mean/var", fill = "", color = "") 

int_sensitivity_plot
ggsave("figures/sensitivity_analysis_stability.pdf", 
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
