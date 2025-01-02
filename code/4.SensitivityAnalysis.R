# Script for sensitivity analysis
library(forecast) # for ARIMA function
library(lmtest) # for ARIMA eigen value
library(ppcor) # for variable correlation
library(tsvr) # for synchrony 
library(tidyverse) # for synchrony 
library(dplyr)
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
load("results/df.stability.summary.csv.gz") 
df.stability.summary.small <- df.stability.summary

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


ggplot(df.sim.sens.std.unscale ,aes(x=as.factor(comp.com)))+
  geom_bar(stat="count")

# only select the community with noise impacting the abundance, to avoid small variance which skew stability 
df.sim.sens.std.unscale <- df.sim.sens.std.unscale %>%
  dplyr::filter(external_factor =="Noisy change" &
  #dplyr::filter(external_factor =="No external factor" &
                comp.com =="two-species community") %>%
  dplyr::filter(stability < 10 & stability > 0)
  
ggplot(df.sim.sens.std.unscale ,aes(x=alpha.0.inter, y=alpha.0.intra,
                                    color=stability))+
  geom_point()  +
  scale_color_gradient() +
  facet_grid(.~function.name) +
  geom_abline(slope=1,intercept=0)

str(df.sim.sens.std.unscale)
ggplot(df.sim.sens.std.unscale, aes(y= stability, 
                                    x=as.factor(function.int),
                                    color=function.int ))+ 
  geom_boxplot()  


stand.variable <- c("alpha.hat.intra","alpha.hat.inter",
                    "alpha.0.intra","alpha.0.inter",
                    "C.intra","C.inter",
                    "N.opt.intra","N.opt.inter")
names(df.sim.sens.std.unscale)
df.sim.sens.std <- df.sim.sens.std.unscale
str(df.sim.sens.std.unscale)
for( n in stand.variable){
  
  df.sim.sens.std[,n] <- scale(df.sim.sens.std[,n])[,1]
}
#df.sim.sens.std[,stand.variable] <- lapply(df.sim.sens.std.unscale[,stand.variable],
                                     # scale)
str(df.sim.sens.std)

ggplot(df.sim.sens.std, aes(y = log(stability), 
                                    x=as.factor(function.int),
                                    color=function.int ))+ 
  geom_boxplot()  

###########################################################################################################
# 4. Sensitivity analysis
##########################################################################################################
#---- Sensitivity analysis -----
stability.distribution <- ggarrange(
  ggplot(df.sim.sens.std, aes(x=log(stability))) + 
    geom_density() + labs(title="Skewed stability")+
    theme_bw(),
  ggplot(df.sim.sens.std, aes(x=log(stability))) + 
    geom_density() +xlim(-7,7) +
    labs(title="Stability to have a normal distribution") +
    theme_bw()
)
stability.distribution 
ggsave("figures/stability.distribution.pdf", 
       stability.distribution,
       width = 10, height = 8)


str(df.sim.sens.std)
df.sim.sens.std <- df.sim.sens.std %>%
  mutate(log.stability = log(stability), 
         function.int = as.factor(function.int)) %>%
  dplyr::filter(!is.infinite(log.stability)) %>%
  dplyr::filter(!is.na(log.stability)) %>%
  dplyr::filter(log.stability < 10 ) %>%
  dplyr::filter(log.stability > -10)
str(df.sim.sens.std)

#df.sim.sens.std <- df.sim.sens.std[!is.infinite(df.sim.sens.std$stability),]
ggplot(df.sim.sens.std, aes(y = log.stability, 
                            x=as.factor(function.int),
                            color=function.int ))+ 
  geom_boxplot() 
#---- run glms for each function with specific parameters----
model.0 <- lm(stability ~ as.factor(function.int),
              df.sim.sens.std)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.0)
par(mfrow=c(1,1))

#view(df.sim.sens.std[which(df.sim.sens.std$function.int==1),])
model.1 <- glm(formula(paste0("log.stability ~ ",
                             paste(c("alpha.0.intra","alpha.0.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==1),],
              family=gaussian)

df.glm <- summary(model.1)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column(var="parameters") %>%
  mutate(funct.form ="1.constant")

pdf(file = "figures/glm.funct1.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.1)
mtext(paste0("function 1 \n stability ~ ",
                paste(c("alpha.0.intra","alpha.0.inter"),
                      collapse = "+")),
      side = 3, line = - 2.5, outer = TRUE)
par(mfrow=c(1,1))
dev.off()


model.2 <- glm(formula(paste0("log.stability~ ",
                             paste(c("alpha.0.intra","alpha.0.inter",
                                     "alpha.hat.intra","alpha.hat.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==2),],
              family=gaussian)

df.glm <- summary(model.2)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column(var="parameters") %>%
  mutate(funct.form ="2.linear") %>%
  bind_rows(df.glm)

pdf(file = "figures/glm.funct2.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.2)
mtext(paste0("function 2 \n stability ~ ",
             paste(c("alpha.0.intra","alpha.0.inter",
                     "alpha.hat.intra","alpha.hat.inter"),
                   collapse = "+")),
      side = 3, line = - 2.5, outer = TRUE)
par(mfrow=c(1,1))
dev.off()


model.3 <- glm(formula(paste0("log.stability ~ ",
                             paste(c("alpha.0.intra","alpha.0.inter",
                                     "alpha.hat.intra","alpha.hat.inter",
                                     "C.intra","C.inter",
                                     "N.opt.intra","N.opt.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==3),],
              family=gaussian)
summary(model.3)
df.glm <- summary(model.3)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column(var="parameters") %>%
  mutate(funct.form ="3.exp") %>%
  bind_rows(df.glm)

pdf(file = "figures/glm.funct3.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.3)
mtext(paste0("function 3 \n stability ~ ",
             paste(c("alpha.0.intra","alpha.0.inter",
                     "alpha.hat.intra","alpha.hat.inter",
                     "C.intra","C.inter",
                     "N.opt.intra","N.opt.inter"),
                   collapse = "+")),
      side = 3, line = - 2.5, outer = TRUE)
par(mfrow=c(1,1))
dev.off()

model.4 <- glm(formula(paste0("log.stability~ ",
                             paste(c("alpha.0.intra","alpha.0.inter",
                                     "alpha.hat.intra","alpha.hat.inter",
                                     "C.intra","C.inter",
                                     "N.opt.intra","N.opt.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==4),],
              family=gaussian)
summary(model.4)
df.glm <- summary(model.4)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column(var="parameters") %>%
  mutate(funct.form ="4.sigmoid") %>%
  bind_rows(df.glm)
write.csv(df.glm, 
          "results/df.glm.parameter.csv")

pdf(file = "figures/glm.funct4.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(model.4)
mtext(paste0("function 4 \n stability ~ ",
             paste(c("alpha.0.intra","alpha.0.inter",
                     "alpha.hat.intra","alpha.hat.inter",
                     "C.intra","C.inter",
                     "N.opt.intra","N.opt.inter"),
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
  geom_point(aes(y=log.stability),alpha=0.1,size=0.5) +
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
  mutate(function.name = "Sigmoid") %>% 
  bind_rows(
    #(tidy(model.0) %>%
     #  mutate(term = "function",
     #         function.int = c(1,2,3,4))),
    (tidy(model.1) %>%
       mutate(function.name = "Traditional constant")),
    (tidy(model.3) %>%
       mutate(function.name = "Exp")),
    (tidy(model.2) %>%
       mutate(function.name = "Linear"))
  ) %>%
  filter(term != "(Intercept)") %>%
  mutate(function.name  = factor(function.name,
                                 levels=c("Traditional constant","Linear","Exp","Sigmoid"))) %>%
  mutate(term = factor(sens_out$term, 
                        levels = 
                          c("alpha.0.intra","alpha.0.inter",
                            "alpha.hat.intra","alpha.hat.inter",
                            "C.intra","C.inter",
                            "N.opt.intra","N.opt.inter")))
        
levels(sens_out$term)

int_sensitivity_plot <-  ggplot() +
  geom_vline(xintercept=0,color="black",size=0.75) +
  geom_bar(data=sens_out,
           aes(x = estimate,
               color = as.factor(function.name), 
               fill = as.factor(function.name), y = as.factor(term)),
           stat = "identity",width =0.6,
           position = position_dodge(), alpha = 0.7) +
  geom_errorbar(data=sens_out,
                aes(x = estimate, xmin = (estimate - std.error), 
                    xmax = (estimate + std.error),
                    color = as.factor(function.name), 
                    y = as.factor(term)),
                position = position_dodge(), 
                width = 0.6, show.legend = FALSE) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  scale_y_discrete(labels =c(expression(paste(alpha["0,ii"])),
                             expression(paste(alpha["0,ij"])),
                               expression(paste(widehat(alpha["ii"]))),
                                 expression(paste(widehat(alpha["ij"]))),
                             expression(paste(C["ii"])),
                             expression(paste(C["ij"])),
                             expression(paste(N["o,ii"])),
                             expression(paste(N["o,ij"])))) +
  scale_x_continuous(breaks=c(-0.4,0.,.4)) + 
  #xaxis
  annotate('text', x = -0.56, y = 0.3, 
           label = "less stable",size=5) +
  annotate("segment", x = -0.5, y = 0.2, 
           xend = -0.7, yend = 0.2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate('text', x = 0.57, y = 0.3, 
           label = "more stable",size=5) +
  annotate("segment", x = 0.5, y = 0.2, 
           xend = 0.7, yend = 0.2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  # alpha zero
  geom_path(data=bracket(0.7,0.02,1.5,2),
            aes(x=x,y=y),color="black",size=1) + # x, width, y, height
  annotate('text', x = 0.83, y = 1.5, 
           label = "Higher value \nmeans more \nfacilitative\ninteractions",size=5) +
  geom_path(data=bracket(0.7,0.02,4.5,3),
            aes(x=x,y=y),color="black",size=1) + # x, width, y, height
  annotate('text', x = 0.83, y = 4.5, 
           label = "Higher value \nmeans less \nfacilitation\n at low density",size=5) +
  geom_path(data=bracket(0.7,0.02,7.5,2),
            aes(x=x,y=y),color="black",size=1) + #
  annotate('text', x = 0.83, y = 7.5, 
           label = "Higher value \nmeans more \nfacilitation\n at low density",size=5) +
  theme_clean() +
  coord_cartesian(ylim = c(0.5, 8.5),xlim = c(-0.7, 0.7), clip = "off",expand=FALSE) + 
  theme(panel.grid.major.y = element_line(color=col_grid ,size = 0.75,
                                          linetype = 2),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(0.1,4,0.1,0.1), "cm"),
        legend.position = "bottom", #c(.87,.90),
        legend.background = element_rect(color = "white"),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, units = "in"),
        title = element_text(size = 16),
        axis.text.y = element_text(size = 16, angle = 0, 
                                   hjust = 1, margin = margin(t =10)),
        
        axis.text.x = element_text(size = 16, angle = 0, 
                                   hjust = 1, margin = margin(t =10)),
        axis.title.x = element_text(size = 18)) +
  #facet_wrap(~comp, nrow = 2) +
  labs(#title="Effect size of parameters for each function on the stability \nof one species in 2-species community",
       y = "", x = "Effect of parameter strength on stability", fill = "", color = "") 

int_sensitivity_plot
col_grid <- scales::alpha("grey92", .6)
bracket <- function(x, width, y, height){
  data.frame(
    x=c(0,1,1,2,1,1,0)/2*(width) + x,
    y=(c(0,1,4,5,6,9,10)/10-0.5)*(height) + y
  )
}

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
