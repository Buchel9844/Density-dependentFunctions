##########################################################################################################
# 1. Compute slopes of growth over time
##########################################################################################################

load("results/df.sim.csv.gz")

df.glm_all <- NULL

for(i in 1:nsims){
  for( function.int in 1:4){
    print(paste0("int ", i,"for funct ",function.int))
    
    # Slope of Ni with init.cond= c(1,Nj*)
    df.Ni <- df.sim[which(df.sim$sim.i == i &
                                  df.sim$function.int == function.int &
                                  df.sim$invader =="Ni"),]
    if(sum(is.na(df.Ni$Ni))>1 |df.Ni$Nj[1] == 0|df.Ni$Ni[2] >1000 | df.Ni$Nj[2] >1000){
      df.Ni.glm <-  data.frame(com.comp.prob.Ni =F, 
                               sim = i,
                               function.int = function.int,
                               intercept.Ni =NA,
                               time.Ni = NA,
                               com.comp.Ni = 0)
    }else{
    
    df.Ni.glm <- as.data.frame(summary(glm(formula = Ni ~ time, data = df.Ni,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      dplyr::select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Ni")) %>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Ni = mean(intercept.Ni + 100*time.Ni),
             com.comp.prob.Ni = case_when(com.comp.Ni>0~T,
                                          com.comp.Ni<0~F),
             sim = i,
             function.int = function.int)
    }
    # Slope of Nj with init.cond= c(Ni*,1)
    df.Nj <- df.sim[which(df.sim$sim.i == i &
                                  df.sim$function.int == function.int &
                                  df.sim$invader =="Nj"),]
    
    if(sum(is.na(df.Nj$Nj))>1 |df.Nj$Ni[1] == 0|df.Nj$Nj[2] >1000 | df.Nj$Ni[2] >1000){
      df.Nj.glm <-  data.frame(com.comp.prob.Nj =F, 
                               sim = i,
                               function.int = function.int,
                               intercept.Nj =NA,
                               time.Nj = NA,
                               com.comp.Nj = 0)
    }else{
    df.Nj.glm <- as.data.frame(summary(glm(formula = Nj ~ time, data = df.Nj,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      dplyr::select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Nj"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Nj = mean(intercept.Nj + 100*time.Nj),
             com.comp.prob.Nj = case_when(com.comp.Nj>0~T,
                                          com.comp.Nj<0~F),
             sim = i,
             function.int = function.int)
    }
    # Slope of Nj with init.cond= c(1,1)
    df.Ni.both <- df.sim[which(df.sim$sim.i == i &
                            df.sim$function.int == function.int &
                            df.sim$invader =="both"),]
    
    if(sum(is.na(df.Ni.both$Ni))>1 |df.Ni.both$Nj[1] == 0|df.Ni.both$Ni[2] >1000 | df.Ni.both$Nj[2] >1000){
      df.Ni.both.glm <-  data.frame(com.comp.prob.Ni.both =F, 
                               sim = i,
                               function.int = function.int,
                               intercept.Ni.both =NA,
                               time.Ni.both = NA,
                               com.comp.Ni.both = 0)
    }else{
    df.Ni.both.glm <- as.data.frame(summary(glm(formula = Ni ~ time, data = df.Ni,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      dplyr::select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Ni.both"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Ni.both = mean(intercept.Ni.both + 100*time.Ni.both),
             com.comp.prob.Ni.both = case_when(com.comp.Ni.both >0 ~ T,
                                               com.comp.Ni.both < 0 ~ F),
             sim = i,
             function.int = function.int)
    }
    # Slope of Nj with init.cond= c(1,1)
  
    df.Nj.both <- df.sim[which(df.sim$sim.i == i &
                            df.sim$function.int == function.int &
                            df.sim$invader =="both"),]
    if(sum(is.na(df.Nj.both$Nj))>1 |df.Nj.both$Ni[1] == 0|df.Nj.both$Nj[2] >1000 | df.Nj.both$Ni[2] >1000){
      df.Nj.both.glm <-  data.frame(com.comp.prob.Nj.both =F, 
                               sim = i,
                               function.int = function.int,
                               intercept.Nj.both =NA,
                               time.Nj.both = NA,
                               com.comp.Nj.both = 0)
    }else{
    df.Nj.both.glm <- as.data.frame(summary(glm(formula = Nj ~ time, data = df.Nj,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      dplyr::select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Nj.both"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Nj.both = mean(intercept.Nj.both + 100*time.Nj.both),
             com.comp.prob.Nj.both = case_when(com.comp.Nj.both>0~T,
                                          com.comp.Nj.both<0~F),
             sim = i,
             function.int = function.int)
    
    }
    df.glm <- full_join(full_join(df.Ni.glm,df.Nj.glm),
                         full_join(df.Ni.both.glm,df.Nj.both.glm)) %>%
      mutate(com.comp.SLOPE = case_when(mean(com.comp.Ni,com.comp.Ni.both) > 0 &
                                         mean(com.comp.Nj,com.comp.Nj.both) > 0 ~ "j_i",
                                       mean(com.comp.Ni,com.comp.Ni.both) <= 0 & 
                                         mean(com.comp.Nj,com.comp.Nj.both) > 0 ~ "j",
                                       mean(com.comp.Ni,com.comp.Ni.both) > 0 & 
                                         mean(com.comp.Nj,com.comp.Nj.both) <= 0 ~ "i",
                                       mean(com.comp.Ni,com.comp.Ni.both) <= 0 & 
                                         mean(com.comp.Nj,com.comp.Nj.both) <= 0 ~ "0"))
    
    df.glm_all <- bind_rows(df.glm_all,df.glm)
  }
}

df.glm_all[which(is.na(df.glm_all$com.comp.SLOPE)),"com.comp.coex"] <- 0
head(df.glm_all)
str(df.glm_all)

save(df.glm_all,
     file="results/df.glm_all.csv.gz")

load("results/df.glm_all.csv.gz")


plot_com.comp.SLOPE <- ggplot(df.glm_all,aes(x=com.comp.SLOPE,
                      group=as.factor(function.int),
                      fill=as.factor(function.int)))+ 
  stat_count()+ theme_bw()  +
  scale_fill_colorblind() + 
  labs(title ="Community composition based on mean slope")
plot_com.comp.SLOPE
 
ggsave(plot_com.comp.SLOPE,
       file = "figures/com.comp.SLOPE.pdf")

##########################################################################################################
# 2. Compute minimal growth rate
##########################################################################################################

load("results/df.sim.csv.gz")

closest<-function(xv,sv){ # function to find the value in a vector closest to a value 
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]
}

df.GRWL <- NULL


for(i in 1:nsims){
  for( function.int in 1:4){
    print(paste0("int ", i,"for funct ",function.int))
    
    # Slope of Ni with init.cond= c(1,Nj*)
    df.GRWL_n <- df.sim[which(df.sim$sim.i == i &
                            df.sim$function.int == function.int),] %>%
      select(-invader) %>%
      mutate_if(is.character,as.numeric)
    
    lower.quantile_i = quantile(df.GRWL_n$Ni, 0.05, na.rm=T)

    GRWL_i <- df.GRWL_n$dNi[which(df.GRWL_n$Ni == closest(df.GRWL_n$Ni,lower.quantile_i))]
    GRWL_i <- min(GRWL_i[which(GRWL_i>=0)])
    
    if(is.na(GRWL_i)){
      GRWL_i = 0
    }
    
    lower.quantile_j = quantile(df.GRWL_n$Nj, 0.05, na.rm=T)
    GRWL_j <- df.GRWL_n$dNj[which(df.GRWL_n$Nj ==closest(df.GRWL_n$Nj,lower.quantile_j))]
    GRWL_j <- min(GRWL_j[which(GRWL_j>=0)])
    
    if(is.na(GRWL_j)){
      GRWL_j = 0
    }
    
    df.GRWL_n <- df.GRWL_n[1,] %>%
      mutate(GRWL_i=GRWL_i,
             GRWL_j=GRWL_j,
             com.comp.GRWL = case_when(GRWL_i > 0 & GRWL_j > 0 ~"j_i",
                                       GRWL_i > 0 & GRWL_j <= 0 ~"i",
                                       GRWL_i <= 0 & GRWL_j > 0 ~"j",
                                       GRWL_i <= 0 & GRWL_j <= 0 ~"0")
               
             )
    df.GRWL <- bind_rows(df.GRWL,df.GRWL_n)
    }
}
  
plot_com.comp.GRWL <- ggplot(df.GRWL,aes(x=com.comp.GRWL,
                                        group=as.factor(function.int),
                                        fill=as.factor(function.int)))+ 
  stat_count()+ theme_bw()  +
  scale_fill_colorblind() +
  labs(title ="Community composition based on GRWlow")
plot_com.comp.GRWL

ggsave(plot_com.comp.GRWL,
       file = "figures/com.comp.GRWL.pdf")

#join both methods of community composition
df.GRWL$sim <- df.GRWL$sim.i
plot_com.comp <- full_join(df.GRWL,df.glm_all) %>%
  gather(com.comp.GRWL,com.comp.coex, key="com.comp.int",value="com.comp.id") %>%
  ggplot(com.comp,aes(x=com.comp.id, group=as.factor(function.int),
                    fill=as.factor(function.int)))+ 
  stat_count() +
  scale_fill_colorblind() +
  facet_grid(.~com.comp.int) +
  theme_bw()

ggsave(plot_com.comp,
       file = "figures/com.comp.pdf")

plot_heatmap_ainit_aslope <- full_join(df.GRWL,df.glm_all) %>%
  gather(com.comp.GRWL,com.comp.SLOPE, key="com.comp.int",value="com.comp.id") %>%
  gather(a_initial.i.i,a_initial.i.j,a_initial.j.i,a_initial.j.j, 
         key="parameter_a_initial",value="a_initial") %>%
  gather(a_slope.i.i,a_slope.i.j,a_slope.j.i,a_slope.j.j, 
         key="parameter_a_slope",value="a_slope") %>%
  mutate(a_slope = case_when(function.int==1~0,
                             T~a_slope)) %>%
  ggplot(aes(x=a_slope,y=a_initial,
                      color=as.factor(com.comp.id)))+ 
  geom_point(alpha=0.5) +
  scale_color_manual(values=safe_colorblind_palette) +
  facet_grid(com.comp.int~function.int, scales="free") +
  theme_bw()
plot_heatmap_ainit_aslope 

ggsave(plot_heatmap_ainit_aslope ,
       file = "figures/plot_heatmap_ainit_aslope.pdf")


plot_heatmap_c_aslope <- full_join(df.GRWL,df.glm_all) %>%
  gather(com.comp.GRWL,com.comp.SLOPE, key="com.comp.int",value="com.comp.id") %>%
  gather(c.i.i,c.i.j,c.j.i,c.j.j, 
         key="parameter_c",value="c") %>%
  gather(a_slope.i.i,a_slope.i.j,a_slope.j.i,a_slope.j.j, 
         key="parameter_a_slope",value="a_slope") %>%
  mutate(a_slope = case_when(function.int==1~0,
                             T~a_slope))%>%
  mutate(c = case_when(function.int==1~0,
                             T~c)) %>%
  ggplot(aes(x=a_slope,y=c,
             color=as.factor(com.comp.id)))+ 
  geom_point(alpha=0.5) +
  scale_color_manual(values=safe_colorblind_palette) +
  facet_grid(com.comp.int~function.int, scales="free") +
  theme_bw()
plot_heatmap_c_aslope

plot_heatmap_N_aslope <- full_join(df.GRWL,df.glm_all) %>%
  gather(com.comp.GRWL,com.comp.SLOPE, key="com.comp.int",value="com.comp.id") %>%
  gather(Nmax.i.i,Nmax.i.j,Nmax.j.i,Nmax.j.j, 
         key="parameter_Nmax",value="Nmax") %>%
  gather(a_slope.i.i,a_slope.i.j,a_slope.j.i,a_slope.j.j, 
         key="parameter_a_slope",value="a_slope") %>%
  mutate(a_slope = case_when(function.int==1~0,
                             T~a_slope))%>%
  mutate(Nmax = case_when(function.int==1~0,
                       T~Nmax)) %>%
  ggplot(aes(x=a_slope,y=Nmax,
             color=as.factor(com.comp.id)))+ 
  geom_point(alpha=0.5) +
  scale_color_manual(values=safe_colorblind_palette) +
  facet_grid(com.comp.int~function.int, scales="free") +
  theme_bw()
plot_heatmap_N_aslope

plot_heatmap_all<- full_join(df.GRWL,df.glm_all) %>%
  gather(com.comp.GRWL,com.comp.SLOPE, key="com.comp.int",value="com.comp.id") %>%
  gather(Nmax.i.i,Nmax.i.j,Nmax.j.i,Nmax.j.j, 
         key="parameter_Nmax",value="Nmax") %>%
  gather(a_initial.i.i,a_initial.i.j,a_initial.j.i,a_initial.j.j, 
         key="parameter_a_initial",value="a_initial") %>%
  gather(c.i.i,c.i.j,c.j.i,c.j.j, 
         key="parameter_c",value="c") %>%
  gather(a_slope.i.i,a_slope.i.j,a_slope.j.i,a_slope.j.j, 
         key="parameter_a_slope",value="a_slope") %>%
  mutate(a_slope = case_when(function.int==1~0,
                             T~a_slope))%>%
  mutate(Nmax = case_when(function.int==1~0,
                          T~Nmax)) %>%
  mutate(c = case_when(function.int==1~0,
                          T~c)) %>%
  gather(a_initial, a_slope, Nmax, c, 
         key="parameter",value="value") %>%
  ggplot(aes(x=value,y=com.comp.int,
             color=as.factor(com.comp.id)))+ 
  geom_point(alpha=0.5) +
  scale_color_manual(values=safe_colorblind_palette) +
  facet_wrap(parameter~function.int, 
             scales= "free_x") +
  theme_bw()
plot_heatmap_all


ggsave(plot_heatmap_all,
       file = "figures/plot_heatmap_all.pdf")


##########################################################################################################
# 3. Standardization
##########################################################################################################
df.sim$sim <- df.sim$sim.i

df.sim.std <- left_join(df.sim,df.glm_all)
                        
# check that no competition outcome is NA
df.glm_all[which(df.glm_all$time==1 & df.glm_all$function.int==4 & is.na(df.glm_all$com.comp.coex)),]
df.sim.std[which(df.sim.std$time==1 & df.sim.std$function.int==4 & is.na(df.sim.std$com.comp.coex)),]

stand.variable <- c(paste0("Nmax",c(".i.j",".i.i",".j.i",".j.j")),
                    paste0("c",c(".i.j",".i.i",".j.i",".j.j")),
                    paste0("a_initial",c(".i.j",".i.i",".j.i",".j.j")),
                    paste0("a_slope",c(".i.j",".i.i",".j.i",".j.j")))

df.sim.std[,stand.variable] <- lapply(df.sim.std[,stand.variable],
                                               scale)

##########################################################################################################
# 5. Sensitivity analysis
##########################################################################################################

# Incorporate interaction terms
df.sim.std <-  df.sim.std %>%
  mutate( com.comp.coex.int = case_when(com.comp.coex =="j_i" ~1,
                                        TRUE ~ 0))

model.0 <- glm( com.comp.coex.int ~ as.factor(function.int), 
               data = df.sim.std,
               family = "binomial")

model.1 <- glm(formula(paste0("com.comp.coex.int  ~ ",
                              paste(paste0("a_initial",c(".i.j",".i.i",".j.i",".j.j")),collapse = "+"))), 
               data = df.sim.std[which(df.sim.std$function.int==1),],
               family = "binomial")

model.2 <- glm(formula(paste0("com.comp.coex.int   ~ ",
                              paste(c(paste0("Nmax",c(".i.j",".i.i",".j.i",".j.j")),
                                      paste0("a_initial",c(".i.j",".i.i",".j.i",".j.j")),
                                      paste0("a_slope",c(".i.j",".i.i",".j.i",".j.j"))),collapse = "+"))), 
               data = df.sim.std[which(df.sim.std$function.int==2),],
               family = "binomial")

model.3 <- glm(formula(paste0("com.comp.coex.int  ~ ",
                              paste(stand.variable,collapse = "+"))), 
               data = df.sim.std[which(df.sim.std$function.int==3),],
               family = "binomial")

model.4 <- glm(formula(paste0("com.comp.coex.int  ~ ",
                              paste(stand.variable,collapse = "+"))), 
               data = df.sim.std[which(df.sim.std$function.int==4),],
               family = "binomial")

# build dataframe to plot
ggplot( df.sim.std, aes(x=com.comp.coex.int, 
                                 fill=as.factor(function.int)))+
  stat_count() +theme_bw()


sens_out <- tidy(model.2) %>% 
  mutate(function.int = 2) %>% 
  bind_rows(
    #(tidy(model.0) %>%
       #mutate(term = "function",
       #       function.int = c(1,2,3,4))),
    #(tidy(model.1) %>%
    #   mutate(function.int = 1)),
    (tidy(model.3) %>%
       mutate(function.int = 3)),
    (tidy(model.4) %>%
       mutate(function.int = 4))
  ) %>%
  filter(term != "(Intercept)")




sens_out$term <- factor(sens_out$term, 
                        levels = 
                          c("function",
                            paste0("a_initial",c(".i.i",".j.j",".i.j",".j.i")),
                            paste0("a_slope",c(".i.i",".j.j",".i.j",".j.i")),
                            paste0("Nmax",c(".i.i",".j.j",".i.j",".j.i")),
                            paste0("c",c(".i.i",".j.j",".i.j",".j.i"))
                          )
)
# change to have the focal species first
x_axis_labs <- c("function",
                 paste0("a_initial",c(".i.i",".j.j",".j.i",".i.j")),
                 paste0("a_slope",c(".i.i",".j.j",".j.i",".i.j")),
                 paste0("Nmax",c(".i.i",".j.j",".j.i",".i.j")),
                 paste0("c",c(".i.i",".j.j",".j.i",".i.j"))
)

my_cols <- qualitative_hcl(n =4)

int_sensitivity_plot <- sens_out %>% 
  
  ggplot(aes(y = estimate, ymin = (estimate - std.error), ymax = (estimate + std.error),
             color = as.factor(function.int), 
             fill = as.factor(function.int), x = term)) +
  #geom_vline(xintercept = c(1:9+0.5), alpha = 0.5, color = "gray50", size = 0.2) +
  geom_bar(stat = "identity",
           position = position_dodge(), alpha = 0.7) +
  geom_errorbar(position = position_dodge(.9), width = 0.5, show.legend = FALSE) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = c(.07,.92),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.key.size = unit(0.1, units = "in"),
        axis.text.x = element_text(size = 12, angle = 45, 
                                   hjust = 1, margin = margin(t =10))) +
  #facet_wrap(~comp, nrow = 2) +
  labs(x = "", y = "Effect Size on \n probability of coexistence", fill = "function", color = "function") +
  scale_x_discrete(labels = x_axis_labs)

int_sensitivity_plot
ggsave("figures/sensitivity_analysis.pdf", width = 10, height = 8)



# plot probability of each function to have coexistence

df <- df.sim.std

fits <- predict(model.0, newdata = df, se.fit = TRUE)

# Get odds from modrl
df$prediction <- exp(fits$fit) 
df$upper <- exp(fits$fit + 1.96 * fits$se.fit)
df$lower <- exp(fits$fit - 1.96 * fits$se.fit)

# Convert odds to probabilities
df$prediction <- df$prediction / (1 + df$prediction)
df$upper <- df$upper / (1 + df$upper)
df$lower <- df$lower / (1 + df$lower)

# Plot probabilities

prediction.Slopes <- ggplot(df, aes(function.int, prediction)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.25, size = 1, position = position_dodge(width = 0.4)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.4)) +
  theme_light(base_size = 16) +
  scale_y_continuous(name = "Probability of coexistence", limits = c(0, 1),
                     labels = scales::percent)
prediction.Slopes
ggsave(prediction.Slopes,
       "figures/prediction.Slopes.pdf")
