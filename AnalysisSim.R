##########################################################################################################
# Compute slopes of growth over time
##########################################################################################################

load("results/df.sim.csv.gz")

df.glm_all <- NULL

for(i in 1:nsims){
  for( function.int in 1:4){
    print(paste0("int ", i,"for funct ",function.int))
    
    # Slope of Ni with init.cond= c(1,Nj*)
    df.Ni <- df.sim[which(df.sim$sim == i &
                                  df.sim$function.int == function.int &
                                  df.sim$invader =="Ni"),]
    if(sum(is.na(df.Ni$Ni))>1) next
    if(df.Ni$Nj[1] == 0) next
    if(df.Ni$Ni[2] >1000 | df.Ni$Nj[2] >1000) next
    
    df.Ni.glm <- as.data.frame(summary(glm(formula = Ni ~ time, data = df.Ni,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Ni"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Ni = intercept.Ni + 100*time.Ni,
             com.comp.prob.Ni = case_when(com.comp.Ni>0~T,
                                          com.comp.Ni<0~F),
             sim = i,
             function.int = function.int)
    
    # Slope of Nj with init.cond= c(Ni*,1)
    df.Nj <- df.sim[which(df.sim$sim == i &
                                  df.sim$function.int == function.int &
                                  df.sim$invader =="Nj"),]
    if(sum(is.na(df.Nj$Nj))>1) next
    if(df.Nj$Ni[1] ==0) next
    if(df.Nj$Ni[2] >1000 | df.Nj$Nj[2] >1000) next
    df.Nj.glm <- as.data.frame(summary(glm(formula = Nj ~ time, data = df.Nj,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Nj"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Nj = intercept.Nj + 100*time.Nj,
             com.comp.prob.Nj = case_when(com.comp.Nj>0~T,
                                          com.comp.Nj<0~F),
             sim = i,
             function.int = function.int)
    
    # Slope of Nj with init.cond= c(1,1)
    df.Ni.both <- df.sim[which(df.sim$sim == i &
                            df.sim$function.int == function.int &
                            df.sim$invader =="both"),]
    if(sum(is.na(df.Ni$Ni))>1) next
    if(df.Ni$Nj[1] ==0) next
    if(df.Ni$Ni[2] >1000 | df.Ni$Nj[2] >1000) next
    df.Ni.both.glm <- as.data.frame(summary(glm(formula = Ni ~ time, data = df.Ni,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Ni"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Ni.both = intercept.Ni + 100*time.Ni,
             com.comp.prob.Ni.both = case_when(com.comp.Ni>0~T,
                                          com.comp.Ni<0~F),
             sim = i,
             function.int = function.int)
    
    # Slope of Nj with init.cond= c(1,1)
    
    df.Nj.both <- df.sim[which(df.sim$sim == i &
                            df.sim$function.int == function.int &
                            df.sim$invader =="both"),]
    if(sum(is.na(df.Nj$Nj))>1) next
    if(df.Nj$Ni[1] ==0) next
    if(df.Nj$Ni[2] >1000 | df.Nj$Nj[2] >1000) next
    df.Nj.both.glm <- as.data.frame(summary(glm(formula = Nj ~ time, data = df.Nj,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Nj"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Nj.both = intercept.Nj + 100*time.Nj,
             com.comp.prob.Nj.both = case_when(com.comp.Nj>0~T,
                                          com.comp.Nj<0~F),
             sim = i,
             function.int = function.int)
    
    
    df.glm <- right_join(df.Ni.glm,df.Nj.glm,
                         df.Ni.both.glm,df.Nj.both.glm) %>%
      mutate(com.comp.coex = case_when(mean(com.comp.Ni,com.comp.Ni.both) > 0 &
                                         mean(com.comp.Nj,com.comp.Nj.both) > 0 ~ "j_i",
                                       mean(com.comp.Ni,com.comp.Ni.both) <= 0 & 
                                         mean(com.comp.Nj,com.comp.Nj.both) > 0 ~ "j",
                                       mean(com.comp.Ni,com.comp.Ni.both) > 0 & 
                                         mean(com.comp.Nj,com.comp.Nj.both) <= 0 ~ "i",
                                       mean(com.comp.Ni,com.comp.Ni.both) <= 0 & 
                                         mean(com.comp.Nj,com.comp.Nj.both) <= 0 ~ "-1"))
    
    df.glm_all <- bind_rows(df.glm_all,df.glm )
  }
}
head(df.glm_all)
str(df.glm_all)

ggplot(df.glm_all,aes(x=time.Ni, 
                      y=time.Nj,color=as.factor(function.int)))+ geom_point(alpha=0.5) + theme_bw() 

ggsave(last_plot(),
       file = "figures/com.comp.SLOPE.pdf")

##########################################################################################################
# 2. Standardization
##########################################################################################################
df.sim.std <- full_join(df.sim,df.glm_all)

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


model.0 <- glm(com.comp.coex ~ as.factor(function.int), 
               data = df.sim.std[which(df.sim.std$prob.coex >= 0),],
               family = "binomial")

model.1 <- glm(formula(paste0("com.comp.coex  ~ ",
                              paste(paste0("a_initial",c(".i.j",".i.i",".j.i",".j.j")),collapse = "+"))), 
               data = subset(df.sim.std,function.int==1 & prob.coex >= 0),
               family = "binomial")

model.2 <- glm(formula(paste0("com.comp.coex  ~ ",
                              paste(c(paste0("Nmax",c(".i.j",".i.i",".j.i",".j.j")),
                                      paste0("a_initial",c(".i.j",".i.i",".j.i",".j.j")),
                                      paste0("a_slope",c(".i.j",".i.i",".j.i",".j.j"))),collapse = "+"))), 
               data = subset(df.sim.std,function.int==2 & prob.coex >= 0),
               family = "binomial")

model.3 <- glm(formula(paste0("com.comp.coex  ~ ",
                              paste(stand.variable,collapse = "+"))), 
               data = subset(df.sim.std,function.int==3 & prob.coex >= 0),
               family = "binomial")

model.4 <- glm(formula(paste0("com.comp.coex  ~ ",
                              paste(stand.variable,collapse = "+"))), 
               data = subset(df.sim.std,function.int==4 & prob.coex >= 0),
               family = "binomial")

# build dataframe to plot
ggplot( df.sim.std, aes(x=com.comp.coex, 
                                 fill=as.factor(function.int)))+
  stat_count() +theme_bw()


sens_out <- tidy(model.1) %>% 
  mutate(function.int = 1) %>% 
  bind_rows(
    (tidy(model.0) %>%
       mutate(term = "function",
              function.int = c(1,2,3,4))),
    (tidy(model.2) %>%
       mutate(function.int = 2)),
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

ggplot(df, aes(function.int, prediction)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.25, size = 1, position = position_dodge(width = 0.4)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.4)) +
  theme_light(base_size = 16) +
  scale_y_continuous(name = "Probability of coexistence", limits = c(0, 1),
                     labels = scales::percent)

