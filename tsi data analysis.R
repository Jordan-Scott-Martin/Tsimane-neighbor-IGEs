#prep workspace############################################

library(shinystan)
library(cmdstanr)
library(rstan)
library(ggplot2)
library(tidybayes)
library(cowplot)
library(posterior)

# setwd("...")
stan.df = readRDS("tsi_fert_data.RDS")

#standardize variables
std.v = function(x) {return( (x - mean(na.omit(x)))/sd(na.omit(x)) ) }
stan.df$age = std.v(stan.df$age)
stan.df$r = std.v(stan.df$mean_r)
stan.df$yearobs = std.v(stan.df$yearobs) 

#no NAs in Stan, missing values for neighbors set to -99
#these arbitrary values are ignored during model estimation
stan.df$r[is.na(stan.df$r)] = -99
stan.df$mean_r[is.na(stan.df$mean_r)] = -99

#Estimate models in Stan##################################

#set path for cmdstan
# set_cmdstan_path("...")

#unadjusted heritability model (no spouse or neighbor effects)
mod_h2 = cmdstan_model(stan_file = "fertility unadj model.stan",
                        stanc_options = list("O1"))


if (!file.exists("fit_tsi_h2.RDS")) {

  fit = mod_h2$sample(
    data = stan.df,
    iter_sampling = 1000,
    iter_warmup = 1500,
    init = 1e-4,
    chains = 6,  
    parallel_chains = 6,
    adapt_delta = 0.99,
    refresh = 20,
    seed = 9)

  fit$save_object(file = "fit_tsi_h2.RDS")

}

fit_h2 = readRDS("fit_tsi_h2.RDS")

#calculate h2 unadjusted for social effects
post = posterior::as_draws_rvars(fit_h2$draws())
median(post$h2); quantile(post$h2, c(0.05, 0.95))

#average IGE model
mod_avg = cmdstan_model(stan_file = "fertility ige model_neighbors_avg.stan",
                       stanc_options = list("O1"))

if (!file.exists("fit_tsi_avg.RDS")) {

  fit = mod_avg$sample(
    data = stan.df,
    iter_sampling = 1000,
    iter_warmup = 1500,
    init = 1e-4,
    chains = 6,  
    parallel_chains = 6,
    adapt_delta = 0.99,
    refresh = 20,
    seed = 9) 


  fit$save_object(file = "fit_tsi_avg.RDS")

}

fit_avg = readRDS("fit_tsi_avg.RDS")

#extract posterior mean and SD
#priors are important for identification and
#effective sampling of fluctuating selection models
post = posterior::as_draws_rvars(fit_avg$draws())
stan.df$mean_bsw = mean(draws_of(post$b_SW))
stan.df$sd_bsw = sd(draws_of(post$b_SW))
stan.df$mean_h2 = mean(draws_of(post$h2))
stan.df$sd_h2 = sd(draws_of(post$h2))

#community effect model
mod_c = cmdstan_model(stan_file = "fertility ige model_neighbors_c.stan",
                     stanc_options = list("O1"))

if (!file.exists("fit_tsi_c.RDS")) {
  fit = mod_c$sample(
    data = stan.df,
    iter_sampling = 1000,
    iter_warmup = 1500,
    init = 1e-4,
    chains = 6,  
    parallel_chains = 6,
    adapt_delta = 0.99,
    refresh = 20) 

  fit$save_object(file = "fit_tsi_c.RDS")
}

fit_c = readRDS("fit_tsi_c.RDS")

#neighborhood effect (frequency- and density-dependent) model
mod_fd = cmdstan_model(stan_file = "fertility ige model_neighbors_fd.stan",
                       stanc_options = list("O1"))

if (!file.exists("fit_tsi_fd.RDS")) {
  fit = mod_fd$sample(
    data = stan.df,
    iter_sampling = 1000,
    iter_warmup = 1500,
    init = 1e-4,
    chains = 6,  
    parallel_chains = 6,
    adapt_delta = 0.99,
    refresh = 20,
    seed = 9) 

  fit$save_object(file = "fit_tsi_fd.RDS")
}

fit_fd = readRDS("fit_tsi_fd.RDS")

#summarize and plot results###############################

fit = readRDS("fit_tsi_avg.RDS")
fit2 = readRDS("fit_tsi_c.RDS")
fit3 = readRDS("fit_tsi_fd.RDS")
post = posterior::as_draws_rvars(fit$draws())
post2 = posterior::as_draws_rvars(fit2$draws())
post3 = posterior::as_draws_rvars(fit3$draws())

#plot environmental effects####
sd_a = draws_of(sqrt(post$b_2^2 + post$b_3^2))
sd_M = draws_of(post$sd_M)
sd_Y = draws_of(post$sd_Y)
sd_C = draws_of(post$alpha)
sd_F = draws_of(post$sd_F)
sd_sig = draws_of(post$sd_E)
de = cbind(sd_a,sd_M,sd_Y,sd_C,sd_F,sd_sig)
colnames(de) = c("Age", "Maternal", "Birth years", "Spatial (community)", "Father / spouse", "Residual")
del = reshape2::melt(de)

dep = 
  ggplot(data = del, aes(x = value, group = Var2)) +
  stat_density(aes(y = after_stat(scaled), fill = Var2, color = Var2), alpha = 0.5, linewidth = 1.25)+
  scale_x_continuous(expand = c(0,0))+
  labs(x = "\n Standard deviation", y = "")+
  facet_wrap(.~ Var2)+
  coord_cartesian(xlim=c(0,1.1))+
  theme(legend.title = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color=NA, linewidth=1, linetype="solid"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(aes(color = "none", fill = "none"))
ggsave("de_plot.png", dep, width = 6, height = 4, units = "in")

#DGEs and heritability####
V_D = post$sd_G^2 #mother col 1, father col 2
median(V_D); quantile(V_D, c(0.05, 0.95))
median(post$h2); quantile(post$h2, c(0.05, 0.95))

#average IGE effects
beta_SW = post$b_SW #social selection
mean_n = mean(stan.df$social_n[stan.df$social_n>0])
V_I = V_D * beta_SW^2 * mean_n^2 #total IGE
V_Im = V_D * beta_SW^2 #marginal IGE
cov_DI = (V_D * beta_SW) #cov DGE-IGE

median(V_I); quantile(V_I, c(0.05, 0.95))
median(V_I - V_D); quantile(V_I - V_D, c(0.05, 0.95)); sum(V_I-V_D>0)/length(V_I)
median(V_Im); quantile(V_Im, c(0.05, 0.95))
median(V_Im - V_D); quantile(V_Im - V_D, c(0.05, 0.95))
median(cov_DI); quantile(cov_DI, c(0.05, 0.95)); sum(cov_DI>0)/length(cov_DI)
median(post$b_SW); quantile(post$b_SW, c(0.05, 0.95))

#total evolvability
mean_r = mean(stan.df$mean_r[stan.df$social_n>0]) 
mean_n = mean(stan.df$social_n[stan.df$social_n>0])
e_w =  ((V_D + mean_n * cov_DI) + #total IGE
          mean_r * (V_I + mean_n * cov_DI)) /post$W_0^2
e_0 = V_D / post$W_0^2 #non-social evolvability

median(e_w); quantile(e_w, c(0.05, 0.95))
median(e_0); quantile(e_0, c(0.05, 0.95))
diff_e = e_w / e_0
median(diff_e); quantile(diff_e, c(0.05, 0.95)); sum(diff_e>1)/length(diff_e)

#inclusive fitness evolvability
e_wIF =  (V_D + mean_r * (V_I + 2*mean_n * cov_DI)) /post$W_0^2
diff_es = e_wIF / e_0
delta_eIF = e_w - e_wIF

median(e_wIF); quantile(e_wIF, c(0.05, 0.95))
median(diff_es); quantile(diff_es, c(0.05, 0.95))
median(delta_eIF); quantile(delta_eIF, c(0.05, 0.95))

#fluctuating selection
median(post2$sd_Cb^2 * post2$sd_G^2); quantile(post2$sd_Cb^2 * post2$sd_G^2,c(0.05,0.95))
cov_c = list()
for(i in 1:nrow(draws_of(post2$W_Cb))){
  cov_c[[i]] = draws_of(post2$sd_G)[i]^2 * draws_of(post2$W_Cb)[i,]
}
cov_cv = do.call(rbind.data.frame, cov_c)
colnames(cov_cv) = 1:82
range(apply(cov_cv, 2, median))

median(post3$b_D); quantile(post3$b_D, c(0.05,0.95))
sum(post3$b_D<0)/length(post3$b_D) 
median(post3$b_I); quantile(post3$b_I, c(0.05,0.95))
sum(post3$b_I)/length(post3$b_I)
median(post3$b_Ir); quantile(post3$b_Ir, c(0.05,0.95))
sum(post3$b_Ir>0)/length(post3$b_Ir)

#plot (co)variance####
df = data.frame(delta = c(draws_of(V_D), draws_of(V_Im), draws_of(V_I), draws_of(cov_DI)), 
                mod = rep(c("a","b","c", "d"), each = length(draws_of(V_D))))
var.p = 
      ggplot(data = df,
      aes(x = delta, y = after_stat(scaled), group = mod, color = mod, fill = mod)) +
      geom_density(linewidth = 1.25, alpha = 0.75) + facet_wrap(. ~ mod, scales = "free")+
      scale_fill_manual(values = c("#58e0ae","#5dabc9","#5d5dc9","#cf8ade"))+
      scale_color_manual(values = c("#42856c","#427487", "#494287", "#764082"))+
      labs(x = "\n Quantitative genetic effects on fertility")+
      theme(legend.title = element_text(face = "bold"),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.line.x = element_line(),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"))+
      guides(aes(color = "none", fill = "none"))+
      geom_vline(data=subset(df, mod=="a"), aes(xintercept = median(V_D)), lwd = 0.5, lty = "dashed")+
      geom_vline(data=subset(df, mod=="b"), aes(xintercept = median(V_Im)), lwd = 0.5, lty = "dashed")+
      geom_vline(data=subset(df, mod=="c"), aes(xintercept = median(V_I)), lwd = 0.5, lty = "dashed")+
      geom_vline(data=subset(df, mod=="d"), aes(xintercept = median(cov_DI)),lwd = 0.5, lty = "dashed")

ggsave("var_plot.png", var.p, width = 5, height = 5, units = "in")


#plot social selection####
df = data.frame(delta = draws_of(beta_SW))

beta_sd = 
      ggplot(data = df,
      aes(x = delta, y = after_stat(scaled))) +
      geom_density(size = 1.25, alpha = 0.75, color = "darkblue", fill = "lightblue") + 
      labs(x = "\n Selection")+
      theme(legend.title = element_text(face = "bold"),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.line.x = element_line(),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"))+
      guides(aes(color = "none", fill = "none"))+
      geom_vline(data=df, aes(xintercept = median(beta_SW)), lwd = 0.5, lty = "dashed")

ggsave("beta_plot.png", beta_sd, width = 4, height = 3, units = "in")


#plot evolvability####
df = data.frame(delta = c(draws_of(e_0), draws_of(e_wIF), draws_of(e_w)), 
                type = rep(c("a","b","c"), each = length(draws_of(e_0))))

evolv = 
      ggplot(data = df,
      aes(x = delta, y = after_stat(scaled), group = type)) +
      geom_density(size = 1.25, alpha = 0.55, aes(color = type, fill = type))+
      scale_fill_manual(values = c("#58e0ae","#5dabc9","#5dc9c7"))+
      scale_color_manual(values = c("#42856c","#427487","#428784"))+
      labs(x = "")+
      coord_cartesian(xlim=c(0,0.039))+
      facet_wrap(.~ type, ncol = 1)+
      theme(legend.title = element_text(face = "bold"),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.line.x = element_line(),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #strip.text.x = element_blank(),
        panel.spacing = unit(2, "lines"))+
      guides(aes(color = "none", fill = "none"))+
      geom_vline(data=subset(df, type=="a"), aes(xintercept = median(e_0)), lwd = 0.5, lty = "dashed")+
      geom_vline(data=subset(df, type=="b"), aes(xintercept = median(e_wIF)), lwd = 0.5, lty = "dashed")+
      geom_vline(data=subset(df, type=="c"), aes(xintercept = median(e_w)), lwd = 0.5, lty = "dashed")

ggsave("evolv_plot.png", evolv, width = 3, height = 4, units = "in")


#plot community fluctuations####
cov_c = list()
for(i in 1:nrow(draws_of(post2$W_Cb))){
  cov_c[[i]] = draws_of(post2$sd_G)[i]^2 * draws_of(post2$W_Cb)[i,]
}
cov_cv = do.call(rbind.data.frame, cov_c)
colnames(cov_cv) = 1:stan.df$C

levels=10
wcbl = reshape2::melt(cov_cv)
com.p =
ggplot(wcbl, aes(x = value, y = variable))+
  geom_vline(xintercept = 0, lty = "dashed", linewidth = 1, color = "black")+
  geom_vline(xintercept = median(post2$b_SW_0 * post2$sd_G^2), alpha = 0.25, lty = "solid", linewidth = 1, color = "maroon4")+
  stat_pointinterval(linewidth = 7, .width = ppoints(levels), alpha=1/levels, color = "violetred2")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_discrete(expand = c(0.02,0.02))+
  xlab("cov(Wd, Wi)")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12,face="bold"),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.line = element_line(linewidth = 1),
        panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = "none", color = "none")
ggsave("com_plot.png", com.p, width = 3, height = 8, units = "in")

#density-dependence####
dp1 = 
  ggplot(data = data.frame(draws_of(post3$b_D)), aes(x = draws_of(post3$b_D))) +
  stat_density(aes(y = after_stat(scaled)), size = 1.25, color = "springgreen3", fill = "lightgreen")+
  geom_vline(xintercept = median(draws_of(post3$b_D)), lty = "dashed")+
  scale_x_continuous(expand = c(0,0))+
  labs(x = "\n beta_D", y = "\n cov(Wd,Wi)")+
  coord_cartesian(xlim=c(-0.3,0.3))+
  theme(legend.title = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color=NA, linewidth=1, linetype="solid"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(aes(color = "none", fill = "none"))
ggsave("bd_plot.png", dp1, width = 4, height = 4, units = "in")

seq = 1:21
beta_d = list()
for(i in 1:length(draws_of(post3$b_D))){
  beta_d[[i]] = (draws_of(post3$b_SW_0)[i] + draws_of(post3$b_D)[i] * seq) * draws_of(post3$sd_G)[i]^2
  }
bd = do.call("rbind", beta_d)
bd = reshape2::melt(bd)
levels=10
dp2 = 
      ggplot(data = bd, aes(x = Var2, y = value)) +
      stat_lineribbon(size=2, .width = ppoints(levels), alpha=1/levels, 
                      color = "springgreen3", fill = "springgreen" )+
      geom_hline(yintercept = 0, lty = "dashed")+
      scale_x_continuous(expand = c(0,0))+
      labs(x = "\n # neighbors", y = "\n cov(Wd,Wi)")+
      #coord_cartesian(ylim=c(-2.25,4.25))+
      theme(legend.title = element_text(face = "bold"),
      panel.border=element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"))+
      guides(aes(color = "none", fill = "none"))
ggsave("bd2_plot.png", dp2, width = 4, height = 4, units = "in")

#frequency-dependence####
dff = data.frame(bi = draws_of(post3$b_I), bir = draws_of(post3$b_Ir))
dff = reshape2::melt(dff)
dfm = data.frame(bi = median(post3$b_I), bir = median(post3$b_Ir))
dfm = reshape2::melt(dfm)

df1 = 
  ggplot(data = dff, aes(x = value, group = variable, color = variable, fill = variable)) +
  stat_density(aes(y = after_stat(scaled)), size = 1.25)+
  scale_x_continuous(expand = c(0,0))+
  scale_color_manual(values = c("royalblue4","darkmagenta"))+
  scale_fill_manual(values = c("royalblue","orchid"))+
  facet_wrap(.~ variable)+
    geom_vline(data = dfm, mapping = aes(xintercept = value), lty = "dashed")+
  labs(x = "\n beta_I", y = "\n cov(Wd,Wi)")+
  coord_cartesian(xlim=c(-3,4.5))+
  theme(legend.title = element_text(face = "bold"),
        panel.border=element_rect(fill=NA,color=NA, linewidth=1, linetype="solid"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(aes(color = "none", fill = "none"))
ggsave("bf_plot.png", df1, width = 4, height = 4, units = "in")

seq = seq(-1, 1, by = 0.05)
beta_i = list()
beta_ir = list()
for(i in 1:length(draws_of(post3$b_I))){
  beta_i[[i]] = (draws_of(post$b_SW)[i] + draws_of(post3$b_I)[i] * seq + draws_of(post3$b_Ir)[i] * seq * 0) * draws_of(post3$sd_G)[i]^2
  beta_ir[[i]] = (draws_of(post$b_SW)[i] + draws_of(post3$b_I)[i] * seq + draws_of(post3$b_Ir)[i] * seq * 1.79) * draws_of(post3$sd_G)[i]^2
}
bi = data.frame(do.call("rbind", beta_i))
bir = data.frame(do.call("rbind", beta_ir))
colnames(bi) = seq
colnames(bir) = seq
bi = reshape2::melt(bi)
bir = reshape2::melt(bir)
bi$variable = as.numeric(as.character(bi$variable))
bir$variable = as.numeric(as.character(bir$variable))
bi$type = "cousin"
bir$type = "sister"
bic = rbind(bi,bir)

df2 = 
      ggplot(data = bic,
      aes(x = variable, y = value, color = type, fill = type)) +
        stat_lineribbon(size=2, .width = ppoints(levels), alpha=1/levels)+
        geom_hline(yintercept = 0, lty = "dashed")+
        labs(x = "\n focal * neighbor DGE", y = "\n cov(Wd,Wi)")+
        scale_color_manual(values = c("royalblue4","darkmagenta"))+
        scale_fill_manual(values = c("royalblue","orchid"))+
        scale_x_continuous(expand=c(0,0))+
        coord_cartesian(xlim=c(-1,1))+
        theme(legend.title = element_text(face = "bold"),
              panel.border=element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
              axis.title.y = element_text(face = "bold"),
              axis.title.x = element_text(face = "bold"),
              axis.line.y = element_line(),
              axis.line.x = element_line(),
              panel.background= element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              panel.spacing = unit(2, "lines"))+
        guides(aes(color = "none", fill = "none"))
ggsave("bf2_plot.png", df2, width = 4, height = 4, units = "in")

#combine####
dpc = plot_grid(dp1,dp2, ncol = 1)
dfc = plot_grid(df1,df2, ncol = 1)
dc = plot_grid(com.p,dpc,dfc, nrow = 1, ncol = 3)
save_plot("dpc_plot.png",dc, base_width = 10, base_height = 5, units = "in")

