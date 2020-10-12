## Load libraries..
library(lme4)

## Define transparency function for plotting..
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
  }

## Load data..
xdata = read.table("blue_tit_data_updated_2020-04-18.csv", sep = ",", head = T, stringsAsFactors= F, na.strings= ".") # data table is available at 

{## Prepare data for "Tarsus" and "Re-capture" models..

  ## Exclude rows with NAs in any of the relevant variables:
  xx = complete.cases(xdata[, c("day_14_tarsus_length", "net_rearing_manipulation", "rear_d0_rear_nest_brood_size", "hatch_year", "rear_nest_OH", "rear_area", "home_or_away", "chick_sex_molec", "hatch_Area", "hatch_mom_Ring", "hatch_nest_breed_ID", "rear_nest_breed_ID", "genetic_dad_ring_.WP_or_EP.")])
  nrow(xdata)
  ydata=droplevels(xdata[xx, ])
  nrow(ydata)

  ## Standardize covariates..
  ydata$z_man = as.vector(scale(ydata$net_rearing_manipulation))
  ydata$z_bs0 = as.vector(scale(ydata$rear_d0_rear_nest_brood_size)) # It is now the number of siblings in the *rearing* nest..
  ydata$z_nh = as.vector(scale(ydata$rear_nest_OH)) # NOTE that rear_nest_OH and hatch_nest_OH are essentially the same..

  ## Manually dummy code and center factor(s):
  ydata$f_year = as.factor(ydata$hatch_year)
  ydata$f_year.2002 = as.numeric(ydata$f_year=="2002")
  ydata$f_year.2003 = as.numeric(ydata$f_year=="2003")
  ydata$f_year.2002 = ydata$f_year.2002-mean(ydata$f_year.2002)
  ydata$f_year.2003 = ydata$f_year.2003-mean(ydata$f_year.2003)
  ydata$f_sexM = ydata$chick_sex_molec - 1
  ydata$f_sexM = ydata$f_sexM - mean(ydata$f_sexM)
  ydata$f_ha = ydata$home_or_away - 1
  ydata$f_ha = ydata$f_ha - mean(ydata$f_ha)
  }

{## Check which random effects are needed..
source("Diagnostic_fcns.r")

feretab = fe.re.tab(
  fe.model = "day_14_tarsus_length ~ 

    I(z_man^2):z_bs0 + I(z_man^2):f_year + I(z_man^2):f_sexM + 
    z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + f_sexM:z_man + 
    I(z_man^2) + 
    z_man + z_bs0 + f_ha + 
    
    z_nh:f_year + I(z_nh^2):f_year + 
    f_year + z_nh + f_sexM + I(z_nh^2) ",
    
  re = c("hatch_nest_breed_ID", "rear_nest_breed_ID", "hatch_mom_Ring", "genetic_dad_ring_.WP_or_EP.", "rear_area"),

  data = ydata
  )
  }

{## Fit "Tarsus model", starting from a full model comprising ALL possible random slopes according to fe.re.tab (only TEST variables..)
  ## We have kept all slopes for main terms where a consistent fraction of the levels of the RE (at least ca. 1/5) contained at least 2 unique values (levels) or at least 2 x 2 unique values for interactions..
  
  tars_full_k = lmer(day_14_tarsus_length ~ 

    I(z_man^2):z_bs0 + I(z_man^2):f_year + I(z_man^2):f_sexM + # 2-way interactions of z_man^2 (test)
    z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + f_sexM:z_man + f_sexM:z_bs0 + # 2-way interactions (test)
    I(z_man^2) + # squared term for z_man (test)
    z_man + z_bs0 + f_ha + # other terms (test)
    
    z_nh:f_year + I(z_nh^2):f_year + # 2-way interactions (control)
    f_year + z_nh + f_sexM + I(z_nh^2) + # other terms (control)
    
    # Random intercepts..
    (1|hatch_nest_breed_ID) + # random intercept of hatching nest
    (1|rear_nest_breed_ID) + # random intercept of rearing nest
    (1|genetic_dad_ring_.WP_or_EP.) + # random intercept of dad (NOTE that we are ignoring the rearing parents, as they are already in the id of rearing nest, and there are very few cases of mother or fathers rearing multiple nests..)
    (1|rear_area) + # random intercept for rearing area (where the nest is)  
    
    # Random slopes..
    # Main terms, test
    (0 + I(z_man^2)|hatch_nest_breed_ID) +
    (0 + I(z_man^2)|genetic_dad_ring_.WP_or_EP.) +
    (0 + I(z_man^2)|rear_area) +
    (0 + z_man|hatch_nest_breed_ID) +
    (0 + z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_man|rear_area) +
    (0 + z_bs0|hatch_nest_breed_ID) +
    (0 + z_bs0|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0|rear_area) +
    
    # Interactions, test
    (0 + z_bs0:z_man|hatch_nest_breed_ID) +
    (0 + z_bs0:z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0:z_man|rear_area) +  
    (0 + I(z_man^2):(f_year.2002+f_year.2003)||rear_area) +  
    (0 + z_man:(f_year.2002+f_year.2003)||rear_area) +
    (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) +
    (0 + f_sexM:z_man|hatch_nest_breed_ID) +
    (0 + f_sexM:z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + f_sexM:z_man|rear_area),

    data = ydata, control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  tars_null_k = lmer(day_14_tarsus_length ~ 
        
    z_nh:f_year + I(z_nh^2):f_year + # 2-way interactions (control)
    f_year + z_nh + f_sexM + I(z_nh^2) + # other terms (control)
    
    # Random intercepts..
    (1|hatch_nest_breed_ID) + # random intercept of hatching nest
    (1|rear_nest_breed_ID) + # random intercept of rearing nest
    (1|genetic_dad_ring_.WP_or_EP.) + # random intercept of dad (NOTE that we are ignoring the rearing parents, as they are already in the id of rearing nest, and there are very few cases of mother or fathers rearing multiple nests..)
    (1|rear_area) + # random intercept for rearing area (where the nest is)  
    
    # Random slopes..
    # Main terms, test
    (0 + I(z_man^2)|hatch_nest_breed_ID) +
    (0 + I(z_man^2)|genetic_dad_ring_.WP_or_EP.) +
    (0 + I(z_man^2)|rear_area) +
    (0 + z_man|hatch_nest_breed_ID) +
    (0 + z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_man|rear_area) +
    (0 + z_bs0|hatch_nest_breed_ID) +
    (0 + z_bs0|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0|rear_area) +
    
    # Interactions, test
    (0 + z_bs0:z_man|hatch_nest_breed_ID) +
    (0 + z_bs0:z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0:z_man|rear_area) +  
    (0 + I(z_man^2):(f_year.2002+f_year.2003)||rear_area) +  
    (0 + z_man:(f_year.2002+f_year.2003)||rear_area) +
    (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) +
    (0 + f_sexM:z_man|hatch_nest_breed_ID) +
    (0 + f_sexM:z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + f_sexM:z_man|rear_area),

    data = ydata, control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
  
  anova(tars_full_k, tars_null_k, test = "Chisq")
  #             npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)   
  # tars_null_k   36 2873.8 3072.8 -1400.9   2801.8                       
  # tars_full_k   51 2866.9 3148.7 -1382.5   2764.9 36.94 15   0.001291 **

  drop1(tars_full_k, test = "Chisq")
  #                   npar    AIC    LRT Pr(Chi)
  # <none>                 2866.9               
  # f_ha                 1 2866.2 1.3335  0.2482
  # I(z_man^2):z_bs0     1 2864.9 0.0354  0.8507
  # I(z_man^2):f_year    2 2866.1 3.1909  0.2028
  # I(z_man^2):f_sexM    1 2864.9 0.0016  0.9685
  # z_bs0:z_man          1 2864.9 0.0186  0.8914
  # f_year:z_man         2 2863.6 0.6727  0.7144
  # z_bs0:f_year         2 2866.4 3.4775  0.1757
  # f_sexM:z_man         1 2865.0 0.0466  0.8291
  # z_bs0:f_sexM         1 2866.3 1.3921  0.2380
  # f_year:z_nh          2 2866.8 3.9037  0.1420
  # f_year:I(z_nh^2)     2 2863.5 0.6308  0.7295

  tars_red1_k = lmer(day_14_tarsus_length ~ 

    # z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + f_sexM:z_man + f_sexM:z_bs0 + # 2-way interactions (test)
    I(z_man^2) + # squared term for z_man (test)
    z_man + z_bs0 + f_ha + # other terms (test)
    
    # z_nh:f_year + I(z_nh^2):f_year + # 2-way interactions (control)
    f_year + z_nh + f_sexM + I(z_nh^2) + # other terms (control)
    
    # Random intercepts..
    (1|hatch_nest_breed_ID) + # random intercept of hatching nest
    (1|rear_nest_breed_ID) + # random intercept of rearing nest
    (1|genetic_dad_ring_.WP_or_EP.) + # random intercept of dad (NOTE that we are ignoring the rearing parents, as they are already in the id of rearing nest, and there are very few cases of mother or fathers rearing multiple nests..)
    (1|rear_area) + # random intercept for rearing area (where the nest is)  
    
    # Random slopes..
    # Main terms, test
    (0 + I(z_man^2)|hatch_nest_breed_ID) +
    (0 + I(z_man^2)|genetic_dad_ring_.WP_or_EP.) +
    (0 + I(z_man^2)|rear_area) +
    (0 + z_man|hatch_nest_breed_ID) +
    (0 + z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_man|rear_area) +
    (0 + z_bs0|hatch_nest_breed_ID) +
    (0 + z_bs0|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0|rear_area),
    
    # Interactions, test
    # (0 + z_bs0:z_man|hatch_nest_breed_ID) +
    # (0 + z_bs0:z_man|genetic_dad_ring_.WP_or_EP.) +
    # (0 + z_bs0:z_man|rear_area) +  
    # (0 + I(z_man^2):(f_year.2002+f_year.2003)||rear_area) +  
    # (0 + z_man:(f_year.2002+f_year.2003)||rear_area) +
    # (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) +
    # (0 + f_sexM:z_man|hatch_nest_breed_ID) +
    # (0 + f_sexM:z_man|genetic_dad_ring_.WP_or_EP.) +
    # (0 + f_sexM:z_man|rear_area),
    
    data = ydata, control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  drop1(tars_red1_k, test = "Chisq")
  #            npar    AIC    LRT   Pr(Chi)    
  # <none>          2836.8                     
  # I(z_man^2)    1 2835.4   0.57  0.449784    
  # z_man         1 2851.6  16.79 4.173e-05 ***
  # z_bs0         1 2842.6   7.72  0.005455 ** 
  # f_ha          1 2836.2   1.35  0.245127    
  # f_year        2 2846.1  13.30  0.001293 ** 
  # z_nh          1 2835.6   0.79  0.375340    
  # f_sexM        1 3238.2 403.32 < 2.2e-16 ***
  # I(z_nh^2)     1 2837.7   2.83  0.092495 . 
  
  tars_red2_k = lmer(day_14_tarsus_length ~ 

    # z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + f_sexM:z_man + f_sexM:z_bs0 + # 2-way interactions (test)
    # I(z_man^2) + # squared term for z_man (test)
    z_man + z_bs0 + f_ha + # other terms (test)
    
    # z_nh:f_year + I(z_nh^2):f_year + # 2-way interactions (control)
    f_year + z_nh + f_sexM + I(z_nh^2) + # other terms (control)
    
    # Random intercepts..
    (1|hatch_nest_breed_ID) + # random intercept of hatching nest
    (1|rear_nest_breed_ID) + # random intercept of rearing nest
    (1|genetic_dad_ring_.WP_or_EP.) + # random intercept of dad (NOTE that we are ignoring the rearing parents, as they are already in the id of rearing nest, and there are very few cases of mother or fathers rearing multiple nests..)
    (1|rear_area) + # random intercept for rearing area (where the nest is)  
    
    # Random slopes..
    # Main terms, test
    # (0 + I(z_man^2)|hatch_nest_breed_ID) +
    # (0 + I(z_man^2)|genetic_dad_ring_.WP_or_EP.) +
    # (0 + I(z_man^2)|rear_area) +
    (0 + z_man|hatch_nest_breed_ID) +
    (0 + z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_man|rear_area) +
    (0 + z_bs0|hatch_nest_breed_ID) +
    (0 + z_bs0|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0|rear_area),
    
    # Interactions, test
    # (0 + z_bs0:z_man|hatch_nest_breed_ID) +
    # (0 + z_bs0:z_man|genetic_dad_ring_.WP_or_EP.) +
    # (0 + z_bs0:z_man|rear_area) +  
    # (0 + I(z_man^2):(f_year.2002+f_year.2003)||rear_area) +  
    # (0 + z_man:(f_year.2002+f_year.2003)||rear_area) +
    # (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) +
    # (0 + f_sexM:z_man|hatch_nest_breed_ID) +
    # (0 + f_sexM:z_man|genetic_dad_ring_.WP_or_EP.) +
    # (0 + f_sexM:z_man|rear_area),

    data = ydata, control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  drop1(tars_red2_k, test = "Chisq")
  #           npar    AIC    LRT   Pr(Chi)    
  # <none>         2829.4                     
  # z_man        1 2847.9  20.53 5.857e-06 ***
  # z_bs0        1 2834.6   7.18  0.007389 ** 
  # f_ha         1 2829.0   1.61  0.204436    
  # f_year       2 2838.3  12.86  0.001614 ** 
  # z_nh         1 2828.3   0.85  0.355340    
  # f_sexM       1 3230.8 403.38 < 2.2e-16 ***
  # I(z_nh^2)    1 2829.8   2.43  0.118977    
  
  
  tars_red3_k = lmer(day_14_tarsus_length ~ 

    # z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + f_sexM:z_man + f_sexM:z_bs0 + # 2-way interactions (test)
    # I(z_man^2) + # squared term for z_man (test)
    z_man + z_bs0 + f_ha + # other terms (test)
    
    # z_nh:f_year + I(z_nh^2):f_year + # 2-way interactions (control)
    f_year + z_nh + f_sexM + 
    # I(z_nh^2) + # other terms (control)
    
    # Random intercepts..
    (1|hatch_nest_breed_ID) + # random intercept of hatching nest
    (1|rear_nest_breed_ID) + # random intercept of rearing nest
    (1|genetic_dad_ring_.WP_or_EP.) + # random intercept of dad (NOTE that we are ignoring the rearing parents, as they are already in the id of rearing nest, and there are very few cases of mother or fathers rearing multiple nests..)
    (1|rear_area) + # random intercept for rearing area (where the nest is)  
    
    # Random slopes..
    # Main terms, test
    # (0 + I(z_man^2)|hatch_nest_breed_ID) +
    # (0 + I(z_man^2)|genetic_dad_ring_.WP_or_EP.) +
    # (0 + I(z_man^2)|rear_area) +
    (0 + z_man|hatch_nest_breed_ID) +
    (0 + z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_man|rear_area) +
    (0 + z_bs0|hatch_nest_breed_ID) +
    (0 + z_bs0|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0|rear_area),
    
    # Interactions, test
    # (0 + z_bs0:z_man|hatch_nest_breed_ID) +
    # (0 + z_bs0:z_man|genetic_dad_ring_.WP_or_EP.) +
    # (0 + z_bs0:z_man|rear_area) +  
    # (0 + I(z_man^2):(f_year.2002+f_year.2003)||rear_area) +  
    # (0 + z_man:(f_year.2002+f_year.2003)||rear_area) +
    # (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) +
    # (0 + f_sexM:z_man|hatch_nest_breed_ID) +
    # (0 + f_sexM:z_man|genetic_dad_ring_.WP_or_EP.) +
    # (0 + f_sexM:z_man|rear_area),

    data = ydata, control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
  
  drop1(tars_red3_k, test = "Chisq")
  #        npar    AIC    LRT   Pr(Chi)    
  # <none>      2829.8                     
  # z_man     1 2848.8  20.94 4.743e-06 ***
  # z_bs0     1 2834.6   6.77  0.009295 ** 
  # f_ha      1 2829.4   1.60  0.206252    
  # f_year    2 2836.4  10.52  0.005204 ** 
  # z_nh      1 2829.0   1.18  0.278199    
  # f_sexM    1 3230.4 402.54 < 2.2e-16 ***  
  
  summary(tars_red3_k)
  # Fixed effects:
  #             Estimate Std. Error t value
  # (Intercept) 16.90410    0.09251 182.731
  # z_man       -0.18592    0.01891  -9.834
  # z_bs0       -0.07944    0.02909  -2.731
  # f_ha         0.03449    0.02725   1.266
  # f_year2002  -0.29929    0.13341  -2.243
  # f_year2003  -0.32264    0.10113  -3.190
  # z_nh        -0.06019    0.05758  -1.045
  # f_sexM       0.48029    0.02244  21.404

  }

{## Plot "Tarsus model"..

  pdf(file = "tarsus_length.pdf", width = 10.5, height = 5.5)
  par(mfrow = c(1,2))

  fe = fixef(tars_red3_k)

  tt = table(ydata$rear_d0_rear_nest_brood_size)
  sd_bs = sd(ydata$rear_d0_rear_nest_brood_size)
  m_bs = mean(ydata$rear_d0_rear_nest_brood_size)

  modal_bs = as.numeric(names(tt)[which(tt == max(tt))])
  max_bs = max(ydata$rear_d0_rear_nest_brood_size)
  min_bs = min(ydata$rear_d0_rear_nest_brood_size)
  
  sd_man = sd(ydata$net_rearing_manipulation)
  m_man = mean(ydata$net_rearing_manipulation)

  year = "2003"
  
  # Females..
  plot(0,0, type = "n", xaxt = "n", xlim = (c(-5,5) - m_man)/sd_man, ylim = c(14, 18), xlab = "net manipulation", ylab = "tarsus length (mm)", main = "a. Females")
  
  abline(v = (0 - m_man)/sd_man)
  
  axis(1, at =  (-10:10 - m_man)/sd_man, labels = -10:10)
  
  legend("bottomleft", col = rgb((modal_bs - c(-6,-3,0,3,6) - min_bs)/(max_bs - min_bs), 1 - (modal_bs - c(-6,-3,0,3,6) - min_bs)/(max_bs - min_bs), 0), lwd = 2, legend = modal_bs - c(-6,-3,0,3,6), title = "BS_0", bty = "n")

  for(d in c(-6,-3,0,3,6)) {

      x = modal_bs + d

      x_z = (x - m_bs)/sd_bs 
      
      sel = ydata[which(ydata$rear_d0_rear_nest_brood_size >= x - 1 & ydata$rear_d0_rear_nest_brood_size <= x + 1 & ydata$chick_sex_molec == 1),]
      
      newdata_1 = data.frame(
        z_man = seq((min(sel$net_rearing_manipulation) - m_man)/sd_man, (max(sel$net_rearing_manipulation) - m_man)/sd_man, l=100), z_bs0 = mean(sel$z_bs0), f_ha = as.numeric(names(table(ydata$f_ha))[1]), f_year = year, z_nh = 0, f_sexM = unique(sel$f_sexM)
        )

      pred_1 = predict(tars_red3_k, newdata = newdata_1, re.form = NA)
      
      aa = aggregate(sel$day_14_tarsus_length, by = list(sel$z_man), FUN = mean)
      nn = aa$Group.1
      mm = aa$x
      cc = as.vector(table(sel$z_man))
      ss = aggregate(sel$day_14_tarsus_length, by = list(sel$z_man), FUN = sd)$x
      ll = mm - ss
      uu = mm + ss
      
      points(nn, mm, cex = sqrt(cc)/3, col = add.alpha(rgb((x-min_bs)/(max_bs-min_bs),1-(x-min_bs)/(max_bs-min_bs),0,1), 0.5), pch = 19)
      
      segments(x0 = nn, x1 = nn, y0 = ll, y1 = uu, col = rgb((x-min_bs)/(max_bs-min_bs),1-(x-min_bs)/(max_bs-min_bs),0,1))
            
      lines(newdata_1$z_man, pred_1, col = rgb((x-min_bs)/(max_bs-min_bs),1-(x-min_bs)/(max_bs-min_bs),0), lwd = 2)
      }

  # Males..
  plot(0,0, type = "n", xaxt = "n", xlim = (c(-5,5) - m_man)/sd_man, ylim = c(14, 18), xlab = "net manipulation", ylab = "tarsus length (mm)", main = "b. Males")
 
  abline(v = (0 - m_man)/sd_man)
  
  axis(1, at =  (-10:10 - m_man)/sd_man, labels = -10:10)
  
  legend("bottomleft", col = rgb((modal_bs - c(-6,-3,0,3,6) - min_bs)/(max_bs - min_bs), 1 - (modal_bs - c(-6,-3,0,3,6) - min_bs)/(max_bs - min_bs), 0), lwd = 2, legend = modal_bs - c(-6,-3,0,3,6), title = "BS_0", bty = "n")

  for(d in c(-6,-3,0,3,6)) {

    x = modal_bs + d

    x_z = (x - m_bs)/sd_bs 
    
    sel = ydata[which(ydata$rear_d0_rear_nest_brood_size >= x - 1 & ydata$rear_d0_rear_nest_brood_size <= x + 1 & ydata$chick_sex_molec == 2),]

    newdata_1 = data.frame(
      z_man = seq((min(sel$net_rearing_manipulation) - m_man)/sd_man, (max(sel$net_rearing_manipulation) - m_man)/sd_man, l=100), z_bs0 = mean(sel$z_bs0), f_ha = as.numeric(names(table(ydata$f_ha))[1]), f_year = year, z_nh = 0, f_sexM = unique(sel$f_sexM)
      )

    pred_1 = predict(tars_red2_k, newdata = newdata_1, re.form = NA)
    
    aa = aggregate(sel$day_14_tarsus_length, by = list(sel$z_man), FUN = mean)
    nn = aa$Group.1
    mm = aa$x
    cc = as.vector(table(sel$z_man))
    ss = aggregate(sel$day_14_tarsus_length, by = list(sel$z_man), FUN = sd)$x
    ll = mm - ss
    uu = mm + ss
    
    points(nn, mm, cex = sqrt(cc)/3, col = add.alpha(rgb((x-min_bs)/(max_bs-min_bs),1-(x-min_bs)/(max_bs-min_bs),0,1), 0.5), pch = 19)
    
    segments(x0 = nn, x1 = nn, y0 = ll, y1 = uu, col = rgb((x-min_bs)/(max_bs-min_bs),1-(x-min_bs)/(max_bs-min_bs),0,1))
        
    lines(newdata_1$z_man, pred_1, col = rgb((x-min_bs)/(max_bs-min_bs),1-(x-min_bs)/(max_bs-min_bs),0), lwd = 2)
    }

  dev.off()
  
  }
    
{## Fit "Re-capture" model..
  ## This is a full model comprising ALL possible random slopes according to fe.re.tab (only TEST variables..)
  ## We have kept all slopes for main terms where a consistent fraction of the levels of the RE (at least ca. 1/5) contained at least 2 unique values (levels) or at least 2 x 2 unique values for interactions..
  
  surv_full_k = glmer(chick_survival_to_first_breed_season ~ 

    z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + f_sexM:z_man + f_sexM:z_bs0 + # 2-way interactions (test)
    z_man + z_bs0 + f_ha + # other terms (test)
    
    z_nh:f_year + I(z_nh^2):f_year + # 2-way interactions (control)
    f_year + z_nh + f_sexM + I(z_nh^2) + # other terms (control)
    
    # Random intercepts..
    (1|hatch_nest_breed_ID) + # random intercept of hatching nest
    (1|rear_nest_breed_ID) + # random intercept of rearing nest
    (1|genetic_dad_ring_.WP_or_EP.) + # random intercept of dad (NOTE that we are ignoring the rearing parents, as they are already in the id of rearing nest, and there are very few cases of mother or fathers rearing multiple nests..)
    (1|rear_area) + # random intercept for rearing area (where the nest is)  
    
    # Random slopes..
    # Main terms, test
    (0 + z_man|hatch_nest_breed_ID) +
    (0 + z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_man|rear_area) +
    (0 + z_bs0|hatch_nest_breed_ID) +
    (0 + z_bs0|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0|rear_area) +
    
    # Interactions, test
    (0 + z_bs0:z_man|hatch_nest_breed_ID) +
    (0 + z_bs0:z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0:z_man|rear_area) +  
    (0 + z_man:(f_year.2002+f_year.2003)||rear_area) +
    (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) +
    (0 + f_sexM:z_man|hatch_nest_breed_ID) +
    (0 + f_sexM:z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + f_sexM:z_man|rear_area),

    data = ydata, family = "binomial", control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

    
  surv_null_k = glmer(chick_survival_to_first_breed_season ~ 
    
    z_nh:f_year + I(z_nh^2):f_year + # 2-way interactions (control)
    f_year + z_nh + f_sexM + I(z_nh^2) + # other terms (control)
    
    # Random intercepts..
    (1|hatch_nest_breed_ID) + # random intercept of hatching nest
    (1|rear_nest_breed_ID) + # random intercept of rearing nest
    (1|genetic_dad_ring_.WP_or_EP.) + # random intercept of dad (NOTE that we are ignoring the rearing parents, as they are already in the id of rearing nest, and there are very few cases of mother or fathers rearing multiple nests..)
    (1|rear_area) + # random intercept for rearing area (where the nest is)  
    
    # Random slopes..
    # Main terms, test
    (0 + z_man|hatch_nest_breed_ID) +
    (0 + z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_man|rear_area) +
    (0 + z_bs0|hatch_nest_breed_ID) +
    (0 + z_bs0|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0|rear_area) +
    
    # Interactions, test
    (0 + z_bs0:z_man|hatch_nest_breed_ID) +
    (0 + z_bs0:z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + z_bs0:z_man|rear_area) +  
    (0 + z_man:(f_year.2002+f_year.2003)||rear_area) +
    (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) +
    (0 + f_sexM:z_man|hatch_nest_breed_ID) +
    (0 + f_sexM:z_man|genetic_dad_ring_.WP_or_EP.) +
    (0 + f_sexM:z_man|rear_area),
    
    data = ydata, family = "binomial", control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  anova(surv_full_k, surv_null_k, test = "Chisq")  
  #             npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # surv_null_k   30 880.93 1046.7 -410.47   820.93                     
  # surv_full_k   40 885.21 1106.2 -402.60   805.21 15.726 10     0.1078

  ## Not significant!
  }

{## Prepare data for "Fledging" model..
  # reduce data to one data point per nest..
  ndata = xdata[which(duplicated(xdata$rear_nest_breed_ID) == F),]
  
  # prepare variables..
  # n. fledged chicks per nest
  ndata$fled_y = ndata$number_chicks_fledged_from_rear_nest
  # n. unfledged chicks per nest
  ndata$fled_n = ndata$rear_d0_rear_nest_brood_size + ndata$net_rearing_manipulation - ndata$number_chicks_fledged_from_rear_nest
  
  # make hatching year a factor..
  ndata$f_year = as.factor(ndata$hatch_year)

  ## Exclude rows with NAs in any of the relevant variables:
  xx = complete.cases(ndata[, c("fled_y", "fled_n", "net_rearing_manipulation", "rear_d0_rear_nest_brood_size", "f_year", "rear_nest_OH", "rear_area")])
  ndata=droplevels(ndata[xx, ])
  ## now ndata is complete wrt the relevant variables
  ## only *now* standardize covariates...
  ## (note that I put all the relevant varriables into ndata)
  ## (note also that I combine scale with as.vector which is because some modelling functions
  ## can't cope with the attributes that scale appends to the resulting object)
  ndata$z_man = as.vector(scale(ndata$net_rearing_manipulation))
  ndata$z_bs0 = as.vector(scale(ndata$rear_d0_rear_nest_brood_size))
  ndata$z_nh = as.vector(scale(ndata$rear_nest_OH))
  ## ... and manually dummy code and center the factor:
  ndata$f_year.2002=as.numeric(ndata$f_year=="2002")
  ndata$f_year.2003=as.numeric(ndata$f_year=="2003")
  ndata$f_year.2002=ndata$f_year.2002-mean(ndata$f_year.2002)
  ndata$f_year.2003=ndata$f_year.2003-mean(ndata$f_year.2003)
  
  ndata$nest.id=as.factor(1:nrow(ndata))

  }

{## Fit "Fledging model"..

  ## FULL MODEL
  caz_full = glmer(cbind(fled_y, fled_n) ~ 
    # z_man:z_bs0:f_year + # 3-way interaction
    z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + # 2-way interactions 
    z_nh:f_year + I(z_nh^2):f_year + # more 2-way interactions
    z_man + z_bs0 + f_year + z_nh + I(z_nh^2) + # other terms
    (1|nest.id) + # random effect of nest
    (1|rear_area) + # random intercept for rearing area (where the nest is)
    # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
    (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
    (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
    data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  ## Compare to null model..
  caz_null = glmer(cbind(fled_y, fled_n) ~ 
    z_nh:f_year + I(z_nh^2):f_year + # more 2-way interactions
    f_year + z_nh + I(z_nh^2) + # other terms
    (1|nest.id) + # random effect of nest
    (1|rear_area) + # random intercept for rearing area (where the nest is)
    # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
    (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
    (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
    data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  anova(caz_full, caz_null, test = "chisq")
  #          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # caz_null   20 1119.6 1199.2 -539.80   1079.6                        
  # caz_full   27 1110.2 1217.6 -528.08   1056.2 23.444  7   0.001426 **
  ## The test predictors *as a whole* have a significant effect on the response..

  ## Look at the significance of the highest-order interaction..
  drop1(caz_full, test = "Chisq")
  #                  npar    AIC    LRT Pr(Chi)  
  # <none>                1110.2                 
  # z_man:z_bs0         1 1109.5 1.3214 0.25034  
  # z_man:f_year        2 1107.3 1.2003 0.54872  
  # z_bs0:f_year        2 1110.2 4.0136 0.13442  
  # f_year:z_nh         2 1113.2 7.1006 0.02872 *
  # f_year:I(z_nh^2)    2 1109.0 2.8162 0.24460  

  ## Reduced_model_1
  caz_r1 = glmer(cbind(fled_y, fled_n) ~ 
    # z_man:z_bs0:f_year + # 3-way interaction
    # z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + # 2-way interactions 
    z_nh:f_year + 
    # I(z_nh^2):f_year + # more 2-way interactions
    z_man + z_bs0 + f_year + z_nh + I(z_nh^2) + # other terms
    (1|nest.id) + # random effect of nest
    (1|rear_area) + # random intercept for rearing area (where the nest is)
    # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
    # (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
    (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
    data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  drop1(caz_r1, test = "Chisq")
  #             npar    AIC     LRT   Pr(Chi)    
  # <none>           1098.3                      
  # z_man          1 1109.3 12.9681 0.0003168 ***
  # z_bs0          1 1100.0  3.6875 0.0548230 .  
  # I(z_nh^2)      1 1102.9  6.5728 0.0103548 *  
  # z_nh:f_year    2 1098.8  4.5166 0.1045286 

  caz_r2 = glmer(cbind(fled_y, fled_n) ~ 
    # z_man:z_bs0:f_year + # 3-way interaction
    # z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + # 2-way interactions 
    # z_nh:f_year + 
    # I(z_nh^2):f_year + # more 2-way interactions
    z_man + z_bs0 + f_year + z_nh + I(z_nh^2) + # other terms
    (1|nest.id) + # random effect of nest
    (1|rear_area) + # random intercept for rearing area (where the nest is)
    # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
    # (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
    (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
    data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  drop1(caz_r2, test = "Chisq")
  #           npar    AIC     LRT   Pr(Chi)    
  # <none>         1098.8                      
  # z_man        1 1109.6 12.7121 0.0003633 ***
  # z_bs0        1 1099.9  3.0260 0.0819388 .  
  # f_year       2 1109.2 14.3106 0.0007807 ***
  # z_nh         1 1101.8  4.9087 0.0267217 *  
  # I(z_nh^2)    1 1114.8 18.0115 2.196e-05 ***

  summary(caz_r2)
  # Fixed effects:
  #             Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)   4.9197     0.4385  11.219  < 2e-16 ***
  # z_man        -0.7109     0.1419  -5.010 5.45e-07 ***
  # z_bs0        -0.2667     0.1546  -1.726   0.0844 .  
  # f_year2002   -1.1669     0.6122  -1.906   0.0566 .  
  # f_year2003   -2.3163     0.5269  -4.396 1.10e-05 ***
  # z_nh         -0.5568     0.2506  -2.221   0.0263 *  
  # I(z_nh^2)    -0.5835     0.1372  -4.251 2.13e-05 ***

  }

{## Plot "Fledging model"..
  pdf(file = "fledging.pdf", width = 10.5, height = 5.5)
  par(mfrow = c(1,2))
  
  fe = fixef(caz_r2)

  tt = table(ndata$rear_d0_rear_nest_brood_size)
  sd_bs = sd(ndata$rear_d0_rear_nest_brood_size)
  m_bs = mean(ndata$rear_d0_rear_nest_brood_size)

  modal_bs = as.numeric(names(tt)[which(tt == max(tt))])
  max_bs = max(ndata$rear_d0_rear_nest_brood_size)
  min_bs = min(ndata$rear_d0_rear_nest_brood_size)
  
  sd_man = sd(ndata$net_rearing_manipulation)
  m_man = mean(ndata$net_rearing_manipulation)

  x = modal_bs
  
  ## Plot proportion of fledged chicks..

  plot(0,0, type = "n", xaxt = "n", xlim = (c(-5,5) - m_man)/sd_man, ylim = c(0,1), xlab = "net manipulation", ylab = "pr. chicks fledging", main = "a. Proportion of fledged chicks")
    
  abline(v = (0 - m_man)/sd_man)
    
  axis(1, at =  (-10:10 - m_man)/sd_man, labels = -10:10)
  
  cols = brewer.pal(length(levels(ndata$f_year)), name = "Set1")
  
  legend("bottomleft", col = c(cols), lwd = 2, legend = levels(ndata$f_year), title = "year", bty = "n")
    
  for(i in 1:length(levels(ndata$f_year))){
    
    year = levels(ndata$f_year)[i]
    
    sel = ndata[which(ndata$f_year == year),]

    newdata_1 = data.frame(
      z_man = seq((min(sel$net_rearing_manipulation) - m_man)/sd_man, (max(sel$net_rearing_manipulation) - m_man)/sd_man, l = 100), z_bs0 = x_z, f_year = year, z_nh = mean(sel$z_nh))
      
    pred_1 = predict(caz_r2, newdata = newdata_1, re.form = NA, type = "response")
    
    aa = aggregate(sel$fled_y/(sel$fled_y + sel$fled_n), by = list(sel$z_man), FUN = mean)
    nn = aa$Group.1
    mm = aa$x
    cc = as.vector(table(sel$z_man))
    
    points(nn, mm, cex = sqrt(cc), col = add.alpha(cols[i], 0.5), pch = 19)
              
    lines(newdata_1$z_man, pred_1, lty = 1, lwd = 2, col = cols[i])
  
    }

  ## Plot absolute number of fledged chicks..
            
  plot(0,0, type = "n", xaxt = "n", xlim = (c(-5,5) - m_man)/sd_man, ylim = c(0,15), xlab = "net manipulation", ylab = "n. chicks fledging", main = "b. Number of fledged chicks")
  
  abline(v = (0 - m_man)/sd_man)
  
  axis(1, at =  (-10:10 - m_man)/sd_man, labels = -10:10)
  
  cols = brewer.pal(length(levels(ndata$f_year)), name = "Set1")
  
  legend("bottomleft", col = c(cols), lwd = 2, legend = levels(ndata$f_year), title = "year", bty = "n")
    
  for(i in 1:length(levels(ndata$f_year))){
    
    year = levels(ndata$f_year)[i]
    
    sel = ndata[which(ndata$f_year == year),]

    newdata_1 = data.frame(
      z_man = seq((min(sel$net_rearing_manipulation) - m_man)/sd_man, (max(sel$net_rearing_manipulation) - m_man)/sd_man, l = 100), z_bs0 = x_z, f_year = year, z_nh = mean(sel$z_nh))
      
    pred_1 = predict(caz_r2, newdata = newdata_1, re.form = NA, type = "response")
        
    aa = aggregate(sel$fled_y/(sel$fled_y + sel$fled_n), by = list(sel$z_man), FUN = mean)
    nn = aa$Group.1
    mm = aa$x
    cc = as.vector(table(sel$z_man))
    
    points(nn, mm*(x + nn*sd_man+m_man), cex = sqrt(cc), col = add.alpha(cols[i], 0.5), pch = 19)
              
    lines(newdata_1$z_man, pred_1*(x + newdata_1$z_man*sd_man+m_man), lty = 1, lwd = 2, col = cols[i])
    
    }

  dev.off()

  }

{## Fit "Re-capture per nest model"..

  # Start from ndata and attach chicks recaptured or not recaptured (surv_y, surv_n)..
  ttt = table(xdata$rear_nest_breed_ID, xdata$chick_survival_to_first_breed_season)
  surv_y = ttt[,2]
  surv_n = ttt[,1]

  ndata$surv_y = surv_y[match(ndata$rear_nest_breed_ID, as.integer(rownames(ttt)))]
  ndata$surv_n = surv_n[match(ndata$rear_nest_breed_ID, as.integer(rownames(ttt)))]

  
  ## FULL MODEL
  caz_s_full = glmer(cbind(surv_y, surv_n) ~ 
    # z_man:z_bs0:f_year + # 3-way interaction
    z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + # 2-way interactions 
    z_nh:f_year + I(z_nh^2):f_year + # more 2-way interactions
    z_man + z_bs0 + f_year + z_nh + I(z_nh^2) + # other terms
    (1|nest.id) + # random effect of nest
    (1|rear_area) + # random intercept for rearing area (where the nest is)
    # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
    (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
    (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
    data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

#   ## Mi è venuto il dubbio che l'interazione 3-way tra anno, bs_0 e manipolazione nel GLMM "survival per nest" potesse essere significativa (alla fine non la avevamo testata..). Non lo è. Comunque, se guardi il grafico "survival_per_nest.pdf" sembra proprio che, mentre nel 2001 e 2003 quelli con tante uova hanno beneficiato della rimozione, questo non è in realtà successo nel 2002. Probabilmente vale la pena di dirlo nella discussione. Qualcosa tipo: "Nonostante l'interazione tra anno e manipolazione (né l'interazione a tre vie tra anno, BS_0 e manipolazione) non siano risultate statisticamente significative nei nostri dati, l'eslporazione visiva del plot in Fig. N suggerisce che nell'anno 2002, che è stato il più produttivo, anche le coppie con un numero molto elevato di nidi possano avere beneficiato dall'incremento sperimentale del numero di uova, portando alla sopravvivenza all'anno successivo un numero maggiore di piccoli. Probabilmente, la conferma statistica di questa ipotesi necessiterebbe di una quantità di dati ancora maggiore. Se confermata, comunque, questa osservazione implicherebbe che la clutch size delle cinciarelle sia "sottodimensionata" per annate favorevoli come sembrerebbe essere stato, localmente, il 2002." 
#     
#   caz_s_full_3 = glmer(cbind(surv_y, surv_n) ~ 
#     z_man:z_bs0:f_year + # 3-way interaction
#     z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + # 2-way interactions 
#     z_nh:f_year + I(z_nh^2):f_year + # more 2-way interactions
#     z_man + z_bs0 + f_year + z_nh + I(z_nh^2) + # other terms
#     (1|nest.id) + # random effect of nest
#     (1|rear_area) + # random intercept for rearing area (where the nest is)
#     # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
#     (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
#     (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
#     data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
# 
#   drop1(caz_s_full_3, test = "Chisq")
#   #                      npar    AIC    LRT Pr(Chi)
#   # <none>                  609.93               
#   # f_year:z_nh           2 610.43 4.5011  0.1053
#   # f_year:I(z_nh^2)      2 607.84 1.9084  0.3851
#   # z_man:z_bs0:f_year    2 606.20 0.2757  0.8712

    
  ## Compare to null model..
  caz_s_null = glmer(cbind(surv_y, surv_n) ~ 
    z_nh:f_year + I(z_nh^2):f_year + # more 2-way interactions
    f_year + z_nh + I(z_nh^2) + # other terms
    (1|nest.id) + # random effect of nest
    (1|rear_area) + # random intercept for rearing area (where the nest is)
    # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
    (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
    (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
    data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  anova(caz_s_full, caz_s_null, test = "chisq")
  #            npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
  # caz_s_null   20 611.14 690.72 -285.57   571.14                        
  # caz_s_full   27 606.20 713.63 -276.10   552.20 18.939  7    0.00838 **

  drop1(caz_s_full, test = "Chisq")
  #                    npar    AIC    LRT Pr(Chi)
  # <none>                606.20               
  # z_man:z_bs0         1 604.47 0.2626  0.6083
  # z_man:f_year        2 605.08 2.8809  0.2368
  # z_bs0:f_year        2 604.92 2.7196  0.2567
  # f_year:z_nh         2 606.60 4.3992  0.1109
  # f_year:I(z_nh^2)    2 604.10 1.8989  0.3869

  ## Reduced_model_1
  caz_s_r1 = glmer(cbind(surv_y, surv_n) ~ 
    # z_man:z_bs0:f_year + # 3-way interaction
    # z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + # 2-way interactions 
    # z_nh:f_year + 
    # I(z_nh^2):f_year + # more 2-way interactions
    z_man + z_bs0 + f_year + z_nh + I(z_nh^2) + # other terms
    (1|nest.id) + # random effect of nest
    (1|rear_area) + # random intercept for rearing area (where the nest is)
    # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
    # (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
    (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
    data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  drop1(caz_s_r1, test = "Chisq")
  #           npar    AIC     LRT   Pr(Chi)    
  # <none>         594.82                      
  # z_man        1 603.09 10.2648 0.0013559 ** 
  # z_bs0        1 593.76  0.9318 0.3343844    
  # f_year       2 606.68 15.8519 0.0003612 ***
  # z_nh         1 617.57 24.7464 6.539e-07 ***
  # I(z_nh^2)    1 594.20  1.3747 0.2410064    

  ## Reduced_model_2
  caz_s_r2 = glmer(cbind(surv_y, surv_n) ~ 
    # z_man:z_bs0:f_year + # 3-way interaction
    # z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + # 2-way interactions 
    # z_nh:f_year + 
    # I(z_nh^2):f_year + # more 2-way interactions
    z_man + z_bs0 + f_year + z_nh + # other terms
    (1|nest.id) + # random effect of nest
    (1|rear_area) + # random intercept for rearing area (where the nest is)
    # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
    # (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
    (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
    data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

  drop1(caz_s_r2, test = "Chisq")
  #          npar    AIC    LRT   Pr(Chi)    
  # <none>      594.20                     
  # z_man     1 602.82 10.619 0.0011196 ** 
  # z_bs0     1 593.05  0.851 0.3562817    
  # f_year    2 604.68 14.481 0.0007169 ***
  # z_nh      1 615.58 23.383 1.328e-06 ***
  
  summary(caz_s_r2)
  # Random effects:
  #  Groups      Name        Variance Std.Dev.
  #  nest.id     (Intercept) 0.23825  0.4881  
  #  rear_area   (Intercept) 0.02378  0.1542  
  #  rear_area.1 z_man       0.00000  0.0000  
  #  rear_area.2 z_bs0       0.00000  0.0000  
  #  rear_area.3 f_year.2002 0.19502  0.4416  
  #  rear_area.4 f_year.2003 0.00000  0.0000  
  # Number of obs: 395, groups:  nest.id, 395; rear_area, 9
  # 
  # Fixed effects:
  #             Estimate Std. Error z value Pr(>|z|)    
  # (Intercept) -2.64270    0.26754  -9.878  < 2e-16 ***
  # z_man       -0.32403    0.08427  -3.845 0.000120 ***
  # z_bs0       -0.10388    0.11233  -0.925 0.355105    
  # f_year2002  -0.88176    0.46836  -1.883 0.059751 .  
  # f_year2003  -1.40658    0.37287  -3.772 0.000162 ***
  # z_nh        -0.96438    0.20380  -4.732 2.22e-06 ***
  # ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Correlation of Fixed Effects:
  #            (Intr) z_man  z_bs0  f_2002 f_2003
  # z_man       0.033                            
  # z_bs0      -0.158  0.133                     
  # f_year2002 -0.796  0.022  0.018              
  # f_year2003 -0.711  0.041  0.208  0.705       
  # z_nh       -0.596  0.049  0.158  0.812  0.670
  # convergence code: 0
  # boundary (singular) fit: see ?isSingular

  ## 
  
  save.image("all_models.RData")
  }

{## Plot "Survival per nest model"..
  pdf(file = "survival_per_nest.pdf", width = 10.5, height = 5.5)
  par(mfrow = c(1,2))
  
  fe = fixef(caz_s_r2)

  tt = table(ndata$rear_d0_rear_nest_brood_size)
  sd_bs = sd(ndata$rear_d0_rear_nest_brood_size)
  m_bs = mean(ndata$rear_d0_rear_nest_brood_size)

  modal_bs = as.numeric(names(tt)[which(tt == max(tt))])
  max_bs = max(ndata$rear_d0_rear_nest_brood_size)
  min_bs = min(ndata$rear_d0_rear_nest_brood_size)
  
  sd_man = sd(ndata$net_rearing_manipulation)
  m_man = mean(ndata$net_rearing_manipulation)

  x = modal_bs
  
  ## Plot proportion of fledged chicks..

  plot(0,0, type = "n", xaxt = "n", xlim = (c(-5,5) - m_man)/sd_man, ylim = c(0,0.2), xlab = "net manipulation", ylab = "pr. chicks surviving", main = "a. Proportion of chicks surviving at 1 year")
    
  abline(v = (0 - m_man)/sd_man)
    
  axis(1, at =  (-10:10 - m_man)/sd_man, labels = -10:10)
  
  cols = brewer.pal(length(levels(ndata$f_year)), name = "Set1")
  
  legend("topright", col = c(cols), lwd = 2, legend = levels(ndata$f_year), title = "year", bty = "n")
    
  for(i in 1:length(levels(ndata$f_year))){
    
    year = levels(ndata$f_year)[i]
    
    sel = ndata[which(ndata$f_year == year),]

    newdata_1 = data.frame(
      z_man = seq((min(sel$net_rearing_manipulation) - m_man)/sd_man, (max(sel$net_rearing_manipulation) - m_man)/sd_man, l = 100), z_bs0 = x_z, f_year = year, z_nh = mean(sel$z_nh))
      
    pred_1 = predict(caz_s_r2, newdata = newdata_1, re.form = NA, type = "response")
    
    aa = aggregate(sel$surv_y/(sel$surv_y + sel$surv_n), by = list(sel$z_man), FUN = mean)
    nn = aa$Group.1
    mm = aa$x
    cc = as.vector(table(sel$z_man))
    
    points(nn, mm, cex = sqrt(cc), col = add.alpha(cols[i], 0.5), pch = 19)
              
    lines(newdata_1$z_man, pred_1, lty = 1, lwd = 2, col = cols[i])
  
    }
  }
  
## Reduced_model_3 (for plotting only!)
  caz_s_r3 = glmer(cbind(surv_y, surv_n) ~ 
    # z_man:z_bs0:f_year + # 3-way interaction
    # z_man:z_bs0 + z_man:f_year + f_year:z_bs0 + # 2-way interactions 
    # z_nh:f_year + 
    # I(z_nh^2):f_year + # more 2-way interactions
    z_man + f_year + z_nh + # other terms
    (1|nest.id) + # random effect of nest
    (1|rear_area) + # random intercept for rearing area (where the nest is)
    # (0 + z_man:z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 3-way interaction in area
    # (0 + z_man:z_bs0|rear_area) + (0 + z_man:(f_year.2002+f_year.2003)||rear_area) + (0 + z_bs0:(f_year.2002+f_year.2003)||rear_area) + # RS for 2-way interaction(s) in area (I did not include date in year as I assume that the weather/season progression is the same everywhere)
    (0 + z_man|rear_area) + (0 + z_bs0|rear_area) + (0 + (f_year.2002+f_year.2003)||rear_area), # RS for other terms
    data = ndata, family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
   
  {## Plot absolute number of fledged chicks..

    pdf(file = "survival_per_nest.pdf", height = 3, width = 8)

    par(mfrow = c(1,3))

    years =  c("2001", "2002", "2003")
    colx = c("Reds", "Blues", "Greens")

    for(j in 1:3){

      plot(0,0, type = "n", xaxt = "n", xlim = (c(-5,5) - m_man)/sd_man, ylim = c(0,2), xlab = "net manipulation", ylab = "n. chicks survived", main = years[j])
      
      abline(v = (0 - m_man)/sd_man)
      
      axis(1, at = (-10:10 - m_man)/sd_man, labels = -10:10, cex.axis = 0.8)
      
      clutch_sizes = seq(4, 13, by = 3)
      
      year = years[j]
      
      cols = brewer.pal(length(clutch_sizes)+1, name = colx[j])
      
      legend("topleft", col = cols[2:length(cols)], lwd = 2, legend = clutch_sizes, title = "BS_0", bty = "n", cex = 0.7)
        
      for(i in 1:length(clutch_sizes)){
        
        x = clutch_sizes[i]
        
        sel = ndata[which(ndata$f_year == year & (ndata$rear_d0_rear_nest_brood_size > x - 2 & ndata$rear_d0_rear_nest_brood_size <= x + 1)),]

        newdata_1 = data.frame(
          z_man = seq((min(sel$net_rearing_manipulation) - m_man)/sd_man, (max(sel$net_rearing_manipulation) - m_man)/sd_man, l = 100), z_bs0 = x_z, f_year = year, z_nh = mean(sel$z_nh))
          
        pred_1 = predict(caz_s_r3, newdata = newdata_1, re.form = NA, type = "response")
            
        aa = aggregate(sel$surv_y/(sel$surv_y + sel$surv_n), by = list(sel$z_man), FUN = mean)
        nn = aa$Group.1
        mm = aa$x
        cc = as.vector(table(sel$z_man))
        
        points(nn, mm*(x + nn*sd_man+m_man), cex = sqrt(cc), col = add.alpha(cols[i+1], 0.9), pch = 19)
                  
        lines(newdata_1$z_man, pred_1*(x + newdata_1$z_man*sd_man+m_man), lty = 1, lwd = 2, col = cols[i+1])
        }
      }

    dev.off()
    }

## Plot number of chicks per nest..
pdf(file = "chicks_per_nest.pdf", height = 6.5, width = 2.75)
cols = brewer.pal(length(levels(ndata$f_year)), name = "Set1")
par(mfrow = c(3,1))
for(i in 1:length(levels(f_year))) {
  year = levels(f_year)[i]
  hist(ndata$rear_d0_rear_nest_brood_size[which(ndata$f_year == year)], xlim = c(1,16), ylim = c(0,40), main = year, breaks = seq(0,18,by = 1), col = add.alpha(cols[i],0.75), xlab = "BS_0", xaxt = "n")
  axis(side = 1, at = seq(2, 16, by = 2)-0.5, labels = seq(2, 16, by = 2))
  }
dev.off()






































