#--------------------------------------#
# Decomposing the NH Black-White racial disparity on the absolute scale 
# using g-estimation of a structural nested mean model
# 
# Updated: 03/01/2023
# K. Labgold
# To view details and output of code, please download the associated html file
#--------------------------------------#


# Estimate Total Absolute Racial Disparity ----

unadj.G.imp <- glm(formula = defG.flag ~ raceeth.1_0,
                   data = model.dat2,
                   family = binomial(link = "identity"))

te <- round(unadj.G.imp$coefficients[2]*10000, 1)

df <- confint(unadj.G.imp)

te.lci <- data.frame(round(df*10000, 1))[2,1]
te.uci <- data.frame(round(df*10000, 1))[2,2]

out.unadj <- cbind(te, te.lci, te.uci)

kableExtra::kable(out.unadj)

# Decomposition Method ----
## Stage 1 ----
### Step 1: Obtain predicted mediator values and create mediator residuals. ----
### HDP ~ Race + C_HDP.SMM
### Mediator model: # m ~ x + c_xy + c_my

mM_HDP <- gest.HDP ~ raceeth.1_0 + 
              rural.1_0 + insure_imp + age.NOT20to34.1_0 + ct.proportion.poverty_imp +
              multiplegest_imp + hdd.gestdiabetes.flag +
              hdd.renal.flag + hdd.prediabetes.flag + hdd.obesity.flag

### Model 8
mod_8.HDP <- glm(formula = mM_HDP,
                 family = binomial(link = "logit"),
                 data = data2) 

### Pull the predicted probabilities for the mediator
data2$pred_med.HDP <- mod_8.HDP$fitted.values 

head(data2$pred_med.HDP)
range(data2$pred_med.HDP)

### Calculate the mediator residual: observed mediator value - predicted probability
data2 <- data2 %>% 
            mutate(mediator_resid.HDP.v2 = gest.HDP - pred_med.HDP)

head(data2$mediator_resid.HDP.v2)
range(data2$mediator_resid.HDP.v2)
### Note: Be sure to explore propensity score overlap between two racial groups.



### Step 2: Regress outcome against predicted mediator residuals, ----
### interaction between exposure and mediator residuals, and mediator-outcome confounders 
### using ordinary least squares to obtain mediator effect estimates.
### SMM ~ Race + HDP.resid + Race:HDP.resid + C_HDP.SMM
### Outcome model: # y ~ x + m.residual + x:m.residual + c_xy + c_my
mY <- defG.flag ~ raceeth.1_0 + mediator_resid.HDP.v2 + raceeth.1_0:mediator_resid.HDP.v2 + 
  rural.1_0 + insure_imp + age.NOT20to34.1_0 + ct.proportion.poverty_imp +
  multiplegest_imp + hdd.gestdiabetes.flag +
  hdd.renal.flag + hdd.prediabetes.flag + hdd.obesity.flag


mod_9 <- glm(formula = mY, 
             family = gaussian(link = "identity"),
             data = data2)

## Intermediate Transformation Between Stages ----
## Create transformed outcome by removing estimated mediator effect from observed outcome.
data2 <- data2 %>% 
          mutate(defG_trans = defG.flag - 
                   (mod_9$coefficients["mediator_resid.HDP.v2"])*gest.HDP - # Mediator
                   (mod_9$coefficients["raceeth.1_0:mediator_resid.HDP.v2"])*gest.HDP*raceeth.1_0) # Exposure x Mediator

head(data2$defG_trans)
range(data2$defG_trans)

## Estimate race-specific risk with mediator effect removed from observed outcome.
cdm.cases = sum(data2$defG_trans)
cdm.cases_nhb = sum(data2$defG_trans[data2$raceeth.1_0 == 1])
cdm.cases_nhw = sum(data2$defG_trans[data2$raceeth.1_0 == 0])

cdm.nhb.risk = round((cdm.cases_nhb/sum(data2$raceeth.1_0 == 1))*10000, 1)
cdm.nhw.risk = round((cdm.cases_nhw/sum(data2$raceeth.1_0 == 0))*10000, 1)

## Stage 2 ----
### Step 1 Obtain predicted exposure values and create exposure residuals. ----
### Null model because there are no hypothesized confounders of Race-SMM (C_xy). mX is included for code completeness.
### Exposure model: x ~ c_xy
mX <- raceeth.1_0 ~ NULL  

mod_10 <- glm(formula = mX, # exposure
              family = binomial(link = "logit"),
              data = data2)

### Pull predicted exposure values
data2$pred_exp <- mod_10$fitted.values 

table(data2$pred_exp)

### Calculate exposure residuals
data2 <- data2 %>% 
          mutate(exp_resid.v2 = raceeth.1_0 - pred_exp)

table(data2$exp_resid.v2)

### Step 2: Regress transformed outcome against predicted exposure residuals and ----
### exposure-outcome confounders using ordinary least squares to estimate counterfactual disparity measure.

mod_11 <- glm(formula = defG_trans ~ exp_resid.v2,
              family = gaussian(link = "identity"),
              data = data2)

### Output:
psi = mod_11$coefficients[2]*10000
stderr = sqrt(diag(vcov(mod_11)))[2]*10000
lci = round(psi - 1.96*stderr, 1)
uci = round(psi + 1.96*stderr, 1)
psi.r = round(psi, 1)

out = cbind(psi.r, lci, uci)

kableExtra::kable(out)
### Interpretation (out): Estimated CDM(No HDP): Excess 41.1 SMM events per 10,000 delivery hospitalizations among 
### NH Black women compared to NH White women.

nhb.risk = round((sum(model.dat2$defG.flag[model.dat2$raceeth.1_0 == 1])/sum(model.dat2$raceeth.1_0 == 1))*10000, 1)
nhw.risk = round((sum(model.dat2$defG.flag[model.dat2$raceeth.1_0 == 0])/sum(model.dat2$raceeth.1_0 == 0))*10000, 1)

out1 = cbind(nhb.risk, nhw.risk)
out2 = cbind(cdm.nhb.risk, cdm.nhw.risk)
out3 = round(rbind(out1, out2), 1)
names(out3) <- c("NH Black Rate per 10,000", "NH White Rate per 10,000")
row.names(out3) <- c("Observed Rate", "Counterfactual Rate w/ HDP Effect Removed from Outcome")

kableExtra::kable(out3)

### Interpretation (out3): Both NH Black and NH White rates decreased when effect of HDP was 
### removed from outcome, but reduction was of greater magnitude for NH Black rate.