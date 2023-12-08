P8108 Survival Analysis Heart Failure
================
Yi Huang
2023-11-10

## Preparation

install and read packages, read functions

``` r
#----------------------------------------------------------------
# CLEAR ENVIRONMENT
#----------------------------------------------------------------
rm(list = ls())

#----------------------------------------------------------------
# INSTALL PACKAGES
#----------------------------------------------------------------
packages <- c("dplyr", "tidyverse", "readxl", 
              "survival", "KMsurv", 
              # ggsurvplot()
              "survminer",
              # stepwiseCox()
              "StepReg")


# Install missing packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], dependencies = TRUE)
}

# Load packages invisibly
invisible(lapply(packages, library, character.only = TRUE))

# Remove variables associated with package installation
rm(packages, installed_packages)

# Read function
source("shared_code/stepwiseCox.R", local = TRUE)
# source("shared_code/stepwiseCox.R", local = TRUE)$value
```

## Load data

``` r
# read data
data <- read_csv("data/heart_failure.csv")

# # data cleaning
data <- data |>
  arrange(TIME) |>
  janitor::clean_names() |>
  mutate(gender = factor(gender),
         smoking = factor(smoking),
         diabetes = factor(diabetes),
         bp = factor(bp),
         # event = factor(event),
         anaemia = factor(anaemia)) |>
  rename(pletelets = pletelets)
```

## Cox model

**Full model**

``` r
# use survival to fit a cox model
cox_model <- coxph(Surv(time, event) ~ gender + smoking + diabetes + bp + anaemia + age + ejection_fraction + sodium + creatinine + pletelets + cpk, data = data)
summary(cox_model)
```

    ## Call:
    ## coxph(formula = Surv(time, event) ~ gender + smoking + diabetes + 
    ##     bp + anaemia + age + ejection_fraction + sodium + creatinine + 
    ##     pletelets + cpk, data = data)
    ## 
    ##   n= 299, number of events= 96 
    ## 
    ##                         coef  exp(coef)   se(coef)      z Pr(>|z|)    
    ## gender1           -2.375e-01  7.886e-01  2.516e-01 -0.944   0.3452    
    ## smoking1           1.289e-01  1.138e+00  2.512e-01  0.513   0.6078    
    ## diabetes1          1.399e-01  1.150e+00  2.231e-01  0.627   0.5307    
    ## bp1                4.757e-01  1.609e+00  2.162e-01  2.201   0.0278 *  
    ## anaemia1           4.601e-01  1.584e+00  2.168e-01  2.122   0.0338 *  
    ## age                4.641e-02  1.048e+00  9.324e-03  4.977 6.45e-07 ***
    ## ejection_fraction -4.894e-02  9.522e-01  1.048e-02 -4.672 2.98e-06 ***
    ## sodium            -4.419e-02  9.568e-01  2.327e-02 -1.899   0.0575 .  
    ## creatinine         3.210e-01  1.379e+00  7.017e-02  4.575 4.76e-06 ***
    ## pletelets         -4.635e-07  1.000e+00  1.126e-06 -0.412   0.6806    
    ## cpk                2.207e-04  1.000e+00  9.919e-05  2.225   0.0260 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                   exp(coef) exp(-coef) lower .95 upper .95
    ## gender1              0.7886     1.2681    0.4816     1.291
    ## smoking1             1.1376     0.8790    0.6953     1.861
    ## diabetes1            1.1501     0.8695    0.7427     1.781
    ## bp1                  1.6092     0.6214    1.0534     2.458
    ## anaemia1             1.5843     0.6312    1.0358     2.423
    ## age                  1.0475     0.9547    1.0285     1.067
    ## ejection_fraction    0.9522     1.0502    0.9329     0.972
    ## sodium               0.9568     1.0452    0.9141     1.001
    ## creatinine           1.3786     0.7254    1.2014     1.582
    ## pletelets            1.0000     1.0000    1.0000     1.000
    ## cpk                  1.0002     0.9998    1.0000     1.000
    ## 
    ## Concordance= 0.741  (se = 0.027 )
    ## Likelihood ratio test= 81.95  on 11 df,   p=6e-13
    ## Wald test            = 87.27  on 11 df,   p=6e-14
    ## Score (logrank) test = 88.39  on 11 df,   p=3e-14

``` r
# # stepwiseCox example
# lung <- survival::lung
# 
# my.data <- na.omit(lung)
# my.data$status
# 
# my.data$status1 <- ifelse(my.data$status==2,1,0)
# 
# data <- my.data
# 
# formula = Surv(time, status1) ~ . - status 
# 
# stepwiseCox(formula=formula,
#             data=my.data,
#             selection="bidirection",
#             select="HQ",
#             method="efron")
```

## Variable Selection

``` r
stepwise_data <- read_csv("data/heart_failure.csv")
## data cleaning for stepwiseCox, variables need to be all numeric
stepwise_data <- stepwise_data |>
  arrange(TIME) |>
  janitor::clean_names()

# Variable selection using stepwise Cox model using Sl
stepwise_model1 <- stepwiseCox(Surv(time, event) ~ gender + smoking + diabetes + bp + anaemia + age + ejection_fraction + sodium + creatinine + pletelets + cpk, data = stepwise_data)

stepwise_model1
```

    ##           Table 1. Summary of Parameters          
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##            Paramters                  Value       
    ## ——————————————————————————————————————————————————
    ## Response Variable              Surv(time, event)   
    ## Included Variable              NULL                
    ## Selection Method               forward             
    ## Select Criterion               SL                  
    ## Entry Significance Level(sle)  0.15                
    ## Method                         efron               
    ## Multicollinearity Terms        NULL                
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                        Table 2. Variables Type                                       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##    class                                            variable                                         
    ## —————————————————————————————————————————————————————————————————————————————————————————————————————
    ## nmatrix.2  Surv(time, event)                                                                          
    ## numeric    gender smoking diabetes bp anaemia age ejection_fraction sodium creatinine pletelets cpk   
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                        Table 3. Process of Selection                        
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##  Step    EnteredEffect    RemovedEffect  DF  NumberIn           SL          
    ## ————————————————————————————————————————————————————————————————————————————
    ## 1     age                               1   1         1.23929981722045e-06   
    ## 2     ejection_fraction                 1   2         6.41568078848927e-07   
    ## 3     creatinine                        1   3         1.97703471022321e-05   
    ## 4     bp                                1   4         0.0285149082859129     
    ## 5     anaemia                           1   5         0.112154426743682      
    ## 6     cpk                               1   6         0.074196169932428      
    ## 7     sodium                            1   7         0.0590078707141594     
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                 Table 4. Selected Varaibles                                
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##  variables1     variables2      variables3  variables4  variables5  variables6  variables7 
    ## ———————————————————————————————————————————————————————————————————————————————————————————
    ## age         ejection_fraction  creatinine  bp          anaemia     cpk         sodium       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                       Table 5. Coefficients of the Selected Variables                                      
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##      Variable               coef              exp(coef)            se(coef)                z                Pr(>|z|)       
    ## ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
    ## age                0.0435743706459777    1.04453767436502   0.0088312390906717    4.93411742096357   8.05139942240973e-07   
    ## ejection_fraction  -0.0474729978018283   0.953636223059251  0.0102739896095381    -4.6206974706063   3.82452124142e-06      
    ## creatinine         0.313898473194111     1.36875076459962   0.0689530824782774    4.5523486682847    5.305031910756e-06     
    ## bp                 0.496532093434627     1.64301356200148   0.213688570525934     2.32362494733599   0.0201456041559075     
    ## anaemia            0.446002041151676     1.56205465492097   0.214992549688069     2.07449998522636   0.0380329001283438     
    ## cpk                0.000210061077954651  1.00021008314233   9.82507258683582e-05  2.13801044316051   0.0325158956606839     
    ## sodium             -0.0456907483452765   0.955337356171037  0.0233585301861612    -1.95606264525779  0.0504577748881677     
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗

``` r
# 7 variables are selected: age, ejection_fraction, creatinine, bp, anaemia, cpk, sodium   

# # Variable selection using stepwise Cox model using AIC
# stepwise_model2 <- stepwiseCox(Surv(time, event) ~ gender + smoking + diabetes + bp + anaemia + age + ejection_fraction + sodium + creatinine + pletelets + cpk, select = "AIC", data = stepwise_data)
# stepwise_model2
# same as SL, 7 variables: age, ejection_fraction, creatinine, bp, anaemia, cpk, sodium   

stepwise_model3 <- stepwiseCox(Surv(time, event) ~ gender + smoking + diabetes + bp + anaemia + age + ejection_fraction + sodium + creatinine + pletelets + cpk, select = "AICc", data = stepwise_data)
stepwise_model3
```

    ##        Table 1. Summary of Parameters       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##         Paramters               Value       
    ## ————————————————————————————————————————————
    ## Response Variable        Surv(time, event)   
    ## Included Variable        NULL                
    ## Selection Method         forward             
    ## Select Criterion         AICc                
    ## Method                   efron               
    ## Multicollinearity Terms  NULL                
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                        Table 2. Variables Type                                       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##    class                                            variable                                         
    ## —————————————————————————————————————————————————————————————————————————————————————————————————————
    ## nmatrix.2  Surv(time, event)                                                                          
    ## numeric    gender smoking diabetes bp anaemia age ejection_fraction sodium creatinine pletelets cpk   
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                      Table 3. Process of Selection                      
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##  Step    EnteredEffect    RemovedEffect  DF  NumberIn        AICc       
    ## ————————————————————————————————————————————————————————————————————————
    ## 1     age                               1   1         994.937602832907   
    ## 2     ejection_fraction                 1   2         970.240960412782   
    ## 3     creatinine                        1   3         952.161508339985   
    ## 4     bp                                1   4         947.543535730665   
    ## 5     anaemia                           1   5         945.247044012968   
    ## 6     cpk                               1   6         942.336538050341   
    ## 7     sodium                            1   7         939.100392487579   
    ## 8     gender                            1   8         938.816678046656   
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                       Table 4. Selected Varaibles                                      
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##  variables1     variables2      variables3  variables4  variables5  variables6  variables7  variables8 
    ## ———————————————————————————————————————————————————————————————————————————————————————————————————————
    ## age         ejection_fraction  creatinine  bp          anaemia     cpk         sodium      gender       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                       Table 5. Coefficients of the Selected Variables                                       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##      Variable               coef              exp(coef)            se(coef)                z                 Pr(>|z|)       
    ## ————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
    ## age                0.0443478623916164    1.04534592818268   0.00886373828574522   5.00329104514934    5.63597512435079e-07   
    ## ejection_fraction  -0.0489345323951675   0.952243468757005  0.0104472640821526    -4.68395668094232   2.81389502185893e-06   
    ## creatinine         0.314161877214093     1.36911134654065   0.0687004808063813    4.57292108478099    4.80971325244378e-06   
    ## bp                 0.47498672100622      1.60799284481703   0.215345509645011     2.20569596175568    0.0274052948613976     
    ## anaemia            0.454708120224179     1.57571339703387   0.215535324007572     2.10966866947622    0.0348869049533196     
    ## cpk                0.000220294490030447  1.00022031875664   9.91832292359842e-05  2.22108608206641    0.0263451315283099     
    ## sodium             -0.0466035038416729   0.954465764583797  0.0233169534882998    -1.99869609316921   0.0456412459841202     
    ## gender             -0.183201088691266    0.832600712314579  0.222766772632833     -0.822389652307891  0.410855166710747      
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗

``` r
# AICc 8 variables: age, ejection_fraction, creatinine, bp, anaemia, cpk, sodium, gender   
```

**consider transformation on left-skewed predictors**

``` r
stepwise_data <- stepwise_data %>%
  mutate(logcpk = log(cpk+1),
         logcre=log(creatinine+1))

# Variable selection using stepwise Cox model using SL
stepwise_model4 <- stepwiseCox(Surv(time, event) ~ gender + smoking + diabetes + bp + 
                                anaemia + age + ejection_fraction + sodium + logcre + 
                                pletelets + logcpk, data = stepwise_data)
stepwise_model4
```

    ##           Table 1. Summary of Parameters          
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##            Paramters                  Value       
    ## ——————————————————————————————————————————————————
    ## Response Variable              Surv(time, event)   
    ## Included Variable              NULL                
    ## Selection Method               forward             
    ## Select Criterion               SL                  
    ## Entry Significance Level(sle)  0.15                
    ## Method                         efron               
    ## Multicollinearity Terms        NULL                
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                       Table 2. Variables Type                                       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##    class                                           variable                                         
    ## ————————————————————————————————————————————————————————————————————————————————————————————————————
    ## nmatrix.2  Surv(time, event)                                                                         
    ## numeric    gender smoking diabetes bp anaemia age ejection_fraction sodium logcre pletelets logcpk   
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                        Table 3. Process of Selection                        
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##  Step    EnteredEffect    RemovedEffect  DF  NumberIn           SL          
    ## ————————————————————————————————————————————————————————————————————————————
    ## 1     logcre                            1   1         1.6871984448572e-07    
    ## 2     age                               1   2         9.78512797059471e-06   
    ## 3     ejection_fraction                 1   3         4.14836454652977e-06   
    ## 4     bp                                1   4         0.0163833264326175     
    ## 5     anaemia                           1   5         0.0689106560670863     
    ## 6     sodium                            1   6         0.137463071634242      
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                           Table 4. Selected Varaibles                          
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##  variables1  variables2     variables3      variables4  variables5  variables6 
    ## ———————————————————————————————————————————————————————————————————————————————
    ## logcre      age         ejection_fraction  bp          anaemia     sodium       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                      Table 5. Coefficients of the Selected Variables                                     
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##      Variable              coef              exp(coef)           se(coef)                z                Pr(>|z|)       
    ## —————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
    ## logcre             1.30763101581221     3.69740423676436   0.293778646211065    4.45107577653123   8.54411832970726e-06   
    ## age                0.0416644317595133   1.04254457519706   0.00904436952880791  4.60667066143248   4.09167266089582e-06   
    ## ejection_fraction  -0.0429559987808823  0.957953540263324  0.0101497273424933   -4.23223179612333  2.31383757283858e-05   
    ## bp                 0.50744191505745     1.66103668263713   0.21181836595948     2.39564644339916   0.0165910851849217     
    ## anaemia            0.414916745814877    1.51424466823816   0.209598507238112    1.97957872545111   0.0477508854851059     
    ## sodium             -0.0366399594020293  0.96402316034887   0.0240870312129267   -1.52114883225484  0.128222492476848      
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗

``` r
# 6 avriables: logcre, age, ejection_fraction, bp, anaemia, sodium  

# Variable selection using stepwise Cox model using SL
stepwise_model5 <- stepwiseCox(Surv(time, event) ~ gender + smoking + diabetes + bp + 
                                anaemia + age + ejection_fraction + sodium + logcre + 
                                pletelets + logcpk, select = "AIC", data = stepwise_data)
stepwise_model5
```

    ##        Table 1. Summary of Parameters       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##         Paramters               Value       
    ## ————————————————————————————————————————————
    ## Response Variable        Surv(time, event)   
    ## Included Variable        NULL                
    ## Selection Method         forward             
    ## Select Criterion         AIC                 
    ## Method                   efron               
    ## Multicollinearity Terms  NULL                
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                       Table 2. Variables Type                                       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##    class                                           variable                                         
    ## ————————————————————————————————————————————————————————————————————————————————————————————————————
    ## nmatrix.2  Surv(time, event)                                                                         
    ## numeric    gender smoking diabetes bp anaemia age ejection_fraction sodium logcre pletelets logcpk   
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                      Table 3. Process of Selection                      
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##  Step    EnteredEffect    RemovedEffect  DF  NumberIn        AIC        
    ## ————————————————————————————————————————————————————————————————————————
    ## 1     logcre                            1   1         993.048326864199   
    ## 2     age                               1   2         975.495411177225   
    ## 3     ejection_fraction                 1   3         956.300370091503   
    ## 4     bp                                1   4         952.539110881546   
    ## 5     anaemia                           1   5         951.230331424551   
    ## 6     sodium                            1   6         951.024200746901   
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                           Table 4. Selected Varaibles                          
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##  variables1  variables2     variables3      variables4  variables5  variables6 
    ## ———————————————————————————————————————————————————————————————————————————————
    ## logcre      age         ejection_fraction  bp          anaemia     sodium       
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ## 
    ##                                      Table 5. Coefficients of the Selected Variables                                     
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
    ##      Variable              coef              exp(coef)           se(coef)                z                Pr(>|z|)       
    ## —————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
    ## logcre             1.30763101581221     3.69740423676436   0.293778646211065    4.45107577653123   8.54411832970726e-06   
    ## age                0.0416644317595133   1.04254457519706   0.00904436952880791  4.60667066143248   4.09167266089582e-06   
    ## ejection_fraction  -0.0429559987808823  0.957953540263324  0.0101497273424933   -4.23223179612333  2.31383757283858e-05   
    ## bp                 0.50744191505745     1.66103668263713   0.21181836595948     2.39564644339916   0.0165910851849217     
    ## anaemia            0.414916745814877    1.51424466823816   0.209598507238112    1.97957872545111   0.0477508854851059     
    ## sodium             -0.0366399594020293  0.96402316034887   0.0240870312129267   -1.52114883225484  0.128222492476848      
    ## ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗

``` r
# same as SL 6 avriables: logcre, age, ejection_fraction, bp, anaemia, sodium  
```

## Check assumptions for Cox model

**full model**

``` r
# Check assumptions with cox.zph
# full model
cox_zph <- cox.zph(cox_model)
plot(cox_zph) # Residual plots
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-7.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-8.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-9.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-10.png)<!-- -->

``` r
# Plot survival curves
ggsurvplot(survfit(cox_model), data = data, conf.int = TRUE)
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-5-11.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-5-12.png)<!-- -->

**Refit model with variables selected from stepwiseCox**

``` r
# refit a model with 7 selected variables 
# bp, anaemia, age, ejection_fraction, sodium, creatinine, cpk
step_model <- coxph(Surv(time, event) ~ bp + anaemia + age + ejection_fraction + sodium + 
                      creatinine + cpk, data = stepwise_data)
summary(step_model)
```

    ## Call:
    ## coxph(formula = Surv(time, event) ~ bp + anaemia + age + ejection_fraction + 
    ##     sodium + creatinine + cpk, data = stepwise_data)
    ## 
    ##   n= 299, number of events= 96 
    ## 
    ##                         coef  exp(coef)   se(coef)      z Pr(>|z|)    
    ## bp                 4.965e-01  1.643e+00  2.137e-01  2.324   0.0201 *  
    ## anaemia            4.460e-01  1.562e+00  2.150e-01  2.074   0.0380 *  
    ## age                4.357e-02  1.045e+00  8.831e-03  4.934 8.05e-07 ***
    ## ejection_fraction -4.747e-02  9.536e-01  1.027e-02 -4.621 3.82e-06 ***
    ## sodium            -4.569e-02  9.553e-01  2.336e-02 -1.956   0.0505 .  
    ## creatinine         3.139e-01  1.369e+00  6.895e-02  4.552 5.31e-06 ***
    ## cpk                2.101e-04  1.000e+00  9.825e-05  2.138   0.0325 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                   exp(coef) exp(-coef) lower .95 upper .95
    ## bp                   1.6430     0.6086    1.0808     2.498
    ## anaemia              1.5621     0.6402    1.0249     2.381
    ## age                  1.0445     0.9574    1.0266     1.063
    ## ejection_fraction    0.9536     1.0486    0.9346     0.973
    ## sodium               0.9553     1.0468    0.9126     1.000
    ## creatinine           1.3688     0.7306    1.1957     1.567
    ## cpk                  1.0002     0.9998    1.0000     1.000
    ## 
    ## Concordance= 0.738  (se = 0.027 )
    ## Likelihood ratio test= 80.58  on 7 df,   p=1e-14
    ## Wald test            = 88.43  on 7 df,   p=3e-16
    ## Score (logrank) test = 87.66  on 7 df,   p=4e-16

``` r
# Check assumptions with model obtained from stepwiseCox
cox_step <- cox.zph(step_model)
plot(cox_step) # Residual plots
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-6-6.png)<!-- -->

``` r
# Plot survival curves
ggsurvplot(survfit(step_model), data = data, conf.int = TRUE)
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-6-7.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-6-8.png)<!-- -->

**Refit model with variables selected from stepwiseCox** log transform
on left-skewed predictors

``` r
# refit a model with 6 selected variables 
# logcre,age, ejection_fraction, bp, anaemia, sodium
step_model4 <- coxph(Surv(time, event) ~ log(creatinine+1)+age + ejection_fraction + bp +
                       anaemia + sodium, data = stepwise_data)
summary(step_model4)
```

    ## Call:
    ## coxph(formula = Surv(time, event) ~ log(creatinine + 1) + age + 
    ##     ejection_fraction + bp + anaemia + sodium, data = stepwise_data)
    ## 
    ##   n= 299, number of events= 96 
    ## 
    ##                          coef exp(coef)  se(coef)      z Pr(>|z|)    
    ## log(creatinine + 1)  1.307631  3.697404  0.293779  4.451 8.54e-06 ***
    ## age                  0.041664  1.042545  0.009044  4.607 4.09e-06 ***
    ## ejection_fraction   -0.042956  0.957954  0.010150 -4.232 2.31e-05 ***
    ## bp                   0.507442  1.661037  0.211818  2.396   0.0166 *  
    ## anaemia              0.414917  1.514245  0.209599  1.980   0.0478 *  
    ## sodium              -0.036640  0.964023  0.024087 -1.521   0.1282    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                     exp(coef) exp(-coef) lower .95 upper .95
    ## log(creatinine + 1)     3.697     0.2705    2.0789    6.5760
    ## age                     1.043     0.9592    1.0242    1.0612
    ## ejection_fraction       0.958     1.0439    0.9391    0.9772
    ## bp                      1.661     0.6020    1.0967    2.5158
    ## anaemia                 1.514     0.6604    1.0041    2.2835
    ## sodium                  0.964     1.0373    0.9196    1.0106
    ## 
    ## Concordance= 0.734  (se = 0.027 )
    ## Likelihood ratio test= 79.39  on 6 df,   p=5e-15
    ## Wald test            = 85.86  on 6 df,   p=<2e-16
    ## Score (logrank) test = 85.53  on 6 df,   p=3e-16

``` r
# Check assumptions with model obtained from stepwiseCox
cox_step4 <- cox.zph(step_model4)
plot(cox_step4) # Residual plots
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-7-5.png)<!-- -->

``` r
# Plot survival curves
ggsurvplot(survfit(step_model4), data = data, conf.int = TRUE)
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-7-6.png)<!-- -->![](Cox_Model_files/figure-gfm/unnamed-chunk-7-7.png)<!-- -->

**Schoenfeld residuals**

``` r
colon_coxph <- coxph(Surv(time, event)~ logcre+age + ejection_fraction + bp +
                       anaemia + sodium, data = stepwise_data)
ggcoxzph(cox.zph(colon_coxph), var = c("logcre"), df = 2, nsmo = 1000)
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggcoxzph(cox.zph(colon_coxph), var = c("age"), df = 2, nsmo = 1000)
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
ggcoxzph(cox.zph(colon_coxph), var = c("ejection_fraction"), df = 2, nsmo = 1000)
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
ggcoxzph(cox.zph(colon_coxph), var = c("bp"), df = 2, nsmo = 1000)
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

``` r
ggcoxzph(cox.zph(colon_coxph), var = c("anaemia"), df = 2, nsmo = 1000)
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->

``` r
ggcoxzph(cox.zph(colon_coxph), var = c("sodium"), df = 2, nsmo = 1000)
```

![](Cox_Model_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->