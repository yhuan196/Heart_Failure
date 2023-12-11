Paramatric Model
================
Yi Huang
2023-12-08

## Parametric models

Weibull model selection:  
Df AIC  
- \<none\>  1273.2  
- sodium 1 1273.5  
- bp 1 1277.8  
- logcre 1 1282.9  
- ef_cat 2 1294.7  
- age 1 1298.1  
  
  
Gompertz model selection:  
DF AIC  
- \<none\>  1273.6  
- anaemia 1 1273.6  
- sodium 1 1274.5  
- bp 1 1278.2  
- logcre 1 1283.6  
- ef_cat 2 1294.8  
- age 1 1298.9  

## Weibull

``` r
fit_weibull = phreg(Surv(time, event) ~ gender + smoking + diabetes + bp + anaemia + 
                  age + ef_cat+ sodium + pletelets + logcre + logcpk, 
                data = model_data, dist = "weibull")
summary(fit_weibull)
```

    ## Covariate             Mean       Coef     Rel.Risk   S.E.    LR p
    ## gender                                                      0.4785 
    ##                0      0.356     0         1 (reference)
    ##                1      0.644    -0.177     0.838     0.250
    ## smoking                                                     0.5645 
    ##                0      0.685     0         1 (reference)
    ##                1      0.315     0.147     1.158     0.254
    ## diabetes                                                    0.3344 
    ##                0      0.572     0         1 (reference)
    ##                1      0.428     0.217     1.243     0.224
    ## bp                                                          0.0151 
    ##                0      0.705     0         1 (reference)
    ##                1      0.295     0.534     1.705     0.215
    ## anaemia                                                     0.1541 
    ##                0      0.610     0         1 (reference)
    ##                1      0.390     0.306     1.358     0.214
    ## age                  59.251     0.051     1.052     0.010   0.0000 
    ## ef_cat                                                      0.0000 
    ##              Low      0.270     0         1 (reference)
    ##           Medium      0.481    -1.165     0.312     0.255
    ##             High      0.249    -0.973     0.378     0.280
    ## sodium              136.855    -0.041     0.960     0.023   0.0898 
    ## pletelets        263968.682    -0.000     1.000     0.000   0.5375 
    ## logcre                0.792     1.170     3.221     0.297   0.0003 
    ## logcpk                5.732     0.120     1.127     0.101   0.2359 
    ## 
    ## Events                    96 
    ## Total time at risk         38948 
    ## Max. log. likelihood      -626.14 
    ## LR test statistic         88.60 
    ## Degrees of freedom        12 
    ## Overall p-value           9.22595e-14

``` r
# AIC Both direction
library(MASS)
stepAIC(fit_weibull)
```

    ## Start:  AIC=1280.28
    ## Surv(time, event) ~ gender + smoking + diabetes + bp + anaemia + 
    ##     age + ef_cat + sodium + pletelets + logcre + logcpk
    ## 
    ##             Df    AIC
    ## - smoking    1 1278.6
    ## - pletelets  1 1278.7
    ## - gender     1 1278.8
    ## - diabetes   1 1279.2
    ## - logcpk     1 1279.7
    ## <none>         1280.3
    ## - anaemia    1 1280.3
    ## - sodium     1 1281.2
    ## - bp         1 1284.2
    ## - logcre     1 1291.1
    ## - ef_cat     2 1302.0
    ## - age        1 1305.7
    ## 
    ## Step:  AIC=1278.61
    ## Surv(time, event) ~ gender + diabetes + bp + anaemia + age + 
    ##     ef_cat + sodium + pletelets + logcre + logcpk
    ## 
    ##             Df    AIC
    ## - gender     1 1276.9
    ## - pletelets  1 1276.9
    ## - diabetes   1 1277.5
    ## - logcpk     1 1277.9
    ## - anaemia    1 1278.5
    ## <none>         1278.6
    ## - sodium     1 1279.3
    ## - bp         1 1282.7
    ## - logcre     1 1289.2
    ## - ef_cat     2 1300.8
    ## - age        1 1303.7
    ## 
    ## Step:  AIC=1276.88
    ## Surv(time, event) ~ diabetes + bp + anaemia + age + ef_cat + 
    ##     sodium + pletelets + logcre + logcpk
    ## 
    ##             Df    AIC
    ## - pletelets  1 1275.1
    ## - diabetes   1 1275.8
    ## - logcpk     1 1276.1
    ## - anaemia    1 1276.7
    ## <none>         1276.9
    ## - sodium     1 1277.6
    ## - bp         1 1281.3
    ## - logcre     1 1287.8
    ## - ef_cat     2 1298.8
    ## - age        1 1301.7
    ## 
    ## Step:  AIC=1275.1
    ## Surv(time, event) ~ diabetes + bp + anaemia + age + ef_cat + 
    ##     sodium + logcre + logcpk
    ## 
    ##            Df    AIC
    ## - diabetes  1 1274.0
    ## - logcpk    1 1274.3
    ## - anaemia   1 1275.0
    ## <none>        1275.1
    ## - sodium    1 1275.9
    ## - bp        1 1279.4
    ## - logcre    1 1286.5
    ## - ef_cat    2 1296.8
    ## - age       1 1299.8
    ## 
    ## Step:  AIC=1274.01
    ## Surv(time, event) ~ bp + anaemia + age + ef_cat + sodium + logcre + 
    ##     logcpk
    ## 
    ##           Df    AIC
    ## - logcpk   1 1273.2
    ## <none>       1274.0
    ## - anaemia  1 1274.1
    ## - sodium   1 1275.2
    ## - bp       1 1278.3
    ## - logcre   1 1284.6
    ## - ef_cat   2 1295.4
    ## - age      1 1297.8
    ## 
    ## Step:  AIC=1273.25
    ## Surv(time, event) ~ bp + anaemia + age + ef_cat + sodium + logcre
    ## 
    ##           Df    AIC
    ## - anaemia  1 1273.2
    ## <none>       1273.2
    ## - sodium   1 1274.1
    ## - bp       1 1277.5
    ## - logcre   1 1283.1
    ## - ef_cat   2 1294.0
    ## - age      1 1297.0
    ## 
    ## Step:  AIC=1273.16
    ## Surv(time, event) ~ bp + age + ef_cat + sodium + logcre
    ## 
    ##          Df    AIC
    ## <none>      1273.2
    ## - sodium  1 1273.5
    ## - bp      1 1277.8
    ## - logcre  1 1282.9
    ## - ef_cat  2 1294.7
    ## - age     1 1298.1

    ## Call:
    ## phreg(formula = Surv(time, event) ~ bp + age + ef_cat + sodium + 
    ##     logcre, data = model_data, dist = "weibull")
    ## 
    ## Covariate          W.mean      Coef Exp(Coef)  se(Coef)    Wald p
    ## bp 
    ##                0    0.705     0         1           (reference)
    ##                1    0.295     0.559     1.749     0.213     0.009 
    ## age                59.251     0.048     1.049     0.009     0.000 
    ## ef_cat 
    ##              Low    0.270     0         1           (reference)
    ##           Medium    0.481    -1.145     0.318     0.249     0.000 
    ##             High    0.249    -0.978     0.376     0.278     0.000 
    ## sodium            136.855    -0.038     0.963     0.024     0.117 
    ## logcre              0.792     1.110     3.034     0.292     0.000 
    ## 
    ## log(scale)                    4.361               3.552     0.220 
    ## log(shape)                   -0.061               0.089     0.491 
    ## 
    ## Events                    96 
    ## Total time at risk         38948 
    ## Max. log. likelihood      -628.58 
    ## LR test statistic         83.72 
    ## Degrees of freedom        6 
    ## Overall p-value           5.55112e-16

``` r
final_weibull <- phreg(formula = Surv(time, event) ~ bp + age + ef_catMedium + 
                          ef_catHigh + sodium + logcre, data = stepwise_data, dist = "weibull")

final_weibull
```

    ## Call:
    ## phreg(formula = Surv(time, event) ~ bp + age + ef_catMedium + 
    ##     ef_catHigh + sodium + logcre, data = stepwise_data, dist = "weibull")
    ## 
    ## Covariate          W.mean      Coef Exp(Coef)  se(Coef)    Wald p
    ## bp                  0.295     0.559     1.749     0.213     0.009 
    ## age                59.251     0.048     1.049     0.009     0.000 
    ## ef_catMedium        0.481    -1.145     0.318     0.249     0.000 
    ## ef_catHigh          0.249    -0.978     0.376     0.278     0.000 
    ## sodium            136.855    -0.038     0.963     0.024     0.117 
    ## logcre              0.792     1.110     3.034     0.292     0.000 
    ## 
    ## log(scale)                    4.361               3.552     0.220 
    ## log(shape)                   -0.061               0.089     0.491 
    ## 
    ## Events                    96 
    ## Total time at risk         38948 
    ## Max. log. likelihood      -628.58 
    ## LR test statistic         83.72 
    ## Degrees of freedom        6 
    ## Overall p-value           5.55112e-16

``` r
1 - pchisq((0.477/0.212876134)^2, df = 1)
```

    ## [1] 0.02504294

``` r
# Extract relevant information
weibull_summary <- data.frame(
  Variable= c(names(final_weibull$coefficients)),
  Coef = round(as.vector(final_weibull$coefficients),4),
  `Exp(Coef)` = round(exp(as.vector(final_weibull$coefficients)),4),
  `se(Coef)` = round(sqrt(as.vector(diag(final_weibull[["var"]]))),4),
  `Wald p` = round(1-pchisq((as.vector(final_weibull$coefficients)/sqrt(as.vector(diag(final_weibull[["var"]]))))^2,1),4))

# Create table using kable
kable(as.data.frame(weibull_summary), caption = "Summary of Gompertz PH Model Fitting", digits = 4) 
```

<table>
<caption>
Summary of Gompertz PH Model Fitting
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
Coef
</th>
<th style="text-align:right;">
Exp.Coef.
</th>
<th style="text-align:right;">
se.Coef.
</th>
<th style="text-align:right;">
Wald.p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
bp
</td>
<td style="text-align:right;">
0.5591
</td>
<td style="text-align:right;">
1.7492
</td>
<td style="text-align:right;">
0.2132
</td>
<td style="text-align:right;">
0.0087
</td>
</tr>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
0.0480
</td>
<td style="text-align:right;">
1.0492
</td>
<td style="text-align:right;">
0.0092
</td>
<td style="text-align:right;">
0.0000
</td>
</tr>
<tr>
<td style="text-align:left;">
ef_catMedium
</td>
<td style="text-align:right;">
-1.1446
</td>
<td style="text-align:right;">
0.3183
</td>
<td style="text-align:right;">
0.2492
</td>
<td style="text-align:right;">
0.0000
</td>
</tr>
<tr>
<td style="text-align:left;">
ef_catHigh
</td>
<td style="text-align:right;">
-0.9780
</td>
<td style="text-align:right;">
0.3761
</td>
<td style="text-align:right;">
0.2780
</td>
<td style="text-align:right;">
0.0004
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium
</td>
<td style="text-align:right;">
-0.0376
</td>
<td style="text-align:right;">
0.9631
</td>
<td style="text-align:right;">
0.0240
</td>
<td style="text-align:right;">
0.1174
</td>
</tr>
<tr>
<td style="text-align:left;">
logcre
</td>
<td style="text-align:right;">
1.1099
</td>
<td style="text-align:right;">
3.0342
</td>
<td style="text-align:right;">
0.2924
</td>
<td style="text-align:right;">
0.0001
</td>
</tr>
<tr>
<td style="text-align:left;">
log(scale)
</td>
<td style="text-align:right;">
4.3611
</td>
<td style="text-align:right;">
78.3415
</td>
<td style="text-align:right;">
3.5518
</td>
<td style="text-align:right;">
0.2195
</td>
</tr>
<tr>
<td style="text-align:left;">
log(shape)
</td>
<td style="text-align:right;">
-0.0613
</td>
<td style="text-align:right;">
0.9405
</td>
<td style="text-align:right;">
0.0890
</td>
<td style="text-align:right;">
0.4911
</td>
</tr>
</tbody>
</table>

``` r
# |> kable_styling(latex_options = c("striped", "hold_position"))
```

## Gompertz

``` r
fit_gompertz = phreg(Surv(time, event) ~ gender + smoking + diabetes + bp + anaemia + 
                  age + ef_cat+ sodium + pletelets + logcre + logcpk, 
                data = model_data, dist = "gompertz")
fit_gompertz
```

    ## Call:
    ## phreg(formula = Surv(time, event) ~ gender + smoking + diabetes + 
    ##     bp + anaemia + age + ef_cat + sodium + pletelets + logcre + 
    ##     logcpk, data = model_data, dist = "gompertz")
    ## 
    ## Covariate          W.mean      Coef Exp(Coef)  se(Coef)    Wald p
    ## gender 
    ##                0    0.356     0         1           (reference)
    ##                1    0.644    -0.175     0.839 NA        NA        
    ## smoking 
    ##                0    0.685     0         1           (reference)
    ##                1    0.315     0.148     1.159 NA        NA        
    ## diabetes 
    ##                0    0.572     0         1           (reference)
    ##                1    0.428     0.220     1.247 NA        NA        
    ## bp 
    ##                0    0.705     0         1           (reference)
    ##                1    0.295     0.549     1.732 NA        NA        
    ## anaemia 
    ##                0    0.610     0         1           (reference)
    ##                1    0.390     0.315     1.370 NA        NA        
    ## age                59.251     0.052     1.053 NA        NA        
    ## ef_cat 
    ##              Low    0.270     0         1           (reference)
    ##           Medium    0.481    -1.174     0.309 NA        NA        
    ##             High    0.249    -0.979     0.376 NA        NA        
    ## sodium            136.855    -0.041     0.960 NA        NA        
    ## pletelets        263968.682    -0.000     1.000 NA        NA        
    ## logcre              0.792     1.178     3.247 NA        NA        
    ## logcpk              5.732     0.120     1.128 NA        NA        
    ## 
    ## log(scale)                   16.674           NA        NA        
    ## log(shape)                   11.772           NA        NA        
    ## 
    ## Events                    96 
    ## Total time at risk         38948 
    ## Max. log. likelihood      -626.29 
    ## LR test statistic         92.50 
    ## Degrees of freedom        12 
    ## Overall p-value           1.62093e-14

``` r
# AIC Both direction
stepAIC(fit_gompertz, direction = "backward")
```

    ## Start:  AIC=1280.58
    ## Surv(time, event) ~ gender + smoking + diabetes + bp + anaemia + 
    ##     age + ef_cat + sodium + pletelets + logcre + logcpk
    ## 
    ##             Df    AIC
    ## - smoking    1 1278.9
    ## - pletelets  1 1279.0
    ## - gender     1 1279.1
    ## - diabetes   1 1279.5
    ## - logcpk     1 1280.0
    ## <none>         1280.6
    ## - anaemia    1 1280.7
    ## - sodium     1 1281.5
    ## - bp         1 1284.9
    ## - logcre     1 1291.5
    ## - ef_cat     2 1302.6
    ## - age        1 1307.6
    ## 
    ## Step:  AIC=1278.92
    ## Surv(time, event) ~ gender + diabetes + bp + anaemia + age + 
    ##     ef_cat + sodium + pletelets + logcre + logcpk
    ## 
    ##             Df    AIC
    ## - gender     1 1277.2
    ## - pletelets  1 1277.2
    ## - diabetes   1 1277.8
    ## - logcpk     1 1278.2
    ## - anaemia    1 1278.9
    ## <none>         1278.9
    ## - sodium     1 1279.7
    ## - bp         1 1283.5
    ## - logcre     1 1289.7
    ## - ef_cat     2 1301.5
    ## - age        1 1305.6
    ## 
    ## Step:  AIC=1277.17
    ## Surv(time, event) ~ diabetes + bp + anaemia + age + ef_cat + 
    ##     sodium + pletelets + logcre + logcpk
    ## 
    ##             Df    AIC
    ## - pletelets  1 1275.4
    ## - diabetes   1 1276.2
    ## - logcpk     1 1276.4
    ## - anaemia    1 1277.2
    ## <none>         1277.2
    ## - sodium     1 1277.9
    ## - bp         1 1282.0
    ## - logcre     1 1288.3
    ## - ef_cat     2 1299.5
    ## - age        1 1303.6
    ## 
    ## Step:  AIC=1275.39
    ## Surv(time, event) ~ diabetes + bp + anaemia + age + ef_cat + 
    ##     sodium + logcre + logcpk
    ## 
    ##            Df    AIC
    ## - diabetes  1 1274.3
    ## - logcpk    1 1274.6
    ## - anaemia   1 1275.4
    ## <none>        1275.4
    ## - sodium    1 1276.2
    ## - bp        1 1280.1
    ## - logcre    1 1287.0
    ## - ef_cat    2 1297.6
    ## - age       1 1301.6
    ## 
    ## Step:  AIC=1274.32
    ## Surv(time, event) ~ bp + anaemia + age + ef_cat + sodium + logcre + 
    ##     logcpk
    ## 
    ##           Df    AIC
    ## - logcpk   1 1273.6
    ## <none>       1274.3
    ## - anaemia  1 1274.5
    ## - sodium   1 1275.5
    ## - bp       1 1279.0
    ## - logcre   1 1285.1
    ## - ef_cat   2 1296.2
    ## - age      1 1299.7
    ## 
    ## Step:  AIC=1273.58
    ## Surv(time, event) ~ bp + anaemia + age + ef_cat + sodium + logcre
    ## 
    ##           Df    AIC
    ## <none>       1273.6
    ## - anaemia  1 1273.6
    ## - sodium   1 1274.5
    ## - bp       1 1278.2
    ## - logcre   1 1283.6
    ## - ef_cat   2 1294.8
    ## - age      1 1298.9

    ## Call:
    ## phreg(formula = Surv(time, event) ~ bp + anaemia + age + ef_cat + 
    ##     sodium + logcre, data = model_data, dist = "gompertz")
    ## 
    ## Covariate          W.mean      Coef Exp(Coef)  se(Coef)    Wald p
    ## bp 
    ##                0    0.705     0         1           (reference)
    ##                1    0.295     0.558     1.747     0.212     0.009 
    ## anaemia 
    ##                0    0.610     0         1           (reference)
    ##                1    0.390     0.303     1.354     0.210     0.149 
    ## age                59.251     0.048     1.049     0.009     0.000 
    ## ef_cat 
    ##              Low    0.270     0         1           (reference)
    ##           Medium    0.481    -1.116     0.328     0.250     0.000 
    ##             High    0.249    -1.001     0.368     0.278     0.000 
    ## sodium            136.855    -0.042     0.959     0.024     0.080 
    ## logcre              0.792     1.110     3.035     0.289     0.000 
    ## 
    ## log(scale)                   16.661             163.817     0.919 
    ## log(shape)                   12.781             163.852     0.938 
    ## 
    ## Events                    96 
    ## Total time at risk         38948 
    ## Max. log. likelihood      -627.79 
    ## LR test statistic         89.50 
    ## Degrees of freedom        7 
    ## Overall p-value           1.11022e-16

``` r
final_gompertz <- phreg(formula = Surv(time, event) ~ bp + anaemia + age + ef_catMedium + 
                          ef_catHigh + sodium + logcre, data = stepwise_data, dist = "gompertz")
```

``` r
# Extract relevant information
gompertz_summary <- data.frame(
  Variable= c(names(final_gompertz$coefficients)),
  Coef = round(as.vector(final_gompertz$coefficients),4),
  `Exp(Coef)` = round(exp(as.vector(final_gompertz$coefficients)),4),
  `se(Coef)` = round(sqrt(as.vector(diag(final_gompertz[["var"]]))),4),
  `Wald p` = round(1-pchisq((as.vector(final_gompertz$coefficients)/sqrt(as.vector(diag(final_gompertz[["var"]]))))^2,1),4))

# Create table using kable
kable(as.data.frame(weibull_summary), caption = "Summary of Weibul PH Model Fitting", digits = 4) 
```

<table>
<caption>
Summary of Weibul PH Model Fitting
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
Coef
</th>
<th style="text-align:right;">
Exp.Coef.
</th>
<th style="text-align:right;">
se.Coef.
</th>
<th style="text-align:right;">
Wald.p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
bp
</td>
<td style="text-align:right;">
0.5591
</td>
<td style="text-align:right;">
1.7492
</td>
<td style="text-align:right;">
0.2132
</td>
<td style="text-align:right;">
0.0087
</td>
</tr>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
0.0480
</td>
<td style="text-align:right;">
1.0492
</td>
<td style="text-align:right;">
0.0092
</td>
<td style="text-align:right;">
0.0000
</td>
</tr>
<tr>
<td style="text-align:left;">
ef_catMedium
</td>
<td style="text-align:right;">
-1.1446
</td>
<td style="text-align:right;">
0.3183
</td>
<td style="text-align:right;">
0.2492
</td>
<td style="text-align:right;">
0.0000
</td>
</tr>
<tr>
<td style="text-align:left;">
ef_catHigh
</td>
<td style="text-align:right;">
-0.9780
</td>
<td style="text-align:right;">
0.3761
</td>
<td style="text-align:right;">
0.2780
</td>
<td style="text-align:right;">
0.0004
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium
</td>
<td style="text-align:right;">
-0.0376
</td>
<td style="text-align:right;">
0.9631
</td>
<td style="text-align:right;">
0.0240
</td>
<td style="text-align:right;">
0.1174
</td>
</tr>
<tr>
<td style="text-align:left;">
logcre
</td>
<td style="text-align:right;">
1.1099
</td>
<td style="text-align:right;">
3.0342
</td>
<td style="text-align:right;">
0.2924
</td>
<td style="text-align:right;">
0.0001
</td>
</tr>
<tr>
<td style="text-align:left;">
log(scale)
</td>
<td style="text-align:right;">
4.3611
</td>
<td style="text-align:right;">
78.3415
</td>
<td style="text-align:right;">
3.5518
</td>
<td style="text-align:right;">
0.2195
</td>
</tr>
<tr>
<td style="text-align:left;">
log(shape)
</td>
<td style="text-align:right;">
-0.0613
</td>
<td style="text-align:right;">
0.9405
</td>
<td style="text-align:right;">
0.0890
</td>
<td style="text-align:right;">
0.4911
</td>
</tr>
</tbody>
</table>

``` r
# |> kable_styling(latex_options = c("striped", "hold_position"))
```
