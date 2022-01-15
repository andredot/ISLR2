ISLR Ch.7
================

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code.

### ex. 1

### ex. 2

### ex. 3

### ex. 4

### ex. 5

**Consider two curves, ˆg1 and ˆg2, defined by \[equation on p.323\]. where g(m) represents the mth derivative of g.**

**(a) As λ →∞, will ˆg1 or ˆg2 have the smaller training RSS?** g1, because since it has less costraints on the shape of the fitting line, it an provide a better fit on training data

**(b) As λ →∞, will ˆg1 or ˆg2 have the smaller test RSS?** Unsure, it really depends on the data. In a setting where more flexibility is required, g1 will probably behave better.

**(c) For λ = 0, will ˆg1 or ˆg2 have the smaller training and test RSS?** They will have exactly the same training and test error since the formula simplifies to the MSE for linear regression.

### ex. 6

**In this exercise, you will further analyze the Wage data set considered throughout this chapter.**

**(a) Perform polynomial regression to predict wage using age. Use cross-validation to select the optimal degree d for the polynomial. What degree was chosen, and how does this compare to the results of hypothesis testing using ANOVA? Make a plot of the resulting polynomial fit to the data.**

``` r
set.seed(123)
wage <- Wage
cv_error <- rep(NA, 5)

for (i in 1:5) { 
  lm_fit <- glm(wage ∼ poly(age , i),
                data = wage)
  cv_error[i] <- cv.glm(wage , lm_fit , K = 10)$delta[1]
}
cv_error
```

    ## [1] 1676.988 1602.473 1597.036 1594.582 1594.424

``` r
ggplot() +
  geom_line( aes( x = 1:5,
                  y = cv_error))
```

![](chapter_7_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
fit_1 <- lm(wage ∼ age , data = wage)
fit_2 <- lm(wage ∼ poly(age , 2), data = wage)
fit_3 <- lm(wage ∼ poly(age , 3), data = wage)
fit_4 <- lm(wage ∼ poly(age , 4), data = wage)
fit_5 <- lm(wage ∼ poly(age , 5), data = wage)

(anova_fit <- anova(fit_1, fit_2, fit_3, fit_4, fit_5) )
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: wage ~ age
    ## Model 2: wage ~ poly(age, 2)
    ## Model 3: wage ~ poly(age, 3)
    ## Model 4: wage ~ poly(age, 4)
    ## Model 5: wage ~ poly(age, 5)
    ##   Res.Df     RSS Df Sum of Sq        F    Pr(>F)    
    ## 1   2998 5022216                                    
    ## 2   2997 4793430  1    228786 143.5931 < 2.2e-16 ***
    ## 3   2996 4777674  1     15756   9.8888  0.001679 ** 
    ## 4   2995 4771604  1      6070   3.8098  0.051046 .  
    ## 5   2994 4770322  1      1283   0.8050  0.369682    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

As we can see from the CV, little improvement can be seen after the 3rd degree polinomial, and ANOVA shows the same thing, suggesting no better model in the 4th compared to the 3rd.

``` r
wage %>%
  select(wage, age) %>% 
  mutate( wage_pred_3 = predict(fit_3)) %>% 
  ggplot() +
    geom_jitter( aes( x = age,
                     y = wage),
                alpha = 0.2) +
    geom_line( aes( x = age,
                    y = wage_pred_3),
               color = "red")
```

![](chapter_7_files/figure-markdown_github/unnamed-chunk-4-1.png)

**(b) Fit a step function to predict wage using age, and perform crossvalidation to choose the optimal number of cuts. Make a plot of the fit obtained.**

``` r
n <- 20
cut_error <- rep(NA, n)

for (i in 2:n) { 
  lm_fit <- glm(wage ∼ cut(age , i, labels = FALSE),
                data = wage)
  cut_error[i] <- cv.glm(wage , lm_fit , K = 10)$delta[1]
}

#cut_error
ggplot() +
  geom_line( aes( x = 1:n,
                  y = cut_error)) +
  coord_cartesian( ylim  = c(1500,max(cut_error)))
```

    ## Warning: Removed 1 row(s) containing missing values (geom_path).

![](chapter_7_files/figure-markdown_github/unnamed-chunk-5-1.png)

Any break more than 9 does not seem to have any significant improvement, so it can be chosen to plot

``` r
breaks <- 9
lm_fit <- lm( wage ∼ cut(age , breaks),
                data = wage)
summary(lm_fit)
```

    ## 
    ## Call:
    ## lm(formula = wage ~ cut(age, breaks), data = wage)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -98.548 -24.896  -4.754  15.590 201.404 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                   73.344      3.030  24.208  < 2e-16 ***
    ## cut(age, breaks)(24.9,31.8]   25.386      3.619   7.014 2.84e-12 ***
    ## cut(age, breaks)(31.8,38.7]   39.746      3.485  11.406  < 2e-16 ***
    ## cut(age, breaks)(38.7,45.6]   45.289      3.404  13.306  < 2e-16 ***
    ## cut(age, breaks)(45.6,52.4]   43.594      3.443  12.660  < 2e-16 ***
    ## cut(age, breaks)(52.4,59.3]   46.060      3.640  12.655  < 2e-16 ***
    ## cut(age, breaks)(59.3,66.2]   45.605      4.391  10.385  < 2e-16 ***
    ## cut(age, breaks)(66.2,73.1]   26.230      7.173   3.657  0.00026 ***
    ## cut(age, breaks)(73.1,80.1]   21.217     11.522   1.841  0.06566 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 40.08 on 2991 degrees of freedom
    ## Multiple R-squared:  0.0799, Adjusted R-squared:  0.07744 
    ## F-statistic: 32.47 on 8 and 2991 DF,  p-value: < 2.2e-16

``` r
## one way to extract the breakpoints
labs <- levels(cut(wage[,"age"],breaks))
levels <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
                upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ),
                intercept = as.numeric(lm_fit[["coefficients"]]))

# add 0 intercept
levels[2:9,"intercept"] <- levels[2:breaks,"intercept"] + levels[1,"intercept"]

ggplot() +
  geom_jitter( data = wage,
              aes( x = age,
                   y = wage),
              alpha = 0.2) +
  geom_segment( aes(x = levels[,"lower"],
                    y = levels[,"intercept"],
                    xend = levels[,"upper"],
                    yend = levels[,"intercept"]),
                color = "red", size = 1) 
```

![](chapter_7_files/figure-markdown_github/unnamed-chunk-6-1.png)

### ex. 7

**The Wage data set contains a number of other features not explored in this chapter, such as marital status (maritl), job class (jobclass), and others. Explore the relationships between some of these other predictors and wage, and use non-linear fitting techniques in order to fit flexible models to the data. Create plots of the results obtained, and write a summary of your findings.**
