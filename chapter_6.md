ISLR Ch.6
================

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code.

### ex. 1

**We perform best subset, forward stepwise, and backward stepwise selection on a single data set. For each approach, we obtain p+1 models, containing 0, 1, 2,. .., p predictors. Explain your answers:**

**(a) Which of the three models with k predictors has the smallest training RSS?** All three models will have the same training RSS, and it will be the lowest in the model will all the p predictors. Since we are talking about training RSS, it is reduced for every predictor that we add since it can also model part of the noise (this clearly deteriorates the test RSS, but the training RSS is always equal or lower). Since the three models are fit on the same dataset (the complete one) they have also the same coefficients.

**(b) Which of the three models with k predictors has the smallest test RSS?** Since they are trained on the same dataset, the smallest test RSS will be in the best subse model (and in the forward or backward ones only if they suggest the same subset of predictors). Since it explores the whole "predictor space" the best subset is definitely slower but assured to find the best one, while the other rely on certain assumptions to reduce the searched combinations.

**(c) True or False:**

**i. The predictors in the k-variable model identified by forward stepwise are a subset of the predictors in the (k+1)-variable model identified by forward stepwise selection.** True, there is no substition of predictors, only addition

\*ii. The predictors in the k-variable model identified by backward stepwise are a subset of the predictors in the (k + 1)variable model identified by backward stepwise selection.\*\* True, only one of them was eliminated in the last step

**iii. The predictors in the k-variable model identified by backward stepwise are a subset of the predictors in the (k + 1)variable model identified by forward stepwise selection.** No, we can not be sure about it (it may be, it may be not)

**iv. The predictors in the k-variable model identified by forward stepwise are a subset of the predictors in the (k+1)-variable model identified by backward stepwise selection.** No, we can not be sure about it (it may be, it may be not)

**v. The predictors in the k-variable model identified by best subset are a subset of the predictors in the (k + 1)-variable model identified by best subset selection.** No, the combination selected may change if the best subset approach is chosen.

### ex. 2

**For parts (a) through (c), indicate which of i. through iv. is correct. Justify your answer.**

**(a) The lasso, relative to least squares, is:**

**i. More flexible and hence will give improved prediction accuracy when its increase in bias is less than its decrease in variance.**

**ii. More flexible and hence will give improved prediction accuracy when its increase in variance is less than its decrease in bias.**

**iii. Less flexible and hence will give improved prediction accuracy when its increase in bias is less than its decrease in variance.**

**iv. Less flexible and hence will give improved prediction accuracy when its increase in variance is less than its decrease in bias.**

The correct option is (iii). Lasso is less flexible, but selecting the proper level of inflexibility will reduce the variance of the model while adding only a moderate quantity of bias, especially when only a few predictors contain most of the signal to correctly fit the model (and the others are shrunk to zero).

**(b) The ridge regression, relative to least squares, is:** (iii) since Ridge regression is very similar (even if no coefficient is brought to zero), where the formula to minimize in order to fit the regression only differ from |betas| (lasso) to (betas)^2 (in the ridge)

**(c) Non-linear methods, relative to least squares, are:** (i) since they are able to better model complex relationships between prediction and response. They are superior when the linear model has a strong bias and a flexible approach can reduce it, with only a modest increase in variance compared to it

### ex. 3

**Suppose we estimate the regression coefficients in a linear regression model by minimizing \[equation on p.283\] for a particular value of s. For parts (a) through (e), indicate which of i. through v. is correct. Justify your answer.**

**i. Increase initially, and then eventually start decreasing in an inverted U shape**

**ii. Decrease initially, and then eventually start increasing in a U shape**

**iii. Steadily increase**

**iv. Steadily decrease**

**v. Remain constant**

**(a) As we increase s from 0, the training RSS will:** (iii) since in the training set the best possible fit is the one of the standard linear regression. Penalizing the coefficients, shrinking them, will increare progressively the training RSS

**(b) As we increase s from 0, the test RSS will:** (ii) since (unless the best coefficients are exactly zero), the initial decrease in the variance provided by the less flexible model compesates the increase in bias introduced by the lower betas. The reduction continues decreases until a point where the bias is growing more rapidly than the decrease in variance and test RSS starts to grow again

**(c) As we increase s from 0, the variance will:** (iv) as the model decrease its flexibility (and the betas are shrunk toward zero)

**(d) As we increase s from 0, the squared bias will:** (iii) increase as the model inflexibility prevent the model from estimating with more precision.

**(e) As we increase s from 0, the irreducible error will:** (v) since it does not depend on the model

### ex. 4

Same answers as in 5 since the different equation does not modify any of the previous considerations.

### ex. 5

**It is well-known that ridge regression tends to give similar coefficient values to correlated variables, whereas the lasso may give quite different coefficient values to correlated variables. We will now explore this property in a very simple setting. Suppose that n = 2, p = 2, x11 = x12, x21 = x22. Furthermore, suppose that y1+y2 = 0 and x11+x21 = 0 and x12+x22 = 0, so that the estimate for the intercept in a least squares, ridge regression, or lasso model is zero: ˆβ0 = 0.**

``` r
   | p1 | p2 | y
-----------------
 1 | x11| x12| y1
 2 |-x11|-x12|-y1
 
 beta0 = 0
```

**(a) Write out the ridge regression optimization problem in this setting.**

``` r
minimize( (y1 - (b1*x11 + b2*x12))^2 + (-y1 - (b1*(-x11) + b2*(-x12) ))^2 )
  with ( |b1|+|b2| )<s
```

**(b) Argue that in this setting, the ridge coefficient estimates satisfy βˆ1 = ˆβ2.**

``` r
# simplify (plase note i forgot x11=x12)
(y1 - b1*x11 - b2*x12)^2 + ( -y1 + b1*-x11 + b2*-x12 )^2

# rewrite 2nd term
(y1 - b1*x11 - b2*x12)^2 + (-(y1 - b1*-x11 - b2*-x12) )^2 

# we can remove the - in our minim. problem
2(y1 - b1*-x11 - b2*-x12)^2

#which is equivalent to minimize
y1 - b1*-x11 - b2*-x12
# or
y1 + b1*x11 + b2*x12

since x11 = x12, the coefficient should be the same
```

**(c) Write out the lasso optimization problem in this setting.**

``` r
   | p1 | p2 | y
-----------------
 1 | x11| x11| y1
 2 |-x11|-x11|-y1
 
 beta0 = 0
```

``` r
minimize( (y1 - (b1*x11 + b2*x12))^2 + (-y1 - (b1*(-x11) + b2*(-x12) ))^2 )
  with (b1+b2)^2 < s
```

**(d) Argue that in this setting, the lasso coefficients βˆ1 and βˆ2 are not unique—in other words, there are many possible solutions to the optimization problem in (c). Describe these solutions**

``` r
# simplify as in previous (b)
y1 + b1*x11 + b2*x12 with (b1+b2)^2 < s

# so we have two possible solution b1 = b2 and b1 = -b2
```

Help (unsure about conclusions)

### ex. 6

**6. We will now explore (6.12) and (6.13) further.**

**(a) Consider (6.12) with p = 1. For some choice of y1 and λ &gt; 0, plot (6.12) as a function of β1. Your plot should confirm that (6.12) is solved by (6.14).**

``` r
beta <- seq(from = -5, to = 5, by = 0.1)
lambda <- 10

# formula with p=1: (y - beta)^2 + lambda*(beta^2)
# which is: y^2 + y*(-2*beta) + (1 + lambda)*beta^2
a <- rep(1., 401)
b <- -2*beta
c <- (1 + lambda)*beta^2

quad <- function(a, b, c)
{
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  if(all(Im(answer) == 0)) answer <- Re(answer)
  if(answer[1] == answer[2]) return(answer[1])
  answer
}
```

``` r
y <- quad(a,b,c)
ggplot() +
  geom_point( aes( x = beta,
                   y = quad(a, b, c)[,1]))
```

**(b) Consider (6.13) with p = 1. For some choice of y1 and λ &gt; 0, plot (6.13) as a function of β1. Your plot should confirm that (6.13) is solved by (6.15)**

Help on the simulation part to prove solution

### ex. 7

Help

### ex. 8

**In this exercise, we will generate simulated data, and will then use this data to perform best subset selection.**

**(a) Use the rnorm() function to generate a predictor X of length n = 100, as well as a noise vector ϵ of length n = 100.**

``` r
set.seed(123)
x <- rnorm( n = 100)
e <- rnorm( n = 100)
```

**(b) Generate a response vector Y of length n = 100 according to the model Y = β0 + β1X + β2X2 + β3X3 + ϵ, where β0, β1, β2, and β3 are constants of your choice.**

``` r
beta0 <- 1.5
beta1 <- 1
beta2 <- 0.05
beta3 <- 0.8

y <- beta0 + beta1*x + beta2*x^2 + beta3*x^3 + e

d <- data.frame(matrix(0, nrow = length(x), ncol = 10))
for (i in 1:length(x)){
  for (j in 1:10) d[i,j] <- x[i]^j
}

d <- add_column(d, y)
d %>%
  gather(-y, key = "var", value = "value") %>% 
  ggplot(aes(x = value,
             y = y)) +
    geom_point() +
    facet_wrap(~ var, scales = "free")
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-11-1.png)

**(c) Use the regsubsets() function to perform best subset selection in order to choose the best model containing the predictors X,X2,. ..,X10.**

``` r
regfit_full <- regsubsets(y ∼ ., 
                          data = d,
                          nvmax = 9)
regfit_full_sum <- summary(regfit_full)
```

**What is the best model obtained according to Cp, BIC, and adjusted R2? Show some plots to provide evidence for your answer, and report the coefficients of the best model obtained. Note you will need to use the data.frame() function to create a single data set containing both X and Y.**

``` r
plot_reg <- function(regfit_sum) {
  par(mfrow = c(2, 2))
  
  plot(regfit_sum$rss , 
       xlab = "Number of Variables", 
       ylab = "RSS",
       type = "l")
  
  plot(regfit_sum$adjr2 , 
       xlab = "Number of Variables", 
       ylab = "Adjusted RSq",
       type = "l")
  points(which.max(regfit_sum$adjr2), 
         regfit_sum$adjr2 [which.max(regfit_sum$adjr2)], 
         col = "red",
         cex = 2, 
         pch = 20)
  
  plot(regfit_sum$cp, 
       xlab = "Number of Variables", 
       ylab = "Cp",
       type = "l")
  points(which.min(regfit_sum$cp),
         regfit_sum$cp[which.min(regfit_sum$cp)], 
         col = "red",
         cex = 2, 
         pch = 20)
  
  plot(regfit_sum$bic , 
       xlab = "Number of Variables", 
       ylab = "BIC",
       type = "l")
  points(which.min(regfit_sum$bic), 
         regfit_sum$bic[which.min(regfit_sum$bic)], 
         col = "red",
         cex = 2, 
         pch = 20)
}

plot_reg(regfit_full_sum)
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-13-1.png)

The best model selected is the one with two variables.

``` r
coef(regfit_full, 2)
```

    ## (Intercept)          X1          X3 
    ##   1.4368231   0.9173009   0.8176441

**(d) Repeat (c), using forward stepwise selection and also using backwards stepwise selection. How does your answer compare to the results in (c)?**

``` r
# forward
regfit_fwd <- regsubsets(y ∼ ., 
                         data = d,
                         nvmax = 9,
                         method = "forward")
regfit_fwd_sum <- summary(regfit_fwd)
plot_reg(regfit_fwd_sum)
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
coef(regfit_fwd, 2)
```

    ## (Intercept)          X1          X3 
    ##   1.4368231   0.9173009   0.8176441

``` r
# backward
regfit_bwd <- regsubsets(y ∼ ., 
                         data = d,
                         nvmax = 9,
                         method = "backward")
regfit_bwd_sum <- summary(regfit_bwd)
plot_reg(regfit_bwd_sum)
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
coef(regfit_fwd, 3)
```

    ## (Intercept)          X1          X2          X3 
    ##  1.47039390  0.92044621 -0.04154309  0.82043630

The best model provided by the forward stepwise selection is the same as the exhaustive one, while the one produced by the backwad is with three variables instea of 2. X1 and X3 are the same of the full model, and X2 is added to them. Clearly, in this last case, even if the X1 and X3 coefficient are very similar,they are not identical.

**(e) Now fit a lasso model to the simulated data, again using X,X2, .. .,X10 as predictors. Use cross-validation to select the optimal value of λ. Create plots of the cross-validation error as a function of λ. Report the resulting coefficient estimates, and discuss the results obtained.**

``` r
# unvalidated model
x <- model.matrix(y ∼ ., d)[, -1] # y is already in y
grid <- 10^seq(from = 10, to = -2, length = 100) 
lasso_mod <- glmnet(x, y, 
                    alpha = 1, 
                    lambda = grid)
```

``` r
#crossvalidated choice of lambda

set.seed(123)
train <- sample(x = 1:dim(x)[1], 
                size = dim(x)[1]/2)
test <- (-train) 
y_test <- y[test]

cv_out <- cv.glmnet( x[train , ], 
                     y[train], 
                     alpha = 1)
( best_lambda <- cv_out[["lambda.min"]] )
```

    ## [1] 0.008677741

``` r
ggplot() +
  geom_pointrange( aes( x = cv_out[["lambda"]],
                        y = cv_out[["cvm"]],
                        ymin = cv_out[["cvlo"]],
                        ymax = cv_out[["cvup"]] )) +
  scale_x_log10() +
  scale_y_log10()
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
# MSE associated with best_lambda
lasso_pred <- predict(lasso_mod , 
                      s = best_lambda , 
                      newx = x[test , ])
mean((lasso_pred - y_test)^2)
```

    ## [1] 0.951781

``` r
# Coef
lasso_coef <- predict(lasso_mod , 
                      type = "coefficients", 
                      s= best_lambda)
lasso_coef[1:11,]
```

    ##   (Intercept)            X1            X2            X3            X4 
    ##  1.502332e+00  9.130216e-01 -1.145606e-01  8.235765e-01  0.000000e+00 
    ##            X5            X6            X7            X8            X9 
    ##  0.000000e+00  0.000000e+00  0.000000e+00  6.487255e-04  0.000000e+00 
    ##           X10 
    ##  8.417303e-05

Best lambda is 0.008, so not very high, but sufficient to shrunk to zero most of the coefficients apart from beta1 and beta3.

**(f) Now generate a response vector Y according to the model Y = β0 + β7X7 + ϵ, and perform best subset selection and the lasso. Discuss the results obtained.**

``` r
beta7 <- 1.5
set.seed(123)
x <- rnorm( n = 100)
y2 <- beta0 + beta7*x^7 + e
d <- d %>% mutate( y = y2)

ggplot(d) + geom_point( aes( x = x,
                             y = y))
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
regfit_full <- regsubsets(y ∼ ., 
                          data = d,
                          nvmax = 9)
regfit_full_sum <- summary(regfit_full)
plot_reg(regfit_full_sum)
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
coef(regfit_full, 1)
```

    ## (Intercept)          X7 
    ##    1.393251    1.499704

``` r
x <- model.matrix(y ∼ ., d)[, -1] # y is already in y
grid <- 10^seq(from = 10, to = -2, length = 100) 
lasso_mod <- glmnet(x, y2, 
                    alpha = 1, 
                    lambda = grid)
#crossvalidated choice of lambda

set.seed(123)
train <- sample(x = 1:dim(x)[1], 
                size = dim(x)[1]/2)
test <- (-train) 
y_test <- y2[test]

cv_out <- cv.glmnet( x[train , ], 
                     y2[train], 
                     alpha = 1)
( best_lambda <- cv_out[["lambda.min"]] )
```

    ## [1] 2.721207

``` r
ggplot() +
  geom_pointrange( aes( x = cv_out[["lambda"]],
                        y = cv_out[["cvm"]],
                        ymin = cv_out[["cvlo"]],
                        ymax = cv_out[["cvup"]] )) +
  scale_x_log10() +
  scale_y_log10()
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-25-1.png)

``` r
# MSE associated with best_lambda
lasso_pred <- predict(lasso_mod , 
                      s = best_lambda , 
                      newx = x[test , ])
mean((lasso_pred - y_test)^2)
```

    ## [1] 4.714608

``` r
# Coef
lasso_coef <- predict(lasso_mod , 
                      type = "coefficients", 
                      s=best_lambda)
lasso_coef[1:11,]
```

    ## (Intercept)          X1          X2          X3          X4          X5 
    ##    1.533337    0.000000    0.000000    0.000000    0.000000    0.000000 
    ##          X6          X7          X8          X9         X10 
    ##    0.000000    1.447825    0.000000    0.000000    0.000000

The best subset selection gave a bizarre result, failing to identify the correct variable, while the lasso shrunk to zero all coefficients apart from the correct beta.

### ex. 9

**In this exercise, we will predict the number of applications received using the other variables in the College data set.**

**(a) Split the data set into a training set and a test set. (b) Fit a linear model using least squares on the training set, and report the test error obtained.**

``` r
set.seed(123)

college <- College
train <- sample( x = 1:dim(college)[1], 
                 size = 0.8*dim(college)[1])
test <- (-train) 
y_test <- college[test, "Apps"]

lm_fit <- lm( Apps ~ .,
              data = college,
              subset = train)
# summary(lm_fit)
y_pred <- predict(lm_fit,
                  interval = "prediction",
                  newdata = college[test,])

errors <- c(lm = mean((y_test - y_pred)^2))
errors["lm"]
```

    ##      lm 
    ## 4610377

**(c) Fit a ridge regression model on the training set, with λ chosen by cross-validation. Report the test error obtained.**

``` r
x <- model.matrix(Apps ∼ .,
                  data = college)[, -1]
grid <- 10^seq(from = 10, to = -2, length = 100) 
ridge_mod <- glmnet(x[train,], 
                    college[train, "Apps"], 
                    alpha = 0, 
                    lambda = grid)

#crossvalidated choice of lambda

cv_out <- cv.glmnet( x[train, ], 
                     college[train, "Apps"], 
                     alpha = 0)
( best_lambda <- cv_out[["lambda.min"]] )
```

    ## [1] 313.5603

``` r
ggplot() +
  geom_pointrange( aes( x = cv_out[["lambda"]],
                        y = cv_out[["cvm"]],
                        ymin = cv_out[["cvlo"]],
                        ymax = cv_out[["cvup"]] )) +
  scale_x_log10() +
  scale_y_log10()
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-28-1.png)

``` r
# MSE associated with best_lambda
ridge_pred <- predict(ridge_mod , 
                      s = best_lambda , 
                      newx = x[test, ])
( errors["ridge"] <- mean((ridge_pred - y_test)^2) )
```

    ## [1] 3938403

``` r
# Coef
ridge_coef <- predict(ridge_mod , 
                      type = "coefficients", 
                      s=best_lambda)
ggplot() +
  geom_point( aes(y = ridge_coef@x,
                 x = ridge_coef@Dimnames[[1]])) +
  coord_flip(ylim = c(-2,+2)) +
  ggtitle( "Ridge MSE = ", mean((ridge_pred - y_test)^2) )
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-29-1.png)

**(d) Fit a lasso model on the training set, with λ chosen by crossvalidation. Report the test error obtained, along with the number of non-zero coefficient estimates.**

``` r
lasso_mod <- glmnet(x[train,], 
                    college[train, "Apps"], 
                    alpha = 1, 
                    lambda = grid)

#crossvalidated choice of lambda

cv_out <- cv.glmnet( x[train, ], 
                     college[train, "Apps"], 
                     alpha = 1)
( best_lambda <- cv_out[["lambda.min"]] )
```

    ## [1] 11.80534

``` r
ggplot() +
  geom_pointrange( aes( x = cv_out[["lambda"]],
                        y = cv_out[["cvm"]],
                        ymin = cv_out[["cvlo"]],
                        ymax = cv_out[["cvup"]] )) +
  scale_x_log10() +
  scale_y_log10()
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-30-1.png)

``` r
# MSE associated with best_lambda
lasso_pred <- predict(lasso_mod , 
                      s = best_lambda , 
                      newx = x[test, ])
( errors["lasso"] <- mean((lasso_pred - y_test)^2) )
```

    ## [1] 2232155

``` r
# Coef
(lasso_coef <- predict(lasso_mod , 
                      type = "coefficients", 
                      s=best_lambda) )
```

    ## 18 x 1 sparse Matrix of class "dgCMatrix"
    ##                        s1
    ## (Intercept) -5.635101e+02
    ## PrivateYes  -5.694126e+02
    ## Accept       1.214179e+00
    ## Enroll       .           
    ## Top10perc    3.670011e+01
    ## Top25perc   -7.221290e+00
    ## F.Undergrad  6.029745e-02
    ## P.Undergrad  .           
    ## Outstate    -4.000401e-02
    ## Room.Board   1.480473e-01
    ## Books        1.038872e-03
    ## Personal     .           
    ## PhD         -4.496929e+00
    ## Terminal    -4.574465e+00
    ## S.F.Ratio    .           
    ## perc.alumni -6.747747e+00
    ## Expend       7.542227e-02
    ## Grad.Rate    8.775119e+00

``` r
ggplot() +
  geom_point( aes(y = lasso_coef@x,
                  x = lasso_coef@i)) +
  coord_flip() +
  ggtitle( "Lasso MSE = ", mean((lasso_pred - y_test)^2) )
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-31-1.png)

**(e) Fit a PCR model on the training set, with M chosen by crossvalidation. Report the test error obtained, along with the value of M selected by cross-validation.**

``` r
pcr_fit <- pcr(Apps ~ .,
               data = college,
               subset = train,
               scale = TRUE,
               validation = "CV")
#summary(pcr_fit)
validationplot(pcr_fit , val.type = "MSEP")
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-32-1.png)

From the plot, 4 components seem to capture most of the information in the dataset

``` r
pcr_pred <- predict(pcr_fit , x[test , ], ncomp = 4) 
( errors["PCAr"] <- mean((pcr_pred - y_test)^2) )
```

    ## [1] 5760389

**(f) Fit a PLS model on the training set, with M chosen by crossvalidation. Report the test error obtained, along with the value of M selected by cross-validation.**

``` r
pls_fit <- plsr(Apps ~ .,
               data = college,
               subset = train,
               scale = TRUE,
               validation = "CV")
#summary(pls_fit)
validationplot(pls_fit , val.type = "MSEP")
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-34-1.png)

After 5 components, it does not seem to be any significant improvement

``` r
pls_pred <- predict(pcr_fit , x[test , ], ncomp = 5) 
( errors["PLSr"] <- mean((pls_pred - y_test)^2) )
```

    ## [1] 5791529

**(g) Comment on the results obtained. How accurately can we predict the number of college applications received? Is there much difference among the test errors resulting from these five approaches?**

``` r
ggplot() +
  geom_col( aes(x = names(errors),
                y = errors)) +
  coord_flip()
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-36-1.png)

``` r
ggplot() +
  geom_density( aes(x = lasso_pred - y_test)) +
  geom_vline(xintercept = 0)
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-37-1.png)

As we can see, most of the prediction are within 1000 for most of the test dataset. The lasso has a substantial advantage over the the other methods, as we can see from the linear model approach.

### ex. 10

**We have seen that as the number of features used in a model increases, the training error will necessarily decrease, but the test error may not. We will now explore this in a simulated data set.**

**(a) Generate a data set with p = 20 features, n =1,000 observations, and an associated quantitative response vector generated according to the model Y = Xβ + ϵ, where β has some elements that are exactly equal to zero.**

``` r
set.seed(123)
p <- 20
n <- 1000
beta <- runif(p, min = -4, max = 4)
beta[5:12] <- 0

data <- data.frame( matrix( data = rnorm(p*n, sd = 4),
                            nrow = 1000,
                            ncol = p) )

e <- rnorm(1000, sd = 2)
y <- rep(0, times = n)
for( i in 1:n) y[i] <- sum(data[i,] * beta) + e[i]

data["y"] <- y
```

**(b) Split your data set into a training set containing 100 observations and a test set containing 900 observations.**

``` r
train <- sample(x = n, size = 100)
test <- (1:n)[-train]
```

**(c) Perform best subset selection on the training set, and plot the training set MSE associated with the best model of each size.**

``` r
regfit_full <- regsubsets(y ∼ ., 
                          data = data[train,],
                          nvmax = 20)
regfit_full_sum <- summary(regfit_full)

train_mat <- model.matrix(y ∼ ., data = data[train, ])

val_errors_train <- rep(NA, p)
for (i in 1:p) { 
  coefi <- coef(regfit_full , id = i)
  pred <- train_mat[, names(coefi)] %*% coefi
  val_errors_train[i] <- mean((data[train, "y"] - pred)^2) 
}

ggplot() +
  geom_line( aes( x = 1:p,
                  y = val_errors_train))
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-40-1.png)

**(d) Plot the test set MSE associated with the best model of each size.**

``` r
test_mat <- model.matrix(y ∼ ., data = data[test, ])

val_errors_test <- rep(NA, p)
for (i in 1:p) { 
  coefi <- coef(regfit_full , id = i)
  pred <- test_mat[, names(coefi)] %*% coefi
  val_errors_test[i] <- mean((data[test, "y"] - pred)^2) 
  }

ggplot() +
  geom_line( aes( x = 1:p,
                  y = val_errors_test)) +
  scale_y_log10()
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-41-1.png)

**(e) For which model size does the test set MSE take on its minimum value?**

``` r
which.min(val_errors_test)
```

    ## [1] 12

**Comment on your results. If it takes on its minimum value for a model containing only an intercept or a model containing all of the features, then play around with the way that you are generating the data in (a) until you come up with a scenario in which the test set MSE is minimized for an intermediate model size.** Since the 8 variable with a beta coefficient set to 0 are noise, the model with 12 variables was expected to be the one with the lowest test MSE, whichis slightly increased with &gt;12 variable models.

**(f) How does the model at which the test set MSE is minimized compare to the true model used to generate the data? Comment on the coefficient values.**

``` r
coef(regfit_full,12)
```

    ## (Intercept)          X1          X2          X3          X4         X13 
    ##  -0.2828308  -1.6689506   2.3304298  -0.8542322   3.0362714   1.4888301 
    ##         X14         X15         X16         X17         X18         X19 
    ##   0.5567702  -3.1288740   3.2763534  -1.9658116  -3.6292481  -1.3683916 
    ##         X20 
    ##   3.7802492

``` r
beta[ beta != 0]
```

    ##  [1] -1.6993798  2.3064411 -0.7281846  3.0641392  1.4205651  0.5810672
    ##  [7] -3.1766025  3.1985998 -2.0312981 -3.6635237 -1.3766342  3.6360292

``` r
ggplot() +
  geom_point( aes(x = names(coef(regfit_full,12)),
                  y = coef(regfit_full,12)),
              alpha = 0.8) +
  geom_point(aes(x = 2:13,
                 y = beta[ beta != 0]),
             colour = "blue") +
  coord_flip()
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-43-1.png)

As we can see, some of the beta coefficients are near each other (X14), while some others are really far away. For sure, all the coefficients that should be zero are zero even in the regfit model. Probably all this distance in the estimates is because of the very small training set (10%) which prevent the models from reaching a higher accuracy in prediction.

\*\*(g) Create a plot displaying of r,where \[formula on p.287\] for a range of values j is the jth coefficient estimate for the best model containing r coefficients. Comment on what you observe. How does this compare to the test MSE plot from (d)?

``` r
beta_pred <- data.frame( matrix(0, nrow = p, ncol = p))

for (i in 1:p) { 
  coefi <- coef(regfit_full , id = i)[-1]
  for(j in names(coefi)) {
    beta_pred[i, j] <- coefi[[j]]
  }}

beta_true <- data.frame( matrix( rep(beta, times = p), 
                                 nrow = p,
                                 byrow = TRUE))

tbse <- rep(NA, times = p) # total betas squadred error
for (i in 1:p) { 
  for( j in names(beta_pred)) {
    beta_pred[i, j] <- (beta_true[i, j] - beta_pred[i, j])^2
  }
  tbse[i] <- sum( beta_pred[i,] )
  }

### Do not use - different formula
### When bij does not exist, it leaves the difference to 0
###
# tbse <- rep(NA, times = p) # total betas squadred error
# for (i in 1:p) { 
#   for( j in names(beta_pred)) {
#     if( beta_pred[i, j] != 0 ) 
#       beta_pred[i, j] <- (beta_true[i, j] - beta_pred[i, j])^2
#   }
#   tbse[i] <- sum( beta_pred[i,] )
#   }
###

tbse_plot <- ggplot() +
  geom_line( aes( x = 1:p,
                  y = sqrt(tbse)))+
  geom_vline( xintercept = which.min(tbse),
              linetype = "dotted")
mse_plot <- ggplot() +
  geom_line( aes( x = 1:p,
                  y = val_errors_test),
             color = "blue") +
  geom_vline( xintercept = which.min(val_errors_test),
              linetype = "dotted") +
  scale_y_log10()
ggarrange(tbse_plot, mse_plot, nrow = 2)
```

![](chapter_6_files/figure-markdown_github/unnamed-chunk-44-1.png)

The curves show a similar behaviour, with a minimum in the same point. I expected the total beta squared difference to have a minimum point at the best model point, since it shows the difference between the real and predicted coefficients.

### ex. 11

\*\*