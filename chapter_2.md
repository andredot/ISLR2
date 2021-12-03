ISLR Ch.2
================

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
execute code within the notebook, the results appear beneath the code.

### ex. 1

**Describe the null hypotheses to which the p-values given in Table 3.4
correspond. Explain what conclusions you can draw based on these
p-values. Your explanation should be phrased in terms of sales, TV,
radio, and newspaper, rather than in terms of the coefficients of the
linear model.**

-   p-value of “intercept” is the probability that if there were no
    relationships in Advertising data, the resulting intercept which
    represent sales (obtained by fitting a linear model with TV, radio
    and newspaper predictors) would be different from 0, when no money
    is spent on predictors
-   p-value of TV/Radio/Newspaper represent the probability that if
    there was no relationship between sales and money spent on
    TV/Radio/Newspaper, the change in sales for each euro spent on
    TV/Radio/Newspaper advertisement would be different from 0
-   in general the null hypothesis is that money spent on TV or Radio or
    Newspapers are not correlated to the average sales of the product
    measured in the database (thus every pattern we see in the dataset
    is random)
-   I’m confident that a LR model represent the reality within an
    acceptable level of accuracy when TV and Radio are used to model
    sales. At the same time, I’m not sure that adding newspaper will
    increase that accuracy, possibly only adding complexity without
    improving results.

### ex. 2

**Carefully explain the differences between the KNN classifier and KNN
regression methods**

-   KNN-c and KNN-r both select the K observation using the same
    definition of “nearest”, thus selecting for a x1,x2…xK point in
    space the K observations that minimize the squared distance from it
-   KNN-c then classifies the observation as the class most represented
    in that set of observations, irrespective of distance from the point
-   KNN-r computes the average of the response on all the K points
    selected, first summing all the y and then dividing them by k
-   Since the average can only make sense in a continuous or discrete
    space (there is no certainty that the predicted result will be one
    of the classification classes), it can not be generalized for cases
    with categorical variables than can not be represented on a x-axis

### ex. 3

**Suppose we have a data set with five predictors, X1 = GPA, X2 = IQ, X3
= Level (1 for College and 0 for High School), X4 = Interaction between
GPA and IQ, and X5 = Interaction between GPA and Level. The response is
starting salary after graduation (in thousands of dollars). Suppose we
use least squares to fit the model, and get βˆ0 = 50, ˆβ1 = 20, ˆβ2
=0.07, ˆβ3 = 35, ˆβ4 =0.01, ˆβ5 = −10**

1.  Which answer is correct, and why?

-   1.  For a fixed value of IQ and GPA, high school graduates earn
        more, on average, than college graduates. **no**

-   2.  For a fixed value of IQ and GPA, college graduates earn more, on
        average, than high school graduates. **no**

-   3.  For a fixed value of IQ and GPA, high school graduates earn
        more, on average, than college graduates provided that the GPA
        is high enough. **yes**

-   4.  For a fixed value of IQ and GPA, college graduates earn more, on
        average, than high school graduates provided that the GPA is
        high enough. **no**

**The effect of beta4 can be ignored, since both GPA and IQ are kept
fixed, and so is the value of their interaction term.**

**The effect of beta5 is 0 when x3 is 0, so the equation for the salary
can be reduced to**

``` r
# BETA = combined effect of Beta0 + Beta1*x1 + Beta2*x2 + beta4*x4
salary_high_school = BETA + beta3*x3 + beta5*x1*x3 
                   = BETA # since x3 = 0
salary_college = BETA + beta3*x3 + beta5*x1*x3 
               = BETA + beta3 + beta5*x1 # since x3 = 1
```

**So, for small values of x1 (GPA), salary_college \> salary_high_school
. But since beta5 is negative, as x1 increases, salary_college will get
progressively smaller**

**(b) Predict the salary of a college graduate with IQ of 110 and a GPA
of 4.0.**

``` r
(salary_student_1 = 50 + 20*4.0 + 0.07*110 + 35*1 + 0.01*4.0*110 - 10*4.0*1)
```

    ## [1] 137.1

3.  True or false: Since the coefficient for the GPA/IQ interaction term
    is very small, there is very little evidence of an interaction
    effect. Justify your answer.

**Short premise: IQ is a random variable with mean 100 and st.dev = 10,
so it is not a good idea to use it in a linear model without
transforming it: the difference between 80 and 90 is not the same
between 90 and 100. Not sure if this also is true for GPA**

**TLDR: False. Full answer: the beta4 coefficient inform us about how
strong is the association between salary and the combination of IQ and
GPA, which is in fact weak (but please note that with GPA>3.5 at least
1/3 of the effect of IQ on salary is given by the interaction term).
But, nothing can be said on the evidence of this interaction, that can
not be appreciated without prediction intervals (which combine both the
confidence interval of our prediction and our uncertainty on the
estimate of the true beta4 value)**

### ex. 4

\*\*I collect a set of data (n = 100 observations) containing a single
predictor and a quantitative response. I then fit a linear regression
model to the data, as well as a separate cubic regression, i.e. Y = β0 +
β1*X + β2*X^2 + β3\*X^3 + ϵ \*\*

1.  Suppose that the true relationship between X and Y is linear, i.e. Y
    = β0 + β1X + ϵ. Consider the training residual sum of squares (RSS)
    for the linear regression, and also the training RSS for the cubic
    regression. Would we expect one to be lower than the other, would we
    expect them to be the same, or is there not enough information to
    tell? Justify your answer

**Adding quadratic and cubic terms introduce a more flexible model
(lower bias, higher variance), which can provide a better fit to the
training data. It depends on how much the error term random fluctuation
is modelled by the additional terms to the Linear Regression (as the
error term approaches zero, the training error of a cubic model will be
more than a linear one, when the true relationship is linear)**

2.  Answer (a) using test rather than training RSS.

**Test RSS will be lower for a linear model than the one of a cubic
model, when the true relationship is linear**

3.  Suppose that the true relationship between X and Y is not linear,
    but we don’t know how far it is from linear. Consider the training
    RSS for the linear regression, and also the training RSS for the
    cubic regression. Would we expect one to be lower than the other,
    would we expect them to be the same, or is there not enough
    information to tell? Justify your answer.

**The more a true relationship is far from a linear one, the more a
cubic model can reduce its bias with only a slight increase of variance.
In training data, a cubic model outperforms the linear one…**

4.  Answer (c) using test rather than training RSS.

**…but in test error it depends on the real case scenario: how much is
the cubic one a better fit to the real function than the linear one? And
how much random noise have been incorporated in the coefficients (due to
high variability of a more complex model)?**

### ex. 5

Help

### ex. 6

Using (3.4), argue that in the case of simple linear regression, the
least squares line always passes through the point (¯x, ¯y).

![Sketch of exercise](img/ch3_ex6_1.jpg)