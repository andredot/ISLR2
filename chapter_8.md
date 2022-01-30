ISLR Ch.8
================

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code.

### ex. 1

Made on paper

### ex. 2

**It is mentioned in Section 8.2.3 that boosting using depth-one trees (or stumps) leads to an additive model: that is, a model of the form **

f(X)= sum\_j( fj(Xj) )

**Explain why this is the case**

Since in boosting, given the current model, we fit a decision tree to the residuals from the model, what happens is that at each iteration of the algorithm we add a binary split (thus adding a 2-step non-continuous function). A model built on addition of step functions is an additive model. So, the model if fit to the residuals, which are calculated form a sum of step functions.

### ex. 3

**Consider the Gini index, classification error, and entropy in a simple classification setting with two classes. Create a single plot that displays each of these quantities as a function of ˆpm1. The x-axis should display ˆpm1, ranging from 0 to 1, and the y-axis should display the value of the Gini index, classification error, and entropy.**

Remember that pmk represents the proportion of training observations in the mth region that are from the kth class. Since we are in 2-class setting, pm2 = 1- pm1.

``` r
pm1 <- seq(from = 0, to = 1, length.out = 1e3)
data <- data.frame(pm1) %>% 
  mutate( pm2 = 1 - pm1,
          cl_err = as.double( map( pm1, ∼ 1 - max(.x, 1 - .x))),
          gini_err = 2 * pm1 * pm2,
          entropy = -1 * ( pm1 * log2(pm1) + pm2 * log2(pm2)) )

ggplot(data,
       aes( x = pm1)) +
  geom_line( aes( y = cl_err),
             color = "blue") +
  geom_line( aes( y = gini_err),
             color = "red") +
  geom_line( aes( y = entropy),
             color = "green")
```

    ## Warning: Removed 2 row(s) containing missing values (geom_path).

![](chapter_8_files/figure-markdown_github/unnamed-chunk-2-1.png)

### ex. 4

made on paper

### ex. 5

**Suppose we produce ten bootstrapped samples from a data set containing red and green classes. We then apply a classification tree to each bootstrapped sample and, for a specific value of X, produce 10 estimates of**

P(Class is Red|X): 0.1, 0.15, 0.2, 0.2, 0.55, 0.6, 0.6, 0.65, 0.7, and 0.75

``` r
pred <- c(0.1, 0.15, 0.2, 0.2, 0.55, 0.6, 0.6, 0.65, 0.7, 0.75)
```

Let's assume decision boundary is at 0.5.

**There are two common ways to combine these results together into a single class prediction. One is the majority vote approach discussed in this chapter. The second approach is to classify based on the average probability. In this example, what is the final classification under each of these two approaches?**

``` r
# majority
( sum( pred > 0.5) )
```

    ## [1] 6

since most of the predictions (6 out of 10) are for red class.

``` r
( mean( pred))
```

    ## [1] 0.45

Since mean of predictions is under the assumed threashold, predicted for average probability is green.

### ex. 6

**Provide a detailed explanation of the algorithm that is used to fit a regression tree.**

Each time we want to split the predictors space (assuming we alredy know the current RSS) we perform the following steps:

-   for each value of each predictor we calculate the resulting RSS of a split made in that point. Formula is standard, to calculate the amount of RSS we imagine to assigas the value for the new half-plane the mean of the predictor points in that region.

-   the next split correspond to the value of the predictor that provides the greatest reduction in the RSS, the new RSS is calculated and we can check if we have reached a stop condition (e.g. each node wih less thant 5 points, max number of leaves,...)

### ex. 7

**In the lab, we applied random forests to the Boston data using mtry = 6 and using ntree = 25 and ntree = 500. Create a plot displaying the test error resulting from random forests on this data set for a more comprehensive range of values for mtry and ntree. You can model your plot after Figure 8.10. Describe the results obtained.**

``` r
boston <- Boston
set.seed (1) 

mtry = 2*(1:6)
ntree = c(10*1:10, 200*1:5)
mse = data.frame( mtry = as.numeric(), 
                  ntree = as.numeric(), 
                  value = as.numeric())

train <- sample (1:nrow(boston), nrow(boston) / 2)
boston_test <- boston[-train , "medv"]

for( i in 1:length(mtry) ) {
  for( j in 1:length(ntree) ){
    rf_boston <- randomForest(medv ∼ ., 
                          data = boston , 
                          subset = train , 
                          mtry = mtry[i],
                          ntree = ntree[j])

     yhat_rf <- predict(rf_boston , 
                        newdata = boston[-train , ])

     mse <- mse %>% 
       add_row( mtry = mtry[i],
                ntree = ntree[j],
                value = mean((yhat_rf - boston_test)^2) )
  }
}

ggplot(mse) +
  geom_contour_filled( aes(x = mtry,
                    y = ntree,
                    z = value)) +
  scale_y_log10()
```

![](chapter_8_files/figure-markdown_github/unnamed-chunk-6-1.png)

As we can see, the MSE ranges between 17 and 27 in our parameter grid, with a test MSE which has a minimum around mtry = 6 and tree = 20.

### ex. 8