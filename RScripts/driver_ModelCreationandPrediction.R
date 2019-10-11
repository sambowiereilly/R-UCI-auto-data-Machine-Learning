# _______________________________________
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# func_ModelFunctions.R
# Custom set of functions to enable systematic candidate and best model selection
# for single x covariate  relationship with y response variable
#
# Functions include (but not limited to);  func_ModelCreation, func_charting_OptimumModelParameters 

# load the glue that holds this whole operation together!
source('func_ModelFunctions.R')

# load source library files
library('splines')
library('psych')
library('car')
library('nlme')
library('mgcv')
library('dplyr')
library('knitr' )
library('ggplot2')
library('gridExtra')

# set seed
set.seed( 123 )

# Declare & Set Global Variables
data = read.csv("auto-data.csv", sep=",", header=T)
n_row <- nrow(data) 

# response column
y_index <- 2
# attribute/covariate columns
x_indexes <- c(4,5,6,7)
#x_covariate <- displacement
#x_covariate <- horsepower    # only horsepower contain NAs (6 NAs in total out of 398) 
#x_covariate <- weight
#x_covariate <- acceleration

n_XCoveriates <- length(x_indexes)
# replace any NA's with mean of column
for(i in 1:n_XCoveriates) {
  data[is.na(data[,x_indexes[i]]), x_indexes[i]] <- round(mean(data[,x_indexes[i]], na.rm = TRUE), digits = 2)
}

# /* ~~~~~~
#  * Please note, within the paper I step through basic data exploration
#  * correlation matrices, box plots etc. 
#  * 
#  * For the code used, please see the bottom of this document
#  * Allows us to keep the cool stuff near the top - easy to hand.


# /* ~~~~~~
#  * Now we can leverage teh fruits of our labour (i.e. 'func_ModelFunctions.R')

#  * To perform full analysis - model tuning and selection for all attributes/covariates
#  * and model types (Polynomial, Bin Smooth & B-Spline) we can simply use the 7 lines of code below

#   Automatic Model & Parameter selection
#   **** for full details on systematic methodology used please see Section 3 of the paper (or check out the 'func_ModelFunctions.R' script!) 
n_row <- nrow(data) 
# mpg
y_index <- 2
# ( displacement, horsepower, weight, acceleration )
x_indexes <- c(4,5,6,7)
# for k-fold cross validation
k <- 5
# +/- tolerance range between train and test prediction accuracy
varianceoffsetallowance <- 0.05
# split data set into train and test using this % proportion
traintestsplit <- 0.7
# finally, run the main function 
# this function calls functions which call functions within the 'func_ModelFunctions.R' script
# from experience, for data set of ~398 rows entire process takes ~ 10 seconds or less
finalModelSelection <- func_OptimalModelSelection( data, x_indexes, y_index, traintestsplit, varianceoffsetallowance, n_row, k)
# view final table 
View(finalModelSelection)


# /* ~~~~~~
#  * In addition to tabular results, we can visualise elements of the journey/methodology too

#  * To view how the parameters/settings for each of the candidate models were selected
#  * simply use the code below - many of the variables mirror those already seen in above example
#   default settings, no need to change these
y_index <- 2
x_indexes <- c(4,5,6,7)

# mix and match these variables - how does increasing/reducing traintestsplit change results
x_index <- x_indexes[2]
chartType <- 'AdjR2'      # options; 'AdjR2', 'RSS'
modelType <- 'BinSmooth'  # options; 'Polynomial', 'BinSmooth', 'BSpline'
charting_showTest <- TRUE
traintestsplit <- 0.7
varianceoffsetallowance <- 0.05
# ^how does varying this offsetallowance (+/- variance range) alter results
# from experience, for data set of ~398 rows entire process takes less than 5 seconds
viewOptimumParameterSelection <- func_charting_OptimumModelParameters(data, x_index, y_index, 
                                                  traintestsplit, varianceoffsetallowance, 
                                                  chartType, modelType, charting_showTest = TRUE )


# /* ~~~~~~
#  * Final charting function allows us to view models over the base data
#  * for example, plot Bin Smooth or B-Spline model over the horsepower x values
# ***~~~~##### PLEASE NOTE FOR THIS EXAMPLE, PLEASE ENSURE YOU HAVE GENERATE 'finalModelSelection' variable above (around lines 75)
#   default settings, no need to change these
y_index <- 2
i_attributes <- c(1, 4, 7, 10) # starting positions of covariates/attributes within finalModelSelectionResults table

# switch this index from 1 through to 4 to see different models for different attributes/covariates
i_attribute <- i_attributes[4]

# Candidate Polynomial, Bin Smooth & B-Spine models are plotted next to each other
# looks pretty good if i'm honest.
par(mfrow=c(1,3))
for(i in i_attribute:(i_attribute+2)){
  modelType <- finalModelSelectionResults[i,5]
  tp_main <- finalModelSelectionResults[i,6]
  tp_Secondary <- finalModelSelectionResults[i,7]
  baseDataSet <- data.frame(y = data[, y_index] , x = data[,finalModelSelectionResults[i,4]])
  baseDataSet[,2] <- scale(baseDataSet[,2])
  train <- func_Sort2DTable(baseDataSet)
  x_ColName <- finalModelSelectionResults[i,3]
  # now call custom charting function
  func_charting_CandidateModels( modelType, tp_main, tp_Secondary, train, x_ColName)
}
par(mfrow=c(1,1))





# /* ~~~~~~
#  * basic data exploration
#  * correlation matrices, box plots etc. 
#  * 
# Lets start by making some pairwise correlations
base_correlationsdata <- data[,c(2,4,5,6,7)]
correlations <- cor( base_correlationsdata )
correlations[ lower.tri( correlations, diag = TRUE ) ] <- NA
correlations <- as.data.frame( as.table( correlations ) )
correlations[,3] <- round(correlations[,3], digits = 2)
# now we are only interested where response mpg is included
yResponse <- 'mpg'
correlations <- correlations %>%
  na.omit( )  %>%
  filter( Var1 == yResponse | Var2 == yResponse )
# final output for report
kable(correlations)

# mean figures for each attribute/covariate
#mean(data$mpg);mean(data$displacement);mean(data$weight);mean(data$acceleration)
round(mean(data$horsepower), digits = 2)

# create table of Column Names & Standard Deviation Results
column_names <- c( 'mpg', 'displacement', 'horsepower', 'weight', 'acceleration')
sd_results <- c(sd(data$mpg), sd(data$displacement), sd(data$horsepower)
                , sd(data$weight), sd(data$acceleration))
sd_results <- round(sd_results, digits = 2)
# final output for report
kable(data.frame( 'ColumnNames' = column_names, 'StandardDeviation' = sd_results))

# Create box plot of columns
distribution_MPG <- ggplot( data = data, aes(y = mpg)) +
  geom_boxplot() +
  ggtitle('Box Plot Distribution:', subtitle = 'MPG') +
  ylab('') +
  xlab('mpg') +
  theme(  axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

distribution_Displacement <- ggplot( data = data, aes(y = displacement)) +
  geom_boxplot() +
  ggtitle('', subtitle = 'Displacement') +
  ylab('') +
  xlab('displacement') +
  theme(  axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

distribution_Horsepower <- ggplot( data = data, aes(y = horsepower)) +
  geom_boxplot() +
  ggtitle('', subtitle = 'Horsepower') +
  ylab('') +
  xlab('horsepower') +
  theme(  axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

distribution_Weight <- ggplot( data = data, aes(y = weight)) +
  geom_boxplot() +
  ggtitle('', subtitle = 'Weight') +
  ylab('') +
  xlab('weight') +
  theme(  axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

distribution_Acceleration <- ggplot( data = data, aes(y = acceleration)) +
  geom_boxplot() +
  ggtitle('', subtitle = 'Acceleration') +
  ylab('') +
  xlab('acceleration') +
  theme(  axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

# Create XY SCatter plot to show pairwise correlation direction
correlation_Displacement <- ggplot( data = data, aes( x = displacement, y = mpg ) ) +
  geom_point( size = 2 ) +
  ggtitle( 'Pair wise relationship:', subtitle = 'Displacement and MPG' ) +
  xlab( 'Displacement' ) +
  ylab( 'MPG' ) +
  geom_smooth(  method = 'lm', formula = y ~ x, se = FALSE) +
  theme_light()

correlation_Horsepower <- ggplot( data = data, aes( x = horsepower, y = mpg ) ) +
  geom_point( size = 2 ) +
  ggtitle( '', subtitle = 'Horsepower and MPG') +
  xlab( 'Horsepower' ) +
  ylab( 'MPG' ) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE) +
  theme_light()

correlation_Weight <- ggplot( data = data, aes( x = weight, y = mpg ) ) +
  geom_point( size = 2 ) +
  ggtitle( '', subtitle = 'Weight and MPG' ) +
  xlab( 'Weight' ) +
  ylab( 'MPG' ) +
  geom_smooth(  method = 'lm', formula = y ~ x, se = FALSE) +
  theme_light()

correlation_Acceleration <- ggplot( data = data, aes( x = acceleration, y = mpg ) ) +
  geom_point( size = 2 ) +
  ggtitle( '', subtitle = 'Acceleration and MPG' ) +
  xlab( 'Acceleration' ) +
  ylab( 'MPG' ) +
  geom_smooth(  method = 'lm', formula = y ~ x, se = FALSE) +
  theme_light()

# arrange in grid for report - Box Plot Distirbutions
grid.arrange( distribution_MPG, distribution_Displacement, 
              distribution_Horsepower, distribution_Weight, 
              distribution_Acceleration, nrow = 1 ) 

# arrange in grid for report - Correlations
grid.arrange( correlation_Displacement , correlation_Horsepower, 
              correlation_Weight, correlation_Acceleration, nrow = 2 ) 



modelTypeList <- c('Polynomial', 'BinSmooth', 'BSpline')
modelTuningParameter1 <- c('nth order', 'bin size', 'quantile knots')
modelExampleValues1 <- c( '1, 2, 3, ... , 13, 14, 15', '20, 25, 30, ..., 240, 245, 250', '1, 2, 3, ... , 13, 14, 15')
modelTuningParameter2 <- c('', '', 'level of degree')
modelExampleValues2 <- c( '', '', '1, 2, 3, ..., 8, 9, 10')

table_TunningParameters <- data.frame( modelTypeList, modelTuningParameter1, modelExampleValues1,
                                       modelTuningParameter2, modelExampleValues2)



