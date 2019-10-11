# _______________________________________
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# func_ModelFunctions.R
# Custom set of functions to enable systematic candidate and best model selection
# for single x covariate  relationship with y response variable
#
# Functions include (but not limited to);  func_ModelCreation, func_charting_OptimumModelParameters 

# load library - to enable bs() function
library(splines)

# /* ~~~~~~
#  *  Custom Model Creation Function
#  *  Goals ~ Generalise Model Creation elements 
#  *  (1 stop shop, reducing code duplication)
#  *
#  *  Models currently supported:   Polynomial Regression(n=1 ~ Linear), Bin Smooth & B-Spline
func_ModelCreation <- function( param_ModelType, 
                                param_y, param_x, 
                                param_pPoly = 1,
                                param_bs_binLength = NULL, param_bs_Knots = NULL,
                                param_bsp_Knots = c(0.5), param_bsp_Degree = 1) {
  # initialise default values
  y_response <- param_y             # pre scaled y and x variables
  x_covariate <- param_x
  n_Rows <- length(x_covariate)     # rows in data set
  p_Poly <- param_pPoly
  m_matrix <- x_covariate           # Model Design Matrix
  q <- 1
  pred_BinSmooth_Table <- NULL
  chart_BinSmooth_Start <- NULL
  chart_BinSmooth_End <- NULL
  chart_BinSmooth_BinLength <- NULL
  predict_newDat <- NULL
  
  # specific logic based on Model Type parsed into function  
  if(param_ModelType == 'BinSmooth'){
    # calculate bespoke Bin Smooth design matrix to enable model creation
    model_BinSmooth <- func_BinSmoothBespoke( x_covariate, param_bs_binLength, param_bs_Knots, n_Rows)
    # update Model Design Matrix
    m_matrix <- model_BinSmooth$bs_matrix 
    chart_BinSmooth_Start <- model_BinSmooth$chart_BinSmooth_Start
    chart_BinSmooth_End <- model_BinSmooth$chart_BinSmooth_End
    chart_BinSmooth_BinLength <- model_BinSmooth$bs_BinLength
    q <- model_BinSmooth$bs_Knots
    # fit BinSmooth specifc model
    base_Model <- lm(y_response ~ 0 + m_matrix)
    # bespoke table to enable custom prediction function (## insert name here!! ##)
    pred_BinSmooth_Table <- data.frame( Coefs = coef(base_Model),
                                        LowerXBounds = model_BinSmooth$bs_PredBasis[,1],
                                        UpperXBounds = model_BinSmooth$bs_PredBasis[,2])
    
    # ***noticed issue with yhat values obtained from lm() above, so we are over-writing the 
    # RSS and Adj R^2 values. Issue links with, lm() can pedict 1 or more values of y for a given value of x..
    # which we cannot unfortunately replicate. Therefore, utilise customer predict function to ensure
    # RSS & Adj R^2 values are like for like when cross fold validation is complete.
    yhat_BinSmooth <- func_ModelPrediction( param_ModelType, y_response , x_covariate, 
                                            n_coeffs = q, NULL, pred_BinSmooth_Table, q)
    
    model_RSS <- yhat_BinSmooth$pred_RSS
    model_AdjRSquared <- yhat_BinSmooth$pred_AdjRSquared
    
  } else {
    # otherwise we can generalise remaining models like so:
    if(param_ModelType == 'BSpline'){
      q <- length(param_bsp_Knots) + param_bsp_Degree + 1
      base_Model <- lm(y_response ~ bs(m_matrix, knots = param_bsp_Knots, degree = param_bsp_Degree))
    } else {
      # for Linear & Polynomial Regression
      q <- p_Poly + 1
      # linear model equiv to nth order polynomial equal to 1
      base_Model <- lm(y_response ~ poly(m_matrix, degree = p_Poly, raw = TRUE))
    }
    #Residual sum of squares
    model_RSS <- deviance(base_Model)
    
    #Adjusted Coefficient of determination: R^2
    #model_AdjRSquared <- summary(base_Model)$adj.r.squared
    model_AdjRSquared <- func_AdjustedRSquared( y_response, model_RSS, n_Rows, q)
    
    #model_MSE <-  mean(sapply( residuals(base_Model) , function(x) { x^2 }))
    #model_AIC <- AIC(base_Model)
  }
  
  # return desired results outlined in assignment file  
  return( list( modelType = param_ModelType,
                modelOutput = base_Model, 
                predhelp_q = q,
                RSS = model_RSS, 
                AdjRSquared = model_AdjRSquared, 
                pred_BinSmooth_Table = pred_BinSmooth_Table,
                chart_BinSmooth_Start = chart_BinSmooth_Start,
                chart_BinSmooth_End = chart_BinSmooth_End,
                chart_BinSmooth_BinLength = chart_BinSmooth_BinLength,
                bsp_Knots = param_bsp_Knots,
                bsp_Degree = param_bsp_Degree) )
}

# /* ~~~~~~
#  * Custom Bin Smooth Function
#  * Goals ~ generate design matrix required for Bin Smooth algorithm
func_BinSmoothBespoke <- function( param_x, 
                                   param_bs_binLength = NULL, 
                                   param_bs_Knots = NULL, n_Rows ) {
  # max knots permitted within model
  # for predictions we require a one to one ratio between x and y 
  max_Knots <- ceiling(length(unique(param_x))*0.95)
  #max_Knots <- length(unique(param_x))
  
  # check for valid inputs ~ if not, simply set knots = 2 and proceed
  if(!is.numeric(param_bs_binLength) && !is.numeric(param_bs_Knots)) {
    bs_Knots <- 2
    bs_BinLength <- ceiling(n_Rows / bs_Knots)
    #print('NO variables assigned')
  } else if (!is.numeric(param_bs_binLength)) {
    bs_Knots <- ifelse(param_bs_Knots < 1, 1, ifelse(param_bs_Knots > max_Knots, max_Knots, param_bs_Knots))
    bs_BinLength <- ceiling(n_Rows / bs_Knots)
    #print('Bin Length NOT Assigned')
  } else {
    bs_BinLength <- ifelse(param_bs_binLength < 1, 1, param_bs_binLength)
    bs_Knots <- ceiling(n_Rows / bs_BinLength)
    # lets check number of Knots does not exceed max possible
    if(bs_Knots > max_Knots){
      #print('hi')
      bs_BinLength <- ceiling(n_Rows / max_Knots)   # if so, rebase bin length accordingly.
      bs_Knots <- ceiling(n_Rows/bs_BinLength)
    }     
    #print('Knots NOT Assigned')
  }
  # initialise Bin Smooth Design Matrix
  bs_matrix <- matrix(0, n_Rows, bs_Knots)
  bs_PredBasis <- matrix(0, bs_Knots, 2)
  
  # geneate design matrix for Bin Smooth Regression  
  for(i in 1:bs_Knots) {
    #bs_start_i <- (bs_BinLength * ifelse(i==1, 0, i-1)) + 1
    bs_start_i <-  ifelse(i==1, 1, (i - 1) * bs_BinLength + 1)
    
    #bs_end_i <- min((bs_BinLength * ifelse(i==1, 1, i)), n_Rows) # min() ensures we do not go out of bounds
    bs_end_i <- min((bs_start_i + bs_BinLength - 1), n_Rows) # min() ensures we do not go out of bounds
    bs_matrix[bs_start_i:bs_end_i , i] <- 1    # change to 1 only when required
    
    bs_PredBasis[i, 1] <- ifelse( i == 1, -Inf, bs_PredBasis[(i-1), 2]) 
    bs_PredBasis[i, 2] <- ifelse( i == bs_Knots, Inf, param_x[bs_end_i])
  }
  # return Bin Smooth Model elements  
  return( list( bs_matrix = bs_matrix, 
                bs_BinLength = bs_BinLength, 
                bs_Knots = bs_Knots, 
                bs_PredBasis = bs_PredBasis,
                chart_BinSmooth_Start = bs_start_i,    # these will be useful for charting later
                chart_BinSmooth_End = bs_end_i) )
}

# /* ~~~~~~
#  *  Custom Model Prediction Function
#  *  Goals ~ Generalise Model Prediction elements 
#  *  (1 stop shop, reducing code duplication)
#  *
#  *  Models currently supported:   Linear, Polynomial & Bin Smooth Regression
func_ModelPrediction <- function( modelType, y_response, new_x, n_coeffs, 
                                  modelOutput = NULL, pred_BinSmooth_Table = NULL, q) {
  # initialise function values
  n_x <- length(new_x)
  yhat <- matrix(0, n_x)
  # unfortunately predict() does not work correctl with BinSmooth model, so quick loop below calculates for us
  # (little bit of a pain but manageable)
  if(modelType == 'BinSmooth') {
    # utilise bespoke Bin Smooth supporting prediction table to enable Bin Smooth predictions
    for ( i in 1:n_x){
      for( j in 1:n_coeffs){
        if( new_x[i] > as.numeric(pred_BinSmooth_Table[j,2]) && new_x[i] <= as.numeric(pred_BinSmooth_Table[j,3])){
          #if( new_x[i] >= as.numeric(pred_BinSmooth_Table[j,2]) && new_x[i] < as.numeric(pred_BinSmooth_Table[j,3])){
          yhat[i] <- as.numeric(pred_BinSmooth_Table[j,1])
          break } } }
  } else {
    # standard predict function for other models
    yhat <- predict(modelOutput, newdata = data.frame(m_matrix = new_x))
  }
  new_RSS <-  sum(sapply((y_response - yhat) , function(x) { x^2 }))
  new_MSE <-  mean(sapply((y_response - yhat) , function(x) { x^2 }))
  
  # so here is the test, what is the adj r^2 for the unseen data
  new_AdjRSquared <- func_AdjustedRSquared( y_response, new_RSS, n_x, q)
  
  # return Model Prediction elements
  return( list( pred_y = yhat, 
                pred_AdjRSquared = new_AdjRSquared, 
                pred_RSS = new_RSS,
                pred_MSE = new_MSE ))
}


# /* ~~~~~~
#  *  Custom Adjusted R Squared function
#  *  Goals ~ Quickly calculate Adjusted R Squared using output from customer Model Creation/Prediction functions above  
#  *
func_AdjustedRSquared <- function( y_Response, model_RSS, n_Rows, q ) {
  #Coefficient of determination: R^2
  #R2 <- 1 - (model_RSS / (t(y_Response)%*%y_Response-(mean(y_Response)**2*n_Rows)))
  R2 <- 1 - (model_RSS / (t(y_Response)%*%y_Response-(mean(y_Response)^2*n_Rows)))  
  #Adjusted Coefficient of determination: R^2
  adjR2 <- 1 - ( (n_Rows-1)/(n_Rows-q) ) * (1-R2)
  # if RSS approaches the total sum of squares it is possible to have negative adjR2
  # this is of no value to is, so overwrite with 0 this is where to occur
  if( adjR2 < 0 ) { adjR2 <- 0  } 
  # return Bin Smooth Model elements  
  return( as.numeric(adjR2) )
}


# /* ~~~~~~
#  *  Custom Model to create candidate models
#  *  Goals ~ Leverage default array of model tuning parameters created in 'func_parameterTuningOptions'
#  *  This function heavily relies on 'func_ModelParamTurning' to achieve generalisable goals - for three models supported.
#  *  (again 1 stop shop, reducing code duplication)
#  *
#  *  Models currently supported:   Polynomial Regression(n=1 ~ Linear), Bin Smooth & B-Spline
func_candidateModels <- function( data, x_indexes, y_index, traintestsplit = 1, 
                                  varianceoffsetallowance, n_modelssupported = 3,
                                  fixedseed = TRUE){
  # allow user to set seed for reproducibility
  if( fixedseed ){set.seed( 123 )}
  samp <- sample(nrow(data), traintestsplit * nrow(data))
  # temporary area to dump results
  temp_Tuning <- list()
  for( i in 1:length(x_indexes)) {
    x_loop <- x_indexes[i]
    x_ColName <- colnames(data)[x_loop]
    # construct x and y data set
    baseDataSet <- data.frame(y = data[, y_index] , x = data[,x_loop])
    baseDataSet[,2] <- scale(baseDataSet[,2]) # Scale covariate
    #mean(baseDataSet[,2]); sd(baseDataSet[,2]) # x_covariateScaled characteristics
    
    # if train test split is set to 1, just use data set as is
    if(traintestsplit == 1){ 
      train <- baseDataSet
      test <- train
    } else {
      # otherwise take train and test samples accordingly
      train <- baseDataSet[samp,]
      test <- baseDataSet[-samp,] 
    }
    # now we can sort the data
    train <- func_Sort2DTable(train)
    test <- func_Sort2DTable(test)
    # update parameterTuning options
    parameterTuning <- func_parameterTuningOptions(nrow(train))
    # temp variables to stoe results
    temp_Results <- list()
    temp_ChartingInfo <- list()
    for (j in 1:n_modelssupported) {
      spt <- func_ModelParamTurning( parameterTuning[[j]], varianceoffsetallowance, train, test )
      temp_Results[[j]] <- spt$return_table  
      temp_ChartingInfo[[j]] <- spt
    }
    table_Results <- bind_rows( temp_Results )
    # append X Covariate details
    table_Results <- cbind( x_loop, table_Results )
    table_Results <- cbind( x_ColName, table_Results )
    # append final result to list
    temp_Tuning[[i]] <- table_Results
  }
  # unlist final results
  table_TuningResults <- bind_rows( temp_Tuning )
  # if only one x index passed, allow function to return chart info
  if( length(x_indexes) > 1) { temp_ChartingInfo <- 'Not Available' }
  # return results to user
  return( list( table_TuningResults = table_TuningResults, 
                charting_CandidateModelInfo = temp_ChartingInfo ))
}

# /* ~~~~~~
#  *  Custom K-Fold Sub Set creation function
#  *  Goals ~ quickly creates k sub sets of data to enable cross fold validation 
#  *  Crucial to enable results to be generalised for 'unseen' data
#  *
func_KFoldSubSets <- function( k, n_row, data, fixedseed = TRUE){
  kFold_DataSetCreation <- data
  kFold_MinKFoldSize <- floor(n_row/(k+1))  
  kFold_Indexes <- list()
  for(i in 1:k) {
    # allow user to set seed for reproducibility
    if( fixedseed ){set.seed( 123 )}
    if(i == k) {
      kFold_Index <- sample(nrow(kFold_DataSetCreation), n_row - ( kFold_MinKFoldSize * k ))
    } else {
      kFold_Index <- sample(nrow(kFold_DataSetCreation), kFold_MinKFoldSize)
    }
    # update remaining data positions
    kFold_DataSetCreation <- kFold_DataSetCreation[-kFold_Index, ]
    
    kFold_Indexes[[i]] <- kFold_Index
  }
  return(kFold_Indexes)
}


# /* ~~~~~~
#  *  Custom Function to sort a 2 column table  
#  *  (ideally table containing y response [,1] and x covariate [,2]
#  *  Sorted data is a crucial requirement for creation of Bin Smooth models.
#  *
func_Sort2DTable <- function( datatable ){
  datatable[,1] <- datatable[,1][order(datatable[,2])]
  datatable[,2] <- sort(datatable[,2])
  return(datatable)
}


# /* ~~~~~~
#  *  Custom Function sets the default list of Model Tuning Parameter options  
#  *  Each model has a set list potential parameters use to create a set of models
#  *  which can later be evaluated (func_OptimumModelTuningValue) 
#  *
func_parameterTuningOptions <- function( train_nrow ){
  # set min and max bin smooth range equal to a proportion of train data set values
  base_pT_BinSmooth <- seq(floor(train_nrow*0.05), ceiling(train_nrow*0.95), 5)
  
  parameterTuning_Poly_nth <- seq(1, 15, 1)
  parameterTuning_BinSmooth_BinSize <- sort(base_pT_BinSmooth, decreasing = TRUE)
  parameterTuning_BSpline_Knots <- seq(1, 15, 1)
  parameterTuning_BSpline_Degree <- seq(1, 10, 1)
  
  parameterTuning <- list(list( modelType = 'Polynomial',
                                parameterTuning = parameterTuning_Poly_nth ),
                          
                          list( modelType = 'BinSmooth',
                                parameterTuning = parameterTuning_BinSmooth_BinSize ),
                          
                          list( modelType = 'BSpline',
                                parameterTuning = list( bsp_Knots = parameterTuning_BSpline_Knots,
                                                        bsp_Degree = parameterTuning_BSpline_Degree) ) )
  return( parameterTuning )
}


# /* ~~~~~~
#  *  Custom Function to create the models using the set of tuning parameters contained   
#  *  within func_parameterTuningOptions function.
#  *  which can later be evaluated (func_OptimumModelTuningValue) 
#  *
func_ModelParamTurning <- function( base_parameterTuning, varianceoffsetallowance, train, test) {
  # set default variables
  base_model <- NULL
  base_modelType <- unlist(base_parameterTuning[1])
  number_MatrixTuningColumns <- 8
  n_main <- 0; n_secondary <- 0; additional_parameterTuning <- 0
  if(base_modelType == 'BSpline') { 
    core_parameterTuning <- unlist(base_parameterTuning[[2]][[1]])
    additional_parameterTuning <- unlist(base_parameterTuning[[2]][[2]],0)  
  } else { core_parameterTuning <- unlist(base_parameterTuning[2]) }
  n_tuningLoops <- length(core_parameterTuning)
  n_tuningDegree <- length(additional_parameterTuning)
  n_maxiterations <- n_tuningLoops * n_tuningDegree
  matrix_tuning <- matrix(0, n_maxiterations, number_MatrixTuningColumns )  
  
  # generate core array of candidate model parameters
  for( n_trace in 1:n_maxiterations){
    if ( base_modelType == 'BSpline' ){ n_secondary <- ifelse( n_secondary == n_tuningDegree, 1, n_secondary + 1 ) 
    } else { n_secondary <- 1 }
    n_main <- ifelse( n_secondary == 1, n_main + 1, n_main ) 
    tp_main <- core_parameterTuning[n_main]
    tp_Secondary <- additional_parameterTuning[n_secondary]
    # create model and predictions
    modelcp <- func_ModelCreationandPrediction( base_modelType, tp_main, tp_Secondary, 
                                                train, test )
    # # store results of interest
    matrix_tuning[n_trace, 1] <- n_trace
    matrix_tuning[n_trace, 2] <- tp_main
    matrix_tuning[n_trace, 3] <- tp_Secondary
    matrix_tuning[n_trace, 4] <- round(modelcp$base_model$AdjRSquared, digits = 2)
    matrix_tuning[n_trace, 5] <- round(modelcp$predictions_model$pred_AdjRSquared, digits = 2)
    # also store Residual Sum of Squares for model train and test data too
    matrix_tuning[n_trace, 7] <- as.numeric(modelcp$base_model$RSS)
    matrix_tuning[n_trace, 8] <- as.numeric(modelcp$predictions_model$pred_RSS)
  }
  # calculate variance between train and test adjusted R squared results
  matrix_tuning[, 6] <- round(abs(matrix_tuning[, 5] - matrix_tuning[, 4]), digits = 2)
  # calculate optimum value based on methodology definition
  tuning_OptimumValue <- func_OptimumModelTuningValue( matrix_tuning, varianceoffsetallowance, base_modelType , number_MatrixTuningColumns)
  # return results to user
  return(  list( return_table = tuning_OptimumValue$return_table, 
                 tuning_optimumParam_index = tuning_OptimumValue$tuning_optimumParam_index,
                 tuning_matrix = matrix_tuning ) )
}


# /* ~~~~~~
#  *  This function combines previous Custom Functions above, mainly func_ModelCreation & func_ModelPrediction
#  *  In practice, methodology used always creates a model and immediately predicts ~ to reduce code duplication
#  *  and improve efficiency best to generalise 
#  *
#  *  Function takes in main tuning parameter (used for every model type), and secondary tuning parameter
#  *  (which is only used for B-Spline models - to represent level of degrees)
func_ModelCreationandPrediction <- function( base_modelType, tp_main, tp_Secondary, train, test ){
  # default values
  bsp_coreKnots <- 0; bsp_degree <- 0 
  # update variables accordingly
  poly_p <- ifelse( base_modelType == 'Polynomial', tp_main, 0)
  bs_binLength <- ifelse( base_modelType == 'BinSmooth', tp_main, 0)
  # calculate Knot Positions & level of Degree for BSpline
  if (base_modelType == 'BSpline') {
    bsp_baseKnots <- 1:tp_main / (tp_main + 1)
    bsp_coreKnots <- quantile(train[,2], probs = bsp_baseKnots)
    bsp_degree <- tp_Secondary }
  # create model from Train Data Data
  base_model <- func_ModelCreation( base_modelType, train[,1], train[,2], 
                                    param_pPoly = poly_p, 
                                    param_bs_binLength = bs_binLength,
                                    param_bsp_Knots = bsp_coreKnots, 
                                    param_bsp_Degree = bsp_degree)  
  # predict results
  predictions_model <- func_ModelPrediction(  base_modelType, 
                                              test[,1] , test[,2], 
                                              length(coef(base_model$modelOutput)), 
                                              base_model$modelOutput, 
                                              base_model$pred_BinSmooth_Table, 
                                              base_model$predhelp_q)
  # return results to user
  return( list( base_model = base_model,
                predictions_model = predictions_model))
}


# /* ~~~~~~
#  *  Custom function takes in a matrix_tuning array 
#  *  In practice, methodology used always creates a model and immediately predicts ~ to reduce code duplication
#  *  and improve efficiency best to generalise 
#  *
func_OptimumModelTuningValue <- function( matrix_tuning, varianceoffsetallowance, base_modelType, number_MatrixTuningColumns){
  # reference to smallest convergence level
  minvalue <- min(matrix_tuning[,6])
  # we want to only focus on those results that fall within smallest convergence level and offset allowance function parameter
  whatwewantbase <- matrix_tuning[ matrix_tuning[, 6] <= (minvalue + varianceoffsetallowance), ]
  # need to order data by second parameter for simpler complexity
  prediction_HighestAccuracywithinTolerance <- whatwewantbase[,5][order(whatwewantbase[,5],decreasing = TRUE)][1]
  # store only results that have the same highest Accuracy (adjusted R squared)  
  candidatResults <- whatwewantbase[ whatwewantbase[,5] == prediction_HighestAccuracywithinTolerance, ]
  # check how many results returned
  number_candidatResults <- length(candidatResults)/number_MatrixTuningColumns
  # if 1, we have our result - if not, order by complexity of model
  #WE WANT THE SIMPLEST MODEL WITH THE HIGHEST Accuracy (akin to low generalisation error)
  if( number_candidatResults == 1){
    result_TuningParam <- candidatResults
    optimumTuningParam_index <- result_TuningParam[1]
  } else {
    if( base_modelType == 'BinSmooth' ) {
      # for bin smooth we want the largest bin length (equivalent to simplest model)
      optimumTuningParam_index <- candidatResults[][order(candidatResults[,2], candidatResults[,4], decreasing = TRUE)][1]
    } else {
      # for b-spline & polynomial we want the smallest values (equivalent to simplest model)
      optimumTuningParam_index <- candidatResults[][order(candidatResults[,2], candidatResults[,3], candidatResults[,6], decreasing = FALSE)][1]
    }
    result_TuningParam <- whatwewantbase[,1:number_MatrixTuningColumns][whatwewantbase[,1]== optimumTuningParam_index]
  }
  # generate data frame format to return
  return_table <- cbind.data.frame( modelType = base_modelType,
                                    tuningParam_Main = result_TuningParam[2], 
                                    tuningParam_Secondary = result_TuningParam[3],
                                    model_AdjRSquared_TrainAccuracy = result_TuningParam[4],
                                    model_AdjRSquared_TestAccuracy = result_TuningParam[5],
                                    model_AccuracyVariance = result_TuningParam[6],
                                    model_RSS_TrainAccuracy = result_TuningParam[7],
                                    model_RSS_TestAccuracy = result_TuningParam[8] )
  #return to user
  return( list( return_table = return_table, 
                tuning_optimumParam_index = optimumTuningParam_index ))
}


# /* ~~~~~~
#  *  Custom Optimal Model Selection function which performs k-fold cross validation
#  *  By utilising the 'func_candidateModels' function to create 12 models (3 for the 4 variables) 
#  *  we can then generate k-fold sub sets using 'func_KFoldSubSets' function
#  *  Now its just a matter of looping through each of the models, perform k-fold iterations 
#  *  made easy using the 'func_ModelCreationandPrediction' function.
#  *
#  *  Store outputs to a table, and measure best model using Mean Square Error (MSE)
func_OptimalModelSelection <- function( data, x_indexes, y_index, traintestsplit, varianceoffsetallowance, n_row, k) {
  
  table_CandidateModels <- func_candidateModels( data, x_indexes, y_index, traintestsplit, varianceoffsetallowance)
  
  table_TuningResults <- table_CandidateModels$table_TuningResults
  
  kFold_Indexes <- func_KFoldSubSets(k, n_row, data )
  
  n_resultstovalidate <- nrow(table_TuningResults)
  generalisation_matrix <- matrix(0, n_resultstovalidate)
  
  for( i in 1:n_resultstovalidate) {
    # default values
    modelType <- table_TuningResults$modelType[i]
    tp_main <- table_TuningResults$tuningParam_Main[i]
    tp_Secondary <- table_TuningResults$tuningParam_Secondary[i]
    bsp_coreKnots <- 0 
    bsp_degree <- 0 
    x_loop <- table_TuningResults$x_loop[i]
    # update variables accordingly
    poly_p <- ifelse( modelType == 'Polynomial', tp_main, 0)
    bs_binLength <- ifelse( modelType == 'BinSmooth', tp_main, 0)
    # create data set
    baseDataSet <- data.frame(y = data[, y_index] , x = data[,x_loop])
    # scale data  
    baseDataSet[,2] <- scale(baseDataSet[,2])
    # now run through k-fold validation for required results  
    kfold_matrix <- matrix(0, k)
    for(j in 1:k) {
      # split into train and test data set 
      k_samp <- kFold_Indexes[[j]]
      train <- baseDataSet[-k_samp,]
      test <- baseDataSet[k_samp,]
      # now we can sort the data
      train <- func_Sort2DTable(train)
      test <- func_Sort2DTable(test)
      # create model and predictions
      modelcp <- func_ModelCreationandPrediction( modelType, tp_main, tp_Secondary, train, test )
      # update k-fold matrix with MSE result from prediction
      kfold_matrix[j] <- modelcp$predictions_model$pred_MSE
    }
    # take the average MSE result of k-fold matrix 
    generalisation_matrix[i] <- mean(kfold_matrix)
  }
  
  # if generalisation score above 100 we dont care - to high!
  generalisation_matrix <- ifelse(generalisation_matrix>100, NA, 
                                  round(generalisation_matrix, digits = 2))
  
  finalModelSelectionResults <- cbind( generalisation_matrix, table_TuningResults )
  
  rank <- matrix(0, n_resultstovalidate)
  rank[order(finalModelSelectionResults$generalisation_matrix)] <- 1:nrow(finalModelSelectionResults)
  
  finalModelSelectionResults <- cbind(rank, finalModelSelectionResults)
  
  return(finalModelSelectionResults)
}


# /* ~~~~~~
#  *  Custom Chart Function - Visualise identification of Optimal Model Parameters
#  *  Building custom chart functions are always 'fun', generalising legend, labels, headers, and colours  
#  *  
#  *  End result here is pretty sweet - hence the use of multiple charts in final paper 
#  *  Flexibility for user to view results by RSS (arent the easiest to interpret - or Adjusted R Squared (personal preference))
#  *  Added boolean 'return_CandidateModelInfo' variable if user would like to see full array of all info used to plot results
func_charting_OptimumModelParameters <- function( data, x_index = 4, y_index = 2, 
                                                 traintestsplit = 1, varianceoffsetallowance = 0.03, 
                                                 chartType = 'AdjR2', modelType = 'Polynomial',
                                                 charting_showTest = TRUE,
                                                 return_CandidateModelInfo = FALSE){
  
  if( chartType == 'RSS' ) { 
    train_index <- 7
    test_index <- 8 
  } else {
    train_index <- 4
    test_index <- 5
  }
  modelType_index <- ifelse( modelType == 'Polynomial', 1, ifelse( modelType == 'BinSmooth', 2, 3))
  chart_xaxis <- ifelse( modelType == 'Polynomial', 'Polynomial Degree', ifelse( modelType == 'BinSmooth', 'Bin Sizes', 'Number of Knot Quantiles'))
  
  
  charting_base <- func_candidateModels(data, c(x_index), y_index, traintestsplit, varianceoffsetallowance)
  
  charting_candidatemodelinfo <- charting_base$charting_CandidateModelInfo[[modelType_index]]
  #charting_candidatemodelinfo$return_table
  
  
  chart_title <- paste( modelType, '-', colnames(data)[x_index], '- Automatic Parameter Tuning' )
  chart_xaxislabel <- paste( modelType, '-', chart_xaxis)
  chart_yaxislabel <- ifelse(chartType == 'RSS', 'Residual Sum of Squares', 'Adjusted R Squared')
  
  if(chartType == 'RSS'){
    minrange <- (min( charting_candidatemodelinfo$tuning_matrix[,train_index], charting_candidatemodelinfo$tuning_matrix[,test_index])*0.95)
    maxrange <- (max( charting_candidatemodelinfo$tuning_matrix[,train_index], charting_candidatemodelinfo$tuning_matrix[,test_index])*1.05 )
  } else {
    minrange <- 0
    maxrange <- 1
  }
  
  plot(charting_candidatemodelinfo$tuning_matrix[,y_index], charting_candidatemodelinfo$tuning_matrix[,train_index], type = 'l', col = 'gray0', lwd = 3,
       ylim = c( minrange, maxrange), xlab = chart_xaxislabel, ylab = chart_yaxislabel,
       main = chart_title)
  
  
  if(chartType == 'RSS'){
    if( charting_showTest ){
      lines(charting_candidatemodelinfo$tuning_matrix[,y_index], charting_candidatemodelinfo$tuning_matrix[,test_index], col = 'blue', lwd = 3)
      legend( x="topright", legend=c("Train", 'Test'),
              col=c("grey0", "blue"), lwd=1, lty=c(1, 1))
    } else {
      legend( x="topright", legend=c("Train"),
              col=c("grey0"), lwd=1, lty=c(1))
    }
  } else {
    criteria_passed <- min(charting_candidatemodelinfo$tuning_matrix[,6]) + varianceoffsetallowance
    params_criteriaMet <- charting_candidatemodelinfo$tuning_matrix[,5]
    params_criteriaMet <- ifelse(charting_candidatemodelinfo$tuning_matrix[,1] == charting_candidatemodelinfo$tuning_optimumParam_index, -1,
                                 ifelse(charting_candidatemodelinfo$tuning_matrix[,6] <= criteria_passed, params_criteriaMet, -1 ))
    
    lines(charting_candidatemodelinfo$tuning_matrix[,y_index], charting_candidatemodelinfo$tuning_matrix[,test_index], col = 'blue', lwd = 3)
    lines(charting_candidatemodelinfo$tuning_matrix[,y_index], charting_candidatemodelinfo$tuning_matrix[,train_index]-criteria_passed, col = 'indianred1', lty = 2, lwd = 1)
    lines(charting_candidatemodelinfo$tuning_matrix[,y_index], charting_candidatemodelinfo$tuning_matrix[,train_index]+criteria_passed, col = 'indianred1', lty = 2, lwd = 1)
    
    points(charting_candidatemodelinfo$tuning_matrix[, y_index], params_criteriaMet, col = 'red', pch = 16, cex=1.5)
    points(charting_candidatemodelinfo$tuning_matrix[charting_candidatemodelinfo$tuning_optimumParam_index,y_index], charting_candidatemodelinfo$tuning_matrix[charting_candidatemodelinfo$tuning_optimumParam_index,test_index], col = 'green', pch = 18, cex=2)
    
    legend( x="topright", 
            legend=c("Train", "Test", 'Upper & Lower Variance Bounds', 'Optimum Value', 'Values within criteria'),
            col=c("grey0", "blue", 'indianred1'), lwd=1, lty=c(1,1, 2, NA, NA), 
            pch=c(NA, NA, NA, NA, NA) )
    legend( x="topright", 
            legend=c("Train", "Test", 'Upper & Lower Variance Bounds', 'Optimum Value', 'Values within criteria'),
            col=c("red", 'green'), lwd=1, lty=c(NA, NA, NA, NA, NA), 
            pch=c(NA, NA, NA, 18, 16) )
  }
  if(return_CandidateModelInfo) { return( charting_candidatemodelinfo ) }
}

# /* ~~~~~~
#  *  Custom Chart Function - Visualise final candidate models 
#  *  Models currently supported:   Polynomial Regression(n=1 ~ Linear), Bin Smooth & B-Spline
#  *  
#  *  This assignment opened the privelage of designing two customer chart functions - cause 1 wasnt enough!
#  *  End result again is pretty sweet, function is generalised to adapt to any of the three support model types
#  *  and works hand in hand with model creation/prediction function above leverging custom model objects returned by those functions
#  *  
#  *  ... okay i did cheat - those custom model objects were specifically tailored to work with functions like these
#  *  cause whats the point in a model if we cant view it?
#  *
#  *  The dynamic legend features works pretty well, however if users want to plot multiple charts
#  *  on one pane, the legend can clutter the chart so adding boolean to remove/show when desired is handy
func_charting_CandidateModels <- function( modelType, tp_main, tp_Secondary, train, xlab = '', showlegend = FALSE) {
  # create model and predictions
  train <- func_Sort2DTable(baseDataSet)
  
  modelcp <- func_ModelCreationandPrediction( modelType, tp_main, tp_Secondary, train, train )
  
  x <- train[,2]
  y <- train[,1]
  
  chart_title <- paste( modelType, 'Model' )
  chart_xaxislabel <- ifelse(xlab != '', paste( 'Scaled', xlab), '')
  
  plot(x,y, main=chart_title, xlab = chart_xaxislabel, ylab = 'MPG', col="blue")
  #pch=20, 
  if ( modelType == 'BinSmooth') {
    # plot bin smooth line
    numberCoefs <- modelcp$base_model$predhelp_q
    coefs_matrix <- coef(modelcp$base_model$modelOutput)
    binlength <- modelcp$base_model$chart_BinSmooth_BinLength
    xstart <- modelcp$base_model$chart_BinSmooth_Start
    xend <- modelcp$base_model$chart_BinSmooth_End
    
    j<-1
    for(i in 1:numberCoefs) {
      if(i>1) lines(c(x[xend],x[xend]), c(as.numeric(coefs_matrix[i-1]), as.numeric(coefs_matrix[i])), col="grey0", lwd=3)
      xstart = j
      if(i>1) lines(c(x[xend],x[xstart]), c(as.numeric(coefs_matrix[i]), as.numeric(coefs_matrix[i])), col="grey0", lwd=3)
      xend = min(j+binlength-1, length(x))
      lines(c(x[xstart],x[xend]), rep(as.numeric(coefs_matrix[i]), 2), col="grey0", lwd=3)
      j<-j+binlength
    }
  } else {
    lines(x, fitted(modelcp$base_model$modelOutput), col='grey0', type='l', lwd=3) 
  }
  
  if(showlegend){
    legend( x="topright", 
            legend=c("Model based on entire data", 'All Data Points'),
            col=c("grey0"), lwd=1, lty=c(1, NA), 
            pch=c(NA, NA) )
    legend( x="topright", 
            legend=c("Model based on entire data", 'All Data Points'),
            col=c( 'blue'), lwd=1, lty=c(NA, NA), 
            pch=c(NA, 1) )
  }
}
