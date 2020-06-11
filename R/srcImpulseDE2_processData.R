### Check and process input

#' Check and process input to runImpulseDE2()
#' 
#' Check validity of input and process count data matrix and annotation
#' into data structures used later in \link{runImpulseDE2}.
#' \link{processData} is structure in the following way:
#' \itemize{
#'    \item Subhelper functions:
#'    \itemize{
#'      \item checkNull() Check whether object was supplied (is not NULL).
#'      \item checkDimMatch() Checks whether dimensions of matrices agree.
#'      \item checkElementMatch() Checks whether vectors are identical.
#'      \item checkNumeric() Checks whether elements are numeric.
#'      \item checkProbability() Checks whether elements are probabilities.
#'      \item checkCounts() Checks whether elements are count data.
#'    }
#'    \item Helper functions:
#'    \itemize{
#'      \item checkData() Check format and presence of input data.
#'      \item nameGenes() Name genes if names are not given.
#'      \item procDESeqAnnotation() Add categorial time variable to 
#'      annotation table.
#'      Add nested batch column if necessary.
#'      Reduce to samples used.
#'      \item reduceCountData() Reduce count data to data which are 
#'      utilised later.
#'    }
#'    \item Script body
#' }
#' 
#' @seealso Called by \link{runImpulseDE2}.
#' 
#' @param matCountData (matrix genes x samples) [Default NULL] 
#' Read count data, unobserved entries are NA.
#' @param dfAnnotation (data frame samples x covariates) 
#' {Sample, Condition, Time (numeric), TimeCateg (str)
#' (and confounding variables if given).}
#' Annotation table with covariates for each sample.
#' @param boolCaseCtrl (bool) 
#' Whether to perform case-control analysis. Does case-only
#' analysis if FALSE.
#' @param boolBeta2 (bool)
#' Whether to model two different slopes for impulse model instead of 
#' assuming onset slope and offset slope are identical.
#' @param vecCovFactor (vector of strings number of categorical covariates)
#' Categorial covariates to adjust for.
#' Names refer to columns in dfAnnotation.
#' @param vecCovContinuous (vector of strings number of continuous covariates)
#' Continuous covariates to adjust for.
#' Names refer to columns in dfAnnotation.

#' @param vecDispersionsExternal (vector length number of
#' genes in matCountData) [Default NULL]
#' Externally generated list of gene-wise dispersion factors
#' which overides DESeq2 generated dispersion factors.
#' @param vecSizeFactorsExternal (vector length number of
#' cells in matCountData) [Default NULL]
#' Externally generated list of size factors which override
#' size factor computation in ImpulseDE2.
#'    
#' @return (list length 5)
#' \itemize{
#'    \item matCountDataProc (matrix genes x samples)
#'    Read count data.
#'    \item dfDESeqAnnotationProc (data frame samples x covariates) 
#'    {Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and confounding variables if given).}
#'    Processed annotation table with covariates for each sample.
#'    \item lsdfCovProc (list of data frames length 3)
#'    One data frame each for "case", "control", and "combined".
#'    Continous covariates are centered and scaled within each group.
#'    Categorical covariates are factored within each group. 
#'    \item vecSizeFactorsExternalProc (numeric vector number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors).
#'    \item vecDispersionsExternalProc (vector number of genes) Gene-wise 
#'    negative binomial dispersion hyper-parameter.
#'    \item strReportProcessing (str) String of stdout of processData().
#' }
#' 
#' @author David Sebastian Fischer
processData <- function(
    dfAnnotation, 
    matCountData,
    boolCaseCtrl, 
    boolBeta2,
    vecCovFactor,
    vecCovContinuous,
    vecDispersionsExternal,
    vecSizeFactorsExternal){
    
    ###############################################################
    # (I) Helper functions
    
    # Check whether object was supplied (is not NULL).
    checkNull <- function(objectInput,strObjectInput){
        if(is.null(objectInput)){
            stop(paste0(strObjectInput,
                         " was not given as input." ))
        }
    }
    # Check whether object does not have NA elements.
    checkNA <- function(objectInput,strObjectInput){
        if(is.na(objectInput)){
            stop(paste0(strObjectInput,
                         " is NA and needs to be specifief." ))
        }
    }
    # Checks whether dimensions of matrices agree.
    checkDimMatch <-function(
        matInput1, matInput2, strMatInput1, strMatInput2){
        if(any(dim(matInput1)!=dim(matInput2))){
            stop(paste0(strMatInput1, 
                         " does not have the dimensions as ", 
                         strMatInput2 , "." ))
        }
    }
    # Checks whether vectors are identical.
    checkElementMatch <- function(vec1, vec2, strVec1, strVec2){
        if(!any(vec1==vec2)){
            stop(paste0(strVec1 ," do not agree with ", 
                         strVec2, "." ))
        }
    }
    # Checks whether elements are numeric
    checkNumeric <- function(matInput, strMatInput){
        if(any(!is.numeric(matInput))){
            stop(paste0(strMatInput, 
                         " contains non-numeric elements. ",
                         "Requires numeric data." ))
        }
    }
    # Checks whether elements are probabilities (in [0,1]).
    checkProbability <- function(matInput, strMatInput){
        checkNumeric(matInput, strMatInput)
        if(any(matInput < 0 | matInput > 1 | is.na(matInput))){
            stop(paste0(strMatInput, 
                         " contains elements outside",
                         " of interval [0,1]." ))
        }
    }
    # Checks whether elements are count data: 
    # non-negative integer finite numeric elements.
    # Note that NA are allowed. 
    # Can be used to check whether element is integer if NA 
    # is checked separately.
    checkCounts <- function(matInput, strMatInput){
        checkNumeric(matInput, strMatInput)
        if(any(matInput[!is.na(matInput)] %% 1 != 0)){
            stop(paste0(strMatInput, 
                         " contains non-integer elements.",
                         " Requires count data." ))
        }
        if(any(!is.finite(matInput[!is.na(matInput)]))){
            stop(paste0(strMatInput, 
                         " contains infinite elements.",
                         " Requires count data." ))
        }
        if(any(matInput[!is.na(matInput)]<0)){
            stop(paste0(strMatInput, 
                         " contains negative elements.",
                         " Requires count data." ))
        }
    }
    
    # Checks whether a model matrix is full rank
    checkFullRank <- function(df){
        matModelMatrix = model.matrix(~., data=df)
        if(ncol(matModelMatrix) != rankMatrix(matModelMatrix)){
            stop(paste0("Matrix with column names {",
                paste0(colnames(df),collapse=','), "} is not full rank.",
                " Check for confounding among these variables."))
        }
    }

    # Check format and presence of input data.
    checkData <- function(
        dfAnnotation, 
        matCountData,
        boolCaseCtrl,
        boolBeta2, 
        vecCovFactor,
        vecCovContinuous,
        vecDispersionsExternal,
        vecSizeFactorsExternal,
        strReportProcessing){
        
        strReportProcessing <- "Processing Details:"
        ### 1. Check that all necessary input was specified
        checkNull(dfAnnotation,"dfAnnotation")
        checkNull(matCountData,"matCountData")
        checkNull(boolCaseCtrl,"boolCaseCtrl")
        checkNull(boolBeta2,"boolBeta2")
        
        
        ### 2. Check annotation table content
        ### a) Check column names
        vecColNamesRequired <- c("Sample","Condition","Time",
                                 vecCovFactor, vecCovContinuous)
        if( !all(vecColNamesRequired %in% 
                 colnames(dfAnnotation)) ){
            stop(paste0(
                "Could not find column ",
                vecColNamesRequired[
                    !(vecColNamesRequired %in% colnames(dfAnnotation))],
                " in annotation table."))
        }
        ### b) Samples
        # Check that sample name do not occur twice
        if(any(duplicated(dfAnnotation$Sample))){
            stop(paste0(
                "[Annotation table] ",
                "Sample names must be unique: Sample(s) ",
                paste0((dfAnnotation$Sample)[
                    duplicated(dfAnnotation$Sample)],collapse=","),
                " is/are duplicated."))
        }
        ### c) Time points
        vecTimepoints <- unique(as.vector( dfAnnotation$Time ))
        # Check that time points are numeric
        checkNumeric(dfAnnotation$Time, "dfAnnotation$Time")
        ### d) Conditions
        lsConditions <- unique( dfAnnotation$Condition )
        # Check that given conditions exisit in annotation table
        if(!("case" %in% lsConditions)){
            stop(paste0(
                "Condition \"case\" does",
                " not occur in annotation table condition column."))
        }
        if(boolCaseCtrl){
            if(!("control" %in% lsConditions)){
                stop(paste0(
                    "Condition \"control\" does",
                    " not occur in annotation table condition column."))
            }
        }
        ### e) Categorial covariates
        if(!is.null(vecCovFactor)){
            # Dummy check: Check that number of batches 
            # for each confounding variable is > 1 and < length(vec)
            for(confounder in vecCovFactor){
                if(length(unique( dfAnnotation[,confounder] ))==1){
                    stop(paste0(
                        "Model matrix based on ",
                        "categorial covariates {", vecCovFactor,
                        "} is not full rank: Only one level",
                        " given for categorical covariate ", confounder, 
                        ". Remove from vecCovFactor or correct",
                        " dfAnnotation."))
                }
                if(length(unique( dfAnnotation[,confounder] )) == length(dfAnnotation[,confounder])){
                    stop(paste0(
                        "vecCovFactor ", confounder, " has as many levels as there are samples.",
                        " If this is supposed to be a continuous covariate,",
                        " define it in 'vecCovContinuous' instead of 'vecCovFactor'."))
                }

            }

        }
        ### f) Continuous covariates
        if(!is.null(vecCovContinuous)){
            # Dummy check: Check that number of unique values 
            # for each covariate is > 1
            # Check numeric
            for(confounder in vecCovContinuous){
                if(length(unique( dfAnnotation[,confounder] ))==1){
                    stop(paste0(
                        "Model matrix based on continous covariates {",
                        paste0(vecCovContinuous, collapse=','),
                        "} is not full rank: Only one unique value",
                        " given for continuous covariate ", confounder, 
                        ". Remove from vecCovContinuous or correct",
                        " dfAnnotation."))
                }
                if(!is.numeric(dfAnnotation[,confounder])){
                    stop(paste0(
                        "vecCovContinuous ", confounder, " is not numeric.",
                        " If this is supposed to be a categorial covariate,",
                        " define it in 'vecCovFactor' instead of 'vecCovContinuous'."))
                }
            }
        }
        ### g) More detailed full rank check
        if(!is.null(vecCovContinuous) | !is.null(vecCovFactor)){
            dfTmp = dfAnnotation[,c(vecCovContinuous, vecCovFactor)]
            # Already verified that vecCovContinuous is numeric
            # Coerce vecCovFactor to character
            for(covFactor in vecCovFactor){
                dfTmp[,covFactor] = paste0('_', dfTmp[,covFactor])
            }
            checkFullRank(dfTmp)
        }
        
        ### 3. Check expression table content
        # Check that all entries in count data table occur in annotation table
        if( any(!(colnames(matCountData) %in% dfAnnotation$Sample)) ){
            warning(paste0("The column(s) ",
                paste0(as.character( colnames(matCountData)[
                    !(colnames(matCountData) %in% dfAnnotation$Sample)] ),
                    collapse=","),
                " in the count data table do(es) not occur",
                " in annotation table and will be ignored.\n"))
        }
        checkNull(rownames(matCountData),"[Rownames of matCountData]")
        checkCounts(matCountData, "matCountData")
        
        ### 4. Check supplied dispersion vector
        if(!is.null(vecDispersionsExternal)){
            # Check that dispersion parameters were named
            if(is.null(names(vecDispersionsExternal))){
                stop(paste0(
                    "vecDispersionsExternal was not named.",
                    " Name according to rownames of matCountData."))
            }
            # Check that one dispersion parameter factors was supplied per gene
            if(any( !(names(vecDispersionsExternal) %in% 
                      rownames(matCountData)) ) |
               any( !(rownames(matCountData) %in% 
                      names(vecDispersionsExternal)) ) ) {
                stop(paste0(
                    "vecDispersionsExternal supplied but names",
                    " do not agree with rownames of matCountData."))
            }
            # Check that dipsersion parameter vector is numeric
            checkNumeric(vecDispersionsExternal, "vecDispersionsExternal")
            # Check that dispersion parameters are positive
            if(any(vecDispersionsExternal <= 0)){
                stop(paste0(
                    "vecDispersionsExternal contains negative",
                    " or zero elements which leads.",
                    "Dispersion parameters must be positive." ))
            }
        }
        
        ### 5. Check supplied size facotrs
        if(!is.null(vecSizeFactorsExternal)){
            # Check that size factors were named
            if(is.null(names(vecSizeFactorsExternal))){
                stop(paste0("vecSizeFactorsExternal was not named.",
                            " Name according to colnames of matCountData."))
            }
            # Check that one size factors was supplied per cell
            if(any( !(names(vecSizeFactorsExternal) %in% 
                      colnames(matCountData)) ) |
               any( !(colnames(matCountData) %in% 
                      names(vecSizeFactorsExternal)) ) ){
                stop(paste0(
                    "vecSizeFactorsExternal supplied but",
                    " names do not agree with colnames of matCountData."))
            }
            # Check that size factors vector is numeric
            checkNumeric(vecSizeFactorsExternal, "vecSizeFactorsExternal")
            # Check that size factors are positive
            if(any(vecSizeFactorsExternal <= 0)){
                stop(paste0(
                    "vecSizeFactorsExternal contains negative ",
                    "or zero elements which leads.",
                    "Size factors must be positive, remove samples if size ",
                    "factor is supposed to be zero." ))
            }
        }
        
        ### Summarise which mode, conditions, samples and
        ### batch were found
        if(boolCaseCtrl) { 
            strReportProcessing <- paste0(
                strReportProcessing,"\n","ImpulseDE2 runs in case-ctrl mode.") 
        } else { 
            strReportProcessing <- paste0(
                strReportProcessing, "\n","ImpulseDE2 runs in case-only mode.") 
        }
        if(boolBeta2) { 
            strReportProcessing <- paste0(
                strReportProcessing,"\n","ImpulseDE2 fits an impulse model with two slopes.") 
        } else { 
            strReportProcessing <- paste0(
                strReportProcessing, "\n","ImpulseDE2 fits an impulse model with a single slope.") 
        }
        strReportProcessing <- paste0(
            strReportProcessing, "\n",
            paste0( "Found time points: ",
                    paste( vecTimepoints, collapse=",") ))
        for(tp in vecTimepoints){
            strReportProcessing <- 
                paste0(strReportProcessing, "\n",
                       paste0(
                           "Case: Found the samples at time point ", 
                           tp,": ", paste0(dfAnnotation[
                               (dfAnnotation$Time %in% tp) &
                                   (dfAnnotation$Condition %in% "case") &
                                   (dfAnnotation$Sample %in% 
                                        colnames(matCountData)),
                               ]$Sample,collapse=","),collapse="," ) )
        }
        if(boolCaseCtrl){
            for(tp in vecTimepoints){
                strReportProcessing <- paste0(
                    strReportProcessing, "\n", paste0(
                        "Control: Found the following samples at time point ", 
                        tp, ":", paste0(dfAnnotation[
                            dfAnnotation$Time %in% tp &
                                dfAnnotation$Condition %in% "control" &
                                dfAnnotation$Sample %in% colnames(matCountData),
                            ]$Sample,collapse=","),collapse="," ) )
            }
        }
        if(!is.null(vecCovFactor)){
            for(confounder in vecCovFactor){
                for(batch in unique( dfAnnotation[,confounder] )){
                    strReportProcessing <- paste0(
                        strReportProcessing, "\n", 
                        "Found the following samples for categorial covariate ", 
                        confounder," and level ", batch, ": ",
                        paste0( dfAnnotation[
                            dfAnnotation[,confounder] %in% batch &
                                dfAnnotation$Sample %in% colnames(matCountData),
                            ]$Sample, collapse=",") )
                }
            }
        }
        if(!is.null(vecCovContinuous)){
            strReportProcessing <- paste0(
                strReportProcessing, "\n",
                "Adjusting for the following continuous covariates: ",
                "{",vecCovContinuous,"}.")
        }
        
        return(strReportProcessing)
    }
    
    # Add categorial time variable to annotation table which
    # differentiates case and control time points (given to DESeq2).
    # Add column with time scalars with underscore prefix.
    # Reduce to samples used.
    procDESeqAnnotation <- function(dfAnnotation,
                               matCountData,
                               boolCaseCtrl,
                               vecCovFactor,
                               vecCovContinuous){

        # Make sure factors are characters
        for(cov in c(vecCovFactor, "Condition", "Sample")){
            if(is.numeric(dfAnnotation[,cov])){
                dfAnnotation[,cov] = paste0('_', dfAnnotation[,cov])
            }
        }
        
        # Make sure vecCovContinuous are numeric 
        for(cov in c(vecCovContinuous, "Time")){
            dfAnnotation[,cov] = as.numeric(dfAnnotation[,cov])
        }
        
        if(!boolCaseCtrl){
            dfAnnotationProc <- 
                dfAnnotation[dfAnnotation$Condition=="case",]
        } else {
            dfAnnotationProc <- 
                dfAnnotation[dfAnnotation$Condition=="case" |
                                 dfAnnotation$Condition=="control",]
        }
        
        # Reduce to samples which occur in count table
        # This allows to represent entire library of samples in annotation
        # file even if not all samples were measured yet. Not measured samples
        # which are not mentioned in count table are ignored. Thereby, the same
        # full annotation table can be used throughout sample collection on 
        # incomplete data sets.
        dfAnnotationProc <- dfAnnotationProc[dfAnnotationProc$Sample %in% 
                                                 colnames(matCountData),]
        
        # Take out columns which are not used
        dfAnnotationProc <- dfAnnotationProc[,c("Sample", "Time", "Condition", 
                                                vecCovFactor, vecCovContinuous)]
        # Add categorial time column for DESeq2
        dfAnnotationProc$TimeCateg <- paste0("_", dfAnnotationProc$Time)
        return(dfAnnotationProc)
    }

    # Split samples into case, control, and combined
    # For each group of samples, center and scale continous covariates;
    # Convert categorical covariates to strings 
    # Return a list of data frames, one data frame per group of samples 
    procCov <- function(dfAnnotation, boolCaseCtrl, vecCovFactor, vecCovContinuous){
        if(is.null(c(vecCovFactor, vecCovContinuous))){
            return(NULL)
        }
        # Get cases, which applies whether or not there are controls 
        dfCase = dfAnnotation[dfAnnotation$Condition == "case", c(vecCovFactor, vecCovContinuous, "Sample")]
        if(boolCaseCtrl){
            dfControl = dfAnnotation[dfAnnotation$Condition == "control", c(vecCovFactor, vecCovContinuous, "Sample")]
            dfCombined = dfAnnotation[dfAnnotation$Condition == "control|case", c(vecCovFactor, vecCovContinuous, "Sample")]
        }else{
            dfControl = NULL
            dfCombined = NULL
        }
        lsdfCov = list(case = dfCase,
            control = dfControl,
            combined = dfCombined)
        # Handle covariates separately in each group of samples
        lsdfCovProc = lapply(lsdfCov, function(df){
            if(is.null(df)){
                return(NULL)
            }
            # Center and scale continuous covariates
            for(cov in vecCovContinuous){
                df[,cov] = scale(df[,cov], center=T, scale=T)
            }
            # Convert categorical covariates to string
            for(cov in vecCovFactor){
                if(!is.character(df[,cov])){
                    df[,cov] = paste0('_', df[,cov])
                }
            }
            rownames(df) = df$Sample
            df$Sample = NULL
            return(df)
            })
        return(lsdfCovProc)
    }

    # Calculate variance inflation factor for all covariates, including Time and Condition (when applicable)
    checkConfounding <- function(dfDESeqAnnotationProc, boolCaseCtrl, 
        vecCovFactor, vecCovContinuous, lsdfCovProc){
            
        if(boolCaseCtrl){
            # check the "combined" data frame 
            df = lsdfCovProc[["combined"]]
            # get condition
            df$Condition = dfDESeqAnnotationProc$Condition[match(dfDESeqAnnotationProc$Sample, rownames(df))]
            # format contrast
            contrast = paste0('tmp ~ ', paste0(c(vecCovFactor, vecCovContinuous, 'Condition', 'Time'), collapse = " + "))
        }else{
            # check the "case" data frame 
            df = lsdfCovProc[["case"]]
            # format contrast
            contrast = paste0('tmp ~ ', paste0(c(vecCovFactor, vecCovContinuous, 'Time'), collapse = " + "))
        }

        # get timepoints 
        df$Time = dfDESeqAnnotationProc$Time[match(dfDESeqAnnotationProc$Sample, rownames(df))]
        # get random outcome
        df$tmp = seq(1:nrow(df))

        # linear regression 
        model = lm(eval(parse(text=contrast)), data=df)
        vif = as.data.frame(ols_vif_tol(model))

        if(any(vif$VIF > 10)){
            collinear_cov = vif[vif$VIF > 10, "Variables"]
            warning(paste0("A variance inflation factor > 10 has been identified ",
                "for the following covariates, which indicates serious collinearity needing correction: {", 
                paste(collinear_cov, collapse=','), "}. We strongly recommend considering the correlation structure ",
                "between covariates before proceeding.\n"))
        }else if(any(vif$VIF > 4)){
            collinear_cov = vif[vif$VIF > 4, "Variables"]
            warning(paste0("A variance inflation factor > 4 has been identified ",
                "for the following covariates, which warrants further investigation: {", 
                paste(collinear_cov, collapse=','), "}.\n"))
        }

    }
    
    # Name genes if names are not given.
    nameGenes <- function(matCountDataProc){
        if(is.null(rownames(matCountDataProc))){
            rownames(matCountDataProc) <- 
                paste0("Gene_", seq(1,nrow(matCountDataProc)))
        }
        return(matCountDataProc)
    }
    
    # Reduce count data to data which are utilised later
    reduceCountData <- function(dfAnnotation, matCountDataProc){
        
        strReportCountRed <- paste0(
            "Input contained ",dim(matCountDataProc)[1]," genes/regions.")
        ### 1. Columns (Conditions):
        # Reduce expression table to columns of considered conditions
        if(boolCaseCtrl){
            vecSampleNames_Case <- as.character(as.vector( 
                dfAnnotation[ 
                    dfAnnotation$Condition %in% "case" &
                        dfAnnotation$Sample %in% colnames(matCountDataProc),
                    ]$Sample ))
            vecSampleNames_Ctrl <- as.character(as.vector( 
                dfAnnotation[ 
                    dfAnnotation$Condition %in% "control" &
                        dfAnnotation$Sample %in% colnames(matCountDataProc),
                    ]$Sample ))
            matCountDataProc <- 
                matCountDataProc[,c(vecSampleNames_Case,vecSampleNames_Ctrl)]
        } else {
            vecSampleNames_Case <- as.character(as.vector( 
                dfAnnotation[ 
                    dfAnnotation$Condition %in% "case" &
                        dfAnnotation$Sample %in% colnames(matCountDataProc),
                    ]$Sample ))
            matCountDataProc <- matCountDataProc[,vecSampleNames_Case]
        }
        # Check that every sample contains at least one observed value (not NA)
        vecboolNASample <- apply(matCountDataProc,2,function(sample)
            all(is.na(sample)) )
        if(any(vecboolNASample)){
            warning(paste0("Sample(s) ",
                paste0(colnames(matCountDataProc)[vecboolNASample], 
                       collapse=","),
                " only contain(s) NA values and will be ",
                "removed from the analysis.\n"))
            matCountDataProc <- matCountDataProc[,!vecboolNASample]
        }
        # Trigger warning if all zero sample is encountered - these are kept!
        vecboolAllZeroSample <- apply(matCountDataProc,2,function(sample) 
            all(sample[!is.na(sample)]==0) )
        if(any(vecboolAllZeroSample)){
            warning(paste0("Sample(s) ",
                paste0(colnames(matCountDataProc)[vecboolAllZeroSample], 
                       collapse=","),
                " only contain(s) zeros (and NAs).",
                " These samples are kept for analysis.\n"))
        }
        
        ### 2. Rows (Genes):
        # Exclude genes with only missing values (NAs)
        vecboolNonZeroGene <- apply(matCountDataProc,1,function(gene)  
            any(!is.na(gene) & gene > 0) )
        matCountDataProc <- matCountDataProc[vecboolNonZeroGene,]
        if(sum(!vecboolNonZeroGene) > 0){
            warning(paste0(sum(!vecboolNonZeroGene), " out of ",
                length(vecboolNonZeroGene),
                " genes do not have obserserved non-zero counts",
                " and are excluded.\n"))
        }
        
        # Sort count matrix column by annotation table
        matCountDataProc <- matCountDataProc[
            ,match(as.vector(dfAnnotation$Sample),colnames(matCountDataProc))]
        
        strReportCountRed <- paste0(strReportCountRed,  
                                    "\n","Selected ",dim(matCountDataProc)[1],
                                    " genes/regions for analysis.")
        return(list( matCountDataProc=matCountDataProc,
                     strReportCountRed=strReportCountRed) )
    }
    
    ###############################################################
    # (II) Main body of function
    
    # Check validity of input
    strReportProcessing <- checkData(
        dfAnnotation=dfAnnotation,
        matCountData=matCountData,
        boolCaseCtrl=boolCaseCtrl,
        boolBeta2=boolBeta2,
        vecCovFactor=vecCovFactor,
        vecCovContinuous=vecCovContinuous,
        vecDispersionsExternal=vecDispersionsExternal,
        vecSizeFactorsExternal=vecSizeFactorsExternal )

    # Process annotation table
    dfDESeqAnnotationProc <- procDESeqAnnotation(dfAnnotation=dfAnnotation,
                                       matCountData=matCountData,
                                       boolCaseCtrl=boolCaseCtrl,
                                       vecCovFactor=vecCovFactor,
                                       vecCovContinuous=vecCovContinuous)

    # Process data frames for model matrices
    lsdfCovProc <- procCov(dfAnnotation=dfAnnotation,
        boolCaseCtrl=boolCaseCtrl,
        vecCovFactor=vecCovFactor,
        vecCovContinuous=vecCovContinuous)

    # Check for confounding and collinearity 
    checkConfounding(
        dfDESeqAnnotationProc=dfDESeqAnnotationProc,
        boolCaseCtrl=boolCaseCtrl,
        vecCovFactor=vecCovFactor,
        vecCovContinuous=vecCovContinuous,
        lsdfCovProc=lsdfCovProc)
    
    # Process raw counts
    matCountDataProc <- nameGenes(matCountDataProc=matCountData)
    lsReduceCounts <- reduceCountData(
        dfAnnotation=dfDESeqAnnotationProc, 
        matCountDataProc=matCountDataProc)
    matCountDataProc <- lsReduceCounts$matCountDataProc
    strReportProcessing <- paste0(strReportProcessing, "\n",
                                  lsReduceCounts$strReportCountRed)
    
    # Reduce externally provided parameters according to reduced data set
    # and reorder according to given data set.
    if(!is.null(vecSizeFactorsExternal)){
        vecSizeFactorsExternalProc <- 
            vecSizeFactorsExternal[colnames(matCountDataProc)]
    } else { vecSizeFactorsExternalProc <-NULL }
    if(!is.null(vecDispersionsExternal)){
        vecDispersionsExternalProc <- 
            vecDispersionsExternal[rownames(matCountDataProc)]
    } else { vecDispersionsExternalProc <-NULL }
    
    return( list(matCountDataProc           = matCountDataProc,
                 dfDESeqAnnotationProc      = dfDESeqAnnotationProc,
                 lsdfCovProc                = lsdfCovProc,
                 vecSizeFactorsExternalProc = vecSizeFactorsExternalProc,
                 vecDispersionsExternalProc = vecDispersionsExternalProc,
                 strReportProcessing        = strReportProcessing ) )
}