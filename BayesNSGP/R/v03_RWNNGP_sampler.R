
#================================================================================================
# Random Walk Metropolis-Hastings MCMC sampler for the Nearest-Neighbor Gaussian Process (NNGP).
# Author: Fabian R. Ketwaroo
# Affiliation: Swiss Ornithological Institute
#================================================================================================

#================================================
# Single Node Reverse Neighbors
#================================================

#' Find reverse neighbors for a single NNGP node
#'
#' Identifies all "children" (reverse neighbors) for a specific target location 
#' \eqn{j} in a Nearest-Neighbor Gaussian Process (NNGP). It searches the neighbor matrix 
#' to find all locations \eqn{k} that condition upon location \eqn{j}.
#'
#' @param target_node Integer. The index of the location \eqn{j} whose reverse 
#'   neighbors you want to find.
#' @param neighbor_idx Integer matrix (\eqn{M \times k}). Neighbor indices for all locations, 
#'   output from \code{computeNeighbors}.
#'
#' @return An integer vector containing the indices of all locations \eqn{k} 
#'   that have \code{target_node} as a neighbor, maintaining the NNGP 
#'   ordering constraint (\eqn{j < k}).
#'
#' @details
#' This function optimizes the reverse lookup for a single node by avoiding the 
#' full \eqn{O(M)} list construction loop. It uses vectorization to scan rows of 
#' the neighbor matrix where the row index \eqn{k} is greater than the target index 
#' \eqn{j}, enforcing the Directed Acyclic Graph (DAG) structure of the NNGP.
#'
#' @examples
#' nn_matrix <- matrix(c(0, 0,
#'                       1, 0,
#'                       1, 2,
#'                       1, 3), nrow = 4, byrow = TRUE)
#' 
#' # Find which nodes have node 1 as a neighbor
#' get_single_reverse_neighbors(target_node = 1, neighbor_idx = nn_matrix)
#'
#' @author Fabian Ketwaroo
#'
#' @export
get_single_reverse_neighbors <- function(target_node, neighbor_idx) {
    
    NN <- neighbor_idx
    M <- nrow(NN)
    
    # Edge case: The last node cannot be a neighbor to any future nodes
    if (target_node >= M) {
        return(integer(0))
    }
    
    # Only search rows after the target_node to enforce j < k (NNGP ordering)
    search_range <- (target_node + 1):M
    
    # Look at the neighbor subsets for those future rows
    sub_NN <- NN[search_range, , drop = FALSE]
    
    # Find which rows contain the target_node
    # rowSums checks across columns for matches
    matches <- rowSums(sub_NN == target_node) > 0
    
    # Return the actual global indices (k) that matched
    return(search_range[matches])
}

#================================================
# NNGP Sampler Control Setup
#================================================

#' Helper to extract NNGP local neighborhood nodes for custom MCMC sampling
#'
#' Pre-calculates the local graph structure, target indices, and specific 
#' coordinate strings within the \code{AD} matrix for a single target node. 
#' This forms the \code{control} list payload needed by \code{sampler_RW_NN_GP}.
#'
#' @param node_id Integer. The index of the specific spatial location currently 
#'   being set up for targeted random-walk sampling (ranges from 1 to \eqn{M}).
#' @param AD Character. The name of the matrix containing the NNGP coefficients 
#'   within the NIMBLE model object (typically passed as a string like \code{"AD"}).
#' @param neighbors.id Integer matrix (\eqn{M \times k}). The forward neighbor 
#'   index tracking matrix where rows indicate the target location.
#' @param Rneighbors.id Integer vector. The reverse neighbor indices (children) 
#'   for \code{node_id}, computed by \code{get_single_reverse_neighbors}.
#' @param N.neighbors Integer vector. A vector of length \eqn{M} storing the 
#'   exact number of active forward neighbors assigned to each location.
#' @param k Integer. The maximum number of neighbors specified in the NNGP configuration.
#'
#' @return A named list containing specific structural elements for the target node:
#' \itemize{
#'   \item \code{update_id}: Integer vector combining the \code{node_id} and its reverse neighbors.
#'   \item \code{Fneighbors.id}: Integer vector containing non-zero forward neighbor indices.
#'   \item \code{AFnodes}: Character vector of parsed string addresses pointing to forward coefficients.
#'   \item \code{ARnodes}: Character vector of parsed string addresses pointing to reverse coefficients.
#'   \item \code{A.neighbors}: Flattened character vector of matrix addresses representing all relevant 
#'         neighborhood coefficients needed to re-evaluate residuals during the loop.
#' }
#'
#' @details
#' This utility is designed to run in an R loop when configuring an MCMC specification 
#' prior to model compilation. By converting matrix indexing operations into explicit 
#' character node paths (e.g., \code{"AD[5,1]"}), it shifts the burden of matrix searching 
#' from the runtime execution loop of the C++ compiled sampler into a one-time R setup cost.
#'
#' It maps out how a change in \code{node_id} will ripple through its forward neighbors 
#' and back through the reverse neighbors whose conditional distributions depend directly 
#' on the target node's value.
#'
#' @author Fabian Ketwaroo
#'
#' @export
RWNNGP_setup <- function(node_id, AD, neighbors.id, Rneighbors.id, N.neighbors, k = k ){
    
    # Forward neighbors
    Fneighbors.id = neighbors.id[node_id, ][neighbors.id[node_id, ] != 0]
    
    Fn <- length(Fneighbors.id)
    
    if(Fn>0){
        
        AFnodes <- paste0(AD,"[", node_id, ",", 1:Fn, "]")
        
    } else  AFnodes <- character(0)
    
    
    J = length(Rneighbors.id)
    if(J>0){
        
        pos <- ARnodes <-  numeric(J) 
        for (j in 1:J) {
            r.id <-  Rneighbors.id[j]
            pos[j] <- which(neighbors.id[r.id,1:k] == node_id) # Position of i in neighbor list of k
            ARnodes[j] <- paste0(AD,"[", r.id, ",", pos[j], "]")
        } 
        
        update_id <- c(node_id, Rneighbors.id)
        
    } else{
        
        pos <- integer(0)
        ARnodes <- character(0)
        update_id <- c(node_id)
        
    } 
    
    
    max_r_neighbors <- max(tabulate(neighbors.id[neighbors.id > 0]))
    N <- max_r_neighbors + 1  # Equivalent to the old length(update_id) 
    S <- N.neighbors[update_id]
    A.n <- matrix(character(0), nrow = N, ncol = k)
    
    NL <- length(update_id)
    
    for (n in 1:NL) {
        
        if(S[n] >0){
            A.n[n, 1:S[n]] <- paste0("AD[", update_id[n], ",", 1:S[n], "]" )
            
        }
        
    }
    
    
    A.neighbors.t <-  unlist(apply(A.n, 1, function(x) x[!is.na(x)]))
    
    out <- list(update_id = update_id, Fneighbors.id = Fneighbors.id, AFnodes = AFnodes, ARnodes = ARnodes, A.neighbors = A.neighbors.t )
    
    return(out) 
}






## create the lists of calcNodes and copyNodes for use in MCMC samplers
mcmc_determineCalcAndCopyNodes <- function(model, target) {
    targetExpanded <- model$expandNodeNames(target)
    modelPredictiveNodes <- model$modelDef$maps$graphID_2_nodeName[model$predictiveNodeIDs]   ## identical to: model$getNodeNames(predictiveOnly = TRUE)
    targetExpandedPPbool <- targetExpanded %in% modelPredictiveNodes
    targetAllPP <- all(targetExpandedPPbool)
    targetAnyPP <- any(targetExpandedPPbool)
    ## if a particular sampler is assigned *jointly to PP and non-PP* nodes, then we're going to bail
    ## out and quit, if the option MCMCusePredictiveDependenciesInCalculations == FALSE.
    ## this is an extreme corner-case, which I think will lead to problems.
    if(targetAnyPP && !targetAllPP && !getNimbleOption('MCMCusePredictiveDependenciesInCalculations'))
        stop('cannot assign samplers jointly to posterior predictive (PP) nodes and non-PP nodes, when MCMCusePredictiveDependenciesInCalculations option is FALSE', call. = FALSE)
    ## if the sampler calling this, itself, is operating exclusively on posterior predictive nodes,
    ## then regardless of how the rest of the model is being sampled (w.r.t. inclusion of posterior predictive nodes),
    ## we'll include 'self' and all stochastic dependencies (the full markov blanket) in the calculations,
    ## which necessarily are taking place entirely within a posterior predictive network of nodes.
    ## this should lead to correct behaviour (consistent samples and joint posteriors) in all cases.
    if(targetAllPP) {
        ## when sampler is operating only on posterior predictive nodes,
        ## then always include all predictive dependencies:
        calcNodes <- model$getDependencies(target, includePredictive = TRUE)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE, includePredictive = TRUE)
        ##calcNodesPPomitted <- character()
        copyNodes <- model$getDependencies(target, self = FALSE)
    } else {
        ## usual case:
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        ##calcNodesPPomitted <- setdiff(model$getDependencies(target, includePredictive = TRUE), calcNodes)
        copyNodes <- calcNodesNoSelf
    }
    isStochCopyNodes <- model$isStoch(copyNodes)
    copyNodesDeterm <- copyNodes[!isStochCopyNodes]
    copyNodesStoch <- copyNodes[isStochCopyNodes]
    ##
    ccList <- list(
        calcNodes = calcNodes,
        calcNodesNoSelf = calcNodesNoSelf,
        ##calcNodesPPomitted = calcNodesPPomitted,
        copyNodesDeterm = copyNodesDeterm,
        copyNodesStoch = copyNodesStoch
    )
    return(ccList)
}





#=======================================================
# NNGP-Specific Random Walk Metropolis-Hastings Sampler
#=======================================================

#' NNGP-Specific Random Walk Metropolis-Hastings Sampler
#'
#' A custom NIMBLE sampler for updating individual spatial random effects under a 
#' Nearest-Neighbor Gaussian Process (NNGP) approximation. This sampler uses a 
#' factorized likelihood approach to perform highly efficient local Metropolis-Hastings sampling
#' with a normal proposal distribution (Metropolis, 1953), implementing the adaptation routine 
#' given in Shaby and Wells (2011).
#'
#' @param model (uncompiled) model on which the MCMC is to be run
#' @param mvSaved \code{modelValues} object to be used to store MCMC samples
#' @param target The node to be sampled (a single spatial random effect scalar).
#' @param control A list of control parameters:
#' \itemize{
#'   \item \code{AD}: Character (default \code{"AD"}). The name of the matrix containing the NNGP coefficients within the NIMBLE model object.
#'   \item \code{neighbors.id}: Character (default \code{"neighbors.id"}). The name of the forward neighbor matrix within the NIMBLE model object.
#'   \item \code{adaptive}: Logical (default \code{TRUE}). Whether to use an adaptive step-size procedure.
#'   \item \code{adaptInterval}: Integer (default \code{200}). Number of iterations between adaptive adjustments.
#'   \item \code{scale}: Numeric (default \code{1}). Initial scale/standard deviation for the random walk proposal.
#' }
#'
#' @details
#' Instead of calculating the full NNGP log-likelihood, which is \eqn{O(M)} 
#' (where \eqn{M} is the total number of spatial locations), this sampler exploits 
#' the Directed Acyclic Graph (DAG) structure of the NNGP to reduce the complexity 
#' of a single-node update to \eqn{O(k^2)}, where \eqn{k} is the number of neighbors. 
#' 
#' During the MCMC configuration stage, the sampler automatically queries the local 
#' graph architecture using \code{get_single_reverse_neighbors} and \code{RWNNGP_setup} 
#' to resolve dependencies internally. When a single node \eqn{w_i} is updated, 
#' only its own conditional density and the conditional densities of its "children" 
#' (the reverse neighbors that depend on it) are affected.
#'
#' The local log-Metropolis-Hastings ratio (\eqn{\log MHR}) is calculated as:
#' \deqn{
#' \log MHR = -\frac{(r_i^{*2} - r_i^2)}{2D_i} + \sum_{j \in \mathcal{R}(i)} -\frac{(r_j^{*2} - r_j^2)}{2D_j}
#' }
#' 
#' where:
#' \itemize{
#'   \item \eqn{r_i} and \eqn{r_i^*} are the current and proposed conditional residuals for the target node.
#'   \item \eqn{\mathcal{R}(i)} is the set of reverse neighbors (indices \eqn{j} such that \eqn{i \in N(j)}).
#'   \item \eqn{r_j} and \eqn{r_j^*} are the current and proposed residuals for those reverse neighbors.
#'   \item \eqn{D} represents the conditional variances.
#' }
#' 
#' The residuals are updated efficiently using the difference \eqn{\delta = w_i^* - w_i}, 
#' such that \eqn{r_j^* = r_j - A_{ji}\delta}, where \eqn{A_{ji}} is the NNGP 
#' regression coefficient.
#' 
#' @references
#' Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., and 
#' Teller, E. (1953). Equation of state calculations by fast computing machines. 
#' *The Journal of Chemical Physics*, 21(6), 1087-1092.
#' 
#' Shaby, B. A. and Wells, M. T. (2011). Exploring an adaptive Metropolis-Hastings 
#' algorithm. *Department of Statistical Science, Cornell University Tech Report*.
#'
#' @author Fabian Ketwaroo
#'
#' @export
sampler_RW_NN_GP <- nimbleFunction(
    name = 'sampler_RW_NN_GP',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        
        ## control list extraction
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scale               <- extractControlElement(control, 'scale',               1)
        
        # NNGP control
        AD                  <- extractControlElement(control, 'AD', "AD" )
        neighbors.id        <- extractControlElement(control, 'neighbors.id', "neighbors.id" )
        
        node_id              <- as.numeric(gsub("[^0-9]", "", target)) # The spatial location index of the target node
        N.neighbors <- apply(neighbors.id, 1, function(x){ sum(x != 0) }) # number of neighbors for each node
        k <- dim(neighbors.id)[2] # number of neighbors
        M <- dim(neighbors.id)[1] # number of spatial locations 
        Dnodes <- paste0(AD,"[", 1:M, ",", k+1, "]") # extract D nodes
        
        Rneighbors.id      <- get_single_reverse_neighbors(node_id, neighbors.id)
        nodes_ext <- RWNNGP_setup(node_id = node_id, AD =  AD, neighbors.id = neighbors.id, Rneighbors.id = Rneighbors.id, N.neighbors = N.neighbors, k= k) 
        Fneighbors.id       <- nodes_ext$Fneighbors.id
        id                  <- nodes_ext$update_id
        AFnodes             <- nodes_ext$AFnodes
        ARnodes             <- nodes_ext$ARnodes
        A.neighbors         <- nodes_ext$A.neighbors
        
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodes
        parentAsScalar <- model$expandNodeNames(model$expandNodeNames(target), returnScalarComponents = TRUE)
        
        
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory      <- c(0, 0)   ## scaleHistory
        acceptanceHistory <- c(0, 0)   ## scaleHistory
        saveMCMChistory <- getNimbleOption('MCMCsaveHistory')
        optimalAR     <- 0.44
        gamma1        <- 0
        ## checks
        if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
        if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
    },
    run = function() {
        
        
        
        D <- values(model, Dnodes)
        
        # Get residuals 
        z.id <- values(model, parentAsScalar)
        
        #id = update_id
        N <- length(id)
        r <- numeric(N)
        S <- N.neighbors[id]
        
        A.N <- values(model, A.neighbors)
        
        
        for (i in 1:N) {
            
            if( S[i] == 0) r[i] <- z.id[id[i]]
            else {
                
                if (N == 1) {
                    r[i] <- z.id[id[i]] - sum(A.N * z.id[neighbors.id[id[i], 1:S[i]]])
                } else {
                    if (i == 1) {
                        r[i] <- z.id[id[i]] - sum(A.N[1:S[1]] * z.id[neighbors.id[id[i], 1:S[i]]])
                    } else {
                        r[i] <- z.id[id[i]] - sum(A.N[(sum(S[1:(i-1)]) + 1):(sum(S[1:i]))] * z.id[neighbors.id[id[i], 1:S[i]]])
                    }
                }
                
            }
            
            
        }
        
        
        
        # Get current value 
        currentValue <- model[[target]]
        
        # Propose new value 
        propValue <- rnorm(1, mean = currentValue,  sd = scale)
        
        # Add propose value to nimble model memory
        model[[target]] <<- propValue
        
        # Difference between proposal and current value 
        delta  <-  propValue -  currentValue 
        
        # forward neighbors residuals 
        Fn <- length(Fneighbors.id)
        
        if( Fn >0 ){
            
            AF <- values(model, AFnodes)  # residual for proposal i from forward neighbors
            
            z.F <- values(model, parentAsScalar[Fneighbors.id] )
            
            r_i_star <-  propValue - sum(AF* z.F )
            
        } else {
            r_i_star <-  propValue
        }
        
        
        # residual for proposal i from reverse neighbors
        J <- length(Rneighbors.id)
        
        
        if(J > 0) {
            
            AR <- values(model, ARnodes) # reverse neighbors residuals of node i
            r_j_star <- numeric(J)
            
            for (j in 1:J) {
                r.id <-  Rneighbors.id[j]
                r_j_star[j] <- r[j+1] - (AR[j]*delta)
            }
            
            # MH ratio: additon of node i and its reverse neighbors
            logMHR <- -(r_i_star^2 - r[1]^2)/(2*D[node_id]) +  sum( -( r_j_star[1:J]^2 - r[2:(J+1)]^2 )/ (2*D[Rneighbors.id[1:J]])  )
            
            
        } else {
            
            # when the node doesn't have any reverse neighbors, the MH ratio just depends on node i
            logMHR <- -(r_i_star^2 - r[1]^2)/(2*D[node_id])
            
        }
        
        
        if(logMHR == -Inf) {
            jump <- FALSE
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
        } else {
            logMHR <- logMHR + model$calculateDiff(calcNodesNoSelf) 
            jump <- decide(logMHR)
            
            if( node_id == (M) ) {
                
                if(jump) {
                    
                    model$calculate(target)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
                } else {
                    
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = FALSE) # back to the current state 
                    
                    # Update log likeihood so it can be saved in NIMBLE memory
                    model$calculate(target) # Loglikeigood at current state
                    nimCopy(from = model, to = mvSaved , row = 1, nodes = target, logProbOnly = TRUE) #  logProb is updated regardless whether node M is accepted or rejected. Only logProb changed not the value
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
                }
                
            } else {
                
                # for the other nodes
                
                if(jump) {
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = FALSE)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
                } else {
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = FALSE)
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
                }
                
            }
            
            
        }
        
        
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(scaleHistory, timesAdapted)                 ## scaleHistory
                    scaleHistory[timesAdapted] <<- scale                ## scaleHistory
                    setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
                    acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
                }
                gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                ## If there are upper and lower bounds, enforce a maximum scale of
                ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
                ## Otherwise, for a poorly-informed posterior,
                ## the scale could grow without bound to try to reduce
                ## acceptance probability.  This creates enormous cost of
                ## reflections.
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        setScale = function(newScale = double()) {
            scale         <<- newScale
            scaleOriginal <<- newScale
        },
        getScaleHistory = function() {       ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(scaleHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
                return(numeric(1, 0))
            }
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(acceptanceHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
                return(numeric(1, 0))
            }
        },
        ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            if(saveMCMChistory) {
                scaleHistory  <<- c(0, 0)    ## scaleHistory
                acceptanceHistory  <<- c(0, 0)
            }
            gamma1 <<- 0
        }
    )
)
