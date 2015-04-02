library(reshape2)
library(plyr)
library(dplyr)

## Extract a sample from an MCMC chain for a variable with the given named indices into a long-format data frame.
## For example, imagine a variable b[i,v] with i in [1..100] and v in [1..3]. An MCMC sample returned from JAGS 
## (for example) would have columns with names like "b[1,1]", "b[2,1]", etc. 
##
## extract_samples(mcmcChain, ~ b[i,v]) would return a data frame with:
##		column ".sample": value in [1..nrow(mcmcChain)]
## 		column "i":		  value in [1..100]
##		column "v":		  value in [1..3]
##      column "b":       value of "b[i,v]" for sample number ".sample" in mcmcChain 
##
## The shorthand ".." can be used to specify one column that should be put into a wide format. For example:
##
## extract_samples(mcmcChain, ~ b[i,..]) would return a data frame with:
##		column ".sample": value in [1..nrow(mcmcChain)]
## 		column "i":		  value in [1..100]
##      column "b1":       value of "b[i,1]" for sample number ".sample" in mcmcChain 
##      column "b2":       value of "b[i,2]" for sample number ".sample" in mcmcChain 
##      column "b3":       value of "b[i,3]" for sample number ".sample" in mcmcChain 
##
## prototypes optionally specifies a list or data.frame where each variable with the same name
## as an index is a prototype for that index -- i.e., its type is the expected type of that index. 
## Those types are used to translate indices back into useful variables (usually factors). The most
## common use of this is to automatically translate indices that correspond to levels of a factor
## back into levels of a factor. Names in prototypes that are not in index_types are ignored.
##
## The simplest use of this is to pass in the data frame from which the original data came. Supported
## data types are factor, ordered, and logical. For example:
##
## 		if data_types$v is a factor, the v column is translated into a factor using factor(v, labels=levels(index_types$v), ordered=is.ordered(index_types$v))
## 		if data_types$v is logical, the v column is translated into a logical using as.logical(v)
##
extract_samples = function(mcmcChain, variable_spec, prototypes=NULL) {
    #first, extract the sample into a data frame
    variable_name = as.character(variable_spec[[2]][[2]])
    index_names = as.character(variable_spec[[2]][-1:-2])
    samples = extract_samples_long_(mcmcChain,
        variable_name, 
        index_names, 
        prototypes)
    
    if (is.null(samples$..)) {
        samples
    } else {
        #if any of the columns were named "..", use it to form a wide version of the data
        samples %>%
            mutate(.. = factor(.., labels=variable_name)) %>%
            dcast(... ~ .., value.var=variable_name)
    }   
}


## Extract a sample from an MCMC chain for a variable with the given named indices into a long-format data frame.
## For example, imagine a variable b[i,v] with i in [1..100] and v in [1..3]. An MCMC sample returned from JAGS 
## (for example) would have columns with names like "b[1,1]", "b[2,1]", etc. 
##
## extract_samples_long_(mcmcChain, "b", c("i", "v")) would return a data frame with:
##		column ".sample": values in [1..nrow(mcmcChain)]
## 		column "i":		  values in [1..100]
##		column "v":		  values in [1..3]
##      column "b":       sample number ".sample" of "b[i,v]" in mcmcChain 
##
## prototypes optionally specifies a list or data.frame where each variable with the same name
## as an index is a prototype for that index -- i.e., its type is the expected type of that index. 
## Those types are used to translate indices back into useful variables (usually factors). Names 
## in prototypes that are not in index_types are ignored.
## The simplest use of this is to pass in the data frame from which the original data came. Supported
## data types are factor, ordered, and logical.
##		if data_types$v is a factor, the v column is translated into a factor using factor(v, labels=levels(index_types$v), ordered=is.ordered(index_types$v))
##		if data_types$v is logical, the v column is translated into a logical using as.logical(v)
extract_samples_long_ = function(mcmcChain, variable_name, index_names, prototypes=NULL) {
    #for each column name in the chain...
    samples = ldply(colnames(mcmcChain), function(colname) {
            #parse it into variable name and indices
            colname_parts = strsplit(colname,"\\,|\\]|\\[")[[1]]
            if (colname_parts[1] == variable_name) {	#this is the variable we want
                #get the values of the indices 
                indices = as.list(as.numeric(colname_parts[-1]))
                names(indices) = index_names
                #get the values of this variable in each sample
                values = list(mcmcChain[,colname])
                names(values) = variable_name
                #put it all together
                data.frame(.sample=1:nrow(mcmcChain), indices, values)
            }
        })
    
    #convert data types back into usable forms
    for (index_name in index_names) {
        if (index_name %in% names(prototypes)) {
            #we have a data type for this index, convert it
            prototype = prototypes[[index_name]]
            if (is.factor(prototype)) {
                samples[[index_name]] = factor(samples[[index_name]], 
                    labels=levels(prototype), ordered=is.ordered(prototype))
            }
            else if (is.logical(prototype)) {
                samples[[index_name]] = as.logical(samples[[index_name]])
            }
        }
    }
    
    samples
}
