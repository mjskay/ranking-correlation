library(plyr)
library(dplyr)
library(tidyr)


#basic data list
data_list = function(...) {
    x = list(...)
    class(x) = c("data_list", "list")
    x
}
c.data_list = function(x, ..., recursive=FALSE) {
    class(x) = class(x)[-which(class(x) == "data_list")]
    x = c(x, ..., recursive=recursive)
    class(x) = c("data_list", class(x))
    x
}
print.data_list = function(x, ...) {
    cat("data_list:\n\n")
    class(x) = class(x)[-which(class(x) == "data_list")]
    print(x, ...) 
}

#conversion into data list from basic data types
as.data_list = function(object, name="", ...) UseMethod("as.data_list")
as.data_list.default = function(object, name="", ...) {
    warning(deparse(name), " has unsupported type ", deparse(class(object)), " and was dropped.")
    data_list()
}
as.data_list.numeric = function(object, name="", ...) {
    data = data_list(object)
    if (name == "") {	#name unspecified, warn
        warning("No name provided for value ", deparse(object, nlines=1))
    }
    names(data) = name
    data
}
as.data_list.logical = function(object, name="", ...) {
    as.data_list(as.numeric(object), name, ...)
}
as.data_list.factor = function(object, name="", .n_name=function(name) paste0("n_", name), ...) {
    data = as.data_list(as.numeric(object), name, .n_name=.n_name, ...)
    if (any(table(object) == 0)) {
        warning("Some levels of factor ", deparse(name), " are unused. This may cause issues if you are using it as an index in a model.")
    }
    data[[.n_name(name)]] = length(levels(object))
    data
}
as.data_list.list = function(object, name="", ...) {
    #go through list and translate variables
    data = data_list()
    for (i in 1:length(object)) {
        data = c(data, as.data_list(object[[i]], names(object)[[i]], ...))
    }
    data
}
as.data_list.data.frame = function(object, name="", .n_name=function(name) paste0("n_", name), ...) {
    #first, translate all variables in the data frame
    data = as.data_list.list(object, name, .n_name, ...)
    #then add "n" column and return final list
    n_name = if (name == "") "n" else .n_name(name)
    data[[n_name]] = nrow(object)
    data
}
as.data_list.data_list = function(object, name="", ...) {
    object
}

## Compose data into a list suitable to be passed into an MCMC sampler (JAGS, BUGS, etc).
##
## Translates each argument into list elements using as.data_list, and then concatenates
## all resulting list elements together.
## Translates a data.frame into a list suitable for use in an MCMC sampler. Does this as follows:
##
## Translates elements as follows:
##
## 		- numerics are included as-is
##		- logicals are translated into numeric using as.numeric
##		- factors are translated into numeric using as.numeric, and an additional
##		  column named .n_name(argument_name) is added with the number of levels in the factor.
##		- lists are translated by translating all elements of the list (recursively)
##		  and adding them to the result.
##		- data.frames are translated by translating every column of the data.frame and
##		  adding them to the result.
##		  A variable named "n" (or .n_name(argument_name) if the data.frame is passed as
##		  a named argument argument_name) is also added with the number of rows in the
##		  data frame.
##		- all other types are dropped (and a warning given)
##
## If you wish to add support for additional types not described above, provide an implementation
## of as.data_list for the type. See as.data_list.numeric, as.data_list.logical, etc for examples.
##
compose_data = function(..., .n_name=function(name) paste0("n_", name)) {
    #translate arguments into a data_list
    objects = list(...)
    if (is.null(names(objects))) {
        #when no named arguments are supplied, we must translate NULL into empty string
        #so that the name argument to as.data_list is still valid
        names(objects) = rep("", length(objects))
    }
    as.data_list(objects, .n_name)
}

## Extract a sample from an MCMC chain for a variable with the given named indices into a long-format data frame.
## For example, imagine a variable b[i,v] with i in [1..100] and v in [1..3]. An MCMC sample returned from JAGS 
## (for example) would have columns with names like "b[1,1]", "b[2,1]", etc. 
##
## extract_samples(mcmcChain, ~ b[i,v]) would return a data frame with:
##		column ".sample": value in [1..nrow(mcmcChain)]
## 		column "i":		  value in [1..20]
##		column "v":		  value in [1..3]
##      column "b":       value of "b[i,v]" for sample number ".sample" in mcmcChain 
##
## The shorthand ".." can be used to specify one column that should be put into a wide format. For example:
##
## extract_samples(mcmcChain, ~ b[i,..]) would return a data frame with:
##		column ".sample": value in [1..nrow(mcmcChain)]
## 		column "i":		  value in [1..20]
##      column "b1":       value of "b[i,1]" for sample number ".sample" in mcmcChain 
##      column "b2":       value of "b[i,2]" for sample number ".sample" in mcmcChain 
##      column "b3":       value of "b[i,3]" for sample number ".sample" in mcmcChain 
##
## prototypes optionally specifies a list or data.frame. Each entry in prototypes with the same name
## as the variable or an index in varible_spec is a used as a prototype for that variable or index -- 
## i.e., its type is taken to be the expected type of that variable or index. 
## Those types are used to translate numeric values of variables back into useful values (usually levels 
## of a factors). 

## The most common use of prototypes is to automatically translate indices that correspond to levels of a factor
## in the original data back into levels of that factor. Names in prototypes that are not in 
## found in variable_spec are ignored.
##
## The usual use of prototypes is to pass in the data frame from which the original data came. 
## Supported types of prototypes are factor, ordered, and logical. For example:
##
## 		if data_types$v is a factor, the v column in the returned samples is translated into 
##      a factor using factor(v, labels=levels(index_types$v), ordered=is.ordered(index_types$v))
##
## 		if data_types$v is a logical, the v column is translated into a logical using as.logical(v)
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
    }
    else {
        #if any column is named "..", use it to form a wide version of the data
        samples %>%
            mutate(.. = factor(.., labels=variable_name)) %>%
            spread_("..", variable_name)
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
    samples = ldply(colnames(mcmcChain), function(colname) {
        #parse column into variable name and indices
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
    for (column_name in c(variable_name, index_names)) {
        if (column_name %in% names(prototypes)) {
            #we have a data type for this index, convert it
            prototype = prototypes[[column_name]]
            if (is.factor(prototype)) {
                samples[[column_name]] = factor(samples[[column_name]], 
                    labels=levels(prototype), ordered=is.ordered(prototype))
            }
            else if (is.logical(prototype)) {
                samples[[column_name]] = as.logical(samples[[column_name]])
            }
        }
    }
    
    samples
}
