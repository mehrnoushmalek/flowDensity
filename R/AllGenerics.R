### Generics ###
if(!isGeneric("flowDensity")){
  setGeneric("flowDensity", function(obj,
                                     channels,
                                     position,
				     singlet.gate,
                                     ...)
    standardGeneric("flowDensity")
  )
}

if(!isGeneric("getflowFrame")){
  setGeneric("getflowFrame", function(obj) standardGeneric("getflowFrame"))
}
