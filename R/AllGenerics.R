### Generics ###

if(!isGeneric("flowDensity")){
  setGeneric("flowDensity", function(obj,
                                     channels,
                                     position,
                                     ...)
    standardGeneric("flowDensity")
  )
}
if(!isGeneric("getflowFrame")){
  setGeneric("getflowFrame", function(obj) standardGeneric("getflowFrame"))
}
