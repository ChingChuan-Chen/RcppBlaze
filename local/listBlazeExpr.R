## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is used to generate needed wrappers for Blaze Expressions
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the 3-Clause BSD License. You should have received
## a copy of 3-Clause BSD License along with RcppBlaze.
## If not, see https://opensource.org/license/BSD-3-Clause.

library(stringr)

blazeExprHeaderFiles <- list.files(
  "inst/include/blaze/math/expressions",
  "^T?[SD](Vec|Mat)([SD](Vec|Mat)|Scalar|Map)\\w*Expr.h",
  full.names = TRUE
)
blazeExprHeaderFiles <- blazeExprHeaderFiles[!str_detect(blazeExprHeaderFiles, "(Inner|Equal)Expr.h$")]

processedExprDefList <- lapply(blazeExprHeaderFiles, function(f){
  templateClassDefs <- str_extract_all(
    paste(readLines(f), collapse = "\n"),
    "template<\\s*[/\\s\\w-\n,]+>[/\\s\\w-\n]+class\\s*T?[SD](Vec|Mat)([SD](Vec|Mat)|Scalar|Map)\\w*Expr"
  )[[1]]
  output <- sapply(templateClassDefs, function(x) sapply(str_split(x, "\n"), function(y) {
    replaceCommentsStartingSpaces <- str_replace(str_replace(y, "^\\s+", ""), "\\s+//[/\\s\\w-]+$", "")
    replaceQuotes <- str_replace(str_replace(replaceCommentsStartingSpaces, "<\\s+", " <"), "\\s+>", "> ")
    str_replace(paste(replaceQuotes, collapse=""), "class\\s+", "class blaze::")
  }), USE.NAMES = FALSE)
  return(output[str_detect(output, str_extract(f, "/([^./]+).h$", 1))])
})

if (any(sapply(processedExprDefList, length) == 0L)) {
  cat("Files not getting blaze Expression definition:\n\t")
  cat(paste(blazeExprHeaderFiles[sapply(processedExprDefList, length) == 0L], collapse="\n\t"), "\n")
  stop("There are files are not getting blaze Expressions")
}

outputList <- lapply(unlist(processedExprDefList), function(x) {
  replaceWrapFunc <- str_replace(x, "class", "SEXP wrap(const")
  templatePars <- str_replace(str_extract_all(replaceWrapFunc, "\\w+[,>]")[[1]], "[,>]", "")

  list(
    forwradDef = sprintf("%s<%s>&);", replaceWrapFunc, paste(templatePars, collapse = ", ")),
    wrapFunc = sprintf("%s<%s>& x) {\n  return RcppBlaze::blaze_wrap(x);\n};\n", replaceWrapFunc, paste(templatePars, collapse = ", "))
  )
})

writeLines(paste(sapply(outputList, `[[`, 1), collapse = "\n"), "local/RcppBlazeForwardExprPart.h")
writeLines(paste(sapply(outputList, `[[`, 2), collapse = "\n"), "local/RcppBlazeWrapExprPart.h")
