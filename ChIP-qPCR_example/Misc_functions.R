#####################################
# Misc helper functions

#############################
# add leading "characters"0" to string
leading.zeros <- function(str, width = 2)
{
  return(out.str = formatC(x = str, width = width, format = "d", flag = "0"))
}

#############################
# replace "." with "_"
replace.with_ <- function(x)
{
  return(gsub(pattern = '\\.', replacement = '_', x = x))
}

#############################
# convert inch to cm
inch2cm <- function(x, inverse = FALSE)
{
  if (inverse == FALSE)
  {
    return(x*2.54)
  }
  else
  {
    return(x/2.54)
  }
}

#############################
# embed fonts into PDF using ghostscript
GSEmbedFonts <- function(input, output, RemoveInputFile = FALSE){
  success <- system2(command = "gs", args = paste0("-q -dNOPAUSE -dBATCH -dPDFSETTINGS=/prepress -sDEVICE=pdfwrite -sOutputFile=", output, " ", input))
  if(RemoveInputFile & success == 0){# success == 0 only if system2 command ran successfully
    file.remove(input)
  }
}


#############################
# convert factor to character in data.frame
factor.to.character <- function(df){
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  return(df)
}

# alias
factorsAsStrings <- function(df){
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  return(df)
}

charactersToNumbers <- function(df){
  i <- sapply(df, is.character)
  df[i] <- lapply(df[i], as.numeric)
  return(df)
}

