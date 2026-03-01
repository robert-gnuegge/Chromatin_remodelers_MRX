# info --------------------------------------------------------------------
# purpose: miscellaneous helper functions
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 12/16/21
# version: 0.1


# convert dash to underscore ----------------------------------------------
# argument: string
# result: string
dash_to_underscore <- function(string){
  gsub(pattern = "-", replacement = "_", x = string)
}


# embed fonts in PDF (using ghostscript) ----------------------------------
# argument: PDF input and output file names (strings)
# result: none
GS_embed_fonts <- function(input, output, remove_input_file = TRUE){
  success <- system2(command = "gs", args = paste0("-q -dNOPAUSE -dBATCH -dPDFSETTINGS=/prepress -sDEVICE=pdfwrite -sOutputFile=", output, " ", input))
  if(remove_input_file & success == 0){  # success == 0 only if system2 command ran successfully
    file.remove(input)
  }
}