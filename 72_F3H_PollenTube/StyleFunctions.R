#install.packages("stringr")
library(stringr)
# install.packages("readxl")
library(readxl)

getSampleExcelSpreadsheet = function(fname="Documentation/muday-144_sample_sheet.xlsx") {
  sample_sheet = read_excel(fname)
  return(sample_sheet)
}

getColorsForSampleCodes = function(sample_names,
                                   sample_sheet=NULL) {
  if (is.null(sample_sheet)) {
    sample_sheet = getSampleExcelSpreadsheet()
  }
  sample_codes = sample_sheet$`Sample Code` 
  sample_colors = sample_sheet$Color
  names(sample_colors) = sample_codes
  to_return = sample_colors[sample_names]
  return(to_return)
}

getSampleGroups = function(sample_names) {
  sample_groups=sub("\\.[789]$",names(sample_names),replacement="")
  return(sample_groups)
}

# example input: "F.28.15.7"
getGenotype = function(sample_name) {
  value=strsplit2(sample_name,"\\.")[1,1]
  return(value)
}

getTreatment=function(sample_name) {
  value=strsplit2(sample_name,"\\.")[1,2]
  ifelse(value=='28',"control","heat stress")
}

getTimePoint=function(sample_name) {
  value = strsplit2(sample_name,"\\.")[1,3]
  return(value)
}

getReplicate=function(sample_name) {
  return (strsplit2(sample_name,"\\.")[1,4])
}

