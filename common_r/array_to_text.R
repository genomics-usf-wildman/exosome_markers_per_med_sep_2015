### $(R) $(ROPTS) -f array_to_text.R --args data_file_name array_name array_file.txt
library("data.table")

arguments <- commandArgs(trailingOnly=TRUE)
data.file.name <- arguments[1]
data.array <- arguments[2]
text.file.name <- arguments[3]

load(data.file.name)
eval(parse(text=c("fwrite(",data.array,",",
             "file=text.file.name,",
             'sep="\\t")')
           )
     )
