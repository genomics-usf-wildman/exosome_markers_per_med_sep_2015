library(data.table)

args <- c("chosen_samples",
          "GTEx_Data_V4_Annotations_SampleAttributesDS.txt",
          "categorized_samples")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
gtex.samples <- fread(args[2])

categorized.samples <- copy(chosen.samples[,list(Sample_Name,SRX)])

gtex.samples <- gtex.samples[,list(SMTSD,SAMPID)]
setnames(gtex.samples,c("SAMPID","SMTSD"),c("SRX","Sample_Name"))
categorized.samples <-
    rbind(categorized.samples,
          gtex.samples)

categorized.samples[,Sample_Group:=Sample_Name]
categories <-
    list(
        ## We're going to ignore differences between fetal days and
        ## sex for gene expression
        combine.fetal.days.and.sex=
             list(match="^(.+),\\s*fetal\\s*(?:day(\\d+))?(?:\\s*[MFU])?$",
                  replace="fetal \\1"),
        ## remove right and left
        ignore.right.left=
            list(match=", (?:right|left)$",
                 replace=""),
        ## remove (heart),
        remove.heart.comma=
            list(match="(?:heart), ",
                 replace=""),
        ## remove (cells|brain) -
        remove.cells.dash=
            list(match="(?:Cells|Brain|Breast|Cervix|Atery|Heart) - ",
                 replace=""),
        swap.dash=
            list(match="^([A-Z]\\w+) - (.+)$",
                 replace="\\1 \\2"),
        ## the one placenta from day 113 we'll just call placenta
        placenta.day.113=
            list(match="^placenta, day113 F$",
                 replace="placenta"),
        ## Uterus term labor and term no labor; combining them as
        ## uterus for the time being
        uterus.term.no.labor=
            list(match="^TNL_rep\\d+$",
                 replace="uterus"),
        uterus.term.labor=
            list(match="^TL_rep\\d+$",
                 replace="uterus"),
        placenta=
            list(match="^[AESX]\\d+_[abc]\\d?(?:\\.\\d+)?$",
                 replace="placenta")
        )
for (category in names(categories)) {
    categorized.samples[grepl(categories[[category]]$match,Sample_Group),
                        Sample_Group:=
                            gsub(categories[[category]]$match,
                                 categories[[category]]$replace,
                                 Sample_Group)]
}

categorized.samples[,Sample_Group:=tolower(Sample_Group)]
         
save(categorized.samples,
     file=args[length(args)])

