#' Automatic Differential Expression Analysis using DESeq2
#'
#' Obtain DESeq2 Normalized results
#'
#' @param sampleTable Filename indicating the sampleTable tab delimited file pointing to HTSeq-count files.
#' @param countTable Filename indicating the raw count table tab delimited file (use for featureCounts output and with colData)
#' @param colData Filename indicating the tab delimited file containing metadata for countTable
#' @param expFilt Numeric indicating the minimum average gene counts to consider a gene as being 'expressed'. If whole number, average gene count across all samples must be greater than this number. Default=0.
#' @param retExplore Boolean indicating if BinfTools::exploreData() function should also be returned. Requires gene symbols, and BinfTools/gpGeneSets packages. Default=F
#' @param retDDS Boolean indicating if the DESeq2 dataset (dds) should be returned for further usage. Defaults to FALSE.
#' @param detID Boolean indicating if ENSEMBL or Flybase gene IDs should be automatically detected for conversion to symbos using BinfTools::getSym(). Defatul=T
#' @return List with DESeq2 analyzed data: results, normalized counts, and conditions
#' @export

autoDE<-function(sampleTable=NULL, countTable=NULL, colData=NULL, expFilt=0, retExplore=F, retDDS=F, detID=T){
  dds<-NULL
  condition<-NULL
  conName<-NULL #contrast name
  if(!is.null(sampleTable)){
    message("sampleTable of HTSeq-count files provided.")
    st<-sampleTable
    if(is.character(sampleTable)){
      message("sampleTable was filename. Reading sampleTable as tab-delimited file.")
      st<-read.delim(sampleTable)
      st[,2]<-as.character(st[,2])
    }
    cols<-colnames(st)[c(3:ncol(st))]
    message(length(cols)," metadata columns for ",nrow(st)," samples found:")
    print(cols)
    if(length(grep("Batch", cols))>0){
      message("Metadata column named 'Batch' found. Using for batch correction in design.")
      cols<-cols[-which(cols %in% "Batch")]
      if(length(cols)>1){
        #With the column 'Batch', we still have more than 2 columns, meaning multivariate analysis
        message("Remaining metadata columns to be combined for multivariate analysis")
        x<-cols[1]
        st$combined<-st[,`x`]
        for(i in 2:length(cols)){
          x<-cols[i]
          st$combined<-paste(st$combined, st[,`x`], sep=".")
        }
        dds<-DESeq2::DESeqDataSetFromHTSeqCount(sampleTable=st, design = ~ Batch + combined)
        condition<-as.character(dds$combined)
        conName<-"combined"
      } else {
        message("One remaining metdata column will be used as design")
        st$condition<-st[,`cols`]
        #print(head(st))
        dds<-DESeq2::DESeqDataSetFromHTSeqCount(sampleTable=st, design = ~ Batch + condition)
        condition<-as.character(dds$condition)
        conName<-"condition"
      }
    } else {
      message("No column named Batch. No batch correction will be performed.")
      if(length(cols)>1){
        #With the column 'Batch', we still have more than 2 columns, meaning multivariate analysis
        message("Multiple metadata columns to be combined for multivariate analysis")
        x<-cols[1]
        st$combined<-st[,`x`]
        for(i in 2:length(cols)){
          x<-cols[i]
          st$combined<-paste(st$combined, st[,`x`], sep=".")
        }
        dds<-DESeq2::DESeqDataSetFromHTSeqCount(sampleTable=st, design = ~ combined)
        condition<-as.character(dds$combined)
        conName<-"combined"
      } else {
        message("One metdata column will be used as design")
        st$condition<-st[,`cols`]
        dds<-DESeq2::DESeqDataSetFromHTSeqCount(sampleTable=st, design = ~ condition)
        condition<-as.character(dds$condition)
        conName<-"condition"
      }
    }
  }
  if(!is.null(countTable) & !is.null(colData)){
    message("countTable and colData provided.")
    ct<-countTable
    cd<-colData
    if(is.character(countTable)){
      message("countTable is filename. Reading countTable from tab-delimited file.")
      ct<-read.delim(countTable)
    }
    ct<-as.matrix(ct)
    if(is.character(colData)){
      message("colData is filename. Reading colData from tab-delimited file.")
      cd<-read.delim(colData, row.names=1)
    }

    cols<-colnames(cd)[c(1:ncol(cd))]
    message(length(cols)," metadata columns for ",nrow(cd)," samples found:")
    print(cols)
    if(length(grep("Batch", cols))>0){
      message("Metadata column named 'Batch' found. Using for batch correction in design.")
      cols<-cols[-which(cols %in% "Batch")]
      if(length(cols)>1){
        #With the column 'Batch', we still have more than 2 columns, meaning multivariate analysis
        message("Remaining metadata columns to be combined for multivariate analysis")
        x<-cols[1]
        cd$combined<-cd[,`x`]
        for(i in 2:length(cols)){
          x<-cols[i]
          cd$combined<-paste(cd$combined, cd[,`x`], sep=".")
        }
        dds<-DESeq2::DESeqDataSetFromMatrix(countData=ct, colData=cd, design = ~ Batch + combined)
        condition<-as.character(dds$combined)
        conName<-"combined"
      } else {
        message("One remaining metdata column will be used as design")
        cd$condition<-cd[,`cols`]
        dds<-DESeq2::DESeqDataSetFromMatrix(countData=ct, colData=cd, design = ~ Batch + condition)
        condition<-as.character(dds$condition)
        conName<-"condition"
      }
    } else {
      message("No column named Batch. No batch correction will be performed.")
      if(length(cols)>1){
        #With the column 'Batch', we still have more than 2 columns, meaning multivariate analysis
        message("Multiple metadata columns to be combined for multivariate analysis")
        x<-cols[1]
        cd$combined<-cd[,`x`]
        for(i in 2:length(cols)){
          x<-cols[i]
          cd$combined<-paste(cd$combined, cd[,`x`], sep=".")
        }
        dds<-DESeq2::DESeqDataSetFromMatrix(countData=ct, colData=cd, design = ~ combined)
        condition<-as.character(dds$combined)
        conName<-"combined"
      } else {
        message("One metdata column will be used as design")
        cd$condition<-cd[,`cols`]
        dds<-DESeq2::DESeqDataSetFromMatrix(countData=ct, colData=cd, design = ~ condition)
        condition<-as.character(dds$condition)
        conName<-"condition"
      }
    }
  }
  message("DESeq2 dataset built!")
  message("Filtering to keep genes with >",expFilt," reads in any sample.")
  dds<-dds[rowMeans(DESeq2::counts(dds))>expFilt,]

  message("Running DESeq2.")
  dds<-DESeq2::DESeq(dds)
  message("Performing PCA.")
  vsd<-DESeq2::varianceStabilizingTransformation(dds)
  print(DESeq2::plotPCA(vsd, intgroup="condition"))
  if(length(grep("Batch", colnames(SummarizedExperiment::colData(dds))))>0){
    vsd$Batch<-factor(vsd$Batch)
    print(DESeq2::plotPCA(vsd, intgroup="Batch"))
    message("Removing batch effect for PCA...")
    SummarizedExperiment::assay(vsd)<-limma::removeBatchEffect(SummarizedExperiment::assay(vsd), vsd$Batch)
    print(DESeq2::plotPCA(vsd, intgroup="condition"))
    print(DESeq2::plotPCA(vsd, intgroup="Batch"))
  }
  res<-NULL
  message(length(unique(condition)), " conditions identified.")
  if(length(unique(condition))==2){
    message("Only one comparison to be performed.")
    message("Assuming first named condition is reference condition: ", unique(condition)[1])
    res<-as.data.frame(DESeq2::results(dds, contrast=c(conName, unique(condition)[2], unique(condition)[1])))
  } else {
    message("Multiple comparisons can be performed.")
    message("Performing pairwise comparisons.")
    res<-list()
    index<-1
    for(i in 1:length(unique(condition))){
      for(k in (i+1):length(unique(condition))){
        if(i<length(unique(condition))){
          compName<-paste(unique(condition)[k], unique(condition)[i], sep="v")
          message("Making results for: ", compName)
          #message("Storing in res[[",index,"]]")
          res[[index]]<-as.data.frame(DESeq2::results(dds, contrast=c(conName, unique(condition)[k], unique(condition)[i])))
          names(res)[index]<-compName
          index<-index+1
        }
      }
    }
  }
  normCounts<-as.data.frame(DESeq2::counts(dds, normalized=T))

  if(isTRUE(retExplore)){
    countList<-NULL
    condList<-NULL
    if(class(res) =="list"){
      countList<-list()
      condList<-list()
      for(i in 1:length(res)){
        countList[[i]]<-normCounts
        condList[[i]]<-condition
        names(countList)[i]<-names(res)[i]
        names(condList)[i]<-names(res)[i]
      }
    } else {
      countList<-normCounts
      condList<-condition
    }
    ID<-FALSE
    type<-NULL
    if(length(grep("ENS", rownames(normCounts)))>0){
      ID<-TRUE
      type<-"Ensembl"
    }
    if(length(grep("FBgn", rownames(normCounts)))>0){
      ID<-TRUE
      type<-"Flybase"
    }
    if(isTRUE(ID)){
      message(type," Gene IDs detected and must be converted to gene symbols for exploreData.")
      opt<-NULL
      if(isTRUE(detID)){
        auto<-F
        if(length(grep("ENSG", rownames(normCounts)))>0){
          message("Human Ensembl IDs detected...")
          opt<-1
          auto<-T
        }
        if(length(grep("ENSMUSG", rownames(normCounts)))>0){
          message("Mouse Ensembl IDs detected...")
          opt<-2
          auot<-T
        }
        if(length(grep("FBgn", rownames(normCounts)))>0){
          message("Flybase gene IDs detected...")
          opt<-3
          auto<-T
        }
        if(isFALSE(auto)){
          message("Cannot detect species automatically...")
          detID<-F
        }
      }
      if(isFALSE(detID)){
        message("Please indicate the species of your samples:")
        print(paste(c(1:3),c("Human","Mouse","Drosophila")))
        opt<-as.numeric(readline("Enter your selection: "))
      }
      options<-c("hsapiens","mmusculus","dmelanogaster")
      targets<-c("HGNC","MGI","FLYBASENAME_GENE")
      res2<-NULL
      if(class(res) =="list"){
        res2<-lapply(res, BinfTools::getSym, obType="res", species=options[opt], target=targets[opt])
        countList<-lapply(countList, BinfTools::getSym, obType="counts", species=options[opt], target=targets[opt])
        res<-lapply(res, BinfTools::getSym, obType="res", species=options[opt], target=targets[opt], addCol=T)
      } else {
        res2<-BinfTools::getSym(res, obType="res", species=options[opt], target=targets[opt])
        countList<-BinfTools::getSym(countList, obType="counts", species=options[opt], target=targets[opt])
        res<-BinfTools::getSym(res, obType="res", species=options[opt], target=targets[opt], addCol=T)
      }
      normCounts<-BinfTools::getSym(normCounts, obType="counts", species=options[opt], target=targets[opt], addCol=T)
    }
    explore<-BinfTools::exploreData(res=res2, counts=countList, cond=condList)
    if(isFALSE(retDDS)){
    return(list(normCounts=normCounts, res=res, condition=condition, explore=explore))
  } else {
    return(list(normCounts=normCounts, res=res, condition=condition, explore=explore, dds=dds))
  }
  } else {
    if(isTRUE(detID)){
      opt<-NULL
      auto<-F
      if(length(grep("ENSG", rownames(normCounts)))>0){
        message("Human Ensembl IDs detected...")
        opt<-1
        auto<-T
      }
      if(length(grep("ENSMUSG", rownames(normCounts)))>0){
        message("Mouse Ensembl IDs detected...")
        opt<-2
        auto<-T
      }
      if(length(grep("FBgn", rownames(normCounts)))>0){
        message("Flybase gene IDs detected...")
        opt<-3
        auto<-T
      }
      if(isFALSE(auto)){
        message("Cannot detect species automatically...")
        detID<-F
      }
      options<-c("hsapiens","mmusculus","dmelanogaster")
      targets<-c("HGNC","MGI","FLYBASENAME_GENE")
      if(class(res) =="list"){
        res<-lapply(res, BinfTools::getSym, obType="res", species=options[opt], target=targets[opt], addCol=T)
      } else {
        res<-BinfTools::getSym(res, obType="res", species=options[opt], target=targets[opt], addCol=T)
      }
      normCounts<-BinfTools::getSym(normCounts, obType="counts", species=options[opt], target=targets[opt], addCol=T)
    }
    if(isFALSE(retDDS)){
    return(list(normCounts=normCounts, res=res, condition=condition))
  } else {
    return(list(normCounts=normCounts, res=res, condition=condition, dds=dds))
  }
  }
}
