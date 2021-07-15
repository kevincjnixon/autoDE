# autoDE
Automated DE analysis by DESeq2

```{r}
#install.packages("devtools")
devtools::install_github("kevincjnixon/AutoDE")
devtools::install_github("kevincjnixon/BinfTools") #Required for converting gene ids to symbols and generating exploreData objects
devtools::install_github("kevincjnixon/gpGeneSets") #Required for exploreData function
library(AutoDE)
```

## Using HTSeq-count and sample Table:
Ensure the directory containing HTSeq-count results files is in your working directory.
Ensure that your sampleTable is a tab-delimited file in your working directory or a data.frame object already loaded in R.
The sampleTable should have the following columns (columns in brackets are optional):
 * sampleName = a unique identifier for each sample
 * fileName = the path pointing to the sample's count file (from the current working directory)
 * sampleCondition = Treatment condition for the sample to be used for comparisons. There must be at least two samples for each condition (replicates) - can be named anything other than 'condition'
 * [Batch] = Indicator for potential batch effects. There must be at least one sample from each condition belonging to the same batch
 * [OtherFactors] = You can include multiple other metadata columns indicating any other treatment conditions (other than sampleCondition). These must have unique column names and there must be at least two replicates in each of these conditions. These will be combined with sampleCondition into their own column named 'condition' to be used for multiple comparisons with DESeq2.
### Usage
 ```{r}
 #Run using tab-delimited sampleTable
 DERes<-autoDE(sampleTable="sampleTable.txt")

 #Run using sampleTable as data.frame already loaded in R
 DERes<-autoDE(sampleTable=sampleTable)
 ```

## Using featureCounts and or raw count matrix with colData:
Ensure the countTable is either a tab-delimited file in the working directory or a data.frame already loaded into R. The countTable should have rows as genes and columns as samples. countTable should contain RAW gene counts.
Ensure the colData is either a tab-delimited file in the working directory or a data.frame already loaded into R.
colData should have rows corresponding (in the same order) to the columns (samples) of countTable and at least one metadata column indicating the sampleCondition (see above - not named 'condition'). If there are batch effects, these should be in a column named 'Batch'. If there are multiple variables, these can also be included in the colData table using unique column names (not 'condition'). If colData is a tab-delimited file, the first column should indicate the unique sample names. If colData is a data.frame object in R, the rownames(colData) should indicate the unique sample names.
### Usage
```{r}
#Run using tab-delimited files
DERes<-autoDE(countTable="countTable.txt", colData="colData.txt")

#Run using loaded data.frame objects:
DEREs<-autoDE(countTable=countTable, colData=colData)
```

## Other arguments
autoDE() has two additional arugments that can be used with either sampleTable or countTable and colData:
 * expFilt - A numeric value indicating the minimum mean gene expression of a gene across all samples in order to be kept for analysis. Defaults to 0 (meaning if the gene is not expressed in any sample, it is not kept for analysis). Using larger values will filter more genes. Note that DESeq2 automatically applies strict filtering on the normalized counts when obtaining differential expression results.
 * retExplore - a Boolean (TRUE/FALSE) indicating if a BinfTools exploreData() object should be returned. This is a shiny app that allows you to browse all results objects generated. See [BinfTools exploreData](https://github.com/kevincjnixon/BinfTools/wiki/Explore-Data) for more information. Note that data must use gene symbols for exploreData. If ensembl gene ids (ENSG, ENSMUS) or flybase gene ids (FBgn) are identified in the data, you will be asked to indicate the species for gene symbol conversion.

### Usage
```{r}
#Using a sampleTable with HTSeq-count, filter genes to have an average raw expression of at least 10 and return the explorData() object
DERes<-autoDE(sampleTable="sampleTable.txt", expFilt=10, retExplore=T)

#Browse the exploreData() object
#Ensure BinfTools and gpGeneSets are loaded
library(BinfTools)
library(gpGeneSets)
DERes$explore #opens shiny app
```

## Results of autoDE()
autoDE() will return the complete DESeq2 analysis results for all possible comparisons (e.g. if you have two conditions only, there will be only one comparison. If you have 3 conditions, there will be 3 comparisons) as a list of results tables. The names within this list will indicate the comparison and the numerator and denominator condition for the log2FoldChange comparison. For example, if you have conditions named 'KO' and 'WT', the results table named 'KOvWT' will indicate that the 'KO' condition is the numerator and 'WT' is the denominator for the log2FoldChange (i.e. positive means higher numerator expression than denominator, negative is vice versa).
```{r}
DERes$res #List containing data.frames of DESeq2 results tables - named for each comparison
```
autoDE() will also return DESeq2 normalized expression as a data.frame. With columns as samples and rows as genes.
```{r}
DERes$normCounts
```

autoDE() will also return a character vector containing the final condition for each sample in the same order as the columns of normCounts. This can be used, along with the other objects, with the BinfTools package functions wherever conditions are required.
```{r}
DERes$cond
```

If retExplore=T, autoDE() will return a Shiny app generated by BinfTools::exploreData() to explore all comparisons and results and perform pahtway enrichment analyses/focused pathway analyses. This can be helpful in determining differentially expressed gene thresholds and quickly identifying any enriched pathways.
```{r}
DERes$explore
```

### Stay tuned for updates!
