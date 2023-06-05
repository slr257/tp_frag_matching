# Script for matching in-silico predicted fragments to measured MS2 data
# Stephanie Rich - Cornell University June 2023

#### 1. Begin by clearing the workspace and installing required R packages --------------------------------------------

# clear workspace
rm(list = ls())

# set working directory
setwd("~Desktop/repos/tp_frag_matching")

# install and load necessary packages
library(mzR)
library(RMassBank)
library(readxl)
library(svMisc)
library(xlsx)

# Note- Some packages may need to be installed from the Bioconductor website
# for example:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("RMassBank")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("mzR")

#### 2. Upload required data into working directory --------------------------------------------

# Create a folder titled "mzXML", add it to the working directory, and input the address here
file.path <- "C:/Users/Desktop/repos/tp_frag_matching"

# Create a folder containing "long reports" containing peak area information in the format exported from XCalibur QuanBrowser
file.path.long <- "C:/Users/Desktop/repos/tp_frag_matching/Long_Reports"

# Add multiple long report names and specific the number of long reports
n <- 7

# Name of long reports
file.name1 <- "EFF_suspectscreening_Cluster1_long"
file.name2 <- "EFF_suspectscreening_Cluster2_long"
file.name3 <- "EFF_suspectscreening_Cluster3_long"
file.name4 <- "EFF_suspectscreening_Cluster4_long"
file.name5 <- "EFF_suspectscreening_Cluster5_long"
file.name6 <- "EFF_suspectscreening_Cluster6_long"
file.name7 <- "EFF_suspectscreening_Cluster7_long"

# Set "Database" file path to excel sheet containing information on in-silico predicted fragments
database.path <- "C:/Users/Desktop/repos/tp_frag_matching/TP_Database_All_Frags.xlsx"

# Load database
database <- as.data.frame(read_excel(database.path, sheet = 1))
names(database)[9:ncol(database)] <- "Fragment"

#### 3. Define functions used for organizing MS2 data and long reports with sample names --------------------------------------------

# function to get MS2 data
getData2 <- function(s){
  peaks <- s@peaksCount
  cols <- c("mz", "intensity", "satellite", "low", "rawOK", "good", "mzCalc", "formula", "dbe", "formulaCount", "dppm", "dppmBest")
  cols.isFilled <- unlist(lapply(cols, function(col) length(slot(s, col)) == peaks))
  cols.filled <- cols[cols.isFilled]
  data <- lapply(cols.filled, function(col) slot(s, col))
  data$stringsAsFactors <- FALSE
  df <- do.call(data.frame,data)
  colnames(df) <- cols.filled
  df
}

filenames <- rep(NA, n)
for(i in 1:n){
  filenames[i] <- paste0(file.path.long, "/", get(paste0("file.name",i)), ".xls")
}

# function to read in long reports and organize them
read_excel_allsheets <- function(filename) {
  sheets <- readxl:::excel_sheets(filename)
  sheets <- sheets[!(sheets %in% c("Component","mdlCalcs"))]
  #x <- lapply(sheets, function(X) xlsx:::read.xlsx(filename, sheetName = X, startRow = 5, colIndex = 1:47))
  x <- lapply(sheets, function(X) readxl:::read_excel(filename, sheet = X, skip = 4))
  x <- lapply(x, as.data.frame)
  x <- lapply(1:length(x), function(X) x[X][[1]][1:(nrow(x[[1]][1])-6),])
  names(x) <- sheets
  x
}

l <- lapply(filenames, function(X) read_excel_allsheets(X))

l2 <- NA
for(i in 1:n){
  temp <- l[[i]]
  l2 <- c(l2, temp)
}
l <- l2[-1]
temp <- NULL
l2 <- NULL


# Now we have the database loaded as a dataframe (database) and the long report as a list (l)

### 4.  Begin setting up loops  for analysis - one for TPs and one for samples --------------------------------------------

# character vector of TP names
names.tp <- names(l)[1:256]

# character vector of sample names
names.sam <- l[1][[1]]$`Sample Name`


# Get sample areas
sam.df <- data.frame(names.tp, row.names = 1)
for(i in names.tp){
  x <-  l[i][[1]]
  x <- subset(x, x$'Sample Type' == "Unknown Sample")
  for(j in x$'Sample Name') sam.df[i, j] <- x[x$'Sample Name'==j, "Area"]
  sam.df[i,] <- as.numeric(sam.df[i,])
  
}

#### 5. Preallocate resulting dataframes (sam.df and frags.df) --------------------------------------------

# fill in this dataframe with T/F of matched fragments!
frags.df <- data.frame(matrix(ncol = length(colnames(sam.df)), nrow = length(rownames(sam.df))))
colnames(frags.df) <- colnames(sam.df)
rownames(frags.df) <- rownames(sam.df)

# dataframe with values of matching fragments
matching.frags.val <- frags.df


# read in all mzXML files in the folder labeled mzXML
# this step is edited to be shorter for testing purposes
filenames <- list.files("mzXML", pattern="*.mzXML", full.names=TRUE)
ldf <- lapply(filenames, openMSfile)

#### 6. Run fragment-matching for-loop --------------------------------------------

# for-loop through tp names from a character vector!
for(x in 1:length(rownames(sam.df))){
  # progress and optimization metrics
  # progress shows value inside parentheses
  progress(x*100/length(rownames(sam.df)))
  start.time <- Sys.time()
  tp <- rownames(sam.df)[x]
  message("Beginning Fragment Search for suspect ", tp)
  
  for(j in 1:length(colnames(sam.df))){
    # pick out a sample by name in quotes
    sam <- colnames(sam.df)[j]
    message("Beginning Fragment Search for sample ", sam)
    
    if(is.na(sam.df[x,j]) == FALSE){ 

      if(is.na(l[tp][[1]]$`Sample Name`[1]) == FALSE){ # if statment to check for blank excel sheets in long reports
        id <- which(sam == l[tp][[1]]$`Sample Name`)[1]
        
        database.sub <- database[database$"Compound Name"==tp,]
        
        # get name of mzXML file name from long report also
        name.file <- l[tp][[1]]$Filename[id]
        
        # access mzXML file from list before the for loop with the correct name
        msRaw <- ldf[[which(paste0("mzXML/",name.file,".mzXML") == filenames)]]
        
        # get mass to find with ppm deviation
        ppm <- 5
        
        mz <- database.sub$`Extracted Mass`
        mz.l <- mz - ppm(mz, ppm, p=TRUE)
        mz.h <- mz + ppm(mz, ppm, p=TRUE)
        
        
        # get retention time from long report (l)
        # the extracted value is the mean RT in minutes for a specific compound, so it is multiplied by 60 to get RT in seconds
        rt <- mean(as.numeric(l[tp][[1]]$RT), na.rm = T)*60
        
        # this is the rt from the long report - "EIC"
        rt.eic <- rt
        
        # upper and lower acceptable RT limits (+/- 30 seconds)
        rt_lim_low <- rt-30
        rt_lim_high <- rt+30
        
        # find MS2 spectra
        # Compound m/z
        fine <- ppm(mz, 5)
        coarse <- 0.5
        
        # this function generates a list with information about the specified peak with rt limits
        # this might take a lot of time
        MSMS <- findMsMsHR.mass(msRaw,mz, 
                                limit.fine = fine, 
                                limit.coarse = coarse,
                                rtLimits = c(rt_lim_low,rt_lim_high),
                                deprofile = NA)
        
        if(MSMS[[1]]@found == TRUE & is.nan(rt.eic) == FALSE){
          
          MSMS.rt <- rep(NA,length(MSMS))
          for(i in 1:length(MSMS)){
            if(MSMS[[1]]@found){
              MSMS.rt[i]<-MSMS[[i]]@parent@rt}else{MSMS.rt<-NA}}
          
          # check polarity of all MSMS scans
          MSMS.pol <- rep(NA,length(MSMS))
          for(i in 1:length(MSMS)){
            if(MSMS[[1]]@found){
              MSMS.pol[i]<-MSMS[[i]]@children@listData[[1]]@polarity}else{MSMS.pol<-NA}}
          
          MSMS.pos <- MSMS[MSMS.pol == 1]
          
          #only include positive MSMS scans
          if(length(MSMS.pos) > 0){MSMS.sel <- which.min(abs(MSMS.rt-rt.eic))
          
          # gets fragment masses and intensities from the peak with the closest RT to the selected RT
          MS2.ms <- getData2(MSMS[[MSMS.sel]]@children@listData[[1]])
          MS2.ms$mz2 <- as.numeric(format(MS2.ms$mz,nsmall=2,digits=1))
          
          # order fragments from largest to smallest intensity
          MS2.ms <- MS2.ms[order(-MS2.ms$intensity),]
          
          # get fragments from database
          frags <- NA
          frags2 <- database.sub[,which(colnames(database.sub) == "Fragment"):ncol(database.sub)]
          for(i in 1:ncol(frags2)){
            temp <- as.numeric(frags2[,i])
            frags <- c(frags, temp)
          }
          
          # get unique fragments from database and round the values to two decimal places
          frags <- unique(na.omit(as.numeric(format(frags[-1],nsmall=2,digits=1))))
          
          # make a data frame of matching fragments
          frags.match <- data.frame(mz = frags[frags %in% MS2.ms$mz2])
          
          # match fragments with MS2.ms
          frags.match3 <- subset(MS2.ms, mz2 %in% frags.match$mz)[,-1]
          frags.match3 <- frags.match3[order(-frags.match3$intensity),]
          message("Matching Frags ", frags.match)
          
          #frags.match2 <- frags.match
          
          message("Most intense frag ", frags.match3[1,2])
          
          
          # boolean statement if there is a matching fragment or not
          frags.match <- ifelse(sum(rowSums(frags.match)) > 0, T, F)
          
          frags.df[x,j] <- frags.match
          }else{frags.df[x,j] <- FALSE}
          
        }else{frags.df[x,j] <- FALSE}
        
        if(frags.df[x,j] == TRUE){matching.frags.val[x,j] <- frags.match3[1,2]}else{matching.frags.val[x,j] <- NA}
        
      }else{frags.df[x,j] <- FALSE}
    }else{frags.df[x,j] <- FALSE}
  } # end sample for-loop (j index)
} # end tp for-loop (x index)


#### 7. Export results as csv files --------------------------------------------

# peak areas of TPs by sample name
write.csv(sam.df, "sample.df.csv")

# TRUE/FALSE of samples with matching fragments
write.csv(frags.df, "frags.df.csv")

# mass of most intense matching fragment
write.csv(matching.frags.val, "matching.fragments.csv")


#### 8. Export results as a formatted Excel Workbook --------------------------------------------

# exporting data
wb <- createWorkbook(type="xlsx")
sheet1 <- createSheet(wb, sheetName = "Peak Area")
sheet2 <- createSheet(wb, sheetName = "Fragments_Boolean")
sheet3 <- createSheet(wb, sheetName = "Diagnostic Fragments")
#sheet4 <- createSheet(wb, sheetName = "KEY")


# Title and sub title styles--------------------
## Define styles
TITLE_STYLE <- CellStyle(wb) + Font(wb,  heightInPoints=16, isBold=TRUE, underline=1) + Alignment(horizontal="ALIGN_CENTER")
SUB_TITLE_STYLE <- CellStyle(wb) + Font(wb,  heightInPoints=14) + Alignment(horizontal="ALIGN_CENTER")
SUB_TITLE_STYLE2 <- CellStyle(wb) + Font(wb,  heightInPoints=14) + Alignment(horizontal="ALIGN_LEFT")


TABLE_FILLED_STYLE   <- CellStyle(wb, fill = Fill(foregroundColor="lightgray"))
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), pen=c("BORDER_THIN", "BORDER_THICK")) 
KEY_STYLE <-            CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=F, horizontal="ALIGN_CENTER")
KEY_DET_STYLE <-        CellStyle(wb, fill = Fill(foregroundColor="lightgreen")) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=F, horizontal="ALIGN_CENTER")
KEY_OVER_STYLE <-       CellStyle(wb, fill = Fill(foregroundColor="lightblue")) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=F, horizontal="ALIGN_CENTER") 
KEY_FRAG_STYLE <-       CellStyle(wb, fill = Fill(foregroundColor="lightyellow")) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=F, horizontal="ALIGN_CENTER") 
#----------------------------------------------------------------------


# cell style for if the peak area is >= 1E5 and there is a diagnostic fragment
DET_STYLE <-          CellStyle(wb, fill = Fill(foregroundColor="lightgreen", backgroundColor="lightgreen"),
                                dataFormat = DataFormat("#,##0.0")) + Font(wb, isBold=TRUE)

# cell style for if fragment is found and peak area is < 1E5

# cell style for if no fragment is found and peak area is > 1E5

# cell style if no fragment found and peak area < 1E5
ND_STYLE <-           CellStyle(wb, fill = Fill(foregroundColor="white", backgroundColor="white"),
                                dataFormat = DataFormat("#,##0.0")) + Font(wb, isBold=TRUE)
NA_STYLE <-           CellStyle(wb) + Font(wb, color="red") + Alignment(horizontal="ALIGN_RIGHT")





addDataFrame(sam.df, sheet1,
             startRow = 3, startColumn = 1,
             showNA = T, characterNA = "NA", row.names = T)
addDataFrame(frags.df, sheet2,
             startRow = 3, startColumn = 1,
             showNA = T, characterNA = "NA", row.names = T)
addDataFrame(matching.frags.val, 
             sheet3, startRow = 3, startColumn = 1, 
             showNA = T, characterNA = "NA", row.names = T)
# addDataFrame(key, sheet4, 
#              startRow = 2, startColumn = 1, row.names = F, col.names = F,
#              colStyle = list("1" = KEY_STYLE))

rows.pa <- getRows(sheet1)
rows.dat  <- getRows(sheet2)
rows.mat <- getRows(sheet3)

for(i in 2:nrow(sam.df)+1){
  for (j in 2:ncol(sam.df)+1){
    cell <- getCells(rows.dat[i], colIndex = j)[[1]]
    val  <- getCells(rows.pa[i], colIndex = j)[[1]]
    mat  <- getCells(rows.mat[i], colIndex = j)[[1]]
    frag <- getCellValue(cell)
    pa <- getCellValue(val)
    mat <- getCellValue(mat)
    
    if(frag == TRUE & is.na(pa) == FALSE){
      style <- "DET"
    }else if(frag == FALSE & is.na(pa) == FALSE){
      style <- "ND"
    } else {style = "NA"}
    
    
    cell2 <- getCells(rows.pa[i], colIndex = j)[[1]]
    if(style == "DET"){
      setCellStyle(cell2, DET_STYLE)
    } else if(style == "ND"){
      setCellStyle(cell2, ND_STYLE)
    } else{setCellStyle(cell2, NA_STYLE)}
    
    
    cell3 <- getCells(rows.mat[i], colIndex = j)[[1]]
    if(style == "DET"){
      setCellStyle(cell3, DET_STYLE)
    } else if(style == "ND"){
      setCellStyle(cell3, ND_STYLE)
    } else{setCellStyle(cell3, NA_STYLE)}
  }
}

## Change column widths
setColumnWidth(sheet1, colIndex = 1, colWidth = 40)
setColumnWidth(sheet1, colIndex = c(2:ncol(sam.df)+1), colWidth = 15)

setColumnWidth(sheet2, colIndex = c(1:9), colWidth = 20)
setColumnWidth(sheet2, colIndex = c(1,6), colWidth = 35)


file.output <- "Fragment_Results_Excel"
saveWorkbook(wb, paste0(file.path,"/",file.output,".xlsx"))

### END ####

