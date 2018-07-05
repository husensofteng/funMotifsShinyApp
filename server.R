# This code is to make ShinyR application for funMotifs database
# By Nour-al-dain Marzouka
# Modified by Husen M. Umer
# Nour.dna.eng@gmail.com

# load the needed packages
library(RPostgreSQL)
library(DT)
library(plotly)
library(rjson)
library(pool)
library(DBI)
library(shiny)
library(data.table)
library(shinyjs)
# repress warnings
#options(warn=-1)
options()

# connect to the database
pool <- dbPool(
  drv = dbDriver("PostgreSQL", max.con = 100),
  dbname = "funmotifsdb",
  host = "130.238.239.17",
  user = "funmotifs",
  password = "ks2018hum"#,
  #idleTimeout = 3600000
)
onStop(function() { poolClose(pool) })

# get tissue names
temp_command <- "Select table_name from INFORMATION_SCHEMA.TABLES where TABLE_TYPE='BASE TABLE';"
tables <- dbGetQuery(pool, temp_command)$table_name
tables <- tables[grep(tables,pattern = "pg_",invert = T)]
tables <- tables[grep(tables,pattern = "sql_",invert = T)]
tables <- tables[grep(tables,pattern = "motif",invert = T)]
tables <- tables[grep(tables,pattern = "cell",invert = T)]

# query for the variants
make_command_SNP <- function(tissue_table, chr_table, start, end, ref, mut){
  # fix the chromosome name
  chr_table <- gsub(x = chr_table, pattern = "chr", replacement = "")
  
  # get the table
  temp_command <- paste0("select concat(chr,':', motifstart, '-', motifend) as motif_pos, name , score , pval , strand , fscore , chromhmm , contactingdomain , dnase__seq , fantom , loopdomain , numothertfbinding , othertfbinding , replidomain , tfbinding , tfexpr , (CASE WHEN (UPPER(posrange * int4range(",start,",",end,")) - LOWER(posrange * int4range(",start,",",end,"))>1) THEN 100 ELSE (CASE when STRAND='-' THEN (motifend-",start,") ELSE (",end,"-motifstart)+1 END) END) as mutposition from chr",chr_table,"motifs,",tissue_table," where chr",chr_table,"motifs.mid=",tissue_table,".mid and posrange && int4range(",start,",",end,",'(]') and tfexpr>0.0;")
  query        <- dbGetQuery(pool, temp_command)
  
  #print("here")
  #print(nrow(query))
  if (nrow(query)==0) {return(data.frame(matrix(NA,nrow=1,ncol=20),stringsAsFactors = F))}
  
  #get the entropy
  query$entropy <- NA
  for (i in 1:nrow(query)) {
    name        <- query[i, "name"]
    mutposition <- query[i, "mutposition"]
    
    if (mutposition != 100){
      temp_command <- paste0("select ((select freq from motifs_pfm where position=",mutposition,
                             " and name='",name,"' and allele='",ref,
                             "') - (select freq from motifs_pfm where position=",mutposition,
                             " and name='",name,"' and allele='",mut,"') ) as entropy;")
      query2        <- dbGetQuery(pool, temp_command)
      query$entropy[i] <- query2[1,1]
    } else {
      query$entropy[i] <- 1
      query$mutposition[i] <- 'multi'
    }
    
  }
  return(query)
}

#make_command_SNP("blood",1,101190407,101190408,"C","A")

#### function to make command for a region
make_command_region <- function(tissue_table, chr_table, start, end){
  # fix the chromosome name
  chr_table <- gsub(x = chr_table, pattern = "chr", replacement = "")
  
  temp_command <- paste0("select  concat(chr,':', motifstart, '-', motifend) as motif_pos, name , score , pval , strand , fscore , chromhmm , contactingdomain , dnase__seq , fantom , loopdomain , numothertfbinding , othertfbinding , replidomain , tfbinding , tfexpr from chr",chr_table,"motifs , ",tissue_table," where chr",chr_table,"motifs.mid=",tissue_table,".mid and posrange && int4range(",start,",",end,",'(]') and tfexpr>0.0;")
  #print(temp_command)
  query        <- dbGetQuery(pool, temp_command)
  
  if (nrow(query)==0) {return(data.frame(matrix(NA,nrow=1,ncol=18),stringsAsFactors = F))}
  return(query)
}

#### function to make command for all tissues
make_command_all_tissues <- function(chr_table, start, end){
  # fix the chromosome name
  chr_table <- gsub(x = chr_table, pattern = "chr", replacement = "")
  
  temp_command <- paste0("select * from chr",chr_table,"motifs,all_tissues where chr",chr_table,
                         "motifs.mid=all_tissues.mid and posrange && int4range(",start,",",end,", '(]');")
  
  query        <- dbGetQuery(pool, temp_command)
  if (nrow(query)==0) {return(data.frame(matrix(NA,nrow=1,ncol=22),stringsAsFactors = F))}
  query        <- query[, !colnames(query) %in% c("mid","posrange")]
  return(query)
}

# function to prepare table from the user text
prepare_table_text <- function(text, sep){
  if (text==""){
    return(NULL)
  } else {
    text <- read.table(text=text,
                       sep = sep, 
                       blank.lines.skip = T, 
                       stringsAsFactors = F)
    return(text)
  }
}

#the limit of the file upload (50MB)
options(shiny.maxRequestSize=50*1024^2)

###########################################################
############### ShinyR server function ####################
###########################################################
shinyServer(function(input, output, session) {
  
  # the tissue list
  output$list_tissues <- renderUI({
    selectInput("tissues", "Select a tissue", 
                sort(tables,decreasing = F), 
                width = 200,
                selected = 'liver', 
                multiple = FALSE)
  })
  
  #create reactive list
  reactive                <- reactiveValues()
  reactive$output_table   <- NULL
  
  # reading the user BED file or the text area
  mydata <- reactive({
    
    inFile <- input$user_file
    
    if (is.null(inFile) & input$text_input_area == ""){
      return(NULL)
    }
    
    if (input$text_input_area != "") {
      tbl <-   prepare_table_text(input$text_input_area, input$sep)
      #validate(
      #  need(ncol(mydata()) > 1, "Please select the right seperator option!")
      #)
    }
    
    if (!is.null(inFile)){
      tbl <- read.csv(inFile$datapath, 
                      header=input$header, 
                      sep=input$sep,  
                      dec = input$dec,stringsAsFactors = F)
    }
    
    return(tbl)
  })

#The following table is shown in the user manual
  column_headers <- data.table(Header = c("motif_pos", "name", "score", "pval", "strand", "fscore", "chromhmm", "contactingdomain", "dnase__seq", "fantom", "loopdomain", "numothertfbinding", "othertfbinding", "replidomain", "tfbinding", "tfexpr"), Description= c("Represents the chromsome number (1-25), start position of the overlapping motif (1-based), and end position of the overlapping motif", "TF name and motif ID", "Motif prediction score calculated by FIMO", "P-value for the motif prediction calculated by FIMO", "DNA strand", "Functionality score of the motif. This is a summarized score to indicate activity of the motif. It is computed using a ligistic regression model for details see the paper", "Chromatin states predicted based on histone modification marks using ChromHMM in the RoadMap Epigenomics project", "HiC contacting domain, the input region is overlapping one of the two sides of an interaction that is identified in HiC experiments. In cases where the tissue type does not have any HiC data, the average of tissue types with HiC contact domain at the loci are given.", "DNase1 hypersensitivity site, the value indicates the signal value as obtained from average DNas1 peaks across relevant data sets", "TSS expression at the loci as obtained from CAGE experiments in the FANTOM5 project ", "Similar to contactingdomain, but the region between the interaction points", "Number of TFs that are bound at the loci (excl. the matching TF)", "List of the bound TFs", "Replication timing category of the loci. ERD: arly replication domain; DTZ: down transition zone; LRD: late replication domain; UTZ: up transition zone", "Matching TF binding. A value above zero indicates evidence of binding. NA refers to the lack of data for the given factor", "Expression level of the TF in the tissue"))


  output$tbl <- renderTable({ head( column_headers, n=20)})
  # Display the user input
   output$show_user_table <- renderTable({
     validate(
      need(ncol(mydata()) > 1, "Please enter or upload your input regions and select the right seperator option!")
     )
     head(mydata(),5)
   })
  
  # print the first 10 rows of the data
  output$output_table <- renderTable({
    if (is.null(input$tissues) | is.null(input$SNP_or_region) | is.null(mydata())) {
      validate(
        need(ncol(mydata()) >0, "Please enter or upload a list of genomic regions or variants in the left panel in 0-based BED format and select a tissue type!")
      )
      return (NULL)
      } else {
        
        validate(
          need(ncol(mydata()) > 1, "Please select the right seperator option!")
        )
        # for all tissues
        if (input$tissues == "all_tissues"){
          make_command_all_tissues
          x <- apply(mydata(),1,function(x){
            make_command_all_tissues(chr_table = gsub(' ','',x[1],fixed = TRUE),
                                     start     = x[2],
                                     end       = x[3])
          })
          
          # merge the output with user columns
          from_user <- mydata()[rep(row.names(mydata()), sapply(x, nrow)), 1:3]
          colnames(from_user) <- c("chr","start","end")
          x <- rbindlist(x)
          x <- cbind(from_user,x)
          reactive$output_table <- x
          return(head(reactive$output_table)[,1:10])
        }
        
        # for regions
        if (input$SNP_or_region == "region") {
          validate(
            need(ncol(mydata()) >= 3, 
                 "Please select the right type of your input data (variants or genomic regions)!")
          )
          
          x <- apply(mydata(),1,function(x){
            make_command_region(tissue_table = input$tissues,chr_table = gsub(' ','',x[1],fixed = TRUE),start = x[2],end = x[3])})
          
          # merge the output with user columns
          from_user <- mydata()[rep(row.names(mydata()), sapply(x, nrow)), 1:3]
          colnames(from_user) <- c("chr","start","end")
          x <- rbindlist(x)
          x <- cbind(from_user,x)
          reactive$output_table <- x
          return(head(reactive$output_table)[,1:10])
        }
        
        # for variants
        if (input$SNP_or_region == "variant") {
          validate(
            need(ncol(mydata()) >=5, 
                 "Please select the right type of your input data (variants or genomic regions)!")
          )
          
          validate(
            need(class(mydata()[1,4])=="character" & class(mydata()[1,5])=="character", 
                 "Column 4 & column 5 should be characters!")
          )
          
          if (ncol(mydata())<5) {
            return(NULL)
          }
          
          x <- apply(mydata(),1,function(x){
            make_command_SNP(tissue_table = input$tissues,
                             chr_table = gsub(' ','',x[1],fixed = TRUE),
                             start = x[2],
                             end = x[3],
                             ref = x[4],
                             mut = x[5])
          })
          
          if (is.null(x)){return(NULL)}
          
          num_row <- sapply(x, nrow)
          num_row[num_row==0] <- 1
          # print(num_row)
          
          # merge the output with user columns
          from_user <- mydata()[rep(row.names(mydata()), num_row), 1:5]
          colnames(from_user) <- c("chr",
                                   "start",
                                   "end",
                                   "ref",
                                   "alt")
          # print(from_user)
          # print(x)
          x <- rbindlist(x)
          
          x <- cbind(from_user,x)
          reactive$output_table <- x
          return(head(reactive$output_table)[,1:12])
        }
      }})
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste("output_motifs.csv", sep = "")
    },
    content = function(file) {
      x <- reactive$output_table    
      # check if the not found regions should be included or not
      if (input$no_NAs) {
        x <- x[!is.na(x[,"fscore"]),]
      }
      write.csv(x, file, row.names = FALSE)
    }
  )

 output$downloadData_logical<- renderUI({
  if(!is.null(mydata()) & !is.null(reactive$output_table)) {
    downloadButton('downloadData', 'Download Output File')
  }
})

 
}) # END of server
