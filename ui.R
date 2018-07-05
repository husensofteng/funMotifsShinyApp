library(shiny)
library(markdown)
library(shinythemes)
library(shinyBS)

shinyUI(
  navbarPage(windowTitle= "funMotifs DB",title=div(a(img(src="funmotifs.png", width="30%", height="50%"), href="http://bioinf.icm.uu.se:3838/funmotifs/"), ""), 
             id="main_page",
             inverse = F, 
             theme = "bootstrap.min.css",#shinytheme("cerulean"),
	     selected="Annotate genomic regions and variants",
             # First panel    
             tabPanel("Annotate genomic regions and variants",
                      fluidPage(titlePanel("Enter or upload genomic regions"),
                                sidebarLayout(
                                  sidebarPanel(
                                    wellPanel(
                                      div(style="display:inline-block;vertical-align:top", 
                                          
                                          radioButtons(inputId = "SNP_or_region",selected = "variant",
                                                       label = "Specify the input type:",
                                                       choices = c("Genomic Regions"="region","Variants"="variant")),
                                          uiOutput('list_tissues')
                                      ),
#                                      tags$hr(), # seperator
                                      
                                      # text input area
                                      div(style="display:inline-block;vertical-align:top",
                                          textAreaInput("text_input_area", "Enter your input here", 
                                                        value = "",
                                                        width = NULL, height = NULL,
                                                        cols = 60, rows = 7, 
							placeholder = "#Sample region input:\nchr1\t936320\t936325\n\n#Sample variant input:\nchr1\t936324\t936324\tA\tG", 
							resize = "both")
                                      ),
                                      
                                      
                                      div(style="display:inline-block;vertical-align:middle",
                                          tags$p(tags$strong("or Upload a BED/CSV file:"), "(sample input files: ", 
						tags$a(href="regions_sample_input.txt", target='blank', 'regions', download = 'regions_sample_input.txt'),' | ',
						tags$a(href="variants_sample_input.txt", target='blank', 'variatns', download = 'variants_sample_input.txt'), ")"),
                                          fileInput('user_file', label=NULL,width = 200, multiple=FALSE,
                                                    accept=c('text/csv', 
                                                             'text/comma-separated- values,text/plain', 
                                                             '.csv'))
					),
					div(style="display:inline-block;vertical-align:top",
                                          radioButtons('sep', 'Separator',inline = T,
                                                       c(Tab='\t',
							 Comma=',',
                                                         Semicolon=';',
                                                         Space=' '),
                                                       '\t'),
					  tags$hr(),
                                          checkboxInput('header', 'Header', FALSE)
                                      ),

                                      # the GO button 
                                      div(style="display:inline-block;vertical-align:top",
                                          submitButton(text= "Submit"),
                                          
                                          tags$hr(), # seperator    
                                          # Download Button
                                          #downloadButton("downloadData", "Download Results"),
                                        uiOutput("downloadData_logical"),  
					checkboxInput('no_NAs', 'Exclude input regions with no matches?', TRUE)
                                      )
)
                                  ) # end of the sidebar
                                  ,
                                  mainPanel(
					bsCollapse(id='outputCollapse', multiple=TRUE,
                                        bsCollapsePanel("A sample from the input data", uiOutput('show_user_table'))),
					#bsCollapsePanel("A sample from the results", uiOutput('output_table'), style="primary")
				
                                    #wellPanel(HTML("<b>A sample from the input data:</b>"),uiOutput('show_user_table'))
                                    
                                    wellPanel(HTML("<h3><b>A sample from the results</b></h3> Download the results for the full set by clicking on the 'Download Output File' button on the left (note that not all rows and columns are shown here)"),tableOutput('output_table'))
                                  ) # end of the mainpanel
                                ) # end of the panel
                      )),
             
             #First panel
              tabPanel("Download",
                       fluidPage(titlePanel(" "),
				fluidRow(
                                   column(3, tags$h3("")),
                                   column(8,offset=0, 
                                    	#tags$strong("Download extracted BED fiels from the funMotifs database"),
					tags$h3("Download extracted fiels from funMotifsDB"),
					tags$hr(),
					includeHTML("download.html"),
					tags$hr()
				))
			    )
			)
,
	tabPanel("Manual", fluidPage(titlePanel("User Guide"),
                                fluidRow(
					column(10, offset=1, 
            					bsCollapse(id='manualCollapse', multiple=TRUE,open='Content', 
					bsCollapsePanel("Content", HTML("
					The current version of the funMotifs database (v1.0) contains annotations for predicted motifs of 510 Transcription Factors (TFs) from the JASPAR2016 vertebrates CORE database, <a href='http://jaspar2016.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt'>click here to see the entire TF list</a>.<br/><br/><h4>Annotations and assays from the following resources were used to annotate the motifs:</h4>ChIP-seq datasets from <a href='https://www.encodeproject.org/'>ENCODE</a>.<br/>DNase1-seq datasets from the ENCODE and <a href='http://www.roadmapepigenomics.org/'>RoadMap Epigenomics</a> projects.<br/>CAGE peaks in promoters and enhancers from the <a href='http://fantom.gsc.riken.jp/5/datafiles/latest'>FANTOM5</a> project.<br/>Chromatin states from RoadMaps<a href='https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/'>15-state core marks model</a>.<br/>Replication domains from <a href='https://doi.org/10.1093/bioinformatics/btv643'>Liu F. et al</a>.<br/>HiC contacting domains from <a href='https://doi.org/10.1016/j.cell.2014.11.021>Rao, S.S., et al</a>.<br/>Gene expression data from <a herf='http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz'>GTEx</a> and ENCODE.<br/>Regulatory elements from MPRAs <a href='https://doi.org/10.1038/nbt.3678'>Ernst, J.. et al.</a>, <a href='https://doi.org/10.1016/j.cell.2018.02.021'>Tewhey, R., et al</a>, and <a href='https://doi.org/10.1016/j.celrep.2016.07.049'>Vockley, C.M., et al</a>.<br/>Further details on the pre-processing of the datasets are listed in their corresponding files on <a href='https://github.com/husensofteng/funMotifs/tree/master/ReadMe'>the project repository in github</a>.

					
					"), style="primary"),
					bsCollapsePanel("What input format is accepted?", tags$p("
					The input should be a table of genomic coordinates with at least three columns: chromosome number, start and end positions. The coordinates should be ", tags$a(href="http://bedtools.readthedocs.io/en/latest/content/overview.html#bed-starts-are-zero-based-and-bed-ends-are-one-based", "BED 0-based"), ".", tags$br(),tags$br(),
                                             "For variants, the input should have at least five columns, first three columns as above, and a 4th column for the refrence allele and a 5th column for the alternative allele. Any additional column from the user will be ignored.",tags$br(), tags$br(), "The columns can be seprated by tab (default), comma, space, or semi-colon. Please select the right separator option on the left side panel.", tags$br(), tags$br(), "Check the Header checkbox if the input contains a header line.
					"), style="success"),
					bsCollapsePanel("What does the output columns represent?", tags$p(tableOutput('tbl'), "*By default input regions or variants that do not overlap any TF motif are not reported in output file. Uncheck the box 'Exclude input regions with no matches?' to retrieve all entries in the input list regardless their overlap with TF motis."), style="success"),
	
					bsCollapsePanel("Batch Analysis", tags$p("We provide extracted files from the database in the download section. For batch analysis, we recommend downloading the funMotifs file for the desired tissue type and use IntersectBed to identify the overlapping motifs.", tags$br(), tags$br(), "In order to search in all motifs regardless of their functionality you can download the", tags$a(href='http://bioinf.icm.uu.se:3838/funmotifs/datafiles.tar.gz', "all_tissues"), "archive file in the download section. The file contains all motifs (~85 million) and the fscore for each tissue type.", tags$br(), tags$br(),"Be aware that, chromosome X,Y and M are represented as 23, 24, and 25, respectively."), style="success")

					),
					tags$hr(),
					wellPanel(
                                      		tags$p("")
                                    ))
                                  )))
		,
#About page
tabPanel("About",
         fluidPage(titlePanel(""),
                   fluidRow(
                     column(10,offset=1, 
                            wellPanel(wellPanel(tags$p(tags$strong("About"), tags$br(), "funMotifs database contains annotated transcription fcator motifs. The aim is to aid researchers in interpreting the noncoding genome. It enables analysis of noncoding variatns and regions. funMotifs is built based on assays and datasets from ENCODE, FANTOM, RoadMap Epigenomics, GTEx and other published works. For futher details about the content and the methods used to build funMotifs please refer to the manual page.")), 
                            wellPanel(tags$p(tags$strong("Availability"), tags$br(), "The motif annotations per tissue type can be downloaded for batch analysis. We also provide the pipeline that is used to build the funMotifs databse. The funMotifs pipeile that is freely available on github allows re-generation of the database using additional datasets:", tags$a(href="https://github.com/husensofteng/funMotifs", "(funMotifs source code)"))),
                            wellPanel(tags$p(tags$strong("Citation"), tags$br(), "Use the following to cite funMotifs: Umer et al. 2018 (in progress).")),
                            wellPanel(tags$p(tags$strong("Contact:"), tags$br(),"For questions and issues regarding funMotifs please contact", tags$a(href="mailto:funmotifs@gmail.com", "funMotifs Support")))
                     )
                   ))
         )
),

		tags$head(tags$style( HTML(' .nav {margin-top:10px;margin-left:300px}')))
  )) # END of UI 
