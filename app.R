## app.R ##
Rlib="/data/boehm/group/shiny_apps/Rlibs3.4.0"

require(shiny,lib.loc=Rlib)
require(shinydashboard,lib.loc=Rlib)
require(rhandsontable,lib.loc=Rlib)

ui <- function(request) {dashboardPage(
    dashboardHeader(title = "Dataset selection"),
    ## Sidebar content
    dashboardSidebar(

        textInput(inputId="group", label="Group", value = "", width = NULL, placeholder = NULL),
        textInput(inputId="owner", label="Project Owner", value = "", width = NULL, placeholder = NULL),
        textInput(inputId="projectid", label="Project ID", value = "", width = NULL, placeholder = NULL),
        actionButton(inputId="selectdataset", label="Select dataset"),
        selectInput(inputId="selectworkflow",label="Select NGS workflow",choices=c("ATAC-seq","ChIP-seq","DNA-mapping","HiC","RNA-seq","WGBS","scRNAseq"))
        
        ),
        
    dashboardBody(
        h2("User input"),
        uiOutput("resultPanels")
               
    )

 )}



server <- function(input, output, session) {
        
       require("yaml",lib.loc=Rlib)



        observeEvent(input$selectdataset, {
        
        if((input$group!="")&(input$owner!="")&(input$projectid!="")){
            inGroup<-isolate(input$group)
            inOwner<-isolate(input$owner)
            inProjectID<-isolate(input$projectid)

            datdir<-system(sprintf("find /data/%s/sequencing_data -name Project_%s_%s_%s -type d",tolower(inGroup),inProjectID,inOwner,inGroup),intern=TRUE)
            #handle single and paired end data
            datPath<-dir(datdir[1],pattern="*.fastq.gz",full.names=TRUE,recursive=TRUE)###allows for only 1 sequencing run!!
            #get rid of optical duplicates
            datPath<-grep("*optical*",datPath,value=TRUE,invert=TRUE)
            if(sum(grepl("*_R2.fastq.gz",datPath))>0){output$PLACEHOLDER<-renderText("Paired reads detected.")
                    output$Read1<-renderText(grep("*_R1.fastq.gz",datPath,value=TRUE))
                    output$Read2<-renderText(grep("*_R2.fastq.gz",datPath,value=TRUE))
                    datshort<-gsub("_R1.fastq.gz","",basename(grep("*_R1.fastq.gz",datPath,value=TRUE)))}
            else{output$PLACEHOLDER<-renderText("Single reads detected.")
                 output$Read1<-renderText(datPath)
                 datshort<-gsub(".fastq.gz","",basename(datPath))}


            ###############initiate reactive table to collect sample information ###############

        values <- reactiveValues()
        DF<-data.frame(SampleID=datshort,Group=(rep("NA",(length(datshort)))),stringsAsFactors = F)
        observe({
          if (!is.null(input$hot)) {
            DF = hot_to_r(input$hot)
          } else {
            if (is.null(values[["DF"]]))
              DF <- DF
            else
              DF <- values[["DF"]]
          }
          values[["DF"]] <- DF
        })

        output$hot <- renderRHandsontable({
          DF <- values[["DF"]]
          if (!is.null(DF))
            rhandsontable(DF, stretchH = "all")
       })

        }
      observeEvent(input$savetable, {
          sampleInfo<-isolate(values[["DF"]])
          sampleInfo<-sampleInfo[!sampleInfo$Group %in% "NA",]
          rownames(sampleInfo)<-sampleInfo$SampleID
          
          ##write table to some file!!!!
          #### report where the file has been written to
          
                  })#end of observe input$savetable
      

        },ignoreInit=TRUE)#end of observe input$selectdataset
  
      ####allow user to select snakemake_workflow version?
   observeEvent(input$selectworkflow, {
     
     inWorkflow<-isolate(input$selectworkflow)
     yamlpath<-file.path("/data/manke/sikora/snakepipes/workflows",inWorkflow,"defaults.yaml") ###update path #,"defaults.yaml"
     paramsYML<-read_yaml(file=yamlpath)
     if (length(paramsYML[["reads"]])>1){
         paramsYML[["reads"]]<-paste(paramsYML[["reads"]],collapse=";")}
     paramsYML[unlist(lapply(paramsYML,is.null))]<-"NA"
     YMLtab<-as.data.frame(do.call(rbind,paramsYML),stringsAsFactors=FALSE)
     colnames(YMLtab)<-"Value"
     output$params<-renderTable(YMLtab,rownames=TRUE)
     
     
     
     
       },ignoreInit=TRUE)#end of observe input$selectworkflow
    
###################################################################################

        output$resultPanels<-renderUI({myTabs<-list(tabPanel(title="Input Data",
                                                      fluidPage(
                                                          fluidRow(
                                                              uiOutput("Read1"),
                                                              uiOutput("Read2")
                                                                   ),
                                                          rHandsontableOutput("hot"),
                                                              actionButton(inputId="savetable",label="Save sample sheet")
                                                            
                                                              )
                                                          ),
                                                    
                                                    tabPanel(title="Workflow parameters",
                                                             fluidPage(
                                                               fluidRow(
                                                                 tableOutput("params")
                                                                       ),
                                                                 actionButton(inputId="saveconfig",label="Save workflow configuration")
                                                               
                                                             )
                                                    )
                                                    
                                                        )

            do.call(tabsetPanel, myTabs)})


}

shinyApp(ui, server)
