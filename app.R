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
        textInput(inputId="pathtodata", label="Data folder", value = "", width = NULL, placeholder = NULL),
        actionButton(inputId="adddataset", label="Add dataset"),
        selectInput(inputId="selectworkflow",label="Select NGS workflow",choices=c("ATAC-seq","ChIP-seq","DNA-mapping","HiC","RNA-seq","WGBS")),
        textInput(inputId="analysistitle", label="Analysis title", value = "", width = NULL, placeholder = NULL)
        
        ),
        
    dashboardBody(
        h2("User input"),
        uiOutput("resultPanels")
               
    )

 )}



server <- function(input, output, session) {
        
       require("yaml",lib.loc=Rlib)
       require("stringi",lib.loc=Rlib)

        values <- reactiveValues()  
        values$datdir<-c()
        observeEvent(input$adddataset, {
      
        
        if((input$group!="")&(input$owner!="")&(input$projectid!="")&(input$pathtodata=="")){
            inGroup<-isolate(input$group)
            inOwner<-isolate(input$owner)
            inProjectID<-isolate(input$projectid)

            values$datdir<-c(values$datdir,system(sprintf("find /data/%s/sequencing_data -name Project_%s_%s_%s -type d",tolower(inGroup),inProjectID,inOwner,inGroup),intern=TRUE)) }
        else if ((input$group=="")&(input$owner=="")&(input$projectid=="")&(input$pathtodata!="")){
            values$datdir<-c(values$datdir,isolate(input$pathtodata))
        }  
            #handle single and paired end data ###allow for multiple sequencing runs
            values$datPath<-dir(values$datdir,pattern="*.fastq.gz",full.names=TRUE,recursive=TRUE)
            #get rid of optical duplicates
            values$datPath<-grep("*optical*",values$datPath,value=TRUE,invert=TRUE)
                                           
        
            if(sum(grepl("*_R2.fastq.gz",values$datPath))>0){output$ispaired<-renderText("Paired reads detected.")
                    output$Read1<-renderText(grep("*_R1.fastq.gz",values$datPath,value=TRUE))
                    output$Read2<-renderText(grep("*_R2.fastq.gz",values$datPath,value=TRUE))
                    values$datshort<-gsub("_R1.fastq.gz","",basename(grep("*_R1.fastq.gz",values$datPath,value=TRUE)))}
            else{output$ispaired<-renderText("Single reads detected.")
                 output$Read1<-renderText(values$datPath)
                 values$datshort<-gsub(".fastq.gz","",basename(values$datPath))}
            if(length(unique(values$datshort))<length(values$datshort)){output$datawarnings<-renderText("Warning! Your dataset contains multiple instances of identically named samples! If these are resequenced samples, consider merging them before continuing with the analysis or exclude them by leaving NAs in the Group column. Otherwise they will be treated as replicates.")}else{output$datawarnings<-renderText("All sample names are unique. Proceed with the analysis.")}
        },ignoreInit=TRUE)#end of observe input$adddataset

            ###############initiate reactive table to collect sample information ###############
            ###use 1 observer to update sample names from added data and onother one to update with manually typed user information
        
        observe({
              input$adddataset
              values$DF<-data.frame(SampleID=values$datshort,Group=(rep("NA",(length(values$datshort)))),stringsAsFactors = F)
        })
            
        observe({
          if(!is.null(input$hot))
            values$DF <- hot_to_r(input$hot)
        })
        
        output$hot <- renderRHandsontable({
          rhandsontable(values$DF, stretchH = "all")
        })

        
      observeEvent(input$savetable, {
          analysisName<-isolate(input$analysistitle)
          ranstring<-stri_rand_strings(n=1,length=8)
          
          sampleInfo<-isolate(values[["DF"]])
          sampleInfo<-sampleInfo[!sampleInfo$Group %in% "NA",]
          rownames(sampleInfo)<-sampleInfo$SampleID
          
          write.table(sampleInfo,file=sprintf("/data/manke/group/shiny/snakepipes_input/%s_%s_sampleSheet.csv",ranstring,analysisName),sep="\t",quote=FALSE)
          output$sIsaved<-renderText("Sample sheet saved.")
          
                  })#end of observe input$savetable
       
       
        
  
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
     
             observeEvent(input$saveconfig, {
               
               
               
               
             },ignoreInit=TRUE)#end of observe input$saveconfig
     
       },ignoreInit=TRUE)#end of observe input$selectworkflow
    
###################################################################################

        output$resultPanels<-renderUI({myTabs<-list(tabPanel(title="Input Data",
                                                      fluidPage(
                                                          box(textOutput("ispaired"),width=12),
                                                          box(textOutput("Read1"),title="Read1",width=12,collapsible=TRUE),
                                                          #fluidRow(box(textOutput("Read2"),title="Read2")),
                                                          rHandsontableOutput("hot"),
                                                          box(textOutput("datawarnings")),
                                                              actionButton(inputId="savetable",label="Save sample sheet",style="background-color: green"),
                                                          box(textOutput("sIsaved"))
                                                            
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
