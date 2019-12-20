## app.R ##
Rlib="/data/manke/sikora/shiny_apps/Rlibs3.5.0_bioc3.7"
.libPaths(Rlib)

library(shiny,lib.loc=Rlib)
library(shinydashboard,lib.loc=Rlib)
library(rhandsontable,lib.loc=Rlib)
library(shinyBS,lib.loc=Rlib)

ui <- function(request) {dashboardPage(
    dashboardHeader(title = "Dataset selection"),
    ## Sidebar content
    dashboardSidebar(

        selectInput(inputId="selectworkflow",label="Select NGS workflow",choices=c("PLEASE SELECT WORKFLOW","ATAC-seq","ChIP-seq","DNA-mapping","HiC","RNA-seq","WGBS")),
        textInput(inputId="analysistitle", label="Analysis title", value = "", width = NULL, placeholder = NULL),
        selectInput(inputId="genome", label="Select organism", choices=c("PLEASE SELECT A GENOME","Zebrafish [zv10]","Fission yeast","Fruitfly [dm6]","Fruitfly [dm3]","Human [hg37]","Human [hg38]","Mouse [mm9]","Mouse [mm10]"), selected = NULL),
        selectInput(inputId="selectformat",label="Select input file format",choices=c("fastq.gz","bam")),
        textInput(inputId="group", label="Group", value = "", width = NULL, placeholder = NULL),
        textInput(inputId="owner", label="Project Owner", value = "", width = NULL, placeholder = NULL),
        textInput(inputId="projectid", label="Project ID", value = "", width = NULL, placeholder = NULL),
        textInput(inputId="pathtodata", label="Data folder", value = "", width = NULL, placeholder = NULL),
        actionButton(inputId="adddataset", label="Add dataset"),
        imageOutput("logo"),
        tags$footer("Copyright 2018 MPI-IE Freiburg Bioinfo Core Unit"),
        bsTooltip(id="group", title="Enter group/department PI name as specified in the sequencing request.", placement = "right", trigger = "hover"),
        bsTooltip(id="owner", title="Enter the name of data owner as specified in the sequencing request.", placement = "right", trigger = "hover"),
        bsTooltip(id="projectid", title="Enter the sequencing project ID/number you have received from the sequencing facility.", placement = "right", trigger = "hover"),
        bsTooltip(id="pathtodata", title="Paste the path to the folder containg the reads you would like to analyze.", placement = "right", trigger = "hover"),
        bsTooltip(id="adddataset", title="Add a set of reads identified using the information above.", placement = "right", trigger = "hover"),
        bsTooltip(id="analysistitle", title="Provide an optional analysis title.", placement = "right", trigger = "hover")
        ),
        
    dashboardBody(
        h2("User input"),
        uiOutput("resultPanels")
               
    )

 )}


server <- function(input, output, session) {

##########load packages  
       library("yaml",lib.loc=Rlib)
       library("stringi",lib.loc=Rlib)
       library("sendmailR",lib.loc=Rlib)

#########define functions  
       sInfoTOyaml<-function(df){
         df2<-df[!df$ChIPgroup %in% "Input",!colnames(df) %in% c("Group","Read1","SampleID","ChIPgroup","Merge")]
         df2$MatchedInput<-as.character(df2$MatchedInput)
         df2$MatchedInput<-paste0("control: ",df2$MatchedInput)
         df2$MarkWidth<-as.character(df2$MarkWidth)
         df2$MarkWidth[!df2$MarkWidth %in% c("Broad","Narrow")]<-"Narrow"
         df2$MarkWidth[grep("Broad",df2$MarkWidth)]<-noquote("broad: True")
         df2$MarkWidth[grep("Narrow",df2$MarkWidth)]<-noquote("broad: False")
         df3<-as.data.frame(t(df2),stringsAsFactors=FALSE)
         l3<-list(df3)
         names(l3)<-"chip_dict"
         l3y<-as.yaml(l3,omap=TRUE)
         l3ymod<-gsub("'","",gsub("- ","",gsub("!omap","",l3y)))
         l4<-yaml.load(l3ymod)
         return(l4)
       }
  

########init reactive values
        values <- reactiveValues()  
        values$datdir<-c()
        values$sInfoDest<-""
        values$chDictDest<-""
        values$genome<-""

####### constant observer on the project input
        observe({
          if ((input$group!="")&(input$owner!="")&(input$projectid!="")&(input$pathtodata!="")){
            showModal(modalDialog(
              title = "REDUNDANT INFORMATION PROVIDED",
              "Please make sure to either fill in the Group-Owner-Project fields AND leave the folder path empty, or the reverse!",
              easyClose = TRUE
            ))} 
           
        })


#######observe add dataset
       observeEvent(input$adddataset, {
        
        dsel<-c("fastq.gz"="Project","bam"="Analysis")
        psel<-c("fastq.gz"="*fastq.gz$","bam"="*.bam$") 
        inFormat<-isolate(input$selectformat)
        
        if((input$group!="")&(input$owner!="")&(input$projectid!="")&(input$pathtodata=="")){
            inGroup<-isolate(input$group)
            inOwner<-isolate(input$owner)
            inProjectID<-isolate(input$projectid)


            values$datdir<-c(values$datdir,system(sprintf("find /data/%s/sequencing_data -path /data/%s/sequencing_data/OxfordNanopore -prune -o -name %s_%s_%s_%s -type d | sort",tolower(gsub("-.+","",inGroup)),tolower(gsub("-.+","",inGroup)),dsel[inFormat],inProjectID,inOwner,inGroup),intern=TRUE)) 
            
        }
            
        else if ((input$group=="")&(input$owner=="")&(input$projectid=="")&(input$pathtodata!="")){
            values$datdir<-c(values$datdir,isolate(input$pathtodata))
                }  
            #handle single and paired end data ###allow for multiple sequencing runs
        
        
            values$datPath<-dir(values$datdir,pattern=psel[inFormat],full.names=TRUE,recursive=TRUE)
            
            if(!isTruthy(values$datdir)){
              showModal(modalDialog(
                title = "DATA NOT FOUND",
                "Please verify the spelling of Group and Owner names (first capital letter of the last name) or the path you provided!",
                easyClose = TRUE
              ))}
            req(values$datdir)
            
            
            if(inFormat %in% "bam"){
              values$ispaired<-FALSE
              values$Read1<-grep(".+.bam$",values$datPath,value=TRUE,perl=TRUE)
              values$datshort<-gsub(".bam","",gsub(".filtered.bam","",gsub(".PCRrm.bam","",basename(values$datPath))))                          
            }
            else{
                #get rid of optical duplicates
                values$datPath<-grep("*optical*",values$datPath,value=TRUE,invert=TRUE)
                if(sum(grepl(".+_R*2.fastq.gz",values$datPath,perl=TRUE))>0){values$ispaired<-TRUE
                        output$ispaired<-renderText("Paired reads detected.")
                        values$Read1<-grep(".+_R*1.fastq.gz",values$datPath,value=TRUE,perl=TRUE)
                        values$Read2<-grep(".+_R*2.fastq.gz",values$datPath,value=TRUE,perl=TRUE)
                        output$Read1<-renderText(values$Read1)
                        output$Read2<-renderText(values$Read2)
                        values$datshort<-gsub("_R*1.fastq.gz","",basename(grep(".+_R*1.fastq.gz",values$Read1,value=TRUE,perl=TRUE)),perl=TRUE)}
                else{values$ispaired<-FALSE
                     output$ispaired<-renderText("Single reads detected.")
                     values$Read1<-grep("*.fastq.gz",values$datPath,value=TRUE)
                     output$Read1<-renderText(values$Read1)
                     values$datshort<-gsub(".fastq.gz","",basename(values$Read1))}
            }
            if(length(unique(values$datshort))<length(values$datshort)){
                values$datwarnings<-"Warning! Your dataset contains multiple instances of identically named samples! If these are resequenced samples, consider merging them before continuing with the analysis or exclude them by leaving NAs in the Group column. Otherwise they will be overwritten."
            }else{values$datwarnings<-"All sample names are unique. Proceed with the analysis."}
            output$datawarnings<-renderText(values$datwarnings)
        },ignoreInit=TRUE)#end of observe input$adddataset
        


