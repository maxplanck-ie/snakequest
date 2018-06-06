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
       require("sendmailR",lib.loc=Rlib)

        values <- reactiveValues()  
        values$datdir<-c()
        observeEvent(input$adddataset, {
      
        
        if((input$group!="")&(input$owner!="")&(input$projectid!="")&(input$pathtodata=="")){
            inGroup<-isolate(input$group)
            inOwner<-isolate(input$owner)
            inProjectID<-isolate(input$projectid)

            values$datdir<-c(values$datdir,system(sprintf("find /data/%s/sequencing_data -name Project_%s_%s_%s -type d | sort",tolower(inGroup),inProjectID,inOwner,inGroup),intern=TRUE)) }
        else if ((input$group=="")&(input$owner=="")&(input$projectid=="")&(input$pathtodata!="")){
            values$datdir<-c(values$datdir,isolate(input$pathtodata))
        }  
            #handle single and paired end data ###allow for multiple sequencing runs
            values$datPath<-dir(values$datdir,pattern="*.fastq.gz",full.names=TRUE,recursive=TRUE)
            #get rid of optical duplicates
            values$datPath<-grep("*optical*",values$datPath,value=TRUE,invert=TRUE)
                                           
        
            if(sum(grepl("*_R2.fastq.gz",values$datPath))>0){values$ispaired<-TRUE
                    output$ispaired<-renderText("Paired reads detected.")
                    values$Read1<-grep("*_R1.fastq.gz",values$datPath,value=TRUE)
                    values$Read2<-grep("*_R2.fastq.gz",values$datPath,value=TRUE)
                    output$Read1<-renderText(values$Read1)
                    output$Read2<-renderText(values$Read2)
                    values$datshort<-gsub("_R1.fastq.gz","",basename(grep("*_R1.fastq.gz",values$datPath,value=TRUE)))}
            else{values$ispaired<-FALSE
                 output$ispaired<-renderText("Single reads detected.")
                 values$Read1<-grep("*_R1.fastq.gz",values$datPath,value=TRUE)
                 output$Read1<-renderText(values$Read1)
                 values$datshort<-gsub(".fastq.gz","",basename(values$datPath))}
            if(length(unique(values$datshort))<length(values$datshort)){output$datawarnings<-renderText("Warning! Your dataset contains multiple instances of identically named samples! If these are resequenced samples, consider merging them before continuing with the analysis or exclude them by leaving NAs in the Group column. Otherwise they will be treated as replicates.")}else{output$datawarnings<-renderText("All sample names are unique. Proceed with the analysis.")}
        },ignoreInit=TRUE)#end of observe input$adddataset

            ###############initiate reactive table to collect sample information ###############
            ###use 1 observer to update sample names from added data and onother one to update with manually typed user information
        
        observe({
              input$adddataset
              values$DF<-data.frame(SampleID=values$datshort,Group=(rep("NA",(length(values$datshort)))),Read1=values$Read1,stringsAsFactors = F)
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
          values$ranstring<-stri_rand_strings(n=1,length=8)
          
          sampleInfo<-isolate(values[["DF"]])
          sampleInfo<-sampleInfo[!sampleInfo$Group %in% "NA",]
          
          #update read1 and read2 stored in reactive values
          if (values$ispaired){
            values$Read1<-sampleInfo$Read1
            values$Read2<-sampleInfo$Read2
            values$Reads<-c(values$Read1,values$Read2)
          } else {
            values$Read1<-sampleInfo$Read1
            values$Reads<-values$Read1
          }
          
          rownames(sampleInfo)<-sampleInfo$SampleID
          values$sInfoDest<-sprintf("/data/manke/group/shiny/snakepipes_input/%s_%s_sampleSheet.tsv",values$ranstring,analysisName)
          
          write.table(sampleInfo,file=values$sInfoDest,sep="\t",quote=FALSE)
          output$sIsaved<-renderText("Sample sheet saved.")
          
                  })#end of observe input$savetable
       
       
   observe({input$selectworkflow
           values$inWorkflow<-input$selectworkflow
           output$configurator<-renderUI({tagList(
             selectInput(inputId="genome", label="Select organism", choices=c("Zebrafish","Fission yeast","Fruitfly","Human","Mouse"), selected = NULL),
             selectInput(inputId="version", label="Select workflow version", choices=c("Placeholder"), selected = "Placeholder")
           )})
           path_to_exec<-paste0("/data/manke/sikora/snakepipes/workflows/",values$inWorkflow,"/",values$inWorkflow)###add version selection
           indir<-sprintf("/data/processing/bioinfo-core/%s/%s_input_reads",inGroup,input$analysistitle,values$ranstring,values$inWorkflow)
           link_cmd<- paste0("ln -t ",indir,' -s ')
           outdir<-sprintf("/data/processing/bioinfo-core/%s/%s_%s_%s_OUT",inGroup,input$analysistitle,values$ranstring,values$inWorkflow)
           genome_sel<-c("Zebrafish"="GRCz10","Fission yeast"="SchizoSPombe_ASM294v2","Fruitfly"="dm6","Human"="hs37d5","Mouse"="mm10")             
           output$from<-renderUI({textInput(inputId="sender",label="Your email address",placeholder="lastname@ie-freiburg.mpg.de")})
           output$freetext<-renderUI({textInput(inputId="comments",label="Your message to the bioinfo facility",placeholder="I can't find my data.",width="600px")})
           output$freetext<-renderUI({textInput(inputId="comments",label="Your message to the bioinfo facility",placeholder="I can't find my data.",width="600px")})
           path_to_DNA_mapping<-paste0("/data/manke/sikora/snakepipes/workflows/DNA-mapping/DNA-mapping")
           
           
           if(values$inWorkflow=="ATAC-seq"){
               
              values$command<-sprintf("mkdir -p %s ; %s ; %s -i %s -o %s %s ; %s -d %s --DOC %s %s ",indir,link_cmd,path_to_DNA_mapping,indir,outdir,genome_sel[input$genome],path_to_exec,outdir,values$sInfoDest,genome_sel[input$genome])
              output$command<-renderText({ values$command })
              
                     }##end of ATACseq
           
           else if(values$inWorkflow=="ChIP-seq"){
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s -i %s -o %s %s ; %s -d %s --DB %s %s",indir,link_cmd,path_to_DNA_mapping,indir,outdir,genome_sel[input$genome],path_to_exec,outdir,values$sInfoDest,genome_sel[input$genome]) ##missing: samples.yaml after genome
             output$command<-renderText({ values$command })

           }##end of ChIP-seq
           
           
           else if(values$inWorkflow=="HiC"){
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s -i %s -o %s %s ; %s -i %s -o %s %s",indir,link_cmd,path_to_DNA_mapping,indir,outdir,genome_sel[input$genome],path_to_exec,indir,outdir,genome_sel[input$genome]) 
             output$command<-renderText({ values$command })
             
           }##end of HiC
           
           
           else if(values$inWorkflow=="DNA-mapping"){
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s -i %s -o %s %s ",indir,link_cmd,path_to_DNA_mapping,indir,outdir,genome_sel[input$genome]) 
             output$command<-renderText({ values$command })
             
           } #end of DNA-mapping
           
           else if(values$inWorkflow=="RNA-seq"){
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s -i %s -o %s --DE %s %s ",indir,link_cmd,path_to_exec,indir,outdir,values$sInfoDest,genome_sel[input$genome]) 
             output$command<-renderText({ values$command })
             
           } #end of RNA-seq
           
           else if(values$inWorkflow=="WGBS"){
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s -ri %s -w %s --sampleInfo %s %s ",indir,link_cmd,path_to_exec,indir,outdir,values$sInfoDest,genome_sel[input$genome]) 
             output$command<-renderText({ values$command })
             
           }
           

           
      })#end of observe input$selectworkflow 
     
             observeEvent(input$savesubmit, {
               cc<-isolate(input$sender)
               from<-sprintf("<sendmailR@@\\%s>", Sys.info()[4])
               to<-"<sikora@ie-freiburg.mpg.de>"  ##change to "bioinfo-core@ie-freiburg.mpg.de"
               subject<-paste0("Analysis request ",isolate(input$analysistitle), "_" ,isolate(values$ranstring))
               msg <- paste0(cc," has requested the following analysis: \n \n", isolate(values$command),  " \n \n User comments: /n /n" ,isolate(input$comments),"\n \n End of message.")
               sendmail(from=sprintf("<%s>",from), to=to, subject=subject, msg=msg,cc=sprintf("<%s>",cc))
               output$eSent<-renderText("Your request has been sent to the MPI-IE Bioinfo facility. A copy was sent to your email address.")

               
             },ignoreInit=TRUE,once=TRUE)#end of observe input$savesubmit
             
 
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
                                                                 uiOutput("configurator"),
                                                                 uiOutput("from"),
                                                                 box(textOutput("command")),
                                                                 uiOutput("freetext")
                                                                       ),
                                                                 actionButton(inputId="savesubmit",label="Save workflow configuration and submit request"),
                                                               box(textOutput("eSent"))

                                                             )
                                                    )
                                                    
                                                        )

            do.call(tabsetPanel, myTabs)})


}

shinyApp(ui, server)
