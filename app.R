## app.R ##
Rlib="/data/boehm/group/shiny_apps/Rlibs3.4.0"

library(shiny,lib.loc=Rlib)
library(shinydashboard,lib.loc=Rlib)
library(rhandsontable,lib.loc=Rlib)



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
        textInput(inputId="analysistitle", label="Analysis title", value = "", width = NULL, placeholder = NULL),
        uiOutput("configurator"),
        imageOutput("logo"),
        tags$footer("Copyright 2018 MPI-IE Freiburg Bioinfo Core Unit")
        ),
        
    dashboardBody(
        h2("User input"),
        uiOutput("resultPanels")
               
    )

 )}



server <- function(input, output, session) {
  
       library("yaml",lib.loc=Rlib)
       library("stringi",lib.loc=Rlib)
       library("sendmailR",lib.loc=Rlib)
  
       sInfoTOyaml<-function(df){
         df2<-df[!df$ChIPgroup %in% "Input",!colnames(df) %in% c("Group","Read1","SampleID","ChIPgroup")]
         df2$MatchedInput<-as.character(df2$MatchedInput)
         df2$MatchedInput<-paste0("control: ",df2$MatchedInput)
         df2$MarkWidth<-as.character(df2$MarkWidth)
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
  

        values <- reactiveValues()  
        values$datdir<-c()
        values$sInfoDest<-""
        values$chDictDest<-""
        observeEvent(input$adddataset, {
      
        
        if((input$group!="")&(input$owner!="")&(input$projectid!="")&(input$pathtodata=="")){
            inGroup<-isolate(input$group)
            inOwner<-isolate(input$owner)
            inProjectID<-isolate(input$projectid)

            values$datdir<-c(values$datdir,system(sprintf("find /data/%s/sequencing_data -name Project_%s_%s_%s -type d | sort",tolower(inGroup),inProjectID,inOwner,inGroup),intern=TRUE)) 
        }
            
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
                 values$Read1<-grep("*.fastq.gz",values$datPath,value=TRUE)
                 output$Read1<-renderText(values$Read1)
                 values$datshort<-gsub(".fastq.gz","",basename(values$datPath))}
            if(length(unique(values$datshort))<length(values$datshort)){
                values$datwarnings<-"Warning! Your dataset contains multiple instances of identically named samples! If these are resequenced samples, consider merging them before continuing with the analysis or exclude them by leaving NAs in the Group column. Otherwise they will be treated as replicates."
            }else{values$datwarnings<-"All sample names are unique. Proceed with the analysis."}
            output$datawarnings<-renderText(values$datwarnings)
        },ignoreInit=TRUE)#end of observe input$adddataset

            ###############initiate reactive table to collect sample information ###############
            ###use 1 observer to update sample names from added data and onother one to update with manually typed user information

        observe({input$selectworkflow
           values$inWorkflow<-input$selectworkflow
           values$inGroup<-input$group
           
           path_to_exec<-paste0("/data/manke/sikora/snakepipes/workflows/",values$inWorkflow,"/",values$inWorkflow)###add version selection
           indir<-sprintf("/data/processing/bioinfo-core/%s/%s_%s_%s_input_reads",values$inGroup,values$analysisName,values$ranstring,values$inWorkflow)
           link_cmd<- paste0("ln -t ",indir,' -s ',paste(values$Reads,collapse=" "))
           
           outdir<-sprintf("/data/processing/bioinfo-core/%s/%s_%s_%s_OUT",values$inGroup,values$analysisName,values$ranstring,values$inWorkflow)
           
           cp_sInfo_cmd<-sprintf("cp -v %s %s",values$sInfoDest,indir)
           values$sInfo_in<-paste0(indir,"/",basename(values$sInfoDest))
           genome_sel<-c("Zebrafish"="GRCz10","Fission yeast"="SchizoSPombe_ASM294v2","Fruitfly"="dm6","Human"="hs37d5","Mouse"="mm10")  
           values$genome<-genome_sel[input$genome]
           output$from<-renderUI({textInput(inputId="sender",label="Your email address",placeholder="lastname@ie-freiburg.mpg.de")})
           output$freetext<-renderUI({textInput(inputId="comments",label="Your message to the bioinfo facility",placeholder="Sample X might be an outlier.",width="600px")})
          
           path_to_DNA_mapping<-paste0("/data/manke/sikora/snakepipes/workflows/DNA-mapping/DNA-mapping")
           
           
           if(values$inWorkflow=="ATAC-seq"){
               
              values$command<-sprintf("mkdir -p %s ; %s ; %s ;%s -i %s -o %s %s ; %s -d %s --DOC %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_DNA_mapping,indir,outdir,values$genome,path_to_exec,outdir,values$sInfo_in,values$genome)
              output$command<-renderText({ values$command })
              
                     }##end of ATACseq
           
           else if(values$inWorkflow=="ChIP-seq"){
             
             cp_chDict_cmd<-sprintf("cp -v %s %s",values$chDictDest,indir)
             values$chDictr_in<-paste0(indir,"/",basename(values$chDictDest))
             
             values$command<-sprintf("mkdir -p %s ; %s  ; %s ; %s ; %s -i %s -o %s %s ; %s -d %s --DB %s %s %s",indir,link_cmd,cp_sInfo_cmd,cp_chDict_cmd,path_to_DNA_mapping,indir,outdir,values$genome,path_to_exec,outdir,values$sInfo_in,values$genome,values$chDictr_in) 
             output$command<-renderText({ values$command })
             
             
           }##end of ChIP-seq
           
           
           else if(values$inWorkflow=="HiC"){
             
             values$command<-sprintf("mkdir -p %s ; %s ;  %s ; %s -i %s -o %s %s ; %s -i %s -o %s %s",indir,link_cmd,cp_sInfo_cmd,path_to_DNA_mapping,indir,outdir,values$genome,path_to_exec,indir,outdir,values$genome) 
             output$command<-renderText({ values$command })
             
           }##end of HiC
           
           
           else if(values$inWorkflow=="DNA-mapping"){
             
             values$command<-sprintf("mkdir -p %s ; %s ;  %s ; %s -i %s -o %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_DNA_mapping,indir,outdir,values$genome) 
             output$command<-renderText({ values$command })
             
           } #end of DNA-mapping
           
           else if(values$inWorkflow=="RNA-seq"){
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s ; %s -i %s -o %s --DE %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_exec,indir,outdir,values$sInfo_in,values$genome) 
             output$command<-renderText({ values$command })
             
           } #end of RNA-seq
           
           else if(values$inWorkflow=="WGBS"){
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s ; %s -ri %s -w %s --sampleInfo %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_exec,indir,outdir,values$sInfo_in,values$genome) 
             output$command<-renderText({ values$command })
             
           } #end of WGBS


        #observe({
        #      input$adddataset
        #      input$selectworkflow
              
              if(values$inWorkflow=="ChIP-seq"){
              values$DF<-data.frame(SampleID=values$datshort,Group=(factor(rep("NA",(length(values$datshort))),levels=c("Control","Treatment"))),Read1=values$Read1,ChIPgroup=factor(rep("NA",(length(values$datshort))),levels=c("ChIP","Input"),ordered=TRUE),MatchedInput=factor(rep("NA",(length(values$datshort))),levels=unique(values$datshort),ordered=TRUE),MarkWidth=factor(rep("NA",(length(values$datshort))),levels=c("Broad","Narrow"),ordered=TRUE),stringsAsFactors = F)} 

              else if(values$inWorkflow=="HiC"){
              values$DF<-data.frame(SampleID=values$datshort,Group=(rep("NA",(length(values$datshort)))),Read1=values$Read1,stringsAsFactors = F)}


               else {

               values$DF<-data.frame(SampleID=values$datshort,Group=(factor(rep("NA",(length(values$datshort))),levels=c("Control","Treatment"))),Read1=values$Read1,stringsAsFactors = F)
               }

       # })
        })#end of observe input$selectworkflow 


        observe({
          if(!is.null(input$hot))
            values$DF <- hot_to_r(input$hot)
        })
        
        output$hot <- renderRHandsontable({
          rhandsontable(values$DF)#, stretchH = "all"
        })

 
      #})#end of observe input$selectworkflow 
     
         
      observeEvent(input$savetable, {
          values$analysisName<-gsub("[^[:alnum:]]", "_", isolate(input$analysistitle))
          values$ranstring<-stri_rand_strings(n=1,length=8)
          
          sampleInfo<-isolate(values$DF)
          ###check for replicates, else issue a warning
          if(sum(is.na(sampleInfo$Group))<length(sampleInfo$Group)){
              sampleInfo$Group<-as.character(sampleInfo$Group)
              sampleInfo<-sampleInfo[!is.na(sampleInfo$Group),]
              rownames(sampleInfo)<-make.names(sampleInfo$SampleID, unique=TRUE)
              sampleInfo<-sampleInfo[order(sampleInfo$Group),]
              gl<-min(ave(1:nrow(sampleInfo),sampleInfo$Group,FUN=length))
              if(gl==1){
                values$datwarnings<-c(values$datwarnings,"SOME SAMPLE GROUPS DON'T HAVE REPLICATES!!!")
              }
              else if (gl==2){
                values$datwarnings<-c(values$datwarnings,"Some of your sample groups have only 2 replicates. This might be suboptimal for some analyses and lead to higher FDR.")
              }
              else if (gl==3){
                values$datwarnings<-c(values$datwarnings,"All sample groups have at least 3 replicates.")
              }
          }else {rownames(sampleInfo)<-make.names(sampleInfo$SampleID, unique=TRUE)}    
          # 
          ###check if ChIPseq is selected, if yes, create a yaml
          if(values$inWorkflow=="ChIP-seq"){
          chip_dict<-sInfoTOyaml(sampleInfo)
          values$chDictDest<-sprintf("/data/manke/group/shiny/snakepipes_input/%s_%s_samples.yaml",values$ranstring,values$analysisName)
          write_yaml(noquote(chip_dict),file=isolate(values$chDictDest))
          
          }
          
          #update read1 and read2 stored in reactive values
          if (values$ispaired){
            vRead1<-sampleInfo$Read1
            vRead2<-gsub("_R1.fastq.gz","_R2.fastq.gz",sampleInfo$Read1)
            values$Reads<-c(vRead1,vRead2)
          } else {
            vRead1<-sampleInfo$Read1
            values$Reads<-vRead1
          }
          
          
          values$sInfoDest<-sprintf("/data/manke/group/shiny/snakepipes_input/%s_%s_sampleSheet.tsv",values$ranstring,values$analysisName)
          
          if(values$inWorkflow!="WGBS"){colnames(sampleInfo)[1:2]<-c("sample","condition")}
          
          write.table(sampleInfo,file=values$sInfoDest,sep="\t",quote=FALSE)
          output$sIsaved<-renderText("Sample sheet saved.")
          
                  })#end of observe input$savetable
      
      output$configurator<-renderUI({tagList(
        selectInput(inputId="genome", label="Select organism", choices=c("Zebrafish","Fission yeast","Fruitfly","Human","Mouse"), selected = NULL)
      )})
       
       
   
             observeEvent(input$savesubmit, {
               bshscript<-sprintf("/data/manke/group/shiny/snakepipes_input/%s_%s_script.sh",values$ranstring,values$analysisName)
               fileConn<-file(bshscript)
               shebang<-"#!/bin/bash"
               exports<-"export PATH=/data/processing/conda/bin:$PATH"
               writeLines(c(shebang,exports,unlist(strsplit(isolate(values$command),split=";"))), fileConn)
               close(fileConn)
               cc<-isolate(input$sender)
               from<-sprintf("<sendmailR@@\\%s>", Sys.info()[4])
               to<-"<sikora@ie-freiburg.mpg.de>"  ##change to "bioinfo-core@ie-freiburg.mpg.de"
               subject<-paste0("Analysis request ",isolate(input$analysistitle), "_" ,isolate(values$ranstring))
               msg <- gsub(";","\n \n",paste0(cc," has requested the following analysis: \n \n", isolate(values$inWorkflow)," \n \n ", values$genome  ," \n \n User comments: \n \n" ,isolate(input$comments),"\n \n Please review the input files and the attached sample sheet before proceeding. Consider treating batch effect if pooling data across sequencing project, especially in case of RNAseq. \n \n End of message. \n \n ",paste(rep("#",times=80),collapse="")," \n \n ", values$command ,"\n \n"))
               if(values$inWorkflow=="ChIP-seq"){
                 sendmail(from=sprintf("<%s>",from), to=to, subject=subject, msg=list(msg,mime_part(isolate(values$sInfoDest)),mime_part(isolate(values$chDictDest)),mime_part(bshscript)),cc=sprintf("<%s>",cc))
               }
               else{
               sendmail(from=sprintf("<%s>",from), to=to, subject=subject, msg=list(msg,mime_part(isolate(values$sInfoDest)),mime_part(bshscript)),cc=sprintf("<%s>",cc))}
               output$eSent<-renderText("Your request has been sent to the MPI-IE Bioinfo facility. A copy was sent to your email address.")

               
             },ignoreInit=TRUE,once=TRUE)#end of observe input$savesubmit
             
 
###################################################################################
             
        output$logo<-renderImage({list(src="/data/manke/sikora/shiny_apps/userIN_to_yaml/MPIIE_logo_sRGB.jpg",width=100,height=100)},deleteFile =FALSE)
             

        output$resultPanels<-renderUI({myTabs<-list(tabPanel(title="Input Data",
                                                      fluidPage(
                                                          box(textOutput("ispaired"),width=12),
                                                          rHandsontableOutput("hot"),
                                                          box(textOutput("datawarnings")),
                                                              actionButton(inputId="savetable",label="Save sample sheet",style="background-color: green"),
                                                          box(textOutput("sIsaved"))
                                                            
                                                              )
                                                          ),
                                                    
                                                    tabPanel(title="User information",
                                                             fluidPage(
                                                               fluidRow(
                                                                 uiOutput("from"),
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
