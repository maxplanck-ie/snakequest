## app.R ##
Rlib="/data/boehm/group/shiny_apps/Rlibs3.4.0"

library(shiny,lib.loc=Rlib)
library(shinydashboard,lib.loc=Rlib)
library(rhandsontable,lib.loc=Rlib)
library(shinyBS,lib.loc=Rlib)


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
        selectInput(inputId="genome", label="Select organism", choices=c("PLEASE SELECT A GENOME","Zebrafish [zv10]","Fission yeast","Fruitfly [dm6]","Fruitfly [dm3]","Human [hg37]","Human [hg38]","Mouse [mm9]","Mouse [mm10]"), selected = NULL),
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
        values$genome<-""
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
                                           
        
            if(sum(grepl(".+_R*2.fastq.gz",values$datPath,perl=TRUE))>0){values$ispaired<-TRUE
                    output$ispaired<-renderText("Paired reads detected.")
                    values$Read1<-grep(".+_R*1.fastq.gz",values$datPath,value=TRUE,perl=TRUE)
                    values$Read2<-grep(".+_R*2.fastq.gz",values$datPath,value=TRUE,perl=TRUE)
                    output$Read1<-renderText(values$Read1)
                    output$Read2<-renderText(values$Read2)
                    values$datshort<-gsub("_R*1.fastq.gz","",basename(grep(".+_R*1.fastq.gz",values$datPath,value=TRUE,perl=TRUE)),perl=TRUE)}
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

        observe({#input$selectworkflow
          #validate(
           # need(input$genome != "PLEASE SELECT A GENOME", "Please select a genome")
          #)
           values$inWorkflow<-input$selectworkflow
           
           path_to_exec<-paste0("module load snakePipes; ",values$inWorkflow)###add version selection
           topdir<-sprintf("/data/processing/bioinfo-core/requests/%s_%s_%s",values$analysisName,values$ranstring,values$inWorkflow)
           indir<-sprintf("%s/fastq",topdir)
           link_cmd<- paste0("ln -t ",indir,' -s ',paste(values$Reads,collapse=" "))
           
           outdir<-sprintf("%s/analysis",topdir)
           
           cp_sInfo_cmd<-sprintf("cp -v %s %s",values$sInfoDest,topdir)
           values$sInfo_in<-paste0(topdir,"/",basename(values$sInfoDest))
           genome_sel<-c("PLEASE SELECT A GENOME"="NONE","Zebrafish [zv10]"="GRCz10","Fission yeast"="SchizoSPombe_ASM294v2","Fruitfly [dm6]"="dm6","Fruitfly [dm3]"="dm3","Human [hg37]"="hs37d5","Human [hg38]"="hg38","Mouse [mm9]"="mm9","Mouse [mm10]"="mm10")  #"PLEASE SELECT A GENOME"="NONE",
           values$genome<-genome_sel[input$genome]
           output$from<-renderUI({textInput(inputId="sender",label="Your email address",placeholder="lastname@ie-freiburg.mpg.de")})
           output$freetext<-renderUI({textInput(inputId="comments",label="Your message to the bioinfo facility",placeholder="Sample X might be an outlier.",width="600px")})
          
           path_to_DNA_mapping<-paste0("module load snakePipes; DNA-mapping")
           
           
           if(values$inWorkflow=="ATAC-seq"){
               
              values$command<-sprintf("mkdir -p %s ; %s ; %s ;%s -i %s -o %s %s ; %s -d %s --DOC %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_DNA_mapping,indir,outdir,values$genome,path_to_exec,outdir,values$sInfo_in,values$genome)
              output$command<-renderText({ values$command })
              
                     }##end of ATACseq
           
           else if(values$inWorkflow=="ChIP-seq"){
             
             cp_chDict_cmd<-sprintf("cp -v %s %s",values$chDictDest,topdir)
             values$chDictr_in<-paste0(topdir,"/",basename(values$chDictDest))
             
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
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s ; %s -i %s -o %s --sampleInfo %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_exec,indir,outdir,values$sInfo_in,values$genome) 
             output$command<-renderText({ values$command })
             
           } #end of WGBS


        #observe({
        #      input$adddataset
        #      input$selectworkflow
              
              if(values$inWorkflow=="ChIP-seq"){
              values$DF<-data.frame(SampleID=values$datshort,Group=(factor(rep("NA",(length(values$datshort))),levels=c("Control","Treatment","NA"))),Read1=values$Read1,ChIPgroup=factor(rep("NA",(length(values$datshort))),levels=c("ChIP","Input"),ordered=TRUE),MatchedInput=factor(rep("NA",(length(values$datshort))),levels=unique(values$datshort),ordered=TRUE),MarkWidth=factor(rep("NA",(length(values$datshort))),levels=c("Broad","Narrow"),ordered=TRUE),stringsAsFactors = F)} 

              else if(values$inWorkflow=="HiC"){
              values$DF<-data.frame(SampleID=values$datshort,Group=(rep("NA",(length(values$datshort)))),Read1=values$Read1,stringsAsFactors = F)}
           
              else if(values$inWorkflow=="WGBS"){
              values$DF<-data.frame(SampleID=values$datshort,Group=(factor(rep("NA",(length(values$datshort))),levels=c("Control","Treatment","WT","Mut","NA"))),PlottingID=values$datshort,Read1=values$Read1,stringsAsFactors = F)}


               else {

               values$DF<-data.frame(SampleID=values$datshort,Group=(factor(rep("NA",(length(values$datshort))),levels=c("Control","Treatment","NA"))),Read1=values$Read1,stringsAsFactors = F)
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
          if(sum(is.na(sampleInfo$Group),sampleInfo$Group %in% "NA")<length(sampleInfo$Group)){
              sampleInfo$Group<-as.character(sampleInfo$Group)
              sampleInfo<-sampleInfo[!is.na(sampleInfo$Group)&!sampleInfo$Group %in% "NA",]
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
          
          if(values$inWorkflow!="WGBS"){colnames(sampleInfo)[1:2]<-c("name","condition")}
          
          write.table(sampleInfo,file=values$sInfoDest,sep="\t",quote=FALSE)
          output$sIsaved<-renderText("Sample sheet saved.")
          
                  })#end of observe input$savetable
      
      
      output$datarequests<- renderUI({tagList(
        checkboxInput(inputId="merge", label="I want to request sample merging.", value = FALSE, width = NULL),
        checkboxInput(inputId="beff", label="I expect batch effect in my data.", value = FALSE, width = NULL)
      )})
   
             observeEvent(input$savesubmit, {
               bshscript<-sprintf("/data/manke/group/shiny/snakepipes_input/%s_%s_script.sh",values$ranstring,values$analysisName)
               fileConn<-file(bshscript)
               shebang<-"#!/bin/bash"
               exports<-"export PATH=/data/processing/conda/bin:$PATH"
               writeLines(c(shebang,exports,unlist(strsplit(isolate(values$command),split=";"))), fileConn)
               close(fileConn)
               merge_request<-ifelse(isolate(input$merge),"I want to request sample merging.","No sample merging is needed.")
               b_eff_request<-ifelse(isolate(input$beff),"I expect batch effect in my data.","No batch effect is expected.")
               cc<-isolate(input$sender)
               from<-sprintf("<sendmailR@%s>", Sys.info()[4])
               to<-"<bioinfo-core@ie-freiburg.mpg.de>"  
               subject<-paste0("Analysis request ",isolate(input$analysistitle), "_" ,isolate(values$ranstring))
               msg <- gsub(";","\n \n",paste0(cc," has requested the following analysis: \n \n Workflow: ", isolate(values$inWorkflow)," \n \n Genome: ", values$genome," \n \n ", merge_request," \n \n ", b_eff_request ," \n \n User comments: \n \n", isolate(input$comments),"\n \n Please review the input files and the attached sample sheet before proceeding. \n \n End of message. \n \n ",paste(rep("#",times=80),collapse="")," \n \n ", values$command ,"\n \n"))
               if(values$inWorkflow=="ChIP-seq"){
                 sendmail(from=sprintf("<%s>",from), to=to, subject=subject, control=list(smtpServer="mail.ie-freiburg.mpg.de"), msg=list(msg,mime_part(isolate(values$sInfoDest)),mime_part(isolate(values$chDictDest)),mime_part(bshscript)),cc=sprintf("<%s>",cc))
               }
               else{
               sendmail(from=sprintf("<%s>",from), to=to, subject=subject, control=list(smtpServer="mail.ie-freiburg.mpg.de"), msg=list(msg,mime_part(isolate(values$sInfoDest)),mime_part(bshscript)),cc=sprintf("<%s>",cc))}
               output$eSent<-renderText("Your request has been sent to the MPI-IE Bioinfo facility. A copy was sent to your email address.")

               
             },ignoreInit=TRUE,once=TRUE)#end of observe input$savesubmit
             
 
###################################################################################
             
        output$logo<-renderImage({list(src="/data/manke/sikora/shiny_apps/userIN_to_yaml/MPIIE_logo_sRGB.jpg",width=100,height=100)},deleteFile =FALSE)
             
        output$walkthrough<-renderText({"<font size=4><ol><li>Select input data. You can do so by providing either a combination of names and project number or by pasting the path to the folder containg your input reads. Click on Add dataset to retrieve the data. Repeat the procedure until you have retrieved all the data you would like to jointly analyze.</li><li>Specify workflow parameters: which kind of analysis should be performed on your data? Which reference genome should be used? Provide an optional title to your analysis.</li><li>Fill in the sample sheet.Provide experimental group information for your samples. For a ChIPseq workflow, provide the information on input-chip matching and mark width.</li><li>Save your sample sheet. This will reset the table to default.</li><li>Fill in your email address, any comments you would like to pass to the bioinformatic facility and check any boxes might be relevant to your data.</li><li>Submit the analysis. Verify the copy of your request in your email box.</li><li>You're done! You will be contacted by the bioinfo facility as soon as your data goes through the requested pipeline.</li></ol></font>"})
        output$workflow<-renderImage({list(src="/data/manke/sikora/shiny_apps/userIN_to_yaml/dev/userIN_to_yaml.workflow.png",width=800,height=600)},deleteFile=FALSE)     

        output$resultPanels<-renderUI({myTabs<-list(tabPanel(title="Walkthrough",
                                                             fluidPage(
                                                               box(htmlOutput("walkthrough"),width=12),
                                                               imageOutput("workflow",inline=TRUE)
                                                               
                                                               
                                                             )
                                                        ),
          
          
                                                    tabPanel(title="Input Data",
                                                      fluidPage(
                                                          box(textOutput("ispaired"),width=12,height=50,title="Read pairing detection"),
                                                          rHandsontableOutput("hot"),
                                                          box(textOutput("datawarnings"),width=4,height=200,title="Data Warnings"),
                                                              actionButton(inputId="savetable",label="Save sample sheet"),
                                                          box(textOutput("sIsaved"),width=4,height=100,title="Table status")
                                                            
                                                              )
                                                          ),
                                                    
                                                    tabPanel(title="User information",
                                                             fluidPage(
                                                               fluidRow(
                                                                 uiOutput("from"),
                                                                 uiOutput("freetext"),
                                                                 uiOutput("datarequests")
                                                                       ),
                                                                 actionButton(inputId="savesubmit",label="Save workflow configuration and submit request"),
                                                               box(textOutput("eSent"),width=4,height=100,title="Request status")

                                                             )
                                                    )
                                                    
                                                        )

            do.call(tabsetPanel, myTabs)})


}

shinyApp(ui, server)
