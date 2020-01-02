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
        

####observe workflow choice and create reactive table

        observe({
        values$inWorkflow<-input$selectworkflow

        if(values$inWorkflow=="ChIP-seq"){
                values$DF<-data.frame(SampleID=values$datshort,Group=(factor(rep("NA",(length(values$datshort))),levels=c("Control","Treatment","NA"))),Merge=factor(rep("NA",(length(values$datshort))),levels=c("NA",unique(values$datshort)),ordered=TRUE),Read1=values$Read1,ChIPgroup=factor(rep("NA",(length(values$datshort))),levels=c("NA","ChIP","Input"),ordered=TRUE),MatchedInput=factor(rep("NA",(length(values$datshort))),levels=c("NA",unique(values$datshort)),ordered=TRUE),MarkWidth=factor(rep("NA",(length(values$datshort))),levels=c("NA","Broad","Narrow"),ordered=TRUE),stringsAsFactors = F)
                output$tabdesc<-renderText({"<font size=4><ul><li>SampleID: automaticaly parsed from read names: do not modify.</li><li>Group: assign samples to Control and Treatment groups. Leave NA for samples you would like to exclude from the analysis.</li><li>Merge: select a Sample ID under which you would like to merge fastq files. Leave NA if not needed.</li><li>Read1: for your information, the identity of the read file. Only 1 of the 2 paired end files will be listed.</li><li>ChIPgroup: define the samples as Input or ChIP.</li><li>MatchedInput: to each ChIP sample, assign the matched Input sample.</li><li>MarkWidth: specify if Narrow or Broad peak calling mode should be used.</li></ul></font>"})
                } 

            else if(values$inWorkflow=="HiC"){
                values$DF<-data.frame(SampleID=values$datshort,Group=(rep("NA",(length(values$datshort)))),Merge=factor(rep("NA",(length(values$datshort))),levels=c("NA",unique(values$datshort)),ordered=TRUE),Read1=values$Read1,stringsAsFactors = F)
                output$tabdesc<-renderText({"<font size=4><ul><li>SampleID: automaticaly parsed from read names: do not modify.</li><li>Group: type in sample group names. Leave NA for samples you would like to exclude from the analysis.</li><li>Merge: select a Sample ID under which you would like to merge fastq files. Leave NA if not needed.</li><li>Read1: for your information, the identity of the read file. Only 1 of the 2 paired end files will be listed.</li></ul></font>"})
                }
           
            else if(values$inWorkflow=="WGBS"){
                values$DF<-data.frame(SampleID=values$datshort,Group=(factor(rep("NA",(length(values$datshort))),levels=c("Control","Treatment","WT","Mut","NA"))),PlottingID=values$datshort,Merge=factor(rep("NA",(length(values$datshort))),levels=c("NA",unique(values$datshort)),ordered=TRUE),Read1=values$Read1,stringsAsFactors = F)
                output$tabdesc<-renderText({"<font size=4><ul><li>SampleID: automaticaly parsed from read names: do not modify.</li><li>Group: assign samples to Control and Treatment or WT and Mut groups. Leave NA for samples you would like to exclude from the analysis.</li><li>PlottingID: type in names to use for plots if deviating from sample names.</li><li>Merge: select a Sample ID under which you would like to merge fastq files. Leave NA if not needed.</li><li>Read1: for your information, the identity of the read file. Only 1 of the 2 paired end files will be listed.</li></ul></font>"})
                }  else {

               values$DF<-data.frame(SampleID=values$datshort,Group=(factor(rep("NA",(length(values$datshort))),levels=c("Control","Treatment","NA"))),Merge=factor(rep("NA",(length(values$datshort))),levels=c("NA",unique(values$datshort)),ordered=TRUE),Read1=values$Read1,stringsAsFactors = F)
               output$tabdesc<-renderText({"<font size=4><ul><li>SampleID: automaticaly parsed from read names: do not modify.</li><li>Group: assign samples to Control and Treatment groups. Leave NA for samples you would like to exclude from the analysis.</li><li>Merge: select a Sample ID under which you would like to merge fastq files. Leave NA if not needed.</li><li>Read1: for your information, the identity of the read file. Only 1 of the 2 paired end files will be listed.</li></ul></font>"})
                   }
              

            #values$DF<-DF
           DF2<-isolate(values$DF)
           values$DF2<-DF2
     })#end of observe input$selectworkflow


            

       observe({
          if(!is.null(input$hot)){
             DF2 <- hot_to_r(input$hot)
             values$DF2<-DF2}
        })
        
        
        output$hot <- renderRHandsontable({
          DF2<-values$DF2
          rhandsontable(DF2)})

##########observe save samplesheet
        observeEvent(input$savetable, {
       
          values$analysisName<-gsub("[^[:alnum:]]", "_", isolate(input$analysistitle))
          values$ranstring<-stri_rand_strings(n=1,length=8)
          
          sampleInfo<-isolate(values$DF2)
          output$hot<-renderRHandsontable({rhandsontable(sampleInfo)})
          ###check for replicates, else issue a warning
          if(sum(is.na(sampleInfo$Group),sampleInfo$Group %in% "NA")<length(sampleInfo$Group)){
              sampleInfo$Group<-as.character(sampleInfo$Group)
              sampleInfo<-sampleInfo[!is.na(sampleInfo$Group)&!sampleInfo$Group %in% "NA",]
              #rownames(sampleInfo)<-make.names(sampleInfo$SampleID, unique=TRUE)
              srep<-table(sampleInfo$SampleID)[unique(sampleInfo$SampleID)]
              if(any(srep>1)){
                  rv<-unlist(sapply(srep,function(X)seq(1,X)))
                  rn<-paste0(sampleInfo$SampleID,"_",rv)
                  rownames(sampleInfo)<-rn} else {rownames(sampleInfo)<-sampleInfo$SampleID}
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
          }else {srep<-table(sampleInfo$SampleID)[unique(sampleInfo$SampleID)]
              if(any(srep>1)){
                rv<-unlist(sapply(srep,function(X)seq(1,X)))
                rn<-paste0(sampleInfo$SampleID,"_",rv)
                rownames(sampleInfo)<-rn} else {rownames(sampleInfo)<-sampleInfo$SampleID}}    
           
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
          
          colnames(sampleInfo)[1:2]<-c("name","condition")
          write.table(sampleInfo,file=values$sInfoDest,sep="\t",quote=FALSE)
          output$sIsaved<-renderText("Sample sheet saved.")
                 
                 

        output$from<-renderUI({textInput(inputId="sender",label="Your email address",placeholder="lastname@ie-freiburg.mpg.de")})
        output$freetext<-renderUI({textInput(inputId="comments",label="Your message to the bioinfo facility",placeholder="Sample X might be an outlier.",width="600px")})

        output$datarequests<- renderUI({
        if(values$inWorkflow %in% c("ATAC-seq","ChIP-seq","DNA-mapping","RNA-seq","WGBS")){
        tagList(
        checkboxInput(inputId="fbam", label="I want to start the analysis from the bam files.", value = FALSE, width = NULL),
        checkboxInput(inputId="merge", label="I want to request sample merging.", value = FALSE, width = NULL),
        checkboxInput(inputId="beff", label="I expect batch effect in my data.", value = FALSE, width = NULL),
        checkboxInput(inputId="nodiff", label="I don't need differential analysis.", value = FALSE, width = NULL),
        checkboxInput(inputId="SE", label="I have single end, NOT paired end data.", value = FALSE, width = NULL))
        }else if (values$inWorkflow %in% "HiC"){
          tagList(checkboxInput(inputId="merge", label="I want to request sample merging.", value = FALSE, width = NULL),
                  checkboxInput(inputId="nodiff", label="I don't need differential analysis.", value = FALSE, width = NULL),
                  checkboxInput(inputId="enz", label="I used DpnII instead of HindIII for restriction digest.", value = FALSE, width = NULL),
                  checkboxInput(inputId="notads", label="I don't need TAD calling.", value = FALSE, width = NULL))
        }
      })
      values$mstr<-reactive({ifelse(input$merge,"--mergeSamples","")})
      values$estr<-reactive({ifelse(input$enz,"--enzyme DpnII","--enzyme HindIII")})


      })#end of observe input$savetable



###########observe submit request
        observeEvent(input$savesubmit, {
               
##########prepare the command
        path_to_exec<-paste0("module load snakePipes; ",values$inWorkflow)###add version selection
           topdir<-sprintf("/data/processing/bioinfo-core/requests/%s_%s_%s",values$analysisName,values$ranstring,values$inWorkflow)
           dsel<-c("fastq.gz"="fastq","bam"="bam")
           indir<-sprintf("%s/%s",topdir,dsel[input$selectformat])
           link_cmd<- paste0("ln -t ",indir,' -s ',paste(values$Reads,collapse=" "))
           
           outdir<-sprintf("%s/analysis",topdir)
           
           cp_sInfo_cmd<-sprintf("cp -v %s %s",values$sInfoDest,topdir)
           values$sInfo_in<-paste0(topdir,"/",basename(values$sInfoDest))
           genome_sel<-c("PLEASE SELECT A GENOME"="NONE","Zebrafish [zv10]"="GRCz10","Fission yeast"="SchizoSPombe_ASM294v2","Fruitfly [dm6]"="dm6","Fruitfly [dm3]"="dm3","Human [hg37]"="hs37d5","Human [hg38]"="hg38","Mouse [mm9]"="mm9","Mouse [mm10]"="mm10")  
           values$genome<-genome_sel[input$genome]
           
          
           path_to_DNA_mapping<-paste0("module load snakePipes; DNA-mapping")
           
           fbam<-c("fastq.gz"="","bam"="--fromBam")
           
           if(values$inWorkflow=="ATAC-seq"){
               
              values$command<-sprintf("mkdir -p %s ; %s ; %s ;%s -i %s -o %s %s ; %s -d %s --sampleSheet %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_DNA_mapping,indir,outdir,values$genome,path_to_exec,outdir,values$sInfo_in,values$genome)
              output$command<-renderText({ values$command })
              
                     }##end of ATACseq
           
           else if(values$inWorkflow=="ChIP-seq"){
             
             cp_chDict_cmd<-sprintf("cp -v %s %s",values$chDictDest,topdir)
             values$chDictr_in<-paste0(topdir,"/",basename(values$chDictDest))
             
             values$command<-sprintf("mkdir -p %s ; %s  ; %s ; %s ; %s -i %s -o %s %s ; %s -d %s --sampleSheet %s %s %s",indir,link_cmd,cp_sInfo_cmd,cp_chDict_cmd,path_to_DNA_mapping,indir,outdir,values$genome,path_to_exec,outdir,values$sInfo_in,values$genome,values$chDictr_in) 
             output$command<-renderText({ values$command })
             
             
           }##end of ChIP-seq
           
           
           else if(values$inWorkflow=="HiC"){
             mstr<-isolate(values$mstr())
             estr<-isolate(values$estr())
             
             values$command<-sprintf("mkdir -p %s ; %s ;  %s ; %s -i %s -o %s %s %s %s",indir,link_cmd,cp_sInfo_cmd,path_to_exec,indir,outdir,estr,mstr,values$genome) 
             output$command<-renderText({ values$command })
             
           }##end of HiC
           
           
           else if(values$inWorkflow=="DNA-mapping"){
             
             values$command<-sprintf("mkdir -p %s ; %s ;  %s ; %s -i %s -o %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_DNA_mapping,indir,outdir,values$genome) 
             output$command<-renderText({ values$command })
             
           } #end of DNA-mapping
           
           else if(values$inWorkflow=="RNA-seq"){
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s ; %s --mode alignment -i %s -o %s --sampleSheet %s %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_exec,indir,outdir,values$sInfo_in,fbam[input$selectformat],values$genome) 
             output$command<-renderText({ values$command })
             
           } #end of RNA-seq
           
           else if(values$inWorkflow=="WGBS"){
             
             values$command<-sprintf("mkdir -p %s ; %s ; %s ; %s -i %s -o %s --sampleSheet %s %s %s ",indir,link_cmd,cp_sInfo_cmd,path_to_exec,indir,outdir,values$sInfo_in,fbam[input$selectformat],values$genome) 
             output$command<-renderText({ values$command })
             
           } #end of WGBS


####write command to script
        bshscript<-sprintf("/data/manke/group/shiny/snakepipes_input/%s_%s_script.sh",values$ranstring,values$analysisName)
        fileConn<-file(bshscript)
        shebang<-"#!/bin/bash"
        #exports<-"export PATH=/data/processing/conda/bin:$PATH"
        exports<-""
        writeLines(c(shebang,exports,unlist(strsplit(isolate(values$command),split=";"))), fileConn)
        close(fileConn)


###############compose email message 
               fbam_request<-ifelse(isolate(input$fbam),"I want to start the analysis from the bam files.","I want to start the analysis from fastq files.")
               merge_request<-ifelse(isolate(input$merge),"I want to request sample merging. Please consider the information I entered in the Merge column of the sample sheet. \n Please update the sample sheet after merging files.","No sample merging is needed.")
               b_eff_request<-ifelse(isolate(input$beff),"I expect batch effect in my data.","No batch effect is expected.")
               nodiff_request<-ifelse(isolate(input$nodiff),"I don't need differential analysis.","I want to request differential analysis.")
               SE_request<-ifelse(isolate(input$SE),"I have single end data.","I have paired end data.")
               enz_request<-ifelse(isolate(input$enz),"I used DpnII.","I used HindIII.")
               notads_request<-ifelse(isolate(input$notads),"I don't need TAD calling.","I need TAD calling.")
               cc<-isolate(input$sender)
               #from<-sprintf("<sendmailR@%s>", Sys.info()[4])
               from<-"sendmailR@ie-freiburg.mpg.de"
               to<-"bioinfo-core@ie-freiburg.mpg.de" 
               #to<-"sikora@ie-freiburg.mpg.de"
               subject<-paste0("Analysis request ",isolate(input$analysistitle), "_" ,isolate(values$ranstring))
               
               if(values$inWorkflow %in% c("ATAC-seq","ChIP-seq","DNA-mapping","RNA-seq","WGBS")){msg <- gsub(";","\n \n",paste0(cc," has requested the following analysis: \n \n Workflow: ", isolate(values$inWorkflow)," \n \n Genome: ", values$genome," \n \n ", fbam_request, " \n \n " ,merge_request," \n \n ", b_eff_request ," \n \n ", nodiff_request, " \n \n ", SE_request, " \n \n User comments: \n \n", isolate(input$comments),"\n \n Please review the input files and the attached sample sheet before proceeding. \n \n End of message. \n \n ",paste(rep("#",times=80),collapse="")," \n \n ", values$command ,"\n \n"))
               }else if(values$inWorkflow %in% "HiC"){msg <- gsub(";","\n \n",paste0(cc," has requested the following analysis: \n \n Workflow: ", isolate(values$inWorkflow)," \n \n Genome: ", values$genome," \n \n ", merge_request," \n \n ", nodiff_request, " \n \n ",enz_request, " \n \n ",notads_request ," \n \n User comments: \n \n", isolate(input$comments),"\n \n Please review the input files and the attached sample sheet before proceeding. \n \n End of message. \n \n ",paste(rep("#",times=80),collapse="")," \n \n ", values$command ,"\n \n"))}
               
               
               if(values$inWorkflow=="ChIP-seq"){
                 sendmail(from=sprintf("%s",from), to=to, subject=subject, control=list(smtpServer="owa.ie-freiburg.mpg.de"), msg=list(msg,mime_part(isolate(values$sInfoDest)),mime_part(isolate(values$chDictDest)),mime_part(bshscript)),cc=sprintf("%s",cc))
               }
               else{
               sendmail(from=sprintf("%s",from), to=to, subject=subject, control=list(smtpServer="owa.ie-freiburg.mpg.de"), msg=list(msg,mime_part(isolate(values$sInfoDest)),mime_part(bshscript)),cc=sprintf("%s",cc))}
               output$eSent<-renderText("Your request has been sent to the MPI-IE Bioinfo facility. A copy was sent to your email address.")

               
             },ignoreInit=TRUE,once=TRUE)#end of observe input$savesubmit



###################################################################################
             
        output$logo<-renderImage({list(src="/data/manke/sikora/shiny_apps/userIN_to_yaml/MPIIE_logo_sRGB.jpg",width=100,height=100)},deleteFile =FALSE)
             
        output$walkthrough<-renderText({"<font size=4><ol><li>Specify workflow parameters: which kind of analysis should be performed on your data? Which reference genome should be used? Provide an optional title to your analysis.</li><li>Select input data. You can do so by providing either a combination of names and project number or by pasting the path to the folder containg your input reads. Click on Add dataset to retrieve the data. The result appears in the Input Data tab. Repeat the procedure until you have retrieved all the data you would like to jointly analyze.</li><li>Fill in the workflow-specific sample sheet. Detailed explanations will be displayed accordingly. Provide new sample names in the Merge column if you wish your reads to be merged.</li><li>Save your sample sheet. This will reset the table to default.</li><li>Fill in your email address, any comments you would like to pass to the bioinformatic facility and check any boxes might be relevant to your data.</li><li>Submit the analysis. Verify the copy of your request in your email box.</li><li>You're done! You will be contacted by the bioinfo facility as soon as your data goes through the requested pipeline.</li></ol></font>"})
        output$workflow<-renderImage({list(src="/data/manke/sikora/shiny_apps/userIN_to_yaml/dev/userIN_to_yaml.workflow.png",width=800,height=600)},deleteFile=FALSE)   
        output$sampleInfoWarning<-renderText({"Please save sample sheet first."})

        output$resultPanels<-renderUI({myTabs<-list(
                                                           
                                                    tabPanel(title="Input Data",
                                                      fluidPage(tags$head(tags$style(HTML(".shiny-output-error-validation {color: green; position: relative; bottom: -300px;}"))),
                                                          box(textOutput("ispaired"),width=12,height=50,title="Read pairing detection"),
                                                          rHandsontableOutput("hot"),
                                                          box(textOutput("datawarnings"),width=4,height=300,title="Data Warnings"),
                                                              actionButton(inputId="savetable",label="Save sample sheet"),
                                                          box(textOutput("sIsaved"),width=4,height=100,title="Table status"),
                                                          box(htmlOutput("tabdesc"),width=12,title="Column description")
                                                            
                                                              )
                                                          ),
                                                    
                                                    tabPanel(title="User information",conditionalPanel(condition= "input.savetable == 1",
                                                             fluidPage(
                                                               fluidRow(
                                                                 uiOutput("from"),
                                                                 uiOutput("freetext"),
                                                                 uiOutput("datarequests")
                                                                       ),
                                                                 actionButton(inputId="savesubmit",label="Save workflow configuration and submit request"),
                                                               box(textOutput("eSent"),width=4,height=100,title="Request status")

                                                             )
                                                    ),
                                                    conditionalPanel(condition= "input.savetable == 0",box(textOutput("sampleInfoWarning")))),
                                                    
                                                    tabPanel(title="Walkthrough",
                                                             fluidPage(
                                                               box(htmlOutput("walkthrough"),width=12),
                                                               imageOutput("workflow",inline=TRUE)
                                                             )
                                                    )
                                                    
                                                        )

            do.call(tabsetPanel, myTabs)
        })

}

shinyApp(ui, server)
