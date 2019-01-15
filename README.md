# userIN_to_yaml
A shiny app that collects user input required for running snakepipes.

User input is used to create a)sample sheet (tsv) containing group information on the samples b)a shell script with the snakemake command c)a samples.yaml file with chip_dict if ChIPseq analysis is requested.
The request is sent to bioinfo core and a copy to the user.


Tab "Input data":
1.User can search for existing *fastq.gz files under /data/group/sequencing_data by providing Group, Project Owner and Sequencing project ID. A simple unix find is used.
  Alternatively, user can paste a path to folder containing *fastq.gz files, e.g. in case of external data use.
  Use action button "Add dataset" to activate file search.
  Files identified this way are listed in the 'Read1' window.
  The action of searching for files can be repeated   - newly added folders are concatenated and files are searched for in all of the provided folders.
  While using one of the search modes, leave the fields belonging to the other search mode empty, or remove any previously entered values.
2.Some warnings are issued if sample names obtained by subsetting the read file names are non-unique: re-sequenced samples might be treated as replicates if overlooked.
3.User is encouraged to provide an analysis title, which is going to be integrated into resulting file and path names.
4.User is expected to edit the sample table, derived from identified read file names. Group information should be filled in for use in differential analyses. Samples with NAs left in group description will be dropped from the analysis. If ChIPseq analysis is desired, user is expected to fill in 3 additional drop-down fields: ChIPgroup (values: ChIP,Input), MatchedInput (values from the SampleID list) as well as MarkWidth (values: Broad,Narrow).
5.By pressing the 'Save sample sheet' button, sample sheet (and samples.yaml in case of ChIP data) is saved to file.
  Additional warnings might be displayed if, based on Group information in the sample table, no replicates are detected.
6.User is expected to select the NGS workflow they would like to have applied on their data from the drop-down menu on the sidebar.

Tab "Workflow parameters":
6.User is expected to select organism and workflow version from the respective drop-down menus.
 An updated shell command is displayed in the text box.
7.User is expected to enter their email address and encouraged to provide any additional information about their project in the respective text fields.
8.By clicking the 'Save workflow configuration and submit request', en email with the copy of the request and the files produced by the app in attachment is sent to the bioinfo core facility, with user in CC.

### this is a test modification
### this is another test modification
