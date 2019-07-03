# ShinyApp

This ShinyApp creates an interactive application that plots exon reads for a group of samples. The input fields are:

- **Gene name** (based on HGNC nomenclature) or **exon information** (in the format of "CNV_chr_start_end"). If the first checkbox is checked, the exon information will be used; if it is unchecked, the gene name will be used.
- **Sample to highlight**: users can input the id of a specific sample, which will be highlighted in the plot.
- **SVD-10**: this checkbox allows users to use datasets that have been processed with the CoNIFER algorithm, which normalizes the data and reduces noise using singular value decomposition (SVD). 10 components are deleted in this application.
- **Display both parents**: highlights the parents of the selected sample.

## Data Processing

The preliminary data processing is executed in process-data.R. It utilizes the CoNIFER algorithm, which normalizes the data and reduces noise using singular value decomposition (SVD). In this application, either none or 10 components are deleted in SVD. The processed dataset is then divided into separate dataframes based on the chromosome.

BioMart is also used to retrieve data on gene names and their respective regions.

## Data Format

The data used in this application are in hdf5 format. Data from different samples are stored in separate hdf5 files but in the same folder named "hdf5". This folder is located in the working directory of the R project.

## Usage

Run process-data.R to process the data. Then run app.R either by clicking on the RunApp button or by calling ```runApp()```. Note that since Shiny is developed by RStudio, the IDE used to run this application must be RStudio.
