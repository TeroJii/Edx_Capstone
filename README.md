# Edx_Capstone

This repository contains the Capstone project for Harvard X's Data Science Professional Certificate:

<https://www.edx.org/professional-certificate/harvardx-data-science>

The project explores the possibility of predicting small molecule retention time in liquid chromatography by using machine learning methods.
The data used in this project is based on the METLIN small molecule dataset published by Xavier Domingo-Almenara et al.
<https://www.nature.com/articles/s41467-019-13680-7>

The data has been expanded by calculating several molecular descriptors, and using them as predictors for retention time.
The code for producing the final dataset for this project is included in CreateDataset.R.

## Files submitted for review

The following files contain the project report as a pdf and Rmd file. The Report.R file contains all the R code related to project analyses

* Report.pdf
* Report.Rmd
* Report.R

The data used for the project is provided in the final_dataset.csv file.

## Supplementary files

The following figure files have been included to help make the pdf-report file smaller in size:

* finalmodelplot.png
* pca_plots.png
* scatterplot.png

The code for producing the plots is included in the Report.Rmd and Report.R files.

Other supplementary files include:

* mybibfile.bib (the references)
* DataSplit.PNG (png-file describing the data resampling methodology)
* Report.hmtl (report in html-format)
* classificationmetrics1.csv (model evaluation metrics saved in csv-format to save time when knitting the report) 
* classificationmetrics1.csv (model evaluation metrics saved in csv-format to save time when knitting the report)
* regressionmetrics.csv (model evaluation metrics saved in csv-format to save time when knitting the report)
