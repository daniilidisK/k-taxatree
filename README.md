# k-taxatree

k-taxatree is a classification workflow written in [R](https://www.r-project.org/), predicting the labels of the first four taxonomic levels (kingdom, phylum, class, order) of metagenomic data with a multi-label Random Forest as the underlying model. The latter accepts as input 6-mer count vectors and as such a method to determine the appropriate k-length was also implemented. 

## Getting started

### Prerequisites

The packages needed to be installed, in order to run the project are:
- from CRAN

```
install.packages(c("parallel", "data.table", "cluster", "Rfast", "plyr", "caret", "stats", "UBL", "splitstackshape", "mlr", "mldr", "dplyr", "hash", "stringr", "randomForestSRC"))
```


### Installing

The project can be downloaded using git:

```
git clone https://github.com/BiodataAnalysisGroup/k-taxatree
```

### Running the project

The project consists of 11 main scripts in the folder [R-scripts](https://github.com/BiodataAnalysisGroup/k-taxatree/tree/main/R-scripts) containing the code required to perform all the steps from selecting the appropriate k-length and constructing the kmer matrix to predicting the labels of the validation subset:
- ```01_k_selection_tool.R```
- ```02_kmer_matrix_creation.R```
- ```03_feature_selection.R```
- ```04_model_hp_optimization.R```
- ```05_final_model.R```
- ```06_unassigned_predictions.R```
- ```07_validation_predictions.R```

In order to run the project, 

The project provides the input datasets and the outputs generated in every step of the workflow. The folder [emp-data](https://github.com/BiodataAnalysisGroup/k-taxatree/tree/main/emp-data) includes the datasets retrived from the Earth Microbiome Project Repository. The [Output]


For more details, please refer to the [wiki](https://github.com/BiodataAnalysisGroup/k-taxatree/wiki).


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
