

### Module I 

This module generates a comprehensive collection of miRNA candidates by aggregating already annotated miRNAs from four plant miRNA databases (i.e., miRBase, PmiREN, sRNAanno, and PsRNA) and predicted miRNAs from user-submitted sRNA-Seq data.

|       Tool       |                          **Input**                           |                            Output                            |                         Applications                         |
| :--------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| *miRNARetrival*  |             Name of species and miRNA databases              | Already annotated miRNAs<br><a href="https://deepngs.nwafu.edu.cn/static/welcome/testData/Test_results/I_Overview.html">Overview</a><br><a href="https://deepngs.nwafu.edu.cn/static/welcome/testData/Test_results/I_Aggregation.html">Aggregation</a> | Aggregate annotated miRNAs provided by different miRNA databases |
|  *miRNAPredict*  | SRA accession numbers or uploaded fastq files<br/><a href="../Test_data/I_miRNAPredict_input.txt">Run accessions</a> | Predicted miRNAs<br/><a href="https://deepngs.nwafu.edu.cn/static/welcome/testData/Test_results/I_Identification.txt">Prediction</a> |              Predict miRNAs from sRNA-Seq data               |
| *genomeRetrival* | Name of species or genome sequences and annotation<br><a href="../Test_data/I_Arabidopsis_trsnsnoRNAs.fa.gz">trsnsnoRNAs.fa.gz</a> | Path of formatted genome sequences and annotation<br/><a href="https://deepngs.nwafu.edu.cn/static/welcome/testData/Test_results/I_Path_of_species.simple">Output</a> |             Get genome sequences and annotation              |
| *miRNATranslate* |          Output from miRNARetrival and miRNAPredict          | miRNA and miRNA precursors with a uniform format<br/><a href="https://deepngs.nwafu.edu.cn/static/welcome/testData/Test_results/I_miRNATranslate_output.txt">Output</a> | Translate annotated and predicted miRNAs into the genomic coordinate system |

### Module II

This module selects a subset of miRNA candidates according to the high-throughput criteria and/or using an ML-based approach.

|    **Tools**     |         **Input**          |                          **Output**                          |           Applications            |
| :--------------: | :------------------------: | :----------------------------------------------------------: | :-------------------------------: |
| *miRNASelection* | Output from miRNATranslate | Selected miRNAs<br><a href="https://deepngs.nwafu.edu.cn/static/welcome/testData/Test_results/II_miRNASelection_output.txt">Output</a> | Select promising miRNA candidates |

### Module III

This module provides the information for all miRNA candidates  for rapid curation of the quality of selected miRNAs.

|    **Tools**     |          **Input**          |                          **Output**                          |          Applications           |
| :--------------: | :-------------------------: | :----------------------------------------------------------: | :-----------------------------: |
| *manualCuration* | Output from MiRNA Selection | Summary and report pages<br/><a href="https://deepngs.nwafu.edu.cn/static/welcome/testData/Test_results/III_manualCuration_output.html">Output</a> | Determine the quality of miRNAs |



> The table above briefly describes each function in this module. The sample data are displayed through a hyperlink. Users can preview the results by clicking it. Due to the limitation of file size, some results may not be displayed properly. For more detailed descriptions about inputs and outputs, please refer to [tutorial for iwa-miRNA](https://github.com/cma2015/iwa-miRNA/tree/master/Tutorials/User_manual.md) or the interpretation of each functional analysis interface in [web server](https://deepngs.nwafu.edu.cn/).

