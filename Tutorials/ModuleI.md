### Module I 

This module generates a comprehensive collection of miRNA candidates by aggregating already annotated miRNAs from four plant miRNA databases (i.e., miRBase, PmiREN, sRNAanno, and PsRNA) and predicted miRNAs from user-submitted sRNA-Seq data.

|  **Tool name**   |                          **Input**                           |                            Output                            |                         Applications                         |
| :--------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| *miRNARetrival*  |             Name of species and miRNA databases              | Already annotated miRNAs<br><a href="../Test_data/miRNARetrieval_output1.html">Overview</a><br><a href="../Test_data/miRNARetrieval_output2.html">Aggregation</a> | Aggregate annotated miRNAs provided by different miRNA databases |
|  *miRNAPredict*  | SRA accessions or uploaded fastq files<br/><a href="../Test_data/miRNAPredict_inputData.tar.gz">collaspeData</a> | Predicted miRNAs<br/><a href="../Test_data/miRNAPredict_output.txt">Prediction</a> |              Predict miRNAs from sRNA-Seq data               |
| *genomeRetrival* |      Name of species or genome sequences and annotation      |      Path of formatted genome sequences and annotation       |             Get genome sequences and annotation              |
| *miRNATranslate* |          Output from miRNARetrival and miRNAPredict          | miRNA and miRNA precursors with a uniformed format<br/><a href="../Test_data/miRTranslate_output.txt">Output</a> | Translate annotated and predicted miRNAs into the genomic coordinate system |

The table above briefly describes each function in this module. The sample data are displayed through a hyperlink. Users can preview the results by clicking it. Due to the limitation of file size, some results cannot be displayed. For more detailed descriptions about inputs and outputs, please refer to the interpretation of each functional analysis interface in web server.