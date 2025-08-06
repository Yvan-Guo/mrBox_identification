# mrBox_identification

Screening of DNA methylation regulatory boxes (mrBox) with convolutional neural networks following data distillation.

# Overview

This code uses a CNN-enhanced framework　that effectively integrated DNA methylation data and genome sequence data to identify DNA methylation regulatory boxes (mrBox). Given the large scale of genomic data and substantial redundancy present, data preparation involved rigorous cleaning and the application of data distillation techniques based on correlation analyses of DNA methylation, which was aimed at improving model efficiency and increasing the biological relevance of the findings.　Then, mrBox are idenitified by a CNN-enhanced MeShClust clustering algorithm.

# Hardware requirement
a standard computer with GPU (tested on RTX4070) and enough RAM to support the in-memory operations.

# OS system
Windows (Windows Server 21H2) or Linux (Ubuntu 20.18)

# Software requirement:
Python v3.9
CUDA v11.7.10
MeshClust v3.0
tensorflow　

# Expected outcome
Recreated mrBox identification in https://doi.org/10.21203/rs.3.rs-7223043/v1
