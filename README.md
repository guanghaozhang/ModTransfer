# ModTransfer
This package is for the unpublished working paper "Transportable Risk Prediction with Heterogeneous Electronic Health Records".

The objective of this package is to address heterogeneity between EHR datasets by estimating mapping matrix. Specifically, Medical codes across electronic health record (EHR) systems and across coding systems are usually heterogeneous. When heterogeneity across different EHR systems is determined by an orthogonal mapping matrix, we proposed estimation methods of the orthogonal mapping matrix such that features from these two heterogeneous EHR systems can be automatically mapped. The estimation methods use sample covariance matrices to obtain the mapping matrix and are therefore scalable in terms of the number of health records. With estimated mapping matrix, one can easily transfer a prediction model trained at source system to the heterogeneous target system.

