[S]
Modern machine learning (ML) methods are becoming increasingly popular for clinical predictions, even though simple logistic regression (LR) often performs better.

[C]
This is partly because clinical data are often small and noisy, and flexible models such as neural networks tend to overfit. In addition, many ML models may be suboptimal because they typically do not incorporate background information, such as directionality (e.g. high levels of a protein are associated with a high risk of cancer), monotonicity (e.g. increasingly higher values of a protein are always associated with the same or a higher risk of cancer), and smoothness (e.g. there is no sharp jump in cancer risk for a small increase in protein level). Logistic regression models on the other hand may not perform well if they do not account for nonlinearities.


The purpose of this project is to develop diagnostic models of early stage cancer from blood-based protein biomarkers.  The dataset contains 50 biomarkers from 680 healthy controls and over 73 cancer cases for several cancers (2-3 cancers will be selected). Data will be divided into a training (80%) and test (20%) set and 3-5 biomarkers selected based on the literature and predictive performance. By ensuring biologically plausible outputs, the proposed models could enhance clinicians' trust in machine learning methods for cancer diagnosis.


[IC]
[Q]
[A]
The models will be developed within a Bayesian framework as this provides the required flexibility to extend the basic LR model with monotonic splines, as well as regularization via hierarchical prior distributions to prevent overfitting. These models will be compared with standard LR models and more flexible ML methods such as random forests. The models will be evaluated on the usual metrics of accuracy, bias, discrimination, and calibration, as well as biological plausibility.
