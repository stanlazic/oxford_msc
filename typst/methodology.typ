= Methodology

== Data

=== Source

=== Cleaning and preprocessing

Removed hemolyzed samples.
Removed Female subjects [Prostate only].
?Removed samples with unknown country.
Removed biomarkers with >20% missing values.
Restricted age range to be similar between Healthy and Cancer subjects.
Centered and scaled data [based on healthy controls] .
Normalize data [based on healthy controls; bestNormalize method].
Grouped similar countries with few samples together.
Imputed missing values [MICE method].
Removed “batch” (age, country) effects [based on healthy controls].


==  Analysis

=== Simulations

=== Additive model

=== Model with interaction effect

== Evaluation metrics

== Software
