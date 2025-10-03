= Introduction


Logistic regression is the most common method for clinical prediction models with a binary outcome [REF]. More recently, a variety of more flexible machine learning or artificial intelligence such as random forests, gradient boosting, and deep learning methods have been tried [REF].


linear nonlinear
interpretable harder to interpret
few parameters many parameters
lower risk of overfitting higher risk of overfitting


Why Bayesian:
regularisation hyperparameters can be learned from the data by placing a prior over them. This removes the need for cross-validation to determine suitable values [REF my comp tox paper]. Furthermore, values for cross-validation are often selected randomly or on a grid, which means that the optimal values may not be found, although they can be closely approximated with enough random values or a dense grid of values, but this can be computationally expensive. Finally, uncertainty in the hyperparameter values is properly accounted for an propagated through the model.

Loss landscape is multi-modal and ML methods that rely on optimization may converge to local optima for parameter values. Bayesian methods can explore the landscape better [REF: Andrew Gordon Wilson papers]

Obtain uncertainty in predictions for individual patients.

Models can be judged by additional criteria (1) CI coverage, (2) CRPS, which gives a more thorough  assessment of performance

== Research question

@Lazic2016
