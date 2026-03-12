# TODO

* implement ptweedie distribution (either in R or Rcpp): use the gamma poisson-structure (to truncate remenrer to use that cdf of gamma < 1, so just cut the num of Poisson distr)

* according to the examples in betanbb, the model is flat on 0: what is wrong?

* revise tests

* add hot_start to staticdistr and empdistr (the argument removes initial zeros)

* is it worh defining and HSP distributional object?

* i don't like how staticdistr has been organised