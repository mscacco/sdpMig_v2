## sdpMig_v2
Optimal migration model on a continuous landscape
This model expands on sdpMig at discrete sites. This new version implements migratory decisions on a continuous gridded landscape, where the decision or staying or migrating is optimised relative to the energy balance (gain vs cost) and the terminal rewards (fitness, i.e. arrival at the breeding ground within a certain time window and with given body conditions).
In this model implementation, gain (energy intake) at every site is determined by the opening day (day of spring start...


## Installation

Development version (package name `sdpMig2`, repository `sdpMig_v2`):

```r
# install.packages("remotes")
remotes::install_github("mscacco/sdpMig_v2")
library(sdpMig2)
```

An example workflow is in `inst/examples/example_modelWorkflow_dynamicGain_v5.R`.
