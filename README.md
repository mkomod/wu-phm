# wu.ss.phm

Fit spike and slab proportional hazards moodels. This is an implementaation of the methods proposed by Ryan Wu, Mihye Ahn and Hojin Yang.

## Installation

```
devtools::install_github("mkomod/wu-phm")
```

## Usage

```
# Y: failure times
# delta: censoring indicator
# X: design matrix

wu.fit(Y, delta, X)
```

### References 

Ryan Wu, Mihye Ahn & Hojin Yang (2021): Spike-and-slab type variableselection in the Cox proportional hazards model for high-dimensional features, Journal of AppliedStatistics, DOI: 10.1080/02664763.2021.1893285
