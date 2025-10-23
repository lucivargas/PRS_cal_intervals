# PRS_calibration_intervals

This repository provides code for calibrating polygenic risk score (PRS) values according to individualized genetic similarity estimates, specifically those derived using the program ADMIXTURE. 

Briefly, in an African American population, the method calibrates each individual’s PRS value by regressing out the genetic similarity to the Yoruba in Ibadan, Nigeria [YRI] reference panel (GSYRI) and dividing the new score by the standard deviation of the PRS values for all participants within the same GSYRI-defined bin (to account for possibility that variance in the PRS is also a function of ancestry proportions). 

The PRS calibration by genetic similarity can be represented as:

$$
PRS_{i,d}=  \frac{(x_i-\mu_i)}{\sigma_d}
$$

Where  
  $i$ represents an individual $i$; <br>
  $d$ is the GSYRI decile grouping for individual $i$; <br>
  $x_i$ is the observed unscaled PRS for person $i$; <br>
  $\mu_i$ is the predicted PRS for person 𝑖 conditional on their GSYRI; <br>
  and $\sigma_d$ is the estimated standard deviation of PRS scores for individuals in GSYRI decile $d$.

This method was used in the following publication: <br>
[Ancestry Calibration of Polygenic Risk Scores Improves Risk Stratification and Effect Estimation in African American Adults. Vargas _et al_., (2025). _medRxiv_ 2025.06.18.25329573; doi: https://doi.org/10.1101/2025.06.18.25329573](https://doi.org/10.1101/2025.06.18.25329573)

