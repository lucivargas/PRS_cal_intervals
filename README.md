# PRS_cal_intervals

This repository provides code for calibrating polygenic risk score (PRS) values according to individualized genetic similarity estimates, specifically those derived using the program ADMIXTURE. 

Briefly, the program calibrates each individual’s PRS value by regressing out the genetic similarity to the Yoruba in Ibadan, Nigeria [YRI] reference panel (GSYRI) and dividing the new score by the standard deviation of the PRS values for all participants within the same GSYRI-defined bin (to account for possibility that variance in the PRS is also a function of ancestry proportions). 

The PRS calibration by genetic similarity can be represented as:

〖PRS〗_(i,d)=  (x_i-μ_i)/σ_d 

Where 𝑖 is  individual 𝑖; 𝑑 is the GSYRI decile grouping for individual 𝑖; 𝒙𝒊 is the observed unscaled PRS for person 𝑖; 𝝁𝒊 is the predicted PRS for person 𝑖 conditional on their GSYRI; and 𝝈𝒅 is the estimated standard deviation of PRS scores for individuals in GSYRI decile 𝒅. 
