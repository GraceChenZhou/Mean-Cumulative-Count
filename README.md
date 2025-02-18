# Mean Cumulative Count (MCC)
MCC is used to reflect mean number of events per person by a given time. 

# Relationship between MCC and cumulative incidence(CumI)

$$MCC(t)=\sum_{p=1}^{m} CumI_p(t)$$, 
where $CumI_p(t)$ represents the cumulative incidence for the pth (p=1,2,\dot, m) occurrence of the event of interest by time t. 

# Practice in R

* Estimation: **MY.MCC** function with the capability to handle left-truncation
* Visulaiztion: at.risk function and mcc.plot function

For details, please refer to the vignette_mcc.html

# Reference

Dong H, Robison LL, Leisenring WM, Martin LJ, Armstrong GT, Yasui Y. Estimating the burden of recurrent events in the presence of competing risks: the method of mean cumulative count. Am J Epidemiol. 2015 Apr 1;181(7):532-40. d. Epub 2015 Feb 17. PMID: 25693770; PMCID: PMC4371763. doi: 10.1093/aje/kwu28

Geskus RB. Cause-specific cumulative incidence estimation and the fine and gray model under both left truncation and right censoring. Biometrics. 2011 Mar;67(1):39-49. doi: 10.1111/j.1541-0420.2010.01420.x. PMID: 20377575.

