# Mean Cumulative Count (MCC)
MCC is used to reflect mean number of events per person by a given time. 

# Relationship between MCC and cumulative incidence(CumI)

$$MCC(t)=\sum_{p=1}^{m} CumI_p(t)$$, 
where $CumI_p(t)$ represents the cumulative incidence for the pth $(p=1,2,\dots,m)$ occurrence of the event of interest by time t. 

# Practice in R

* Estimation: **MY.MCC** function with the capability to handle left-truncation
* Visulaiztion: at.risk function and mcc.plot function

**For details, please refer to the vignette_mcc.html**

# Reference

Dong H, Robison LL, Leisenring WM, Martin LJ, Armstrong GT, Yasui Y. Estimating the burden of recurrent events in the presence of competing risks: the method of mean cumulative count. Am J Epidemiol. 2015 Apr 1;181(7):532-40. d. Epub 2015 Feb 17. PMID: 25693770; PMCID: PMC4371763. doi: 10.1093/aje/kwu28

Geskus RB. Cause-specific cumulative incidence estimation and the fine and gray model under both left truncation and right censoring. Biometrics. 2011 Mar;67(1):39-49. doi: 10.1111/j.1541-0420.2010.01420.x. PMID: 20377575.

Bhakta N, Liu Q, Ness KK, Baassiri M, Eissa H, Yeo F, Chemaitilly W, Ehrhardt MJ, Bass J, Bishop MW, Shelton K, Lu L, Huang S, Li Z, Caron E, Lanctot J, Howell C, Folse T, Joshi V, Green DM, Mulrooney DA, Armstrong GT, Krull KR, Brinkman TM, Khan RB, Srivastava DK, Hudson MM, Yasui Y, Robison LL. The cumulative burden of surviving childhood cancer: an initial report from the St Jude Lifetime Cohort Study (SJLIFE). Lancet. 2017 Dec 9;390(10112):2569-2582. doi: 10.1016/S0140-6736(17)31610-0. Epub 2017 Sep 8. PMID: 28890157; PMCID: PMC5798235.

Mulrooney DA, Hyun G, Ness KK, Bhakta N, Pui CH, Ehrhardt MJ, Krull KR, Crom DB, Chemaitilly W, Srivastava DK, Relling MV, Jeha S, Green DM, Yasui Y, Robison LL, Hudson MM. The changing burden of long-term health outcomes in survivors of childhood acute lymphoblastic leukaemia: a retrospective analysis of the St Jude Lifetime Cohort Study. Lancet Haematol. 2019 Jun;6(6):e306-e316. doi: 10.1016/S2352-3026(19)30050-X. Epub 2019 May 8. PMID: 31078468; PMCID: PMC6756152.



