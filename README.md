# pm_pcap_ed
PM2.5, persistent cold air pools, and respiratory emergency department visits in the western United States

Short-term associations between PM2.5 air pollution, persistent cold air pool (PCAP) conditions, and respiratory emergency department (ED) visits across nine western US cities, 2001–2020.

This repository contains the R code to reproduce the analysis presented in the article:

Riss CR, Krall JR , Darrow LA, Kelly KE, Holmes HA, Strickland MJ. Ambient air pollution, temperature inversions, and respiratory emergency department visits in nine western U.S. Cities. 

This project was funded by grant R01ES032810-10 from the National Institute of Environmental Health Sciences

Study overview
Setting: nine western US cities (Sacramento, Modesto, Fresno, Visalia, Bakersfield, Reno, Las Vegas, Salt Lake City, Provo), 2001–2020

Exposure: daily city-level PM2.5 from harmonized regulatory monitors ([EPA AQS])

Interaction: PCAP events, classified from ERA5 reanalysis using the normalized valley heat deficit (Boomsma & Holmes, 2025), analyzed as a binary PCAP indicator (pcap01)

Outcomes: daily counts of respiratory ED visits, overall and by ICD-9/ICD-10 sub-group ([e.g., all resipiratory, upper and lower airway disease, asthma, bronchitis, pneumonia])

Design: multi-city time-series analysis

Stage 1: city-specific quasi-Poisson GAMs, with linear PM2.5 term interacted with PCAP status and adjusted for [long-term trend, temperature, day of week, holidays, ...]

Stage 2: random-effects multivariate meta-analysis of the city-specific spline coefficients (mixmeta) to reconstruct pooled exposure–response curves, plus metafor::rma() pooling of linear estimates

Complementary analysis: hierarchical GAMs fit with mgcv::bam() using factor-smooth specifications following Pedersen et al. (2019)

R scripts
