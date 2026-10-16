# Psychometric evaluation of the PHQ-9 and a simulated computer adaptive test

**Author:** Jasmine Saxena

## Motivation

As many research assistant job adverts ask for experience in psychometrics and adaptive testing, I wanted to test whether my skills were up to standard as expected in these roles. 

It also linked back to my work at Waterloo Community Counselling (WCC), where I studied the effectiveness of cultural and linguistic client–therapist matching. Matching went beyond demographic characteristics to psychological fit, which is linked to a stronger therapeutic alliance. This relied on clients’ self-report questionnaires. However, we noticed that with longer questionnaires, respondent burden increased and answers became vaguer, a pattern consistent with respondent fatigue and satisficing. I wanted to see whether computer adaptive testing could reach the same measurement precision with fewer items.

## Approach

Using PHQ-9 data from 5,068 US adults (NHANES 2017–2018), I:

- checked reliability (alpha and omega)
- tested whether the PHQ-9 measures one thing (confirmatory factor analysis)
- checked whether it works the same way across sex (measurement invariance)
- fitted an item response theory model to see how informative each item is
- simulated a computer adaptive test to see how many items people would actually need

## Findings 

The PHQ-9 was reliable (omega = 0.93), measured a single factor, and worked the same way for men and women.

The suicidal-thoughts item had the lowest item-total correlation (0.38), which suggested it was a weak item. But its factor loading from the ordinal methods was 0.75, in line with the others. Because 96% of people answered “Not at all”, there was very little variation to correlate, which pulled the correlation down. This showed me why it matters to treat questionnaire items as ordinal rather than continuous.

The PHQ-9 measured precisely for people with at least mild symptoms (a total score of 5 or more), but much less for the 62% of the sample who scored below that range. This showed up in the adaptive test too. For people with mild to severe symptoms, it reached a precise score in about three items instead of nine. For people with few symptoms, it still needed about eight, because the PHQ-9 doesn’t have items aimed at that end of the scale.

Full results are in `output/results.txt`file.

## Limitations

- Generalisability: The sample was US adults from the general population (mean age 50). Arguably, adolescents may interpret items differently, and clinical samples have a much wider spread of severity, so the results may not apply to services like CAMHS. Repeating the analysis with the adolescent version (PHQ-A) in a clinical sample would address this.
- No criterion validity: Without a clinical diagnosis to compare against, I could test how consistently the PHQ-9 measures, but not whether it correctly identifies people with depression. This would need data from a structured diagnostic interview.
- Optimistic CAT performance: The item parameters were estimated and tested on the same respondents, which overstates how well the adaptive test would perform on new people. Calibrating on one half of the sample and testing on the other (cross-validation) would give a more realistic estimate.
- Post-hoc simulation: The CAT re-used answers from the full nine-item form, thereby assuming people answer the same way regardless of test length. It can't show whether a shorter test reduces respondent fatigue, which would need a separate study which gives participants the adaptive version directly.

## Data

Centers for Disease Control and Prevention (CDC), National Center for Health Statistics (NCHS). *National Health and Nutrition Examination Survey Data, 2017–2018.* Hyattsville, MD: U.S. Department of Health and Human Services. [Depression screener documentation](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DPQ_J.htm)

## Instructions

1. Install [R](https://cran.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/).
2. Download this repository and open `analysis.R` in RStudio.
3. Set the working directory to the file's location (Session → Set Working Directory → To Source File Location) and click Source.
4. The script installs any missing packages and downloads the data. Results are saved in the `output` folder.

[I wrote the code with help from an AI assistant (Claude), then ran it, checked the output and interpreted the results myself.]

## Packages

`haven`, `psych`, `lavaan`, `mirt`, `mirtCAT`
