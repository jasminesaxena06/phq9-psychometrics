# Psychometric evaluation of the PHQ-9 and a simulated computer adaptive test

**Author:** Jasmine Saxena

## Aim

This project examines how well the Patient Health Questionnaire-9 (PHQ-9) measures depression in a large population sample. It also tests whether a computer adaptive test (CAT) built from the same items could reach a precise score with fewer questions.

## Data

The data come from the [National Health and Nutrition Examination Survey (NHANES) 2017–2018](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DPQ_J.htm). This is a publicly available US survey that administers the PHQ-9 to adults aged 18 and over. The script downloads the depression screener (`DPQ_J`) and demographics (`DEMO_J`) files automatically. Responses coded "refused" or "don't know" are treated as missing, and only respondents with all nine items answered are analysed.

## Methods

1. **Descriptives:** response distributions for each item and the total score.
2. **Reliability:** Cronbach's alpha, ordinal alpha from polychoric correlations, and corrected item-total correlations.
3. **Confirmatory factor analysis:** a one-factor model for ordinal data (WLSMV estimator, `lavaan`), with ordinal omega calculated from the standardised loadings.
4. **Measurement invariance by sex:** a configural model compared with a model holding loadings and thresholds equal across groups (scaled chi-square difference test and change in CFI).
5. **Item response theory:** Samejima's graded response model (`mirt`), with item fit (S-X²), overall fit (C2), and item and test information curves.
6. **Computer adaptive test simulation:** 500 real response patterns run through a CAT (`mirtCAT`). The CAT chooses the most informative next item and stops when the standard error reaches 0.3, 0.4 or 0.5. CAT scores are compared with full-test scores.

## Results

<!-- Fill this in from output/results.txt after running the script. -->

- **Sample:**
- **Reliability:**
- **Factor structure:**
- **Invariance by sex:**
- **IRT:**
- **CAT simulation:**

## Limitations

- The sample is US adults, so the findings may not generalise to children and young people or to clinical samples.
- The CAT is calibrated and tested on the same data, which makes its performance look somewhat better than it would be in a new sample.
- NHANES survey weights are not applied, so the estimates describe this sample rather than the US population.
- The PHQ-9 has a strong floor effect in general population samples, which limits measurement precision at low severity.

## How to run

1. Install [R](https://cran.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/).
2. Download this repository and open the folder in RStudio.
3. Run `analysis.R`. Missing packages are installed automatically, and the run takes a few minutes.
4. The results are saved in the `output` folder: `results.txt`, the CSV tables, and the figures.

## Packages

`haven`, `psych`, `lavaan`, `mirt`, `mirtCAT`
