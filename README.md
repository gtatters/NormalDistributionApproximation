# Distribution Behaviour & Sampling Inference

An interactive Shiny app for BIOL 3P96 (Biostatistics) at Brock University.

## What this app does

This app lets you draw samples from a wide range of distributions and observe
how the shape of the data, the distribution of paired differences, and the
sampling distribution of the mean all change with sample size and distribution
type. A Q-Q plot with a 95% confidence envelope and a Shapiro-Wilk normality
test are provided to help you judge whether a sample looks normally distributed.

## Distributions available

Normal, t (df = 3), Log-normal, Uniform, Exponential, Chi-square (df = 3),
Gamma (shape = 2), Beta (2,5), Poisson (λ = 3), Binomial (n = 20, p = 0.3).

## Modes

- **Single experiment** — histogram and Q-Q plot for one sample
- **Paired differences (X − Y)** — distribution of differences between two
  independent samples from the same distribution
- **Many experiments (CLT demo)** — sampling distribution of the mean across
  many repeated experiments, with a normal overlay

## Learning goals

- See how skewed or heavy-tailed distributions affect sample histograms and Q-Q plots
- Understand that paired differences are often more symmetric than the raw data
- Watch the Central Limit Theorem in action: means of repeated samples approach
  normality even when the underlying distribution is not normal

## Course context

Developed for BIOL 3P96 — Biostatistics, Brock University.
Built with R and Shiny (base R graphics only).