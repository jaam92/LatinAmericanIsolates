# Understanding the hidden complexity of Latin American population isolates

Data files and code accompanying:

> Mooney, J. A., Huber, C. D., Service, S., Sul, J. H., Marsden, C. D., Zhang, Z., Sabatti, C., Ruiz-Linares, A., Bedoya, G., Costa Rica/Colombia Consortium for Genetic Investigation of Bipolar Endophenotypes, Freimer, N. B., & Lohmueller, K. E. (2018). Understanding the hidden complexity of Latin American population isolates. *American Journal of Human Genetics*, 103(5), 707–726. [Full text on PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC6218714/)

This repository contains the summary files used to compare runs of homozygosity (ROH), identity-by-descent (IBD) sharing, and the burden of deleterious variation in whole-genome sequences from two Latin American population isolates, Antioquia, Colombia (`CO`) and the Central Valley of Costa Rica (`CR`), with populations from the 1000 Genomes Project (build hg19).

## Populations

| Code | Population |
|---|---|
| `CO` | Antioquia, Colombia (isolate) |
| `CR` | Central Valley, Costa Rica (isolate) |
| `CLM` | Colombians from Medellín, Colombia (1000 Genomes) |
| `MEX` / `MXL` | Mexican ancestry from Los Angeles, USA (1000 Genomes; files use `MEX`, plots relabel it `MXL`) |
| `PEL` | Peruvians from Lima, Peru (1000 Genomes) |
| `PUR` | Puerto Ricans from Puerto Rico (1000 Genomes) |
| `CEU` | Utah residents with Northern and Western European ancestry (1000 Genomes) |
| `FIN` | Finnish in Finland (1000 Genomes) |
| `YRI` | Yoruba in Ibadan, Nigeria (1000 Genomes) |

## Data availability

| Resource | Source | Availability |
|---|---|---|
| Colombian (`CO`) and Costa Rican (`CR`) samples | A subset of individuals genotyped on arrays in [PubMed 26712028](https://pubmed.ncbi.nlm.nih.gov/26712028/), whole-genome sequenced as described in [PMC7039961](https://pmc.ncbi.nlm.nih.gov/articles/PMC7039961/) | Summary files only (this repository) |
| Reference populations (`CLM`, `MEX`, `PEL`, `PUR`, `CEU`, `FIN`, `YRI`) | 1000 Genomes Project, hg19/GRCh37 | Publicly downloadable from the [International Genome Sample Resource](https://www.internationalgenome.org/) |

Individual-level genotype data for the Colombian and Costa Rican samples are not publicly available; only the summary files described below are provided. The 1000 Genomes data can be downloaded directly from the project.

## Repository structure

| Path | Contents |
|---|---|
| `Results_SeatSeq_results_JzWholeGenome_v4_Aug17_AllSeqIndivs_Final.txt` | Per-individual counts of variants by annotation (see below). |
| `FinalROH/` | Filtered ROH calls for all sequenced individuals and for the unrelated subset. |
| `IBDSeq/` | IBD segments per population from IBDSeq, plus `PlotIBDSegs_RunScoreTest.R`, which filters segments, computes normalized IBD scores relative to the Finnish, runs the permutation test for IBD enrichment, and plots the results. |

## Summary files

| File | Description |
|---|---|
| `Results_SeatSeq_results_JzWholeGenome_v4_Aug17_AllSeqIndivs_Final.txt` | Counts for each individual of nonsynonymous, synonymous, loss-of-function, putatively neutral, and putatively deleterious sites annotated with SeattleSeq, along with GERP-based scores. Column names give the counting method for each annotation (for example `Count_HomDer_Calls`, `Count_HetNonSyn_Calls`, `Total_DerAlleles_LOF`, `GERPScore_AllDerAlleles`). |
| `FinalROH/TrueROH_AllSeqIndivs_perIndivFINAL_Aug17.txt` | ROH in the N = 449 sequenced individuals remaining after QC: start, end, number of well-covered sites, ROH length, proportion of the ROH covered by sites, individual, chromosome, and population. |
| `FinalROH/TrueROH_Unrelated_perIndivFINAL_Aug17.txt` | ROH in the unrelated subset (N = 30 per population), in the same format. |
| `IBDSeq/<POP>_allchroms_IBDSeq_FormattedRscript.out` | IBD segments output by IBDSeq for each population (nine files) and used as input for IBDNe: start, end, length in cM and bp, SNP count, proportion covered, chromosome, and population. |

## Software

IBDSeq and IBDNe for IBD detection and effective population size inference, SeattleSeq for variant annotation, and GERP scores for conservation. The R script uses ggplot2, data.table, and scales. See the Methods of the paper for versions and parameters.

## Usage

To reproduce the IBD score comparisons, run the script from inside `IBDSeq/` so it can find the segment files:

```bash
cd IBDSeq
Rscript PlotIBDSegs_RunScoreTest.R
```

The script keeps segments between 3 and 20 cM whose SNP coverage falls within one standard deviation of the mean, normalizes summed IBD per population following Nakatsuka et al. (2017), and expresses scores relative to the Finnish.

## Citation

If you use these files, please cite the paper:

```bibtex
@article{mooney2018latinamerican,
  title   = {Understanding the Hidden Complexity of Latin American Population Isolates},
  author  = {Mooney, Jazlyn A. and Huber, Christian D. and Service, Susan and Sul, Jae Hoon and Marsden, Clare D. and Zhang, Zhongyang and Sabatti, Chiara and Ruiz-Linares, Andr{\'e}s and Bedoya, Gabriel and {Costa Rica/Colombia Consortium for Genetic Investigation of Bipolar Endophenotypes} and Freimer, Nelson B. and Lohmueller, Kirk E.},
  journal = {American Journal of Human Genetics},
  volume  = {103},
  number  = {5},
  pages   = {707--726},
  year    = {2018},
  doi     = {10.1016/j.ajhg.2018.09.013}
}
```

## Contact

For questions, please open an issue on this repository or contact the corresponding authors listed in the paper.
