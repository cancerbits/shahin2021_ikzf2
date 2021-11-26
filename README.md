# Supplementary code repository for: Germline biallelic mutation affecting the transcription factor Helios causes pleiotropic defects of immunity

Authors: Tala Shahin,<sup>1,2,3,4</sup> Hye Sun Kuehn,<sup>5</sup> Mohamed R. Shoeb,<sup>1</sup> Lisa M. Gawriyski,<sup>6</sup> Sarah Giuliani,<sup>1,2</sup> Peter Repiscak,<sup>1</sup> Birgit Hoeger,<sup>1,2</sup> Özlem Yüce Petronczki,<sup>1,2</sup> Sevgi Köstel Bal,<sup>1,2</sup> Samaneh Zoghi,<sup>1,2</sup> Jasmin Dmytrus,<sup>1,2</sup> Davide Seruggia,<sup>1</sup> Irinka Castanon,<sup>1</sup> Nima Rezaei,<sup>7,8,9</sup> Markku Varjosalo,<sup>6</sup> Florian Halbritter,<sup>1</sup> Sergio Rosenzweig,<sup>5</sup> Kaan Boztug<sup>1,2,3</sup>

Affiliations:

<sup>1</sup>St. Anna Children’s Cancer Research Institute (CCRI), Vienna, Austria

<sup>2</sup>Ludwig Boltzmann Institute for Rare and Undiagnosed Diseases, Vienna, Austria

<sup>3</sup>CeMM Research Center for Molecular Medicine of the Austrian Academy of Sciences, Vienna, Austria

<sup>4</sup>Department of Pediatrics and Adolescent Medicine, Medical University of Vienna, Vienna, Austria

<sup>5</sup>Immunology Service, Department of Laboratory Medicine, Clinical Center, National Institutes of Health, Bethesda, Maryland, 20814, USA

<sup>6</sup>Institute of Biotechnology, Helsinki Institute of Life Science, Proteomics Unit, University of Helsinki, Helsinki, Finland

<sup>7</sup>Research Center for Immunodeficiencies, Children's Medical Center, Tehran University of Medical Sciences, Tehran, Iran

<sup>8</sup>Department of Immunology, School of Medicine, Tehran University of Medical Sciences, Tehran, Iran

<sup>9</sup>Network of Immunity in Infection, Malignancy and Autoimmunity (NIIMA), Universal Scientific Education and Research Network (USERN), Tehran, Iran

## Abstract:

Helios, a member of the Ikaros family of transcription factors, is predominantly expressed within the T-cell lineage, including developing thymocytes, activated T cells and regulatory T cells (Tregs). Studies in mice have emphasized its role in maintenance of Treg immunosuppressive functions through stabilizing Foxp3 expression and silencing the IL-2 locus. Yet, its contribution to human immune homeostasis and the precise mechanisms by which Helios regulates other T-cell subsets remain unresolved. Here, we studied a patient with recurrent respiratory infections and hypogammaglobulinemia and identified a germline homozygous missense mutation in IKZF2 encoding Helios (p.Ile325Val). We show that HeliosI325V retains DNA-binding and dimerization properties, but loses interaction with several partners, including epigenetic remodelers HDAC1, HDAC3 and the ATAC complex. Single-cell RNA-sequencing of peripheral blood mononuclear cells revealed gene expression signatures indicative of a shift towards pro-inflammatory, effector-like status in patient T cells, while analysis of patient Tregs has shown a dysregulated transcriptome. Intriguingly, patient T cells were shown to be terminally differentiated, with a pronounced defect in proliferation and a marked reduction in IL-2 production upon stimulation, that was rescued by reverting the mutation in Helios back to wild type in patient-derived T cells. Collectively, we identify a novel germline-encoded inborn error of immunity and highlight a role for Helios in conventional T cells, whereby interactions with specific binding partners is necessary to mediate the transcriptional programs that enable T-cell homeostasis in health and disease.

# Repository structure:

* `src/RNA-seq` - data analysis code for bulk- and single-cell RNA-seq

# Analysis workflow:

* RNA-seq: Separate scripts are dedicated to each analysis steps to faciliate understanding and maintaining the workflow. The analysis pipeline is defined in `src/RNA-seq/initialize.Rmd` where scripts are loaded sequentially.  For optimal reproducibility, we used a Docker container `mo-ikzf2-v1.Dockerfile`, which contains R and all dependent libraries preinstalled.

# Data preparation and loading

1. Download and, if needed, preprocess the EGA data into a defined input directory. Metadata table can be found in the Supplementary table S4 Excel worksheet accompanying the manuscript.
2. For RNA-seq analysis, set the path to the input directory at the top of `src/RNA-seq/initialize.Rmd` and use `mo-ikzf2-v1.Dockerfile` to run the analysis workflow.

Input data can be obtained from EGA (link below).

## Links:

* European Genome-Phenome Archive (EGA): <a href="https://ega-archive.org/studies/EGAS00001005675">https://ega-archive.org/studies/EGAS00001005675</a>
* Paper: <a href="https://doi.org/10.1126/sciimmunol.abe3981">https://doi.org/10.1126/sciimmunol.abe3981</a>
