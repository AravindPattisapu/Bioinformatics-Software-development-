##SOFTWARE DEVELOPMENT
Created for the software development group project - QMUL MSc Bioinformatics 

Web application that retrieves SNP information for a genomic region of interest in Homo sapiens chromosome 21 and calculates specific summary statistics for specified populations. Users can search by RS value, gene name or aliases and finally genomic positions. They can select the summary statistics of interest alongside which populations they would like data outputs for. Additionally, for the sliding window statistics calculations, the user can input their window size preference which would output the required results. 
The summary statistics are returned in a table with the option of viewing graphs for statistics versus populations. Moreover, the user is presented with interactive graphs for each statistic, including the sliding window values against the genomic positions. The graphs allow the selection/deselection of populations for comparing results and can be downloaded as a PNG to present key findings. The summary statistics can be downloaded as a CSV file to store for later usage by the user. 
The summary statistics provided by the website are nucleotide diversity, haplotype diversity, Tajima’s D and FST. Information for the British, Finnish, Colombian, Punjabi and Telugu populations are present in the database and is utilised by the web application. The website is simple to navigate through and produces results which can be utilised by users in developing biological research.

Project background:
Population genomics is the science that studies the genetic diversity within and between 
populations. By analysing the similarities and differences of genetic variation between individuals 
from different populations, it is possible to infer historical (both neutral and adaptive) events that 
characterized the evolution of said species. Inferences from population genomic data are widely 
applied in conservation genetics (e.g., molecular monitoring of endangered species), evolutionary 
biology (e.g., demographic reconstruction), and precision medicine, among many other fields.
The discipline has received notable attention in recent years thanks to the technological 
advancements of sequencing genomes at large-scale. The advent of next-generation sequencing 
(NGS) technologies has allowed researchers to obtain vast amounts of genomic data for many 
samples belonging to several populations. After sequencing, mapping, filtering and variant calling, 
data is usually reported in variant call format (VCF) files which contain genotype (or haplotype) 
information for sequenced sampled at each called single nucleotide polymorphism (SNP). In the case 
of human population genomics, the most popular database is The International Genome Sample 
Resource which collects sequencing data from multiple human populations, mostly from the 1000 
Genomes Project data set.
