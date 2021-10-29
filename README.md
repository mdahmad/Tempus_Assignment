# Tempus Bioinformatics Challenge
**Name: Maria Ahmad** <br>
**Email: ahmad.maria98@gmail.com** <br>
**Phone: (678) 800 8135** <br>
**Position applying for: Bioinformatics Analyst** <br>
#### Description
Goal: Prototype a variant annotation tool. Each variant in the VCF must be annotated with the following pieces of information:
1. Type of variation (substitution, insertion, CNV, etc.) and their effect (missense, silent,
intergenic, etc.). If there are multiple effects, annotate with the most deleterious
possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from ExAC API (API documentation is available here:
http://exac.hms.harvard.edu/).
6. Any additional annotations that you feel might be relevant.

#### Annotation workflow
The notebook contains three steps: 1) exploratory analysis, 2) reading and preparing the VCF, and 3) annotating the VCF. 
### Part 0: Input
The notebook takes the VCF file and the VCF output filename prefix as inputs. 
`vcf = "Challenge_data_renamed.vcf" <br>
file_prefix = "Challenge_data"`
### Part 1: Exploratory analysis
This part examines the information in the INFO column, the information in the FORMAT columns, the number of variants in the file, the number of samples. This helps with figuring out what tags will be used to calculate the percentage of reads supporting the variant. 
### Part 2: Read and prepare the VCF
The VCF is read into a pandas dataframe, the VCF comments are save in a list, and the VCF dataframe ID column is populated with IDs in the **CHR-REF-ALT-POS** format.
### Part 3: Annotate the VCF
<p> From the exploratory analysis, the following information is already in the VCF: variant type, depth of sequence coverage at site of variation, and number of reads supporting the variant. The percentage of reads supporting the variant can be calculated using the depth of sequence and the number of reads supporting the variant. The allele frequency and variant impact need to be pulled from the ExAC API output. </p> 
<p> For each variant, the variant information is pulled from the ExAC browser. The allele frequency, variant consequences, variant's most severe consequence, and variant's gene/transcript biological classification is pulled. For each variant, there can be multiple VEP annotations, thus multiple 'most severe' consequence. The most severe consequences from the list is the first one. Source of rankings: https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences. </p> 
