# vcf2avinput

vcf2avinput is a Python script for converting VCF files to ANNOVAR avinput format. It aims to replace convert2annovar in [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) software. The VCF files generated by following callers are supported：GATK4（mutect2, haplotypecaller), DeepVariant, Strelka, Varscan.

## Installation and Dependency

Currently it's only a simple python script so no need to install. Just copy to you folder and run it. Nonetheless it needs [vcfpy](https://pypi.org/project/vcfpy/) package to run.

```bash
pip install vcfpy
```
## Input and output

It takes a left-normalized and siplited VCF file as input and it outputs an avinput format file. After that you can use ANNOVAR to annotate the avinput file.

```bash
table_annovar.pl YOUR_INPUT.avinput [DATA_BASE_FOLDER] -buildver [hg19|hg38] --outfile [OUTPUT_PREFIX] -remove -protocol [DATA_BASE_NAMES] -operation [OPERATION] -nastring . --otherinfo --polish
```
To make a left-normalized and splitted VCF file you can use bcftools as recommended [here](https://annovar.openbioinformatics.org/en/latest/articles/VCF/)

```bash
bcftools norm --force -m-both -o ex1.step1.vcf ex1.vcf.gz  # split
bcftools norm --force -f human_g1k_v37.fasta -o ex1.step2.vcf ex1.step1.vcf  # left normalize
```

## Usage

```python
# generate my_sample.avinput
vcf2avinput.py my_sample.vcf

# only extract sample S1
vcf2avinput.py -s S1 my_sample.vcf

# manually assign type as mutect2
vcf2avinput.py -t mutect my_sample.vcf

# show help
vcf2avinput.py -h
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
