## vcfmerge

_vcfmerge_ - a small Python script that merges somatic SNV/InDel calls from three individual (bgzipped) VCF files:

   - MuTect2 VCF (SNVs + InDels)
   - Strelka2 SNV VCF
   - Strelka2 InDel VCF

Two resulting VCF files are produced:

  - **<tumor_id>_<normal_id>_vcfmerge_all.vcf** - contains all calls, both rejected and PASSed
  - **<tumor_id>_<normal_id>_vcfmerge_somatic.vcf** - contains all somatic calls (PASS by one or both of Strelka2 and MuTect2)

The script produces multiple dedicated VCF INFO tags in the resulting output files to simplify downstream annotation, most importantly:

 - __TDP__ - total sequencing depth of variant site in tumor (i.e. _DP_ in tumor sample, MuTect2 values have priority)
 - __TVAF__ - allelic fraction of alternate allele in tumor (i.e. _AF_ in tumor sample, MuTect2 values have priority)
 - __CDP__ - total sequencing depth of variant site in control (i.e. _DP_ in control sample, MuTect2 values have priority)
 - __CVAF__ - allelic fraction of alternate allele in control (i.e. _AF_ in control sample, MuTect2 values have priority)
 - __CALLERS__ - any of _mutect2_, _strelka2_, or _mutect2,strelka2_ (called by both)
 - __MNV_SUPPORT_STRELKA__ - as Strelka2 does not properly call multinucleotide variants (MNVs or block substitutions), the script gathers consecutive SNVs (all PASS) from Strelka2 calls when they are found as an MNV (PASS) in MuTect2

### Usage

	usage: vcfmerge.py [-h] [--mutect_vcf MUTECT_VCF] [--strelka_snv_vcf STRELKA_SNV_VCF] [--strelka_indel_vcf STRELKA_INDEL_VCF] [--compress] [--force_overwrite] tumor_sample_id control_sample_id output_dir

	Merge somatic calls (SNVs/InDels) from multiple VCF files into a single VCF

	positional arguments:
	  tumor_sample_id       Sample ID for the tumor sample
	  control_sample_id     Sample ID for the control sample
	  output_dir            Directory for output files

	optional arguments:
	  -h, --help            show this help message and exit
	  --mutect_vcf MUTECT_VCF
	                        Bgzipped VCF input file with somatic query variants (SNVs) called with MuTect (version 2.x). (default: None)
	  --strelka_snv_vcf STRELKA_SNV_VCF
	                        Bgzipped VCF input file with somatic query variants (SNVs) called with Strelka (version 2.x). (default: None)
	  --strelka_indel_vcf STRELKA_INDEL_VCF
	                        Bgzipped VCF input file with somatic query variants (InDels) called with Strelka (version 2.x) (default: None)
	  --compress            Compress output VCF with bgzip + tabix (default: False)
	  --force_overwrite     Overwrite existing output files (default: False)

### Notes

- The tumor and control sample identifiers provided as input arguments _must_ match the names provided in the individual VCF files (sample columns)
- Note that the script currently fully ignores multi-allelic sites (i.e. sites with multiple alternate alleles). However, from our experience so far, it seems that limited sites of this nature contain somatic events with a _PASS_ status
