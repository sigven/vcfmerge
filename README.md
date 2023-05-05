## vcfmerge

_vcfmerge_ - a small Python script that merges somatic SNV/InDel calls from three individual VCF files:

   - MuTect2 VCF (SNVs + InDels)
   - Strelka2 SNV VCF
   - Strelka2 InDel VCF

Two resulting VCF files are produced:

  - __<tumor_id>_<normal_id>_all.vcf.gz__ - contains all calls, both rejected and PASSed
  - __<tumor_id>_<normal_id>_somatic.vcf.gz__ - contains all somatic calls (PASS by one or both of Strelka2 and MuTect2)

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

- The tumor and control sample identifiers must match the names provided in the individual VCF files (sample columns)
