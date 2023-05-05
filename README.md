# vcfmerge

* A Python script that merges somatic SNV/InDel calls from three VCF files:
  - MuTect2 VCF (SNVs + InDels)
  - Strelka2 SNV VCF
  - Strelka2 InDel VCF

* Two resulting VCF files are produced
  - __<tumor_id>_<normal_id>_all.vcf.gz__ - contains all calls, both rejected and PASSed
  - __<tumor_id>_<normal_id>_somatic.vcf.gz__ - contains all somatic calls (PASS by one or both of Strelka2 and MuTect2)


