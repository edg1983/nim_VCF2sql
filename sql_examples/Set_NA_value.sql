PRAGMA journal_mode = MEMORY;
PRAGMA synchronous = OFF;
PRAGMA threads = 8;

.timer ON

UPDATE vars SET sup_comphet_dad = NULL WHERE sup_comphet_dad == "";
UPDATE vars SET sup_comphet_mum = NULL WHERE sup_comphet_mum == "";
UPDATE vars SET max_pop_af = NULL WHERE max_pop_af == "";
UPDATE vars SET closest_gene_dist = NULL WHERE closest_gene_dist == "";
UPDATE vars SET CADD_PhredScore = NULL WHERE CADD_PhredScore == "";
UPDATE vars SET CADD_RawScore = NULL WHERE CADD_RawScore == "";
UPDATE vars SET DANN_score = NULL WHERE DANN_score == "";
UPDATE vars SET CCR_pct = NULL WHERE CCR_pct == "";
UPDATE vars SET LoF_tolerance = NULL WHERE LoF_tolerance == "";
UPDATE vars SET FIRE_score = NULL WHERE FIRE_score == "";
UPDATE vars SET LinSight = NULL WHERE LinSight == "";
UPDATE vars SET ReMM_score = NULL WHERE ReMM_score == "";
UPDATE vars SET dbscSNV_ada = NULL WHERE dbscSNV_ada == "";
UPDATE vars SET dbscSNV_rf = NULL WHERE dbscSNV_rf == "";
UPDATE vars SET PhyloP100 = NULL WHERE PhyloP100 == "";
UPDATE vars SET REVEL_score = NULL WHERE REVEL_score == "";
UPDATE vars SET ExPECTO_MAX_VALUE = NULL WHERE ExPECTO_MAX_VALUE == "";
UPDATE vars SET NCBoost_score = NULL WHERE NCBoost_score == "";
UPDATE vars SET dbNSFP_PP2_HDIV = NULL WHERE dbNSFP_PP2_HDIV == "";
UPDATE vars SET dbNSFP_PP2_HVAR = NULL WHERE dbNSFP_PP2_HVAR == "";
UPDATE vars SET SpliceAI_SNP_SpliceAI_max = NULL WHERE SpliceAI_SNP_SpliceAI_max == "";
UPDATE vars SET SpliceAI_SNP_SpliceAI_DP = NULL WHERE SpliceAI_SNP_SpliceAI_DP == "";
UPDATE vars SET SpliceAI_INDEL_SpliceAI_max = NULL WHERE SpliceAI_INDEL_SpliceAI_max == "";
UPDATE vars SET SpliceAI_INDEL_SpliceAI_DP = NULL WHERE SpliceAI_INDEL_SpliceAI_DP == "";
UPDATE vars SET MCAP_score = NULL WHERE MCAP_score == "";
UPDATE vars SET Reg_constraint = NULL WHERE Reg_constraint == "";

UPDATE genes SET pLI_exac = NULL WHERE pLI_exac == "";
UPDATE genes SET pLI_gnomad = NULL WHERE pLI_gnomad == "";
UPDATE genes SET max_transcript_length = NULL WHERE max_transcript_length == "";
UPDATE genes SET GDI_score = NULL WHERE GDI_score == "";
UPDATE genes SET GDI_phred = NULL WHERE GDI_phred == "";
UPDATE genes SET EDS = NULL WHERE EDS == "";
UPDATE genes SET RVIS = NULL WHERE RVIS == "";

UPDATE genes_hpo SET gado_zscore = NULL WHERE gado_zscore == "";
UPDATE genes_hpo SET exomiser_gene_pheno_score = NULL WHERE exomiser_gene_pheno_score == "";
UPDATE genes_hpo SET gado_perc = NULL WHERE gado_perc == "";