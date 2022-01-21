PRAGMA journal_mode = MEMORY;
PRAGMA synchronous = OFF;
PRAGMA threads = 8;

/* denovo vars  
- sup_dnm = tot affected 
- het unaffected = 0
- cohort_het <= 5;
- n_all >= 3 (at least trios) */
DROP TABLE IF EXISTS dnm_DDD_genes;
CREATE TABLE dnm_DDD_genes AS
SELECT DISTINCT v.pid,v.gene,
    g.gado_perc,
    v.chr,v.start,v.ref,v.alt,v.consequence,v.aa_change,
    v.max_pop_af,v.cohort_af,v.cohort_het,v.cohort_hom,
    v.SpliceAI_SNP_SpliceAI_max,v.SpliceAI_INDEL_SpliceAI_max,v.CADD_PhredScore,
    v.ReMM_score, v.CADD_PhredScore,
    v.sup_dnm,v.het_aff,v.tot_aff,f.n_all
FROM rare_vars AS v 
INNER JOIN family AS f ON f.pid = v.pid 
INNER JOIN genes_hpo AS g ON (g.gene = v.gene AND g.pid = v.pid) 
WHERE
    v.gene IN ('AP2S1','ARF1','ARHGAP35','ATP6V0A1','DDX23','DSP','FBXW7','GIGYF1','GNB2','H3F3A','HNRNPD','KCNK3','MFN2','MIB1','MMGT1','MSL2','PRKAR1B','PRPF8','PSMC5','RAB14','SATB1','SPEN','SRRM2','TCF7L2','TFE3','U2AF2','UPF1','ZFHX4') AND
    (sup_dnm = tot_aff OR (het_aff = tot_aff AND tot_aff = 1)) AND 
    (het_unaff = 0 AND hom_unaff = 0) AND
    n_all >= 3 AND
    cohort_het <= 5
    ORDER BY v.pid ASC, g.gado_perc DESC;

/* DUMP RESULT TABLES */
.mode tabs
.header ON

.output dnm_DDD_genes.tsv
SELECT * FROM dnm_DDD_genes;
