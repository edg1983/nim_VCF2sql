PRAGMA journal_mode = MEMORY;
PRAGMA synchronous = OFF;
PRAGMA threads = 8;

/* RECESSIVE SPLICE AFFECTING VARS */
/* homozygous vars 
- spliceAI score >= 0.6 
- sup rec = tot affected 
- hom unaffected = 0 */
DROP TABLE IF EXISTS rare_vars_splice;
CREATE TABLE rare_vars_splice AS
SELECT DISTINCT v.pid,v.gene,
    g.gado_perc,
    v.chr,v.start,v.ref,v.alt,v.consequence,v.aa_change,
    v.max_pop_af,v.cohort_af,v.cohort_het,v.cohort_hom,
    v.SpliceAI_SNP_SpliceAI_max,v.SpliceAI_INDEL_SpliceAI_max,v.CADD_PhredScore,
    v.hom_aff,v.het_aff,v.het_unaff,v.tot_aff,f.n_all
FROM rare_vars AS v 
INNER JOIN family AS f ON f.pid = v.pid 
INNER JOIN genes_hpo AS g ON (g.gene = v.gene AND g.pid = v.pid) 
WHERE 
    consequence NOT IN ("enhancer_variant","silencer_variant","insulator_variant","bivalent_variant") AND
    (
        (SpliceAI_SNP_SpliceAI_max >= 0.6 AND SpliceAI_SNP_SpliceAI_max != "NA") OR 
        (SpliceAI_INDEL_SpliceAI_max >= 0.6 AND SpliceAI_INDEL_SpliceAI_max != "NA")
    ) AND 
    (sup_rec = tot_aff AND hom_unaff = 0)
    ORDER BY v.pid ASC, g.gado_perc DESC;
/* Run Time: real 538.747 user 3.186301 sys 1.032161 */

/* denovo vars  
- spliceAI score >= 0.6 
- sup_dnm = tot affected 
- het unaffected = 0
- cohort_het <= 5;
- n_all >= 3 (at least trios) */
DROP TABLE IF EXISTS rare_vars_splice_dnm;
CREATE TABLE rare_vars_splice_dnm AS
SELECT DISTINCT v.pid,v.gene,
    g.gado_perc,
    v.chr,v.start,v.ref,v.alt,v.consequence,v.aa_change,
    v.max_pop_af,v.cohort_af,v.cohort_het,v.cohort_hom,
    v.SpliceAI_SNP_SpliceAI_max,v.SpliceAI_INDEL_SpliceAI_max,v.CADD_PhredScore,
    v.sup_dnm,v.het_aff,v.het_unaff,v.tot_aff,f.n_all
FROM rare_vars AS v 
INNER JOIN family AS f ON f.pid = v.pid 
INNER JOIN genes_hpo AS g ON (g.gene = v.gene AND g.pid = v.pid) 
WHERE
    consequence NOT IN ("enhancer_variant","silencer_variant","insulator_variant","bivalent_variant") AND
    (
        (SpliceAI_SNP_SpliceAI_max >= 0.6 AND SpliceAI_SNP_SpliceAI_max != "NA") OR 
        (SpliceAI_INDEL_SpliceAI_max >= 0.6 AND SpliceAI_INDEL_SpliceAI_max != "NA")
    ) AND
    (sup_dnm = tot_aff AND het_unaff = 0) AND n_all >= 3 AND
    cohort_het <= 5
    ORDER BY v.pid ASC, g.gado_perc DESC;
/* Run Time: real 42.712 user 3.248706 sys 1.016674 */

/* comphet selection */
/* 1. Select combos where
- at least 1 var has spliceAI >= 0.6 but not coding var
- num_aff = tot_aff */
DROP TABLE IF EXISTS comphet_splice;
CREATE TABLE comphet_splice AS 
SELECT * FROM comphet_rare WHERE 
    num_aff = tot_aff AND 
    (
        v1 IN (SELECT var_dbID FROM rare_vars WHERE
            (
                (SpliceAI_SNP_SpliceAI_max >= 0.6 AND SpliceAI_SNP_SpliceAI_max != "NA") OR 
                (SpliceAI_INDEL_SpliceAI_max >= 0.6 AND SpliceAI_INDEL_SpliceAI_max != "NA")
            ) AND
            consequence IN ("splice_acceptor_variant","splice_donor_variant","5_prime_UTR_variant","promoter_variant","intron_variant","synonymous_variant","non_coding_transcript_exon_variant","splice_region_variant")
        ) OR 
        v2 IN (SELECT var_dbID FROM rare_vars WHERE 
            (
                (SpliceAI_SNP_SpliceAI_max >= 0.6 AND SpliceAI_SNP_SpliceAI_max != "NA") OR 
                (SpliceAI_INDEL_SpliceAI_max >= 0.6 AND SpliceAI_INDEL_SpliceAI_max != "NA")
            ) AND
            consequence IN ("splice_acceptor_variant","splice_donor_variant","5_prime_UTR_variant","promoter_variant","intron_variant","synonymous_variant","non_coding_transcript_exon_variant","splice_region_variant")
        )
    );
/* Run Time: real 840.648 user 6.490805 sys 2.008852 */
CREATE INDEX comphet_splice_v1_idx ON comphet_splice (v1);
CREATE INDEX comphet_splice_v2_idx ON comphet_splice (v2);

/* 2. Further select combos where the other var is coding */
DROP TABLE IF EXISTS comphet_splice_coding;
CREATE TABLE comphet_splice_coding AS 
SELECT * FROM comphet_splice WHERE 
    v1 in (SELECT var_dbID FROM rare_vars WHERE 
        consequence IN ("stop_gained","frameshift_variant","exonic_sv","missense_variant","start_lost","stop_lost")) OR 
    v2 IN (SELECT var_dbID FROM rare_vars WHERE 
        consequence IN ("stop_gained","frameshift_variant","exonic_sv","missense_variant","start_lost","stop_lost"));
/* Run Time: real 66.541 user 21.944704 sys 41.492084 */

/* 3. Annotate vars in comphet */
DROP TABLE IF EXISTS comphet_splice_coding_anno;
CREATE TABLE comphet_splice_coding_anno AS 
SELECT ch.pid, ch.gene, g.gado_perc, ch.var_dbID,
    v.chr,v.start,v.end,v.ref,v.alt,v.consequence,v.aa_change,
    v.SpliceAI_SNP_SpliceAI_max, v.SpliceAI_INDEL_SpliceAI_max,v.CADD_PhredScore,
    v.max_pop_af,v.cohort_af,v.cohort_het,v.cohort_hom,
    ch.num_aff,ch.tot_aff,f.n_all
FROM comphet_splice_coding AS ch
INNER JOIN family AS f ON f.pid = ch.pid 
INNER JOIN vars AS v ON (v.var_dbID = ch.v1 OR v.var_dbID = ch.v2)
INNER JOIN genes_hpo AS g ON (g.gene = ch.gene AND g.pid = ch.pid)
ORDER BY 
    ch.pid ASC,
    g.gado_perc DESC,
    ch.var_dbID ASC;
/* Run Time: real 2.809 user 0.004028 sys 0.006039 */


/* DUMP RESULT TABLES */
.mode tabs
.header ON

.output splice_recessive_vars.tsv
SELECT * FROM rare_vars_splice;
.output splice_dnm_vars.tsv
SELECT * FROM rare_vars_splice_dnm;
.output splice_comphet_vars.tsv
SELECT * FROM comphet_splice_coding_anno;


DROP TABLE IF EXISTS comphet_2splice;
CREATE TABLE comphet_2splice AS 
SELECT * FROM comphet_rare WHERE 
    num_aff = tot_aff AND 
    (
        v1 IN (SELECT var_dbID FROM rare_vars WHERE
            (
                (SpliceAI_SNP_SpliceAI_max >= 0.6 AND SpliceAI_SNP_SpliceAI_max != "NA") OR 
                (SpliceAI_INDEL_SpliceAI_max >= 0.6 AND SpliceAI_INDEL_SpliceAI_max != "NA")
            ) AND
            consequence IN ("splice_acceptor_variant","splice_donor_variant","5_prime_UTR_variant","promoter_variant","intron_variant","synonymous_variant","non_coding_transcript_exon_variant","splice_region_variant")
        ) AND 
        v2 IN (SELECT var_dbID FROM rare_vars WHERE 
            (
                (SpliceAI_SNP_SpliceAI_max >= 0.6 AND SpliceAI_SNP_SpliceAI_max != "NA") OR 
                (SpliceAI_INDEL_SpliceAI_max >= 0.6 AND SpliceAI_INDEL_SpliceAI_max != "NA")
            ) AND
            consequence IN ("splice_acceptor_variant","splice_donor_variant","5_prime_UTR_variant","promoter_variant","intron_variant","synonymous_variant","non_coding_transcript_exon_variant","splice_region_variant")
        )
    );