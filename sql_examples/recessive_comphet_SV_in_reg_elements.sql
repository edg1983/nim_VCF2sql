PRAGMA journal_mode = MEMORY;
PRAGMA synchronous = OFF;
PRAGMA threads = 8;

/* DELETIONS AFFECTING DISTANT REGULATORY ELEMENTS */
/* homozygous vars
- var type is DEL,DUP,INS,INV,DEL:MEI
- sup rec = tot affected 
- hom unaffected = 0 */
DROP TABLE IF EXISTS rare_vars_reg_SV;
CREATE TABLE rare_vars_reg_SV AS
SELECT DISTINCT v.pid,v.gene,
    g.gado_perc,
    v.chr,v.start,v.end,v.ref,v.alt,v.consequence,
    v.max_pop_af,v.cohort_af,v.cohort_het,v.cohort_hom,
    v.hom_aff,v.het_aff,v.het_unaff,v.tot_aff,f.n_all
FROM rare_vars AS v 
INNER JOIN family AS f ON f.pid = v.pid 
INNER JOIN genes_hpo AS g ON (g.gene = v.gene AND g.pid = v.pid) 
WHERE
    var_type IN ("DEL","DUP","INS","INV","DEL:MEI") AND
    consequence IN ("silencer_variant","enhancer_variant","insulator_variant","bivalent_variant") AND
    (sup_rec = tot_aff AND hom_unaff = 0)
    ORDER BY v.pid ASC, g.gado_perc DESC;
/* Run Time: real 38.277 user 10.985971 sys 20.664649 */

/* denovo vars
- var type is DEL,DUP,INS,INV,DEL:MEI
- sup_dnm = tot affected 
- het unaffected = 0
- n_all >= 3 (at least trios) 
- cohort het <= 5 */
DROP TABLE IF EXISTS rare_vars_reg_SV_dnm;
CREATE TABLE rare_vars_reg_SV_dnm AS
SELECT DISTINCT v.pid,v.gene,
    g.gado_perc,
    v.chr,v.start,v.end,v.ref,v.alt,v.consequence,
    v.max_pop_af,v.cohort_af,v.cohort_het,v.cohort_hom,
    v.sup_dnm,v.het_aff,v.het_unaff,v.tot_aff,f.n_all
FROM rare_vars AS v 
INNER JOIN family AS f ON f.pid = v.pid 
INNER JOIN genes_hpo AS g ON (g.gene = v.gene AND g.pid = v.pid) 
WHERE
    var_type IN ("DEL","DUP","INS","INV","DEL:MEI") AND
    consequence IN ("silencer_variant","enhancer_variant","insulator_variant","bivalent_variant") AND
    (sup_dnm = tot_aff AND het_unaff = 0) AND n_all >= 3 AND
    cohort_het <= 5
    ORDER BY v.pid ASC, g.gado_perc DESC;
/* Run Time: real 33.571 user 10.779041 sys 20.519248 */


/* comphet selection */
/* 1. Select combos where at least 1 var is a SV affecting reg elements  */
DROP TABLE IF EXISTS comphet_reg_sv;
CREATE TABLE comphet_reg_sv AS 
SELECT * FROM comphet_rare WHERE 
    num_aff = tot_aff AND
    ( 
        v1 IN (SELECT var_dbID FROM rare_vars WHERE
            var_type IN ("DEL","DUP","INS","INV","DEL:MEI") AND
            consequence IN ("silencer_variant","enhancer_variant","insulator_variant","bivalent_variant")
        ) OR 
        v2 IN (SELECT var_dbID FROM rare_vars WHERE 
            var_type IN ("DEL","DUP","INS","INV","DEL:MEI") AND
            consequence IN ("silencer_variant","enhancer_variant","insulator_variant","bivalent_variant")
        )
    );
/* Run Time: real 95.671 user 22.141626 sys 41.826534 */
CREATE INDEX comphet_regsv_v1_idx ON comphet_reg_sv (v1);
CREATE INDEX comphet_regsv_v2_idx ON comphet_reg_sv (v2);

/* 2. Select combos where the other var is coding */
DROP TABLE IF EXISTS comphet_reg_sv_coding;
CREATE TABLE comphet_reg_sv_coding AS 
SELECT * FROM comphet_reg_sv WHERE 
    v1 in (SELECT var_dbID FROM rare_vars WHERE 
        consequence IN ("stop_gained","frameshift_variant","missense_variant","start_lost","stop_lost","splice_acceptor_variant","splice_donor_variant")) OR 
    v2 IN (SELECT var_dbID FROM rare_vars WHERE 
        consequence IN ("stop_gained","frameshift_variant","missense_variant","start_lost","stop_lost","splice_acceptor_variant","splice_donor_variant"));
/* Run Time: real 66.172 user 22.657623 sys 41.113654 */

/* 3. Annotate vars in comphet */
DROP TABLE IF EXISTS comphet_reg_sv_coding_anno;
CREATE TABLE comphet_reg_sv_coding_anno AS 
SELECT ch.pid, ch.gene, g.gado_perc, ch.var_dbID,
    v.chr,v.start,v.end,v.ref,v.alt,v.consequence,v.aa_change,v.CADD_PhredScore,
    v.max_pop_af,v.cohort_af,v.cohort_het,v.cohort_hom,
    ch.num_aff,ch.tot_aff,f.n_all
FROM comphet_reg_sv_coding AS ch
INNER JOIN family AS f ON f.pid = ch.pid 
INNER JOIN vars AS v ON (v.var_dbID = ch.v1 OR v.var_dbID = ch.v2)
INNER JOIN genes_hpo AS g ON (g.gene = ch.gene AND g.pid = ch.pid)
ORDER BY 
    ch.pid ASC,
    g.gado_perc DESC,
    ch.var_dbID ASC;
/* Run Time: real 6.865 user 0.005728 sys 0.021117 */


/* DUMP RESULT TABLES */
.mode tabs
.header ON
.output regSV_recessive_vars.tsv
SELECT * FROM rare_vars_reg_SV;
.output regSV_dnm_vars.tsv
SELECT * FROM rare_vars_reg_SV_dnm;
.output regSV_comphet_vars.tsv
SELECT * FROM comphet_reg_sv_coding_anno;