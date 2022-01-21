PRAGMA journal_mode = MEMORY;
PRAGMA synchronous = OFF;
PRAGMA threads = 8;

/* MAKE RARE VARS TABLES */
/* Select rare vars pop_af < 0.01 and cohort_af < 0.05 */
CREATE TABLE rare_vars AS SELECT * FROM vars WHERE 
    max_pop_af < 0.01 AND cohort_af < 0.05;
CREATE INDEX rarevars_vardbid_idx ON rare_vars (var_dbID);
CREATE INDEX vars_hetafftotaff_idx ON rare_vars (het_aff, tot_aff);
CREATE INDEX vars_homafftotaff_idx ON rare_vars (hom_aff, tot_aff);
CREATE INDEX vars_spliceAIindel_idx ON rare_vars (SpliceAI_INDEL_SpliceAI_max);
CREATE INDEX vars_spliceAIsnp_idx ON rare_vars (SpliceAI_SNP_SpliceAI_max);
CREATE INDEX vars_dnmtotaff_idx ON rare_vars (sup_dnm, tot_aff);

/* Select comphet combos where both vars are rare */
CREATE TABLE comphet_rare AS SELECT * FROM comphet WHERE 
    v1 IN (SELECT var_dbID FROM rare_vars) AND 
    v2 IN (SELECT var_dbID FROM rare_vars);
CREATE INDEX comphet_rare_v1_idx ON comphet_rare (v1);
CREATE INDEX comphet_rare_v2_idx ON comphet_rare (v2); 
CREATE INDEX comphet_rare_naff_idx ON comphet_rare (num_aff, tot_aff);