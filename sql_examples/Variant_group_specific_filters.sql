SELECT DISTINCT v.pid,v.gene,
g.gado_perc,
v.chr,v.start,v.ref,v.alt,v.var_type,v.consequence,v.aa_change,
v.max_pop_af,v.cohort_af,v.cohort_het,v.cohort_hom,
v.SpliceAI_SNP_SpliceAI_max,v.SpliceAI_INDEL_SpliceAI_max,
v.CADD_PhredScore, v.ReMM_score, v.LinSight,
v.hom_aff,v.tot_aff,f.n_all
FROM rare_vars AS v 
INNER JOIN family AS f ON f.pid = v.pid 
INNER JOIN genes_hpo AS g ON (g.gene = v.gene AND g.pid = v.pid) 
WHERE v.pid LIKE "%VHL%" AND v.gene == "VHL" AND
(
    (consequence IN ("enhancer_variant","silencer_variant","insulator_variant","bivalent_variant") AND
        (ReMM_score >= 0.8 AND ReMM_score != "NA") OR 
        (LinSight >= 0.5 AND LinSight != "NA")
    ) OR
    var_type IN ("DEL","DUP","INS","INV")
    OR
    (consequence IN ("synonymous_variant","intron_variant","splice_region_variant") AND
        (SpliceAI_SNP_SpliceAI_max >= 0.5 AND SpliceAI_SNP_SpliceAI_max != "NA") OR 
        (SpliceAI_INDEL_SpliceAI_max >= 0.5 AND SpliceAI_INDEL_SpliceAI_max != "NA")
    ) OR
    (consequence IN ("missense_variant") AND
        CADD_PhredScore >= 20
    ) OR
    consequence IN ("frameshift_variant","splice_acceptor_variant","splice_donor_variant","stop_gained","start_lost","stop_lost") 
)
AND (het_aff = tot_aff AND het_unaff = 0);