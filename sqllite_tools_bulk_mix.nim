# sqllite_tools
# Author: Edoardo Giacopuzzi
# Process cohort varan output into sqllite DB

import argparse
import times
import math
import strformat
import db_sqlite
import tables
import sequtils
import zip/gzipfiles

const
    interval = 50000
    COMMENT_CHAR = "#"
    SEP = "\t"
    #Following section defines col names from the var2reg idx file
    COL_PID = "pedigree"
    COL_GENEFILE = "gene_file"
    COL_VARFILE = "variant_file"
    COL_COMPHET = "comphet_file"
    COL_PEDFILE = "pedigree_file"
    COL_GADOFILE = "gado_file"
    COL_N_SAMPLES = "n_all_samples"
    COL_N_AFFECTED = "n_affected_samples"
    #Cols in GADO file containing gene symbol and gene rank
    GADO_GENE_COL = 1
    GADO_RANK_COL = 2 
    #Cols idx used to create unique var ID (chr,pos,ref,alt) from vars.tsv.gz
    VAR_COLS_ID = @[3,4,6,7]
    VAR_RECID_COL = 0
    #Following section defines columns to be extracted from genes.tsv.gz
    #1. Cols containing HPO specific scores
    GENE_COLS_HPO = @[1,2,3]
    #2. Cols containing gene specific scores
    GENE_COLS_GENE = concat(@[1],toSeq(7..13))
    #Columns in comphet.tsv.gz containing gene, v1, v2 and num_aff
    COMPHET_COLS = @[1,2,3,4]
    COMPHET_V1_COL = 2 
    COMPHET_V2_COL = 3
    COMPHET_RECID_COL = 0
    #Following section defines the schemas for SQLite tables
    VARSCORES_TABLE_SCHEMA = {
        "var_uniqueID": "VARCHAR",
        "chr": "VARCHAR",
        "start": "INTEGER",
        "end": "INTEGER",
        "ref": "VARCHAR",
        "alt": "VARCHAR",
        "var_type": "VARCHAR",
        "known_ids": "VARCHAR",
        "max_pop_af": "FLOAT",
        "cohort_af": "FLOAT",
        "cohort_calls": "INTEGER",
        "cohort_het": "INTEGER",
        "cohort_hom": "INTEGER",
        "CADD_PhredScore": "FLOAT",
        "CADD_RawScore": "FLOAT",
        "DANN_score": "FLOAT",
        "CCR_pct": "FLOAT",
        "LoF_tolerance": "FLOAT",
        "FIRE_score": "FLOAT",
        "LinSight": "FLOAT",
        "ReMM_score": "FLOAT",
        "dbscSNV_ada": "FLOAT",
        "dbscSNV_rf": "FLOAT",
        "PhyloP100": "FLOAT",
        "REVEL_score": "FLOAT",
        "ExPECTO_MAX_VALUE": "FLOAT",
        "NCBoost_score": "FLOAT",
        "dbNSFP_PP2_HDIV": "FLOAT",
        "dbNSFP_PP2_HVAR": "FLOAT",
        "SpliceAI_SNP_SpliceAI_max": "FLOAT",
        "SpliceAI_SNP_SpliceAI_DP": "FLOAT",
        "SpliceAI_INDEL_SpliceAI_max": "FLOAT",
        "SpliceAI_INDEL_SpliceAI_DP": "FLOAT",
        "MCAP_score": "FLOAT",
        "Reg_constraint": "FLOAT",
        "LowMappability": "INTEGER",
        "LowComplexity": "INTEGER",
        "SegDup": "INTEGER",
        "RefSeqReg": "INTEGER",
        "DNase": "INTEGER",
        "TFBS": "INTEGER",
        "PanelApp": "INTEGER",
        "TopVariableGenes": "INTEGER",
        "UCNE": "INTEGER"
    }.toOrderedTable
    VAR_TABLE_SCHEMA = {
        "var_dbID": "VARCHAR UNIQUE",
        "var_uniqueID": "VARCHAR",
        "pid": "VARCHAR",
        "tot_aff": "INTEGER",
        "rec_id": "VARCHAR",
        "var_id": "VARCHAR",
        "gene": "VARCHAR",
        "aa_change": "VARCHAR",
        "reg_id": "VARCHAR",
        "db_source": "VARCHAR",
        "reg_type": "VARCHAR",
        "closest_gene_dist": "INTEGER",
        "consequence": "VARCHAR",
        "d_score": "FLOAT",
        "a_score": "FLOAT",
        "sup_rec": "INTEGER",
        "sup_dom": "INTEGER",
        "sup_dnm": "INTEGER",
        "sup_comphet_mum": "VARCHAR",
        "sup_comphet_dad": "VARCHAR",
        "num_hq_calls": "INTEGER",
        "hom_aff": "INTEGER",
        "hom_unaff": "INTEGER",
        "het_aff": "INTEGER",
        "het_unaff": "INTEGER"
    }.toOrderedTable
    GENES_TABLE_SCHEMA = {
        "gene": "VARCHAR UNIQUE",
        "pLI_exac": "FLOAT",
        "pLI_gnomad": "FLOAT",
        "max_transcript_length": "INTEGER",
        "GDI_score": "FLOAT",
        "GDI_phred": "FLOAT",
        "EDS": "FLOAT",
        "RVIS": "FLOAT"
    }.toOrderedTable
    GENES_HPO_SCHEMA = {
        "pid": "VARCHAR",
        "gene": "VARCHAR",
        "gado_zscore": "FLOAT",
        "exomiser_gene_pheno_score": "FLOAT",
        "gado_perc": "FLOAT"
    }.toOrderedTable
    COMPHET_TABLE_SCHEMA = {
        "var_dbID": "VARCHAR UNIQUE",
        "pid": "VARCHAR",
        "tot_aff": "INTEGER",
        "gene": "VARCHAR",
        "v1": "VARCHAR",
        "v2": "VARCHAR",
        "num_aff": "INTEGER"
    }.toOrderedTable
    SAMPLES_TABLE_SCHEMA = {
        "pid": "VARCHAR",
        "sample_id": "VARCHAR",
        "father_id": "VARCHAR",
        "mother_id": "VARCHAR",
        "sex": "INTEGER",
        "affected": "INTEGER"
    }.toOrderedTable
    FAMILY_TABLE_SCHEMA = {
        "pid": "VARCHAR UNIQUE",
        "n_all": "INTEGER",
        "n_affected": "INTEGER"
    }.toOrderedTable
    GENOS_TABLE_SCHEMA = {
        "var_uniqueID": "VARCHAR",
        "sampleID": "VARCHAR",
        "pid": "VARCHAR",
        "GT": "INTEGER",
        "GQ": "INTEGER"
    }.toOrderedTable

#get seq of keys from Table
proc `~`[A,B](tb: OrderedTable[A, B]): seq[A] =
    var keys_list: seq[A]
    for k in tb.keys:
        keys_list.add(k)
    return keys_list

#get seq of values from Table
#proc `@`[A,B](tb: OrderedTable[A, B]): seq[B] =
#    var values_list: seq[string]
#    for v in tb.values:
#        values_list.add(v)
#    return values_list

proc elapsedtime(start_time: float32, end_time: float32, item: string = "item"): float {.discardable.} =
    var elapsed_time = end_time - start_time
    var mins = floor(elapsed_time / 60)
    var secs = elapsed_time - (60 * mins)
    var hours = floor(mins / 60)
    mins = mins - (60 * hours)
    echo item, " in ", fmt"{hours:02.0f}", ":", fmt"{mins:02.0f}", ":", fmt"{secs:02.2f}"
    #return(elapsed_time)

proc tokenize(value: string, sep: string): seq[string] {.inline.} =
    var strvalue: string
    strvalue = value
    strvalue.stripLineEnd
    result = strvalue.split(sep)

proc selectIdx[T](s: seq[T], pred: proc(x: T): bool {.closure.}): seq[int] {.inline.} =
    result = newSeq[int]()
    for i in 0 ..< s.len:
        if pred(s[i]):
            result.add(i)

proc subsetSeq(inseq: seq[string], idx: seq[int]): seq[string] =
    var subset_seq: seq[string]
    for i in idx:
        subset_seq.add(inseq[i])
    return subset_seq

proc checkSuccess(success: bool, mode: string, process_name: string = "process") =
    if not success:
        case mode
        of "warn", "warning", "w": 
            echo "WARNING! " & process_name & " failed!"
            discard
        of "q", "quit":
            quit("FATAL! " & process_name & " failed!")
        else: discard

proc makeSQLTable(db: DbConn, table_name: string, table_schema: OrderedTable[system.string, system.string], exists_ok: bool = false, unique_vals: seq[string] = @[]): bool =
    var cols_definition: seq[string]
    var exists_sql = ""
    var unique_val = ""
    if exists_ok: exists_sql = "IF NOT EXISTS "
    if len(unique_vals) > 0: unique_val = ",UNIQUE(" & unique_vals.join(",") & ")"
    for k,v in table_schema.pairs:
        cols_definition.add(k & " " & v)
    var sql_string = cols_definition.join(",")
    let sql_command = "CREATE TABLE " & exists_sql & table_name & "(" & sql_string & unique_val & ")"  
    result = db.tryExec(sql(sql_command))

proc makeInsertSql(table_name: string, n_values: int): string {.inline} =
    let placeholders = repeat("?", n_values)
    result = "INSERT OR IGNORE INTO " & table_name & " VALUES (" & placeholders.join(",") & ")"

proc insertRecord(db: DbConn, sql_statement: SqlQuery, ins_values: seq[string]): bool {.discardable, inline.} =
    #let placeholders = repeat("?", len(ins_values))
    #var sql_command = "INSERT OR IGNORE INTO " & table_name & " VALUES (" & placeholders.join(",") & ")"
    result = db.tryExec(sql_statement, ins_values)

#Read a file with header and yield for each line a table[col_name: value]
#Lines starting with comment_char are skipped
iterator readfile(filename: string, sep: string, comment_char: string): OrderedTable[system.string, system.string] {.closure.} =
    var is_header_line = true
    var col_names: seq[string]
    var tab_out = initOrderedTable[string, string]()

    if filename.endsWith(".gz"):
        let in_stream = newGzFileStream(filename)
        while not in_stream.atEnd():
            var line = in_stream.readLine()
            if not line.startsWith(comment_char):
                var tokens = tokenize(line, sep)
                if is_header_line: 
                    col_names = tokens
                    is_header_line = false
                    continue
                for pairs in zip(col_names, tokens):
                    let (col_name, field) = pairs
                    tab_out[col_name] = field
                yield tab_out

    else:
        let in_stream = open(filename)
        defer: in_stream.close()

        var line: string
        while in_stream.readLine(line):
            if not line.startsWith(comment_char):
                var tokens = tokenize(line, sep)
                if is_header_line: 
                    col_names = tokens
                    is_header_line = false
                    continue
                for pairs in zip(col_names, tokens):
                    let (col_name, field) = pairs
                    tab_out[col_name] = field
                yield tab_out

#Return the header line of a file as seq[string]
#This is the first line seen after skipping lines starting with comment_char
proc readHeader(filename: string, sep: string, comment_char: string): seq[string] =
    var header_line: seq[string]
    if filename.endsWith(".gz"):
        let in_stream = newGzFileStream(filename)
        while not in_stream.atEnd():
            var line = in_stream.readLine()
            if not line.startsWith(comment_char):
                header_line = tokenize(line, sep)
                break
    else:
        let in_stream = open(filename)
        defer: in_stream.close()

        var line: string
        while in_stream.readLine(line):
            if not line.startsWith(comment_char):
                header_line = tokenize(line, sep)
                break
    return header_line

# Parse arguments
var p = newParser("cohort varan to SQLite"):
    option("-i", "--input", help="path to cohort varan index file")
    option("-o", "--outdb", help="output db files")
    flag("-v", "--verbose", help="verbose log")
    flag("-g", "--genos", help="store also genotypes")
var 
   argv = commandLineParams()
   opts = p.parse(argv)

# Quit if help is used
if opts.help:
    quit "END of help message"

let start_time = cpuTime()
let db = open(opts.outdb, "", "", "")
db.exec(sql"PRAGMA synchronous = OFF")
db.exec(sql"PRAGMA journal_mode = MEMORY")

#Make sql insert statements
let
    genes_table_insert = sql(makeInsertSql("genes", len(GENES_TABLE_SCHEMA)))
    vars_table_insert = sql(makeInsertSql("vars", len(VAR_TABLE_SCHEMA)))
    vars_scores_insert = sql(makeInsertSql("vars_scores", len(VARSCORES_TABLE_SCHEMA)))
    comphet_table_insert = sql(makeInsertSql("comphet", len(COMPHET_TABLE_SCHEMA)))
    family_table_insert = sql(makeInsertSql("family", len(FAMILY_TABLE_SCHEMA)))
    samples_table_insert = sql(makeInsertSql("samples", len(SAMPLES_TABLE_SCHEMA)))
    geneshpo_table_insert = sql(makeInsertSql("genes_hpo", len(GENES_HPO_SCHEMA)))
    genos_table_insert = sql(makeInsertSql("genos", len(GENOS_TABLE_SCHEMA)))

echo "Creating SQLite tables"
checkSuccess(
    makeSQLTable(db, "vars", VAR_TABLE_SCHEMA, true), "q", "vars table creation")
checkSuccess(
    makeSQLTable(db, "vars_scores", VARSCORES_TABLE_SCHEMA, true), "q", "scores table creation")
checkSuccess(
    makeSQLTable(db, "genes", GENES_TABLE_SCHEMA, true), 
    "q", "genes table creation")
checkSuccess(
    makeSQLTable(db, "genes_hpo", GENES_HPO_SCHEMA, true, @["pid", "gene"]), 
    "q", "genes_hpo table creation")
checkSuccess(
    makeSQLTable(db, "comphet", COMPHET_TABLE_SCHEMA, true), 
    "q", "comphet table creation")
checkSuccess(
    makeSQLTable(db, "samples", SAMPLES_TABLE_SCHEMA, true, @["pid", "sample_id"]), 
    "q", "samples table creation")
checkSuccess(
    makeSQLTable(db, "family", FAMILY_TABLE_SCHEMA, true), 
    "q", "family table creation")
if opts.genos:
    checkSuccess(
        makeSQLTable(db, "genos", GENOS_TABLE_SCHEMA, true, @["var_uniqueID", "sampleID"]), 
        "q", "genos table creation")

db.exec(sql"BEGIN TRANSACTION")

var 
    n_processed_lines: int
    start_time_item: float
    start_time_pid: float
    is_header_line: bool
    n_cases = 0

for tokens in readfile(opts.input, "\t", "#"):
    start_time_pid = cpuTime()
    var pid = tokens[COL_PID]
    inc n_cases
    echo "### " & $n_cases & ": " & tokens["pedigree"] & " ###"
    
    let family_table_values = @[pid, tokens[COL_N_SAMPLES], tokens[COL_N_AFFECTED]]
    insertRecord(db, family_table_insert, family_table_values)

    #PROCESS GADO FILES
    echo "# Processing GADO score files"
    n_processed_lines = 0
    start_time_item = cpuTime()
    var gado_perc_score = initTable[string, string]()

    is_header_line = true
    var in_stream = newGzFileStream(tokens[COL_GADOFILE])
        
    while not in_stream.atEnd():
        var line = in_stream.readLine()
        if line.startsWith(COMMENT_CHAR): continue
        var line_values = tokenize(line, SEP)
        if is_header_line: 
            is_header_line = false
            continue
        inc n_processed_lines
    
        #made a table with gene_symbol: rank
        gado_perc_score[line_values[GADO_GENE_COL]] = line_values[GADO_RANK_COL] 
    
    in_stream.close()
    elapsedtime(start_time_item, cpuTime(), "\t" & $n_processed_lines & " processed")
    
    echo "\tcalculating GADO percentiles"
    for k, v in gado_perc_score.mpairs:
        let perc_val = 1 - (parseInt(v) / n_processed_lines)
        v = fmt"{perc_val:1.3f}"

    #PROCESS PED FILE
    echo "# Processing PED file"
    n_processed_lines = 0
    start_time_item = cpuTime()

    in_stream = newGzFileStream(tokens[COL_PEDFILE])   
    
    while not in_stream.atEnd():
        var line = in_stream.readLine()
        if line.startsWith(COMMENT_CHAR): continue
        var line_values = tokenize(line, SEP)
        inc n_processed_lines
    
        #add pid and insert values into SQLite
        var samples_table_values = line_values
        
        #Insert records into SQLite tables
        insertRecord(db, samples_table_insert, samples_table_values)
        #if not success: inc n_failed_lines
    
    in_stream.close()
    elapsedtime(start_time_item, cpuTime(), "\t" & $n_processed_lines & " processed")

    #PROCESS GENE FILES
    echo "# Processing gene file"
    n_processed_lines = 0
    start_time_item = cpuTime()
    
    is_header_line = true
    in_stream = newGzFileStream(tokens[COL_GENEFILE])
        
    while not in_stream.atEnd():
        let line = in_stream.readLine()
        if line.startsWith(COMMENT_CHAR): continue
        var line_values = tokenize(line, SEP)
        if is_header_line: 
            is_header_line = false
            continue
        inc n_processed_lines
    
        #subset columns to go into various tables
        var gene_table_values = subsetSeq(line_values, GENE_COLS_GENE)
        var hpo_table_values = subsetSeq(line_values, GENE_COLS_HPO)
        
        #Insert pid where necessary
        hpo_table_values.add(gado_perc_score.getOrDefault(hpo_table_values[0], "NA"))
        hpo_table_values.insert(@[pid], 0)
        
        #Insert records into SQLite tables
        insertRecord(db, genes_table_insert, gene_table_values)
        #if not success: inc n_failed_lines
        insertRecord(db, geneshpo_table_insert, hpo_table_values)
        #if not success: inc n_failed_lines 
    
    in_stream.close()
    elapsedtime(start_time_item, cpuTime(), "\t" & $n_processed_lines & " processed")

    #PROCESS VAR FILES
    echo "# Processing vars file"
    let header = readHeader(tokens[COL_VARFILE], "\t", "#")
    let GT_GQ_cols = selectIdx(header, proc(x: string): bool = x.contains("GT_") or x.contains("GQ_"))
    let vars_table_cols = selectIdx(header, proc(x: string): bool = x in ~VAR_TABLE_SCHEMA)
    let vars_scores_cols = selectIdx(header, proc(x: string): bool = x in ~VARSCORES_TABLE_SCHEMA)

    var
        GT_cols: seq[int]
        GQ_cols: seq[int]
        samples: seq[string]
    if opts.genos:
        GT_cols = selectIdx(header, proc(x: string): bool = x.contains("GT_"))
        GQ_cols = selectIdx(header, proc(x: string): bool = x.contains("GQ_"))
        samples = subsetSeq(header, GT_cols)
        samples.applyIt(it.replace("GT_",""))

    n_processed_lines = 0
    start_time_item = cpuTime()
    is_header_line = true
    in_stream = newGzFileStream(tokens[COL_VARFILE])
    var
        GT_values: seq[string]
        GQ_values: seq[string]
        line: string
        line_values: seq[string]
        vars_values: seq[string]
        vars_scores_values: seq[string]

    while not in_stream.atEnd():
        line = in_stream.readLine()
        if line.startsWith(COMMENT_CHAR): continue
        line_values = tokenize(line, SEP)
        if is_header_line: 
            is_header_line = false
            continue
        inc n_processed_lines
        
        #make unique ids for var
        let var_id = line_values[VAR_COLS_ID[0]] & "_" &
            line_values[VAR_COLS_ID[1]] & "_" &
            line_values[VAR_COLS_ID[2]] & "_" &
            line_values[VAR_COLS_ID[3]]
        let var_dbid = pid & "_" & line_values[VAR_RECID_COL]

        if opts.genos:
            GT_values = subsetSeq(line_values, GT_cols)
            GQ_values = subsetSeq(line_values, GQ_cols)

        #subset columns to go into vars table 
        #line_values.delete(GT_GQ_cols[0], GT_GQ_cols[^1])
        vars_values = subsetSeq(line_values, vars_table_cols)
        vars_scores_values = subsetSeq(line_values, vars_scores_cols)

        #insert var_dbid and pid
        vars_values.insert(@[var_dbid, var_id, pid, tokens[COL_N_AFFECTED]], 0)
        vars_scores_values.insert(@[var_id], 0)

        #Insert records into SQLite table
        insertRecord(db, vars_table_insert, vars_values)
        insertRecord(db, vars_scores_insert, vars_scores_values)

        if opts.genos:
            for i,s in samples:
                insertRecord(db, genos_table_insert, @[var_id, s, pid, GT_values[i], GQ_values[i]])

        if opts.verbose:
            if floorMod(n_processed_lines, interval) == 0:
                elapsedtime(start_time_item, cpuTime(), "\t" & $n_processed_lines & " lines")   
        
    in_stream.close()
    elapsedtime(start_time_item, cpuTime(), "\t" & $n_processed_lines & " processed")

    #PROCESS COMPHET FILES
    echo "# Processing comphet file"
    n_processed_lines = 0
    start_time_item = cpuTime()
    is_header_line = true
    in_stream = newGzFileStream(tokens[COL_COMPHET])

    while not in_stream.atEnd():
        var line = in_stream.readLine()
        if line.startsWith(COMMENT_CHAR): continue
        var line_values = tokenize(line, SEP)
        if is_header_line: 
            is_header_line = false
            continue
        inc n_processed_lines

        #create a unique comphet_id = pid_recid
        let comphet_dbid = pid & "_" & line_values[COMPHET_RECID_COL]
        
        #update recid of single vars to reflect new varids
        line_values[COMPHET_V1_COL] = pid & "_" & line_values[COMPHET_V1_COL]
        line_values[COMPHET_V2_COL] = pid & "_" & line_values[COMPHET_V2_COL]
    
        insertRecord(db, comphet_table_insert, 
            concat(@[comphet_dbid, pid, tokens[COL_N_AFFECTED]], subsetSeq(line_values, COMPHET_COLS)))

        if opts.verbose:
            if floorMod(n_processed_lines, interval) == 0:
                elapsedtime(start_time_item, cpuTime(), "\t" & $n_processed_lines & " lines")   
        
    in_stream.close()
    elapsedtime(start_time_item, cpuTime(), "\t" & $n_processed_lines & " processed")

    elapsedtime(start_time_pid, cpuTime(), "# case completed")

db.exec(sql"END TRANSACTION")

#Create SQLLite indexes
echo "### Creating SQLLite indexes ###"
start_time_item = cpuTime()
var sql_tabs = @["vars","comphet"]
for t in sql_tabs:
    var sql_command = "CREATE INDEX " & t & "_vardbid_idx ON " & t & " (var_dbID)"
    db.exec(sql(sql_command))
    sql_command = "CREATE INDEX " & t & "_pid_idx ON " & t & " (pid)"
    db.exec(sql(sql_command))

sql_tabs = @["genes","genes_hpo"]
for t in sql_tabs:
    let sql_command = "CREATE INDEX " & t & "_gene_idx ON " & t & " (gene)"
    db.exec(sql(sql_command))

db.exec(sql"CREATE INDEX vars_varid_idx ON vars (var_uniqueID)")
db.exec(sql"CREATE INDEX vars_consequence_idx ON vars (consequence)")
db.exec(sql"CREATE INDEX vars_gene_idx ON vars (gene)")
db.exec(sql"CREATE INDEX scores_cohortaf_idx ON vars_scores (cohort_af)")
db.exec(sql"CREATE INDEX scores_maf_idx ON vars_scores (max_pop_af)")
db.exec(sql"CREATE INDEX scores_comboaf_idx ON vars_scores (max_pop_af, cohort_af)")
db.exec(sql"CREATE INDEX scores_cadd_idx ON vars_scores (CADD_PhredScore)")
db.exec(sql"CREATE INDEX scores_remm_idx ON vars_scores (ReMM_score)")
db.exec(sql"CREATE INDEX scores_varid_idx ON vars_scores (var_uniqueID)")
db.exec(sql"CREATE INDEX comphet_v1_idx ON comphet (v1)")
db.exec(sql"CREATE INDEX comphet_v2_idx ON comphet (v2)")
db.exec(sql"CREATE INDEX family_pid_idx ON family (pid)")   
db.exec(sql"CREATE INDEX genes_hpo_pidgene_idx ON genes_hpo (pid, gene)")
db.exec(sql"CREATE INDEX genes_hpo_zscore_idx ON genes_hpo (gado_zscore)")
db.exec(sql"CREATE INDEX genes_hpo_gadoperc_idx ON genes_hpo (gado_perc)")

if opts.genos:
    db.exec(sql"CREATE INDEX genos_varid_idx ON genos (var_uniqueID)")
    db.exec(sql"CREATE INDEX genos_sampleid_idx ON genos (sampleID)")
    db.exec(sql"CREATE INDEX genos_gtgq_idx ON genos (GT, GQ)")
elapsedtime(start_time_item, cpuTime(), "indexes created")

db.close()
elapsedtime(start_time, cpuTime(), "Completed")