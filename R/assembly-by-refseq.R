#' RMetaflow
#' Module: assembly-by-refseq
#' Initiated August 2021
#' Patrick C. Seed
#' Stanley Manne Children's Research Institute
#' Ann and Robert H. Lurie Children's Hospital
#' Northwestern University
#' MIT license
#' 

if ( !commands[STEP][[1]]$completed ){
  
  library(tidyverse)
  library(slurmR)
  library(R.utils)
  library(future.apply)
  library(dendextend)
  
  printGlue("#---- Running {names(commands[STEP][1])}", log = TRUE)
  printGlue("#---- Let's do this!")

  # 01 - Set up step variables ---------------------------------------------------
  
  in.dir = dirNormalize(current.step[[1]]$from)
  out.dir = dirNormalize(names(current.step))
  dirIf(out.dir)
  
  # 02 - Assemble all input files using megahit ----------------------------------
  
  ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
  out.samples = sapply(list.files(out.dir, pattern = current.step[[1]]$out_ext, recursive = TRUE), 
                       function(x) str_split(x, ".contigs")[[1]][1]) %>% as.vector()
  
  if( length(out.samples) > 0 ){
    ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
    remaining.samples = samples[!samples%in%out.samples]
  }else{
    remaining.samples = samples
  }
  
  dirIf(str_glue("{out.dir}/megahit"))
  
  if( length(ls) > 0 ){
    
    fseq = ls[ seq(from = 1, to = length(ls), by = 2)]
    rseq = ls[ seq(from = 2, to = length(ls), by = 2)]
    outseq = str_glue("{out.dir}/megahit/{remaining.samples}_out")
    
    job.list = list(fseq = fseq, rseq=rseq, outseq = outseq,
                    remaining.samples = remaining.samples)
    
    printGlue("#---- Preparing samples for assembly using paired reads", log = TRUE)  
    
    calls <- purrr::pmap(job.list, function(fseq, rseq, outseq, remaining.samples){
      str_glue("rm -rf {outseq} || true && \\
      {assembly.by.refseq} megahit {variables.megahit} --tmp-dir {projectdir}/temp_paired_assemblies -1 {fseq} -2 {rseq} -o {outseq} --out-prefix {remaining.samples} \\
      2>&1 | tee -a log/assembly-by-refseq.log" ) %>% str_squish()
    })
    
        # calls = calls[1] # for test purposes only
    
    dirIf("temp_paired_assemblies")
    
    if( slurm ){
      
    assembly_by_refseq_slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},  
                                            plan = "none", 
                                            overwrite = TRUE, 
                                            mc.cores = threads,
                                            njobs = length(calls),
                                            sbatch_opt = list(partition = partition, account = account, time = "02:00:00"),
                                            preamble = c("module load singularity", "module load R/4.0.3"),
                                            job_name = paste0("Slurm-", projectname, "-assembly-by-refseq"), 
                                            tmp_path = projectdir)
      
      sbatch(assembly_by_refseq_slurm)
      
      wait_slurm(assembly_by_refseq_slurm, freq = 60, timeout = -1)
      
      Slurm_clean(assembly_by_refseq_slurm)
      
    }else{
      lapply(calls, function(x) {print(x[[1]]); system(x[[1]])})  
    }
  }
  
  # 03 - Cleanup Assembly Process ------------------------------------------------------
  
  printGlue("# ---------- Clean up in Aisle 5 ", log = TRUE)
  printGlue("# ---------- Removing intermediate files from megahit", log = TRUE)
  
  
  ls.contigs = list.files(str_glue("{out.dir}/megahit"), pattern = "contigs.fa", recursive = TRUE, 
                          full.names = TRUE) %>%
    str_subset(pattern = "intermediate|final|k[0-9][0-9]", negate = TRUE)
   
  sapply(ls.contigs, function(x) system(str_glue("mv {x} {out.dir}/megahit/{basename(x)}")))
  
  dir.to.remove = list.dirs(str_glue("{out.dir}/megahit"), full.names = TRUE)[-1] %>%
    str_subset(pattern = "intermediate_contigs", negate = TRUE)
  
  sapply(dir.to.remove, function(x) system(str_glue("rm -R {x}")))
  

  # 04 - Pre- and post- pruning stats of contigs and prune small contigs and get stats using prinseq-lite ---------------------
  
  dirIf(str_glue("{out.dir}/prinseq"))
  
  in.contigs = list.files(str_glue("{out.dir}/megahit"), pattern = "contigs.fa", 
                          full.names = TRUE) %>% sort()
  prefixes = basename(in.contigs) %>% str_replace_all(., ".contigs.fa", "")
  out.contigs = str_glue("{out.dir}/prinseq/{prefixes}.prinseq.pruned.contigs")
  pre.stats = str_glue("{out.dir}/prinseq/{prefixes}.pre.prune.stats")
  post.stats = str_glue("{out.dir}/prinseq/{prefixes}.post.prune.stats")

  calls = str_glue("{assembly.by.refseq} prinseq-lite.pl -fasta {in.contigs} -stats_len -stats_info -stats_assembly > \\
    {pre.stats} && {assembly.by.refseq}  prinseq-lite.pl -fasta {in.contigs} -min_len 1000 -out_good {out.contigs} && \\
    {assembly.by.refseq} prinseq-lite.pl -fasta {out.contigs}.fasta -stats_len -stats_info -stats_assembly > {post.stats}") %>% 
      str_squish()
  
  if( slurm ){
    
    assembly_by_refseq_slurm <- Slurm_lapply(calls, function(x) {print(x[1]); system(x[1])},  
                                             plan = "none", 
                                             overwrite = TRUE, 
                                             mc.cores = threads,
                                             njobs = length(calls),
                                             sbatch_opt = list(partition = partition, account = account, time = "00:05:00"),
                                             preamble = c("module load singularity", "module load R/4.0.3"),
                                             job_name = paste0("Slurm-", projectname, "-assembly-by-refseq"), 
                                             tmp_path = projectdir)
    
    sbatch(assembly_by_refseq_slurm)
    
    wait_slurm(assembly_by_refseq_slurm, freq = 60, timeout = -1)
    
    Slurm_clean(assembly_by_refseq_slurm)
    
  }else{
    lapply(calls, function(x) {print(x[1]); system(x[1])})  
  }

  # ------- Bring together stats tables and cleanup
  ls.pre.stats = list.files( path = str_glue("{out.dir}/prinseq"), pattern = "pre.prune.stats", full.names = TRUE) %>% sort()
  prefixes = basename(ls.pre.stats) %>% str_replace_all(., "pre.prune.stats", "")
  job.list = list(samples = prefixes, ls.pre.stats = ls.pre.stats)
  
  df.pre.stats = purrr::pmap(job.list, function(samples, ls.pre.stats){
    temp.df = data.table::fread(ls.pre.stats) %>% as.data.frame()
    temp.df = data.frame(Sample = samples, temp.df, Stat = "Pre")
    return(temp.df)
  })

  df.pre.stats <- do.call(rbind, df.pre.stats)
  
  ls.post.stats = list.files( path = str_glue("{out.dir}/prinseq"), pattern = "post.prune.stats", full.names = TRUE) %>% sort()
  job.list = list(samples = prefixes, ls.post.stats = ls.post.stats)
  
  df.post.stats = purrr::pmap(job.list, function(samples, ls.post.stats){
    temp.df = data.table::fread(ls.post.stats) %>% as.data.frame()
    temp.df = data.frame(Sample = samples, temp.df, Stat = "Post")
    return(temp.df)
  })
  
  df.post.stats <- do.call(rbind, df.post.stats)
  
  df.stats = rbind(df.pre.stats, df.post.stats)
  
  colnames(df.stats) <- c("Sample", "Program", "Measure", "Value", "Stat")
  
  df.stats <- df.stats %>% select(-Program) %>%
    pivot_wider(names_from = c("Measure"), values_from = "Value")
  
  write_tsv(df.stats, file = str_glue("{out.dir}/prinseq/megahit.contigs.summary.stats.tsv"))
  
  sapply(c(ls.pre.stats, ls.post.stats), function(x) system(str_glue("rm {x}")))
  
  # 05 - Clean up contig headers to make more sense (>contig_XXX | sample) ------------

  ls.contigs = list.files(str_glue("{out.dir}/prinseq"), pattern = "pruned.contigs.fasta", full.names = TRUE) %>% sort()
  
  slurm_steps = seq_along(ls.contigs)
  
  if(slurm){
    
    Slurm.fix.contigs <- Slurm_sapply(slurm_steps, 
                                      function(x){
                                        in.contigs = ls.contigs[x]
                                        in.samples = samples[x]
                                        in.fasta = read_lines(in.contigs)
                                        contig.header.positions = grep(">", in.fasta)
                                        new.contig.headers = str_glue(">contig_{contig.header.positions}:{in.samples}")
                                        in.fasta[contig.header.positions] = new.contig.headers
                                        data.table::fwrite(list(in.fasta), file = str_glue("{out.dir}/prinseq/{in.samples}.prinseq.headerrev.contigs.fasta"), sep = " ")
                                      },
                                      export = c("ls.contigs", "samples", "out.dir"),
                                      plan = "none", 
                                      overwrite = TRUE, 
                                      mc.cores = 4,
                                      njobs = length(slurm_steps),
                                      sbatch_opt = list(partition = partition, account = account, time = "00:05:00"),
                                      preamble = c("module load singularity", "module load R/4.0.3"),
                                      job_name = paste0("Slurm-", projectname, "-fixContigs"), 
                                      tmp_path = projectdir)
                                    
    sbatch(Slurm.fix.contigs )
  
    wait_slurm(Slurm.fix.contigs , freq = 60, timeout = -1)
    
    Slurm_clean(Slurm.fix.contigs )
  }else{
    future.apply::future_lapply(as.list(slurm_steps), function(x){
      in.contigs = ls.contigs[x]
      in.samples = basename(ls.contigs[x]) %>% str_replace_all(., ".prinseq.pruned.contigs.fasta", "")
      in.fasta = read_lines(in.contigs)
      contig.header.positions = grep(">", in.fasta)
      new.contig.headers = str_glue(">{in.samples}:contig_{contig.header.positions}")
      in.fasta[contig.header.positions] = new.contig.headers
      data.table::fwrite(list(in.fasta), file = str_glue("{out.dir}/prinseq/{in.samples}.prinseq.headerrev.contigs.fasta"), sep = " ")
    })
  }

  system(str_glue("rm {out.dir}/prinseq/*prinseq.pruned*"))
  
  
  # 06 - Make BAM file output  ---------------------------------
  
  dirIf(str_glue("{out.dir}/metaerg/bams"))
  
  ls.contigs = list.files(str_glue("{out.dir}/prinseq"), pattern = "prinseq.headerrev.contigs", full.names = TRUE) %>% 
    sort()
  prefixes = basename(ls.contigs) %>% str_replace_all(., ".prinseq.headerrev.contigs.fasta", "")
  fseq = list.files(str_glue("qc/deduplicate"), pattern = "_R1.fastq", full.names = TRUE) %>% 
    str_subset(pattern = prefixes) %>%
    sort()
  rseq = list.files(str_glue("qc/deduplicate"), pattern = "_R2.fastq", full.names = TRUE) %>% 
    str_subset(pattern = prefixes) %>% sort()
  out.sam = str_glue("{out.dir}/metaerg/bams/{prefixes}.sam")
  out.bam.unsorted = str_glue("{out.dir}/metaerg/bams/{prefixes}.unsorted.bam")
  out.bam = str_glue("{out.dir}/metaerg/bams/{prefixes}.bam")
  

  calls = str_glue("[ ! -e {out.bam} ] && {assembly.by.refseq} minimap2 -ax sr -K 4G -t {threads} {ls.contigs} {fseq} {rseq} | \\
                   {assembly.by.refseq} samtools view -b - | \\
                   {assembly.by.refseq} samtools sort -@ {threads} -l 4 -o {out.bam} -") %>%
    str_squish()
 
  if(slurm){
    
    Slurm.make.bam <- Slurm_sapply(calls, function(x){print(x[[1]]); system(x[[1]])},
                                 plan = "none", 
                                 overwrite = TRUE, 
                                 mc.cores = threads,
                                 njobs = length(calls),
                                 sbatch_opt = list(partition = partition, account = account, time = "00:30:00"),
                                 preamble = c("module load singularity", "module load R/4.0.3"),
                                 job_name = paste0("Slurm-", projectname, "-bams"), 
                                 tmp_path = projectdir)
    
    sbatch(Slurm.make.bam )
    
    wait_slurm(Slurm.make.bam , freq = 60, timeout = -1)
    
    Slurm_clean(Slurm.make.bam )
  }else{
    #plan(multisession)
    future_lapply(calls, function(x){print(x[1]); system(x[1])})
  }
  
  # 07 - Run CoverM for depth of coverage per contig -------------------
  
  dirIf(str_glue("{out.dir}/metaerg/coverm"))
  
  ls.bams = list.files(str_glue("{out.dir}/metaerg/bams"), pattern = ".bam", full.names = TRUE) %>% 
    sort()
  prefixes = basename(ls.bams) %>% str_replace_all(., ".bam", "")
  coverm.depth.out = str_glue("{out.dir}/metaerg/coverm/{prefixes}.coverm.depth.txt")
  
  if(threads > 12){threads.coverm = 12}else{threads.coverm = threads}
  
  calls = str_glue("[ ! -e {coverm.depth.out} ] && {coverm} coverm contig -t {threads}  -b {ls.bams} \\
                   -m metabat -o {coverm.depth.out} --exclude-supplementary --proper-pairs-only") %>%
    str_squish()
  
  
  if(slurm){
    
    Slurm.coverm <- Slurm_sapply(calls, function(x){print(x[[1]]); system(x[[1]])},
                                      plan = "none", 
                                      overwrite = TRUE, 
                                      mc.cores = threads.coverm,
                                      njobs = length(calls),
                                      sbatch_opt = list(partition = partition, account = account, time = "00:15:00"),
                                      preamble = c("module load singularity", "module load R/4.0.3"),
                                      job_name = paste0("Slurm-", projectname, "-coverm"), 
                                      tmp_path = projectdir)
    
    sbatch(Slurm.coverm )
    
    wait_slurm(Slurm.coverm , freq = 60, timeout = -1)
    
    Slurm_clean(Slurm.coverm )
  }else{
    lapply(calls, function(x){print(x[[1]]); system(x[[1]])})
  }
 
  # 08 - Run METAERG on each sample -----------------------
  
  ls.contigs = list.files(str_glue("{out.dir}/prinseq"), pattern = "prinseq.headerrev.contigs", full.names = TRUE) %>% 
    sort()
  metaerg.out = str_glue("{out.dir}/metaerg/{prefixes}")
  metaerg.zipped.products = str_glue("{out.dir}/metaerg/{prefixes}.tar.gz")
  
  calls = str_glue("[ ! -e {metaerg.zipped.products} ] || {metaerg} metaerg.pl {variables.metaerg} --force --prefix {prefixes} \\
    --depth {coverm.depth.out} --outdir {metaerg.out} {ls.contigs}") %>% 
    str_squish()
  
  calls = str_glue("{metaerg} metaerg.pl {variables.metaerg} --force --prefix {prefixes} \\
    --depth {coverm.depth.out} --outdir {metaerg.out} {ls.contigs}") %>% 
    str_squish()
  
  if(slurm){
    
    Slurm.metaerg <- Slurm_sapply(calls, function(x){print(x[[1]]); system(x[[1]])},
                                 plan = "none", 
                                 overwrite = TRUE, 
                                 mc.cores = 24,
                                 njobs = length(calls),
                                 sbatch_opt = list(partition = partition, account = account, time = "24:00:00",
                                                   mem = "64G"),
                                 preamble = c("module load singularity", "module load R/4.0.3", mem = "64G"),
                                 job_name = paste0("Slurm-", projectname, "-metaerg"), 
                                 tmp_path = projectdir)
    
    sbatch(Slurm.metaerg )
    
    wait_slurm(Slurm.metaerg , freq = 60, timeout = -1)
    
    Slurm_clean(Slurm.metaerg )
  }else{
    lapply(calls, function(x){print(x[[1]]); system(x[[1]])})
  }
  
  
  # 09 - Aggregate the metaerg datatables across samples ---------------
  
  ls.dir = list.dirs(str_glue("{out.dir}/metaerg")) %>% 
    str_subset(pattern = "data") %>%
    sort()
  r = sapply(samples, function(x) grepl(x, ls.dir) %>% sum())
  dir.sample.names <- names(r)[r == 1]
  
  files.to.join = c("cds.gene2casgene.tab.txt", 
                    "cds.gene2kegg.tab.txt",
                    "cds.gene2pfam.tab.txt",
                    "go.cds.profile.tab.txt",
                    "master.tsv.txt",
                    "rRNA.tab.txt",
                    "cds.gene2sprot.tab.txt",
                    "kegg.cds.profile.tab.txt",
                    "master.tsv.txt",
                    "taxon.cds.profile.tab.txt",
                    "cds.gene2ec.tab.txt",
                    "cds.gene2ko.tab.txt",
                    "cds.gene2tigrfam.tab.txt",
                    "ko.cds.profile.tab.txt",
                    "taxon.lsu.profile.tab.txt",
                    "cds.gene2genomedb.tab.txt",
                    "crispr.tab.txt",    
                    "taxon.ssu.profile.tab.txt",
                    "cds.gene2go.tab.txt",
                    "cds.gene2metacyc.tab.txt",
                    "ec.cds.profile.tab.txt",
                    "master.stats.txt",
                    "pfam.cds.profile.tab.txt",
                    "tigrfam.cds.profile.tab.txt"
  )
  
  for(i in files.to.join){
    printGlue("Pulling together sample data tables for {i}")
    saved = FALSE
    ls.files = list.files(str_glue("{out.dir}/metaerg"), pattern = i, recursive = TRUE, full.names = TRUE) %>%
      str_subset(pattern = "all", negate = TRUE) %>% 
      sort()
    joined.table = lapply(ls.files, function(x) {
      dt = data.table::fread(x, check.names = FALSE, fill = TRUE) %>% as.data.frame()
      r = sapply(samples, function(q) grepl(q, x) %>% sum())
      current.sample <- names(r)[r == 1]
      #print(dim(dt))
      if(nrow(dt)>0){
        if( i == "master.stats.txt" & ncol(dt) == 2){
          dt = data.frame(Sample = current.sample, dt)
          return(dt)
        }else{
          if(i != "master.stats.txt"){
            dt = data.frame(Sample = current.sample, dt)
            return(dt)
          }
        }
      }
    })
    
    joined.table <- data.table::rbindlist(joined.table)
    
    if(dim(joined.table)[1] > 0 ){
      if(i == "master.stats.txt"){
        colnames(joined.table) = c("Sample", "CountType", "Count")
        joined.table.wide = joined.table %>%     
          pivot_wider(names_from = "Sample", values_from = "Count")
        
        data.table::fwrite(joined.table.wide, file = 
                             str_glue("{out.dir}/metaerg/all.{i}"))
        saved=TRUE
        
      }
      
      if(i=="ko.cds.profile.tab.txt"){
        colnames(joined.table) = c("Sample", "KO", "Name", "Definition",
                                   "Count", "Count_pct", 
                                   "Abund", "Abund_pct",
                                   "Abund_Sample_Depth", "Abund_Sample_Depth_pct")
        joined.table.wide.pct = joined.table %>% 
          select(Sample, KO, Abund_Sample_Depth_pct) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth_pct")
        
        joined.table.wide.count = joined.table %>% 
          select(Sample, KO, Abund_Sample_Depth) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth")
        
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Percentage.{i}"))
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Counts.{i}"))
        data.table::fwrite(joined.table, file = 
                             str_glue("{out.dir}/metaerg/all.{i}"))
        saved=TRUE
        
      }
      
      if(i == "pfam.cds.profile.tab.txt"){
        colnames(joined.table) = c("Sample", "Accession", "ID", "Definition",
                                   "Count", "Count_pct", 
                                   "Abund", "Abund_pct",
                                   "Abund_Sample_Depth", "Abund_Sample_Depth_pct")
        joined.table.wide.pct = joined.table %>% 
          select(Sample, Accession, Abund_Sample_Depth_pct) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth_pct")
        
        joined.table.wide.count = joined.table %>% 
          select(Sample, Accession, Abund_Sample_Depth) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth")
        
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Percentage.{i}"))
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Counts.{i}"))
        data.table::fwrite(joined.table, file = 
                             str_glue("{out.dir}/metaerg/all.{i}"))
        saved=TRUE
        
      }
      
      if(i == "tigrfam.cds.profile.tab.txt"){
        colnames(joined.table) = c("Sample", "Accession", "Name", "Function",
                                   "Count", "Count_pct", 
                                   "Abund", "Abund_pct",
                                   "Abund_Sample_Depth", "Abund_Sample_Depth_pct")
        joined.table.wide.pct = joined.table %>% 
          select(Sample, Accession, Abund_Sample_Depth_pct) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth_pct")
        
        joined.table.wide.count = joined.table %>% 
          select(Sample, Accession, Abund_Sample_Depth) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth")
        
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Percentage.{i}"))
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Counts.{i}"))
        data.table::fwrite(joined.table, file = 
                             str_glue("{out.dir}/metaerg/all.{i}"))
        saved=TRUE
        
      }
      
      if(i=="ec.cds.profile.tab.txt"){
        colnames(joined.table) = c("Sample", "EnzymeID", "Name",
                                   "Count", "Count_pct", 
                                   "Abund", "Abund_pct",
                                   "Abund_Sample_Depth", "Abund_Sample_Depth_pct")
        joined.table.wide.pct = joined.table %>% 
          select(Sample, EnzymeID, Abund_Sample_Depth_pct) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth_pct")
        
        joined.table.wide.count = joined.table %>% 
          select(Sample, EnzymeID, Abund_Sample_Depth) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth")
        
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Percentage.{i}"))
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Counts.{i}"))
        data.table::fwrite(joined.table, file = 
                             str_glue("{out.dir}/metaerg/all.{i}"))
        
        saved=TRUE
      }
      
      if(i=="kegg.cds.profile.tab.txt"){
        colnames(joined.table) = c("Sample", "PathwayID", "Name", "KOs",
                                   "Total_Family", "Total_Family_Found",
                                   "Count", "Count_pct", 
                                   "Abund", "Abund_pct",
                                   "Abund_Sample_Depth", "Abund_Sample_Depth_pct")
        joined.table.wide.pct = joined.table %>% 
          select(Sample, PathwayID, Abund_Sample_Depth_pct) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth_pct")
        
        joined.table.wide.count = joined.table %>% 
          select(Sample, PathwayID, Abund_Sample_Depth) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth")
        
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Percentage.{i}"))
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Counts.{i}"))
        data.table::fwrite(joined.table, file = 
                             str_glue("{out.dir}/metaerg/all.{i}"))
        
        saved=TRUE
      }
      
      if(i== "taxon.cds.profile.tab.txt" | i == "taxon.lsu.profile.tab.txt" |
         i == "taxon.ssu.profile.tab.txt"){
        colnames(joined.table) = c("Sample", "Taxon", 
                                   "Count", "Count_pct", 
                                   "Abund", "Abund_pct",
                                   "Abund_Sample_Depth", "Abund_Sample_Depth_pct")
        joined.table.wide.pct = joined.table %>% 
          select(Sample, Taxon, Abund_Sample_Depth_pct) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth_pct") %>%
          separate(col="Taxon", 
                   into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
        
        joined.table.wide.count = joined.table %>% 
          select(Sample, Taxon, Abund_Sample_Depth) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth") %>%
          separate(col="Taxon", 
                   into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
        
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Percentage.{i}"))
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Counts.{i}"))
        data.table::fwrite(joined.table, file = 
                             str_glue("{out.dir}/metaerg/all.{i}"))
        
        saved=TRUE
      }
      
      if(i=="go.cds.profile.tab.txt"){
        colnames(joined.table) = c("Sample", "GO_ID", "Name", 
                                   "Count", "Count_pct", 
                                   "Abund", "Abund_pct",
                                   "Abund_Sample_Depth", "Abund_Sample_Depth_pct")
        joined.table.wide.pct = joined.table %>% 
          select(Sample, GO_ID, Abund_Sample_Depth_pct) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth_pct")
        
        joined.table.wide.count = joined.table %>% 
          select(Sample, GO_ID, Abund_Sample_Depth) %>%
          pivot_wider(names_from = "Sample", values_from = "Abund_Sample_Depth")
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Percentage.{i}"))
        data.table::fwrite(joined.table.wide.count, file = 
                             str_glue("{out.dir}/metaerg/all.Counts.{i}"))
        data.table::fwrite(joined.table, file = 
                             str_glue("{out.dir}/metaerg/all.{i}"))
        
        saved=TRUE
      }
      
      if(i == "master.stats.txt"){
        colnames(joined.table) = c("Sample", "Category", "Count")
      }
      
      if(! saved){
        data.table::fwrite(joined.table, file = str_glue("{out.dir}/metaerg/all.{i}"))
      }
    }
  }
  
  # 10 - diamond to assign contigs to progenomes references ---------------
  
  printGlue("Assigning contigs to reference proteins", log = TRUE)
  
  ls.contigs = list.files(str_glue("{out.dir}/prinseq"), pattern = "prinseq.headerrev.contigs.fasta", full.names = TRUE) %>% 
    sort()
  out.samples = str_glue("{basename(ls.contigs) %>% str_replace_all(., '.prinseq.headerrev.contigs.fasta', '')}") %>% 
    c(.)
  
  dirIf(str_glue("{out.dir}/semibin/diamond-maps"))
  
  calls <- str_glue("{assembly.by.refseq} diamond blastx -q {ls.contigs} -d  databases/diamond/progenomes.proteins.dmnd \\
             -o {out.dir}/semibin/diamond-maps/{out.samples}.dmnd.tsv")
  
  if( slurm ) {
    Slurm.diamond <- Slurm_sapply(calls, function(x){print(x[[1]]); system(x[[1]])},
                                  plan = "none", 
                                  overwrite = TRUE, 
                                  mc.cores = threads,
                                  njobs = length(calls),
                                  sbatch_opt = list(partition = partition, account = account, time = "01:00:00"),
                                  preamble = c("module load singularity", "module load R/4.0.3"),
                                  job_name = paste0("Slurm-", projectname, "-diamond"), 
                                  tmp_path = projectdir)
    
    sbatch(Slurm.diamond)
    
    wait_slurm(Slurm.diamond , freq = 60, timeout = -1)
    
    Slurm_clean(Slurm.diamond)
  }else{
    future.apply::future_lapply(calls, function(x){print(x[[1]]); system(x[[1]])})
  }

  # 11 - sort into contigs into bins based on diamond assignments --------------
  
  printGlue("Sorting contigs into bins", log = TRUE)
  
  ls.dmnd <- list.files(str_glue("{out.dir}/semibin/diamond-maps"), pattern = "dmnd.tsv", full.names = TRUE)
  
  if(file.exists(str_glue("{out.dir}/semibin/bins/.diamond.completed.RData"))){
    load(str_glue("{out.dir}/semibin/bins/.diamond.completed.RData"))
    ls.dmnd <- ls.dmnd[!ls.dmnd%in%completed.diamond]
  }else{
    completed.diamond = c()
  }
  
  ncbi_taxonomy_progenomes <- data.table::fread("databases/diamond/proGenomes2.1_specI_lineageNCBI.tab", 
                                                header = FALSE)
  ncbi_taxonomy_progenomes[ncbi_taxonomy_progenomes == ""] <- "99999999_Unclassified"
  
  
  if( slurm ) {
    
    slurm_steps = seq_along(ls.dmnd)
    
    Slurm.make.bins <- Slurm_sapply(slurm_steps, function(x){
      print(str_glue("[{Sys.time()}] Loading diamond file {ls.dmnd[x]}"))
      in.dmnd <- data.table::fread(ls.dmnd[x],header = FALSE)
      
      unique.contigs <- in.dmnd$V1 %>% unique()
      print(str_glue("[{Sys.time()}] Getting best matches for contigs from diamond hits"))
      future::plan("multisession")
      unique.tab <- future.apply::future_lapply(unique.contigs, function(x){
        return(in.dmnd %>% filter(V1 == x ) %>% arrange(-V3, -V7) %>% head(n=1))
      }) 
      unique.tab<- do.call(rbind, unique.tab)
      unique.tab$genome = sapply(unique.tab$V2, function(x) {
        paste(str_split(x, "\\.")[[1]][1:2], collapse = ".")})
      
      print(str_glue("[{Sys.time()}] Matching hits, genomes, and NCBI Tax IDs"))
      m = match(unique.tab$genome, ncbi_taxonomy_progenomes$V1)
      unique.tab$tax_ID = ncbi_taxonomy_progenomes$V6[m] 
      
      print(str_glue("[{Sys.time()}] Getting bins"))
      get.bins <- unique.tab$tax_ID %>% unique() %>%
        str_replace_all(., " |\\.|:|\\[|\\]|\\(|\\)|\\'|\\-", "_") 
      
      print(str_glue("[{Sys.time()}] Making bins in {out.dir}/semibin/bins"))
      future::plan("multisession")
      print(str_glue("Set up a multisession plan ...."))
      make.dir <- lapply(get.bins, function(g) {
        calls = str_glue("mkdir -p {out.dir}/semibin/bins/{g}") 
        system(calls)
      })
      
      file.to.load = str_glue("{out.dir}/prinseq/{basename(ls.dmnd[x]) %>% str_replace_all(., '.dmnd.tsv', '')}.prinseq.headerrev.contigs.fasta")
      print(str_glue("[{Sys.time()}] Loading the contig file {file.to.load}"))
      in.contig.fasta <- Biostrings::readDNAStringSet(file.to.load)
      
      print(str_glue("[{Sys.time()}] Putting contigs into bins (this may take a while)"))
      put.bins <- future.apply::future_lapply(seq(1, nrow(unique.tab), 1), function(r){
        current.contig = unique.tab$V1[r]
        current.bin = unique.tab$tax_ID[r]%>%
          str_replace_all(., " |\\.|:|\\[|\\]|\\(|\\)|\\'|\\-", "_") 
        picked.contig = in.contig.fasta[current.contig]
        Biostrings::writeXStringSet(picked.contig, filepath = str_glue("{out.dir}/semibin/bins/{current.bin}/{names(picked.contig)}.fasta"))
      }, 
      future.packages = c("Biostrings", "tidyverse"),
      future.globals = c("in.contig.fasta", "unique.tab", "out.dir"))
      
      print(str_glue("[{Sys.time()}] Completed binning for {file.to.load}"))
      
      completed.diamond = c(completed.diamond, ls.dmnd[x])
      save(completed.diamond, 
           file = str_glue("{out.dir}/semibin/bins/.diamond.completed.RData"))
    },
    export = c("out.dir", "ls.dmnd", "ncbi_taxonomy_progenomes", "threads", "completed.diamond"),
    plan = "none", 
    overwrite = TRUE, 
    mc.cores = 24,
    njobs = length(slurm_steps),
    sbatch_opt = list(partition = partition, account = account, time = "04:00:00"),
    preamble = c("module load singularity", "module load R/4.0.3"),
    job_name = paste0("Slurm-", projectname, "-makeBins"), 
    tmp_path = projectdir)
    
    sbatch(Slurm.make.bins )
    
    wait_slurm(Slurm.make.bins , freq = 60, timeout = -1)
    
    Slurm_clean(Slurm.make.bins)
    
  }else{
    bin_steps = seq_along(ls.dmnd)
    make.bins <- sapply(bin_steps, function(x){
      print(str_glue("[{Sys.time()}] Loading diamond file {ls.dmnd[x]}"))
      in.dmnd <- data.table::fread(ls.dmnd[x],header = FALSE)
      
      unique.contigs <- in.dmnd$V1 %>% unique()
      print(str_glue("[{Sys.time()}] Getting best matches for contigs from diamond hits"))
      future::plan("multisession")
      unique.tab <- future.apply::future_lapply(unique.contigs, function(x){
        return(in.dmnd %>% filter(V1 == x ) %>% arrange(-V3, -V7) %>% head(n=1))
      }) 
      unique.tab<- do.call(rbind, unique.tab)
      unique.tab$genome = sapply(unique.tab$V2, function(x) {
        paste(str_split(x, "\\.")[[1]][1:2], collapse = ".")})
      
      print(str_glue("[{Sys.time()}] Matching hits, genomes, and NCBI Tax IDs"))
      m = match(unique.tab$genome, ncbi_taxonomy_progenomes$V1)
      unique.tab$tax_ID = ncbi_taxonomy_progenomes$V6[m] 
      
      print(str_glue("[{Sys.time()}] Getting bins"))
      get.bins <- unique.tab$tax_ID %>% unique() %>%
        str_replace_all(., " |\\.|:|\\[|\\]|\\(|\\)|\\'|\\-", "_") 
      
      print(str_glue("[{Sys.time()}] Making bins in {out.dir}/semibin/bins"))
      future::plan("multisession")
      print(str_glue("Set up a multisession plan ...."))
      make.dir <- lapply(get.bins, function(g) {
        calls = str_glue("mkdir -p {out.dir}/semibin/bins/{g}") 
        system(calls)
      })
      
      file.to.load = str_glue("{out.dir}/prinseq/{basename(ls.dmnd[x]) %>% str_replace_all(., '.dmnd.tsv', '')}.prinseq.headerrev.contigs.fasta")
      print(str_glue("[{Sys.time()}] Loading the contig file {file.to.load}"))
      in.contig.fasta <- Biostrings::readDNAStringSet(file.to.load)
      
      print(str_glue("[{Sys.time()}] Putting contigs into bins (this may take a while)"))
      put.bins <- future.apply::future_lapply(seq(1, nrow(unique.tab), 1), function(r){
        current.contig = unique.tab$V1[r]
        current.bin = unique.tab$tax_ID[r]%>%
          str_replace_all(., " |\\.|:|\\[|\\]|\\(|\\)|\\'|\\-", "_") 
        picked.contig = in.contig.fasta[current.contig]
        Biostrings::writeXStringSet(picked.contig, filepath = str_glue("{out.dir}/semibin/bins/{current.bin}/{names(picked.contig)}.fasta"))
      }, 
      future.packages = c("Biostrings", "tidyverse"),
      future.globals = c("in.contig.fasta", "unique.tab", "out.dir"))
      
      print(str_glue("[{Sys.time()}] Completed binning for {file.to.load}"))
      
      completed.diamond = c(completed.diamond, ls.dmnd[x])
      save(completed.diamond, 
           file = str_glue("{out.dir}/semibin/bins/.diamond.completed.RData"))
    })
    
  }
  

# 12 - Super assemble contigs in each bin ----------------------------------

  ls.dirs<- system(str_glue("ls {out.dir}/semibin/bins"), intern = TRUE) # the R version list.dirs is SLOW!
  calls = sapply(ls.dirs, function(x){
    #>>> collect the contigs files in each
    call1 = str_glue("[ -e {out.dir}/semibin/bins/{x}/combined.fasta ] && rm {out.dir}/semibin/bins/{x}/combined.fasta")
    ca112 = str_glue("ls {out.dir}/semibin/bins/{x}/*.fasta | xargs cat >> {out.dir}/semibin/bins/{x}/combined.fasta")
    #>>> flye
    call3 = str_glue("singularity run -B /mnt  -B /run/shm:/run/shm containers/assembly-by-refseq.sif flye \\
    --subassemblies {out.dir}/semibin/bins/{x}/combined.fasta -o {out.dir}/semibin/bins/{x}/flye \\
                     --meta --threads {threads} --scaffold") %>%
      str_squish()
    #>>> raven 
    call4 = str_glue("mkdir -p {out.dir}/semibin/bins/{x}/raven && \\
                     singularity run -B /mnt  -B /run/shm:/run/shm containers/assembly-by-refseq.sif  /NGStools/raven/build/bin/raven \\
                      -t {threads} --weaken \\
                     {out.dir}/semibin/bins/{x}/combined.fasta > {out.dir}/semibin/bins/{x}/raven/raven.fasta ") %>%
      str_squish()
    calls = paste(c(call1, call2, call3, call4), collapse = " && ") %>%
      str_squish()
    return(calls)
  })
  
  
}