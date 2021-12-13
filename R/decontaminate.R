#' RMetaflow
#' Module: decontaminate
#' Initiated March 2021
#' Patrick C. Seed
#' Stanley Manne Children's Research Institute
#' Ann and Robert H. Lurie Children's Hospital
#' Northwestern University
#' MIT license
#' 

if ( !commands[STEP][[1]]$completed ){
  
  library(tidyverse)
  library(slurmR)
  
  printGlue("#---- Running {names(commands[STEP][1])}", log = TRUE)
  printGlue("#---- Gonna wash that junk right out of my {current.step[[1]]$in_ext}", log = TRUE)

  # set up step variables ---------------------------------------------------
  
  in.dir = dirNormalize(current.step[[1]]$from)
  out.dir = dirNormalize(names(current.step[1]))
  dirIf(out.dir)
  
# make combined minimap2 reference ----------------------------------------

  decon.l = length(variables.decontaminate)
  decon.ext = variables.decontaminate[1]
  decon.list = list()
  for(i in seq(from = 2, to = decon.l, by = 1)){
    decon.list[[i-1]] = list.files(variables.decontaminate[i], pattern = decon.ext, full.names = TRUE)  
  }
  decon.list = decon.list %>% unlist()
  
  if( ! file.exists(str_glue("{out.dir}/ref.mmi"))){
    printGlue("#-------- Preparing the indexed reference for decontamination", log = TRUE)
    printGlue("#-------- Good time for a walk or a coffee", log = TRUE)
    calls = str_glue("cat {paste(decon.list, collapse = ' ')} > {out.dir}/ref.fna  && \\
             {assembly.by.refseq} /NGStools/CONCOCT/scripts/cut_up_fasta.py {out.dir}/ref.fna -c 10000 -o 0 \\
                     --merge_last -b {out.dir}/ref.bed > {out.dir}/ref.10K.fna && \\
                     {qc} minimap2 -d {out.dir}/ref.10K.mmi {out.dir}/ref.10K.fna && \\
                     rm {out.dir}/ref.fna") %>% 
      str_squish()
    if ( slurm ){
      decon.slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},  
                                  plan = "none", 
                                  overwrite = TRUE, 
                                  njobs = length(calls),
                                  mc.cores = threads,
                                  sbatch_opt = list(partition = partition, account = account, time = "01:00:00", mem="96G"),
                                  preamble = c("module load singularity", "module load R/4.0.3"),
                                  job_name = paste0("Slurm-", projectname, "-decon-ref"), 
                                  tmp_path = projectdir)
      
      sbatch(decon.slurm)
      
      wait_slurm(decon.slurm, freq = 60, timeout = -1)
      
      Slurm_clean(decon.slurm)
      
    }else{
      lapply(calls, function(x) {print(x[[1]]); system(x[[1]])})
  }
  

# map sample seqs against decon references --------------------------------

  if( file.exists(str_glue("{out.dir}/ref.mmi"))){
    printGlue("#-------- I'm indexed and loaded. Ready to get cleaning!", log = TRUE)
    ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
    out.samples = sapply(list.files(out.dir, pattern = ".decon"), 
                         function(x) str_split(x, ".decon")[[1]][1]) %>% as.vector()
    
    if( length(out.samples) > 0 ){
      ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
      remaining.samples = samples[!samples%in%out.samples]
    }else{
      remaining.samples = samples
    }
    
    if( length(ls) > 0 ){
      printGlue("#---- Preparing samples for decontamination", log = TRUE)  
      split_prefix = paste0(out.dir, "/", "tempsplit.", remaining.samples)
      fseq = ls[ seq(from = 1, to = length(ls), by = 2)]
      rseq = ls[ seq(from = 2, to = length(ls), by = 2)]
      bam.out = paste0(out.dir, "/", remaining.samples, ".sorted.bam")
      fout = paste0(out.dir, "/", remaining.samples, ".decon_R1.fastq.gz")
      rout = paste0(out.dir, "/", remaining.samples, ".decon_R2.fastq.gz")
      
      
      job.list = list(fseq = fseq, rseq = rseq,
                      bam.out = bam.out, split_prefix = split_prefix,
                      fout = fout, rout = rout)
      
      calls <- purrr::pmap(job.list, function(fseq, rseq, bam.out, fout, rout, split_prefix){
        str_glue("[ ! -e {fout} ] && {qc} minimap2 -a -x sr -K 4G -2 -t {threads} \\
        --split-prefix {split_prefix} {out.dir}/ref.10K.fna {fseq} {rseq} | \\
        {qc} samtools view -@ {threads} -b -  | \\
        {qc} samtools view -@ {threads} -u -h -f 12 -F 256 - | \\
        {qc} samtools sort -t {threads} -o {bam.out} - && \\
        {qc} reformat.sh in={bam.out} out=stdout.fq primaryonly | \\
        {qc} reformat.sh in=stdin.fq out1={fout} \\
        out2={rout} interleaved addcolon && \\
        rm {bam.out}") %>% str_squish()}) %>% 
        unlist()
      
      printGlue("#-------- Mapping reads against decontamination reference", log = TRUE)
      
      if( slurm ){
        printGlue("#-------- Slurm and burn, Baby!")
        decon.slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},  
        plan = "none", 
        overwrite = TRUE, 
        mc.cores = 12,
        njobs = length(calls),
        sbatch_opt = list(partition = partition, account = account, time = "06:00:00"),
        preamble = c("module load singularity", "module load R/4.0.3"),
        job_name = paste0("Slurm-", projectname, "-decon"), 
        tmp_path = projectdir)
        
        sbatch(decon.slurm)
        
        wait_slurm(decon.slurm, freq = 60, timeout = -1)
        
        Slurm_clean(decon.slurm)
      }
    }else{
      printGlue("#-------- I asked for a Slurm and all I got was this", log = TRUE)
      printGlue("#-------- was this lousy serial processing system with {threads}", log = TRUE)
      
      future.apply::future_lapply(calls, system)  
    }
  }
  
  printGlue("#-------- Those sequences are now clean as a whistle!", log = TRUE)
  ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
  out.samples = sapply(list.files(out.dir, pattern = ".decon"), 
                       function(x) str_split(x, ".decon")[[1]][1]) %>% as.vector()
  
  if( length(out.samples) > 0 ){
    ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
    if( length(ls) == 0 ){
      commands[[STEP]]$completed = TRUE
      updateYaml()
    }
  }
  }
}


source.complete = TRUE