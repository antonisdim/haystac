# only fetch genomes
rip database_build \
--query_name <str> \
--entrez_query <str>  \
--entrez_email <str> \
--entrez_batchsize <int> \
--entrez_rank <str> \
--custom_seq_file <str> \
--custom_acc_file <str> \
--mem_mb <float> \
--mem_rescale_factor <float> \
--with_entrez_query \
--with_refseq_rep \
--build_db

# Up to the bowtie2 indices for the individual genomes
rip database_build \
--query_name <str> \
--entrez_query <str>  \
--entrez_email <str> \
--entrez_batchsize <int> \
--entrez_rank <str> \
--custom_seq_file <str> \
--custom_acc_file <str> \
--mem_mb <float> \
--mem_rescale_factor <float> \
--with_entrez_query \
--with_refseq_rep \
--index_db

# Up to the bowtie2 filtering aln idx 
rip database_build \
--query_name <str> \
--entrez_query <str>  \
--entrez_email <str> \
--entrez_batchsize <int> \
--entrez_rank <str> \
--custom_seq_file <str> \
--custom_acc_file <str> \
--mem_mb <float> \
--mem_rescale_factor <float> \
--with_entrez_query \
--with_refseq_rep \
--aln_filter_index

# Do everything db related
rip database_build \
--query_name <str> \
--entrez_query <str>  \
--entrez_email <str> \
--entrez_batchsize <int> \
--entrez_rank <str> \
--custom_seq_file <str> \
--custom_acc_file <str> \
--mem_mb <float> \
--mem_rescale_factor <float> \
--with_entrez_query \
--with_refseq_rep \
--complete_db_creation

# Only trim adapters
rip analyse_sample \
--sample_name <str> \
--fastq <str> \
--fastq_R1 <str> \
--fastq_R2 <str> \
--mismatch_probability <float> \
--bowtie2_treads <int> \
--sra <str> \
--input_mode <str> \
--trim_adapters \
--specific_genera <str1 str2 str3 ...> \
--data_preprocess

# Only do everything up to the filtering alignment 
rip analyse_sample \
--sample_name <str> \
--fastq <str> \
--fastq_R1 <str> \
--fastq_R2 <str> \
--mismatch_probability <float> \
--bowtie2_treads <int> \
--sra <str> \
--input_mode <str> \
--trim_adapters \
--specific_genera <str1 str2 str3 ...> \
--aln_filter

# Only do everything up to the metagenomic alignments
rip analyse_sample \
--sample_name <str> \
--fastq <str> \
--fastq_R1 <str> \
--fastq_R2 <str> \
--mismatch_probability <float> \
--bowtie2_treads <int> \
--sra <str> \
--input_mode <str> \
--trim_adapters \
--specific_genera <str1 str2 str3 ...> \
--aln_meta

# Only do everything up to calculating the assignment probabilities
rip analyse_sample \
--sample_name <str> \
--fastq <str> \
--fastq_R1 <str> \
--fastq_R2 <str> \
--mismatch_probability <float> \
--bowtie2_treads <int> \
--sra <str> \
--input_mode <str> \
--trim_adapters \
--specific_genera <str1 str2 str3 ...> \
--probabilities

# Only do everything up to calculating the posterior abundances
rip analyse_sample \
--sample_name <str> \
--fastq <str> \
--fastq_R1 <str> \
--fastq_R2 <str> \
--mismatch_probability <float> \
--bowtie2_treads <int> \
--sra <str> \
--input_mode <str> \
--trim_adapters \
--specific_genera <str1 str2 str3 ...> \
--abundances

# Only do everything up to mapdamage analysis
rip analyse_sample \
--sample_name <str> \
--fastq <str> \
--fastq_R1 <str> \
--fastq_R2 <str> \
--mismatch_probability <float> \
--bowtie2_treads <int> \
--sra <str> \
--input_mode <str> \
--trim_adapters \
--specific_genera <str1 str2 str3 ...> \
--mapdamage
