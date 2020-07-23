# the first time you run... prompt for email address

# include ALL persistent config settings
rip config
  # TODO mandatory
  --email '<str>'
  --cache '/path/to/what/used/to/be/database'  # TODO default to ~/rip/

  # TODO optional
  --batchsize '<int>'
  --bowtie2-scaling  '<int>'  # mem-rescaling
  --bowtie2-treads   '<int>'  # needs to be documented that this effects --mem-scaling
  --mismatch_probability '<float>'

  --etc...

# outputs to > ~/.rip/config.yaml

#########################

# default flags, available for ALL commands
  --cores '<int>'    # default to ALL cores
  --memory '<float>' # default to ALL memory
  --snakemake "--any-smk-flag --etc -blah"  # none of these options are documented, but we can link to the smk CLI in our docs

#########################

rip database --help

# Up to the bowtie2 indices for the individual genomes
rip database
  --mode fetch

  --output '/path/to/what/used/to/be/query_name'
  --refseq-rep              # optional, default False
  --query '<str>'           # optional, default ""
  --query-file '<filename>' # optional, default ""
  --rank '<str>'            # optional, default "species"
  --sequences '<filename>'  # optional, default ""
  --accessions '<filename>' # optional, default ""

# if you run fetch, it must produce a settings file

# Up to the bowtie2 filtering aln idx
rip database
  --mode index
  --output '/path/to/db/name'

#########################

# Do everything db related
rip database
  --mode build              # optional, default "build"

  --output '/path/to/db'
  --rank '<str>'            # optional, default "species"
  --refseq-rep              # optional, default False
  --query '<str>'           # optional, default ""
  --sequences '<filename>'  # optional, default ""
  --accessions '<filename>' # optional, default ""

# must have at least one of (--refseq-rep / --query / -sequences / --accessions)

# if you run build, it must produce a settings file
# if you run build... you must not have first run fetch
# if you run either fetch or build, you cannot run a second time with the same ID

################################################################################

# Only trim adapters
rip sample
  --output '/path/to/sample'

  --fastq '<str>'           # optional, default ""
  --fastq-r1 '<str>'        # optional, default ""
  --fastq-r2 '<str>'        # optional, default ""
  --sra '<str>'             # optional, default ""
  --collapse '<str>'        # optional, default False
  --trim-adapters           # optional, default True
  --adapter-remove-flags    # TODO Antony to add flags for specifying adapter sequences

# must provide either (--fastq XOR (--fastq-r1 AND --fastq-r2) XOR --sra)
# --collapse can only be used with (--fastq-r1 AND --fastq-r2)
# --sra should know if PE or SE (will have to query the metadata)

# TODO do we want to run fastqc here?
# TODO make sure that this produces ALL the sample specific targets

################################################################################

# Only do everything up to the filtering alignment 
rip analyse
  --mode filter  # TODO find better name for this
  --database '/path/to/database'
  --sample '/path/to/sample'
  --output '/path/to/results'

# --database must existing and be complete!
# --sample "preprocess" must be complete (otherwise throw error)

# Only do everything up to the metagenomic alignments
rip analyse
  --mode align
  --database '/path/to/database'
  --sample '/path/to/sample'
  --output '/path/to/results'

# Only do everything up to the metagenomic alignments
rip analyse
  --mode likelihoods
  --database '/path/to/database'
  --sample '/path/to/sample'
  --output '/path/to/results'

#########################

# Only do everything up to calculating the assignment probabilities
rip analyse
  --mode probabilities
  --database '/path/to/database'
  --sample '/path/to/sample'
  --genera '<str1 str2 str3 ...>'  # optional, default ""
  --output '/path/to/results'

# Only do everything up to calculating the posterior abundances
rip analyse
  --mode abundances
  --database '/path/to/database'
  --sample '/path/to/sample'
  --genera '<str1 str2 str3 ...>'  # optional, default ""
  --output '/path/to/results'

# Only do everything up to mapdamage analysis
rip analyse
  --mode mapdamage
  --database '/path/to/database'
  --sample '/path/to/sample'
  --genera '<str1 str2 str3 ...>'  # optional, default ""
  --output '/path/to/results'

