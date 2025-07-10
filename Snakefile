configfile: 'config/config.yaml'

rule all:
    input:
        #FASTQC
        expand('out/fastqc/raw/{sample}/{sample}_L001_R1_001_fastqc.html', sample=config['samples']),
        expand('out/fastqc/raw/{sample}/{sample}_L001_R2_001_fastqc.html', sample=config['samples']),
        expand('out/fastqc/raw/{sample}/{sample}_L001_R1_001_fastqc.zip', sample=config['samples']),
        expand('out/fastqc/raw/{sample}/{sample}_L001_R2_001_fastqc.zip', sample=config['samples']),

        #TRIM_GALORE
        expand("out/trim/{sample}_val_1.fq.gz", sample=config["samples"]),
        expand("out/trim/{sample}_val_2.fq.gz",sample=config["samples"]),

        #FASTQC ON TRIMMED READS
        expand("out/fastqc/trimmed/{sample}/{sample}_val_1_fastqc.html", sample=config["samples"]),
        expand("out/fastqc/trimmed/{sample}/{sample}_val_2_fastqc.html", sample=config["samples"]),
        expand("out/fastqc/trimmed/{sample}/{sample}_val_1_fastqc.zip", sample=config["samples"]),
        expand("out/fastqc/trimmed/{sample}/{sample}_val_2_fastqc.zip", sample=config["samples"]),

        #SPADES
        expand("out/scaffolds/{sample}_scaffolds.fasta", sample=config["samples"]),

        #MLST
        expand("out/mlst/{sample}_mlst.tsv", sample=config["samples"]),

        #PROKKA
        expand("out/prokka/{sample}/{sample}.faa", sample=config["samples"]),
        expand("out/prokka/{sample}/{sample}.fna", sample=config["samples"]),
        expand("out/prokka/{sample}/{sample}.gbk", sample=config["samples"]),
        expand("out/prokka/{sample}/{sample}.gff", sample=config["samples"]),

        #NCBI AMRFINDER PLUS
        expand("out/amrfinder/{sample}.tsv", sample=config["samples"]),

        #PLASMIDFINDER
        expand("out/plasmids/{sample}/results_tab.tsv", sample=config["samples"]),


rule plasmid_finder:
	input:
		fasta = "out/prokka/{sample}/{sample}.fna",
	output:
		out_dir = directory("out/plasmids/{sample}"),
		plasmids = "out/plasmids/{sample}/results_tab.tsv",
	shell:
		"""
		mkdir -p out/plasmids/{wildcards.sample}

	    python3 ~/plasmidfinder/plasmidfinder.py -i {input.fasta} -o {output.out_dir} -x -p ~/plasmidfinder/plasmidfinder_db
		"""
rule amrfinder:
	input:
		fna = "out/prokka/{sample}/{sample}.fna",
		faa = "out/prokka/{sample}/{sample}.faa",
		gff = "out/prokka/{sample}/{sample}.gff",
		mlst = "out/mlst/{sample}_mlst.tsv",
	output:
		amr = "out/amrfinder/{sample}.tsv",
	threads: config["threads"]
	shell:
		"""
		SAMPLE_STRAIN=$(gawk '{{print $2}}' {input.mlst})


        case $SAMPLE_STRAIN in
            ecoli*              )   ORGANISM="Escherichia";;
            klebsiella          )   ORGANISM="Klebsiella_pneumoniae";;
            *                   )   ORGANISM="";;
	    unknown             )   ORGANISM="";;
	    smaltophilia        )   ORGANISM="";;
        esac

		mkdir -p out/amrfinder

		amrfinder -a prokka --gff {input.gff} -n {input.fna} -p {input.faa} -O $ORGANISM --plus --threads 1 > {output.amr}


		"""



rule prokka:
	input:
		scaffold = "out/scaffolds/{sample}_scaffolds.fasta",
		mlst = "out/mlst/{sample}_mlst.tsv",
	output:
		faa = "out/prokka/{sample}/{sample}.faa",
		fna = "out/prokka/{sample}/{sample}.fna",
		gbk = "out/prokka/{sample}/{sample}.gbk",
		gff = "out/prokka/{sample}/{sample}.gff",
		out_dir = directory("out/prokka/{sample}"),
	params:
		kingdom = "Bacteria",
	threads: config["threads"]
	shell:
		"""
		#SAMPLE_STRAIN=$(echo {input.mlst} | gawk '{{print $2}}')
		SAMPLE_STRAIN=$(gawk '{{print $2}}' {input.mlst})


        case $SAMPLE_STRAIN in
	    smaltophilia        )   GENUS="Stenotrophomonas"; SPECIES="Maltophilia";; 
            ecoli*              )   GENUS="Escherichia"; SPECIES="Coli";;
            klebsiella          )   GENUS="Klebsiella"; SPECIES="Pneumoniae";;
            *                   )   GENUS=""; SPECIES="";;
	    unknown             )   GENUS=""; SPECIES="";;
        esac

		mkdir -p out/prokka/{wildcards.sample}
		prokka --kingdom {params.kingdom} --genus $GENUS --species $SPECIES --usegenus --compliant --centre UoN --cpus {threads} --prefix {wildcards.sample} --outdir {output.out_dir} --force {input.scaffold}

		"""


rule mlst:
	input:
		scaffold = "out/scaffolds/{sample}_scaffolds.fasta",
	output:
		file = "out/mlst/{sample}_mlst.tsv",
	threads: config["threads"]
	shell:
		"mlst --threads {threads} --quiet --label {wildcards.sample}_mlst {input.scaffold} > {output.file}"

rule copy_scaffolds:
    input:
        scaffolds = "out/spades/{sample}/scaffolds.fasta",
    output:
        copied_scaffolds = "out/scaffolds/{sample}_scaffolds.fasta",
    shell:
        "cp {input.scaffolds} {output.copied_scaffolds}"


rule spades:
    input:
        read1 = "out/trim/{sample}_val_1.fq.gz",
        read2 = "out/trim/{sample}_val_2.fq.gz",
    output:
        out_dir = directory("out/spades/{sample}"),
        scaffolds = "out/spades/{sample}/scaffolds.fasta",
    threads: config["threads"]
    shell:
        "spades.py --isolate -1 {input.read1} -2 {input.read2} -o {output.out_dir} -t {threads} -m 10"


rule fastqc_trimmed:
	input:
		read1 = "out/trim/{sample}_val_1.fq.gz",
		read2 = "out/trim/{sample}_val_2.fq.gz",
	output:
		html1 = "out/fastqc/trimmed/{sample}/{sample}_val_1_fastqc.html",
		html2 = "out/fastqc/trimmed/{sample}/{sample}_val_2_fastqc.html",
		zip1 = "out/fastqc/trimmed/{sample}/{sample}_val_1_fastqc.zip",
		zip2 = "out/fastqc/trimmed/{sample}/{sample}_val_2_fastqc.zip",
	threads: config["threads"]
	shell:
		"mkdir -p out/fastqc/trimmed/{wildcards.sample} && "
		"fastqc -t {threads} -o out/fastqc/trimmed/{wildcards.sample} {input.read1} {input.read2}"


rule trim_galore:
	input:
		read1 = "data/{sample}_L001_R1_001.fastq.gz",
		read2 = "data/{sample}_L001_R2_001.fastq.gz",
	output:
		trimmed_r1 = "out/trim/{sample}_val_1.fq.gz",
		trimmed_r2 = "out/trim/{sample}_val_2.fq.gz",
	params:
		quality = config["trimming"]["quality"],
		length = config["trimming"]["length"],
	threads: config["threads"]
	shell:
		"mkdir -p out/trim && "
		"trim_galore -q {params.quality} --length {params.length} --trim-n "
		"--basename {wildcards.sample} -o out/trim --paired {input.read1} {input.read2} --cores {threads}"

rule fastqc_reads:
	input:
		read1 = "data/{sample}_L001_R1_001.fastq.gz",
		read2 = "data/{sample}_L001_R2_001.fastq.gz",
	output:
		html1 = "out/fastqc/raw/{sample}/{sample}_L001_R1_001_fastqc.html",
		html2 = "out/fastqc/raw/{sample}/{sample}_L001_R2_001_fastqc.html",
		zip1 = "out/fastqc/raw/{sample}/{sample}_L001_R1_001_fastqc.zip",
		zip2 = "out/fastqc/raw/{sample}/{sample}_L001_R2_001_fastqc.zip"
	threads: config["threads"]
	shell:
		"mkdir -p out/fastqc/raw/{wildcards.sample} && "
		"fastqc -t {threads} -o out/fastqc/raw/{wildcards.sample} {input.read1} {input.read2}"
		
