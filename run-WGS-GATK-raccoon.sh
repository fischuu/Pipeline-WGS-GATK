# Perform the WGS-GATK pipeline on Puhti
# Dataset: GBS data from Raccoon
# Daniel Fischer (daniel.fischer@luke.fi)

# To create samples_raccoon.tsv run:
# ls /scratch/project_2001746/RaccoonGATK/FASTQ/CONCATENATED/ | cut -f1 -d 'R' | sort | uniq | sed 's/.$//' | sed '1 i\sample' > /scratch/project_2001746/RaccoonGATK/samples_raccoon.tsv
# ls /scratch/project_2001746/RaccoonGATK/FASTQ/CONCATENATED/ | cut -f1 -d 'R' | sort | uniq | sed 's/.$//' | sed "s/$/\t1/" | sed "s/$/\tILLUMINA/" | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,"../FASTQ/CONCATENATED/"$1"_R1_001.merged.fastq.gz","../FASTQ/CONCATENATED/"$1"_R2_001.merged.fastq.gz"}' | sed '1 i\sample\tunit\tplatform\tfq1\tfq2' > /scratch/project_2001746/RaccoonGATK/pipeline/units_raccoon.tsv
#$1 is target file


module load bioconda/3
# Uncomment this, once the environemnt is created
source activate /projappl/project_2001746/conda_envs/Genotype
#source activate /projappl/project_2001289/FAANGlncRNA

# snakemake -s Snakefile --configfile /scratch/project_2001746/RaccoonGATK/pipeline/WGS-GATK-pipeline_config_raccoon.yaml --forceall --rulegraph | dot -Tpdf > WGS-GATK-rulegraph.pdf

snakemake -s Snakefile \
          --configfile /scratch/project_2001746/RaccoonGATK/pipeline/WGS-GATK-pipeline_config_raccoon.yaml \
          -j 200 \
          --use-conda \
          --latency-wait 60 \
          --cluster-config WGS-GATK_puhti.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --job-name={cluster.job-name} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory}" $1 
          