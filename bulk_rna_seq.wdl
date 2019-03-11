workflow bulk_rna_seq {
    String experiment_name
    String experiment_type
    String organism
    String input_data_type
    File reference_transcriptome

    String read_pairs_file
    String samples_described_file
    String samples_compared_file

    # run bcl2fastq
    call bcl2fastq_and_define_read_pairs {
        input:
            input_data_type = input_data_type,
            experiment_name = experiment_name,
            experiment_type = experiment_type,
            read_pairs_file = read_pairs_file,
            bcl2fastq_and_define_read_pairs_RAM = bcl2fastq_and_define_read_pairs_RAM,
            bcl2fastq_and_define_read_pairs_disk_space = bcl2fastq_and_define_read_pairs_disk_space
    }

    # build reference transcriptome index for kallisto
    call build_kallisto_index {
        input:
            reference_transcriptome = reference_transcriptome,
            experiment_name = experiment_name,
            experiment_type = experiment_type,
            build_kallisto_index_RAM = build_kallisto_index_RAM,
            build_kallisto_index_disk_space = build_kallisto_index_disk_space
    }

    # get pairs of reads
    Map[String, String] read_pairs = read_map(bcl2fastq_and_define_read_pairs.read_pairs_tsv)

    # for each read pair ...
    scatter (pair in read_pairs) {
        # define files matching read pair
        File zipped_R1 = "${pair.left}"
        File zipped_R2 = "${pair.right}"

        # filter read name to contain only sample information
        String sample_name = sub(sub(basename(zipped_R1), ".fastq.gz", ""), "_R[0-9]", "")

        # unzip each file (if zipped) -- this is necessary for fastqc
        call unzip_file {
            input:
                zipped_file_1 = zipped_R1,
                zipped_file_2 = zipped_R2,
                unzip_file_RAM = unzip_file_RAM,
                unzip_file_disk_space = unzip_file_disk_space
        }

        # perform fastqc
        call perform_fastqc {
            input:
                experiment_name = experiment_name,
                experiment_type = experiment_type,
                file_1 = unzip_file.file_1,
                file_2 = unzip_file.file_2,
                sample_name = sample_name,
                perform_fastqc_RAM = perform_fastqc_RAM,
                perform_fastqc_disk_space = perform_fastqc_disk_space
        }

        # get transcript quantification
        call perform_kallisto_quantification {
            input:
                experiment_name = experiment_name,
                experiment_type = experiment_type,
                transcriptome_index = build_kallisto_index.reference_index,
                file_1 = unzip_file.file_1,
                file_2 = unzip_file.file_2,
                sample_name = sample_name,
                perform_kallisto_quantification_RAM = perform_kallisto_quantification_RAM,
                perform_kallisto_quantification_disk_space = perform_kallisto_quantification_disk_space
        }
    }

    # transform transcript counts to gene counts (to be used in edgeR) and also run differential expression
    # via Sleuth and edgeR
    call counts_and_differential_expression_output {
        input:
            experiment_name = experiment_name,
            experiment_type = experiment_type,
            transcript_counts_tar = perform_kallisto_quantification.transcript_counts_tar,
            organism = organism,
            samples_described_file = samples_described_file,
            samples_compared_file = samples_compared_file,
            counts_and_differential_expression_output_RAM = counts_and_differential_expression_output_RAM,
            counts_and_differential_expression_output_disk_space = counts_and_differential_expression_output_disk_space
    }

    # perform multiqc
    call perform_multiqc {
        input:
            experiment_name = experiment_name,
            experiment_type = experiment_type,
            transcript_counts_tar = perform_kallisto_quantification.transcript_counts_tar,
            perform_multiqc_RAM = perform_multiqc_RAM,
            perform_multiqc_disk_space = perform_multiqc_disk_space
    }
}

##################
## define tasks ##
##################

task bcl2fastq_and_define_read_pairs {
    String input_data_type
    String experiment_name
    String experiment_type
    String read_pairs_file

    command <<<
        set -e

        if [[ ${input_data_type} == 'bcl' ]]; then
            # copy bcl data from cloud bucket
            gsutil -m cp -r gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/bcl ./

            # run bcl2fastq
            mkdir -p ./fastq
            bcl2fastq --runfolder-dir ./bcl \
                      --output-dir ./fastq \
                      --mask-short-adapter-reads 10 \
                      --minimum-trimmed-read-length 10 \
                      --no-lane-splitting

            # copy Reads data back into google bucket
            gsutil -m cp -r ./fastq/* gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/fastq/

            # create read_pairs file and copy to google bucket
            fastq_files=$(find ./fastq -name "*.fastq.gz" ! -name "Undetermined*" -type f -exec basename {} \; | sort) # get list of non-'Undetermined' fastq files and ensure they're paired together via sorting
            echo $fastq_files > read_pairs.tsv # write fastq names to file
            sed -i 's/ /\n/g' read_pairs.tsv # replace every whitespace with a newline character
            sed -i -e "s/^/gs:\/\/genomics_xavier_bucket\/${experiment_type}\/${experiment_name}\/fastq\//g" read_pairs.tsv # add path to google bucket at beginning of each file name
            sed -i -e 'N;s/\n/\t/' read_pairs.tsv # for every line corresponding to 'Read 1', replace newline with tab
            gsutil -m cp read_pairs.tsv ${read_pairs_file} # write read_pairs.tsv file we generated to the specified read_pairs_file location
        else
            # localize read_pairs file
            gsutil -m cp ${read_pairs_file} ./read_pairs.tsv
        fi
    >>>

    output {
        File read_pairs_tsv = "read_pairs.tsv"
    }

    runtime {
        docker: "gcr.io/genomics-xavier/bulk_rna_seq"
        zones: "us-east1-b us-east1-c us-east1-d"
		cpu: 4
  		memory: "50GB"
  		preemptible: 2
  		disks: "local-disk 160 HDD"
    }
}

task build_kallisto_index {
    File reference_transcriptome
    String experiment_name
    String experiment_type

    Int build_kallisto_index_RAM
    Int build_kallisto_index_disk_space

    String reference_transcriptome_basename = sub(basename(reference_transcriptome), ".fa.gz", "")

    command <<<
        set -e

        # list reference indices in google bucket
        indices=$(gsutil ls gs://rnaseq_reference_data/kallisto_indices/)

        # check if reference index already exists in google cloud bucket
        index_exists=0
        for index in $indices; do
            if [[ $index == *"${reference_transcriptome_basename}.idx" ]]; then
                index_exists=1
                break
            fi
        done

        # if reference index already exists, then copy from google bucket -- otherwise, generate the reference index and copy to google bucket
        if [! -z "$indices" ] && [[ $index_exists == 1 ]]; then
            echo "reference index already exists"
            # localize reference index
            gsutil -m cp gs://rnaseq_reference_data/kallisto_indices/${reference_transcriptome_basename}.idx ./
        else
            echo "generating reference index"
            # generate reference index
            kallisto index --index=${reference_transcriptome_basename}.idx \
                           ${reference_transcriptome}

            # copy reference index to google bucket
            gsutil -m cp ${reference_transcriptome_basename}.idx gs://rnaseq_reference_data/kallisto_indices/
        fi
    >>>

    output {
        File reference_index = "${reference_transcriptome_basename}.idx"
    }

    runtime {
        docker: "gcr.io/genomics-xavier/bulk_rna_seq"
        zones: "us-east1-b us-east1-c us-east1-d"
		cpu: 4
  		memory: "16GB"
  		preemptible: 2
  		disks: "local-disk 80 HDD"
    }
}

task unzip_file {
    File zipped_file_1
    File zipped_file_2

    Int unzip_file_RAM
    Int unzip_file_disk_space

    String unzipped_basename_1 = sub(basename(zipped_file_1), ".gz", "")
    String unzipped_basename_2 = sub(basename(zipped_file_2), ".gz", "")

    command <<<
        set -e

        if [[ ${zipped_file_1} == *.gz ]]; then
            # unzip files if they are zipped
            unpigz -c ${zipped_file_1} > ${unzipped_basename_1}
            unpigz -c ${zipped_file_2} > ${unzipped_basename_2}
        fi
    >>>

    output {
        File file_1 = unzipped_basename_1
        File file_2 = unzipped_basename_2
    }

    runtime {
        docker: "gcr.io/genomics-xavier/bulk_rna_seq"
        zones: "us-east1-b us-east1-c us-east1-d"
		cpu: 4
  		memory: "8GB"
  		preemptible: 2
  		disks: "local-disk 80 HDD"
    }
}

task perform_fastqc {
    String experiment_name
    String experiment_type
    File file_1
    File file_2
    String sample_name

    Int perform_fastqc_RAM
    Int perform_fastqc_disk_space

    command <<<
        set -e

        # perform fastqc
        mkdir -p ${sample_name}
        fastqc -f fastq -o ${sample_name} ${file_1} ${file_2}

        # tar fastqc output directory
        tar zcfv fastqc_${sample_name}.tar.gz ${sample_name}

        # copy fastqc output to google bucket
        gsutil -q -m cp -r ${sample_name} gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/fastqc/
    >>>

    output {
        File fastqc_output = "fastqc_${sample_name}.tar.gz"
    }

    runtime {
        docker: "gcr.io/genomics-xavier/bulk_rna_seq"
        zones: "us-east1-b us-east1-c us-east1-d"
        cpu: 4
        memory: "4GB"
        preemptible: 2
        disks: "local-disk 80 HDD"
    }
}

task perform_kallisto_quantification {
    String experiment_name
    String experiment_type
    File transcriptome_index
    File file_1
    File file_2
    String sample_name

    Int perform_kallisto_quantification_RAM
    Int perform_kallisto_quantification_disk_space

    command <<<
        set -e

        # perform quantification
        mkdir -p transcript_counts_${sample_name}
        kallisto quant --index=${transcriptome_index} \
                       --output-dir=transcript_counts_${sample_name} \
                       --bootstrap-samples=100 \
                       --threads=11 \
                       ${file_1} ${file_2} > transcript_counts_${sample_name}/stdout_qc.txt 2>&1

        # tar output directory
        tar zcfv transcript_counts_${sample_name}.tar.gz transcript_counts_${sample_name}

        # copy output to google bucket
        gsutil -q -m cp transcript_counts_${sample_name}/* gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/transcript_counts/${sample_name}/
    >>>

    output {
        File transcript_counts_tar = "transcript_counts_${sample_name}.tar.gz"
    }

    runtime {
        docker: "gcr.io/genomics-xavier/bulk_rna_seq"
        zones: "us-east1-b us-east1-c us-east1-d"
		cpu: 12
  		memory: "8GB"
  		preemptible: 2
  		disks: "local-disk 80 HDD"
    }
}

task counts_and_differential_expression_output {
    String experiment_name
    String experiment_type
    Array[File] transcript_counts_tar
    String organism
    String samples_described_file
    String samples_compared_file

    Int counts_and_differential_expression_output_RAM
    Int counts_and_differential_expression_output_disk_space

    command <<<
        set -e

        # NOTE: not yet sure how to use output from perform_kallisto_quantification directly here,
        #       doing so would probably be a better solution than copying from google bucket

        # copy transcript counts from google bucket
        # TODO: should probably change this to using output from last step
        gsutil -m cp -r gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/transcript_counts ./

        if [ ! -z "${samples_described_file}" ]; then
            # copy samples described and compared files
            gsutil -m cp ${samples_described_file} ./samples_described.tsv
            gsutil -m cp ${samples_compared_file} ./samples_compared.tsv
        else
            touch ./samples_compared.tsv
            touch ./samples_described.tsv
        fi

        # run kallisto_output
        mkdir -p ./kallisto_output
        Rscript /scripts/get_kallisto_counts.R ./transcript_counts ${organism} samples_described.tsv samples_compared.tsv ./kallisto_output
        # tar output
        tar zcfv kallisto_output.tar.gz ./kallisto_output
        # copy output to google bucket
        gsutil -m cp -r ./kallisto_output gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/

        # run edgeR
        mkdir -p ./edgeR_output
        if [ -s samples_described.tsv ] && [ -s samples_compared.tsv ]; then
            Rscript /scripts/run_edgeR.R ./kallisto_output/kallisto_gene_counts_rounded.csv samples_described.tsv samples_compared.tsv ./edgeR_output
        fi

        # tar output
        tar zcfv edgeR_output.tar.gz ./edgeR_output
        # copy output to google bucket
        if [ -s samples_described.tsv ] && [ -s samples_compared.tsv ]; then
            gsutil -m cp -r ./edgeR_output gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/
        fi
    >>>

    output {
        File kallisto_counts = "kallisto_output.tar.gz"
        File edgeR_output = "edgeR_output.tar.gz"
    }

    runtime {
        docker: "gcr.io/genomics-xavier/bulk_rna_seq"
        zones: "us-east1-b us-east1-c us-east1-d"
		cpu: 4
  		memory: "40GB"
  		preemptible: 2
  		disks: "local-disk 100 HDD"
    }
}

task perform_multiqc {
    String experiment_name
    String experiment_type
    Array[File] transcript_counts_tar

    Int perform_multiqc_RAM
    Int perform_multiqc_disk_space

    command <<<
        set -e

        mkdir -p ./multiqc_output

        # copy fastqc info
        gsutil -q -m cp -r gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/fastqc ./multiqc_output

        # copy kallisto info
        gsutil -q -m cp -r gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/transcript_counts ./multiqc_output

        # run multiqc
        cd ./multiqc_output
        multiqc ./
        cd ..

        # push multiqc report to google bucket
        gsutil -q -m cp ./multiqc_output/multiqc_report.html gs://genomics_xavier_bucket/${experiment_type}/${experiment_name}/
    >>>

    output {
        File multiqc_output = "multiqc_output/multiqc_report.html"
    }

    runtime {
        docker: "gcr.io/genomics-xavier/bulk_rna_seq"
        zones: "us-east1-b us-east1-c us-east1-d"
		cpu: 4
  		memory: "8GB"
  		preemptible: 2
  		disks: "local-disk 80 HDD"
    }
}
