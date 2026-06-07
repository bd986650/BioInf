version 1.0

workflow MappingPipeline {

  File ref
  File r1
  File r2

  call FastQC {
    input:
      r1 = r1,
      r2 = r2
  }

  call BWA {
    input:
      ref = ref,
      r1 = r1,
      r2 = r2
  }

  call SamtoolsSort {
    input:
      bam = BWA.bam
  }

  call Flagstat {
    input:
      bam = SamtoolsSort.sorted_bam
  }

  call QC {
    input:
      flagstat_txt = Flagstat.flagstat
  }

  call FreeBayes {
    input:
      ref = ref,
      bam = SamtoolsSort.sorted_bam
  }

  output {
    File vcf = FreeBayes.vcf
    File qc  = QC.qc_result
    File flagstat = Flagstat.flagstat
  }
}

task FastQC {
  File r1
  File r2

  command <<<
    fastqc ~{r1} ~{r2}
  >>>

  output {
    File report1 = glob("*_fastqc.zip")[0]
  }

  runtime {
    docker: "biocontainers/fastqc:v0.11.9_cv8"
  }
}

task BWA {
  File ref
  File r1
  File r2

  command <<<
    bwa mem ~{ref} ~{r1} ~{r2} > out.sam
    samtools view -bS out.sam > out.bam
  >>>

  output {
    File bam = "out.bam"
  }
}

task SamtoolsSort {
  File bam

  command <<<
    samtools sort ~{bam} -o sorted.bam
  >>>

  output {
    File sorted_bam = "sorted.bam"
  }
}

task Flagstat {
  File bam

  command <<<
    samtools flagstat ~{bam} > flagstat.txt
  >>>

  output {
    File flagstat = "flagstat.txt"
  }
}

task QC {
  File flagstat_txt

  command <<<
    RATE=$(grep " mapped (" ~{flagstat_txt} | head -1 | sed -E 's/.*\(([0-9.]+)%.*/\1/')
    THRESHOLD=90

    echo "Mapped: $RATE" > qc.txt

    awk -v r=$RATE 'BEGIN{
      if(r>90) print "OK";
      else print "not OK";
    }' >> qc.txt
  >>>

  output {
    File qc_result = "qc.txt"
  }
}

task FreeBayes {
  File ref
  File bam

  command <<<
    freebayes -f ~{ref} ~{bam} > variants.vcf
  >>>

  output {
    File vcf = "variants.vcf"
  }
}