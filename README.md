# CoverageAnalysis

Analysis of Illumina High troughtput sequencing
Coverage analysis 


# Analysis Illumina sequencing

- zipping all files

      a=0;for folder in $(ls -d */); do cd $folder; for i in $(ls *.fastq); do echo $i;a=$((a + 1));bsub -q normal -L /bin/bash -J zip$a -u ivan.mateusgonzalez@epfl.ch  -N  "gzip $i ";done; cd ..; done

- Trimming reads qual20

      a=0;for folder in $(ls -d */); do cd $folder;a=0;for i in $(ls *_R1.fastq.gz | sort -t'_' -k2); do a=$((a + 1)); b=$(echo $i | cut -d'.' -f1 | cut -d'_' -f1); c=$(echo $i | sed 's/R1/R2/g');echo $i" "$c ; bsub -q normal -L /bin/bash -J TRIMMO$a -u ivan.mateusgonzalez@epfl.ch  -N  "  java -jar /home/imateus2/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $i $c $(echo $i | cut -d'_' -f1,2)_Trimm_Pair.fastq $(echo $i | cut -d'_' -f1,2)_Trimm_Unpair.fastq $(echo $c | cut -d'_' -f1,2)_Trimm_Pair.fastq $(echo $c | cut -d'_' -f1,2)_Trimm_Unpair.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50";  done ; cd ..; done


-  Mapping reads to template

    - Novoalign (option 1) 
        - Create index
        
              /home/imateus2/software/novocraft/novoindex -s 1 -n template_A template_*.nix template_*.fasta  # change A. B. C.

        - Map reads
        
              a=0;for folder in $(ls -d */); do cd $folder; a=0;for i in $(ls *R1*Pair*.fastq | sort -t'_' -k2); do a=$((a + 1)); b=$(echo $i | cut -d'.' -f1 | cut -d'_' -f1); c=$(echo $i | sed 's/R1/R2/g'); echo $i" "$c ;bsub -q normal -L /bin/bash -J novo$a -u ivan.mateusgonzalez@epfl.ch  -N  "/home/imateus2/software/novocraft/novoalign -d/scratch/cluster/monthly/imateus2/V_cholera_noem/EXP2/References/template_C/template_C.nix -f $i $c -o SAM 2> 'stats'$i'.txt' > $i'.sam' ";  done ; cd ..; done

    - Bowtie2 (option 2)  
        - Create index 
                  
                  /home/imateus/software/bowtie2-2.3.2/bowtie2-build /home/imateus/Documents/V_cholerae/test1/ref_VCA0107/reference_template/All_chromosomes-Sa5Y-VCA0107-frt-kan-frt.fasta All_chromosomes-Sa5Y-VCA0107-frt-kan-frt_db

        - Align reads
        
                  # a single file
                  /home/imateus/software/bowtie2-2.3.2/bowtie2 -x /home/imateus/Documents/V_cholerae/test1/ref_VCA0107/reference_template/All_chromosomes-Sa5Y-VCA0107-frt-kan-frt_db -1 4703_S6_R1_Trimm_PairQ1.fastq -2 4703_S6_R2_Trimm_PairQ1.fastq -S 4703_S6_R1_Trimm_PairQ1.fastq_VCA0107.sam --no-unal
        
                  # all files
                  a=0;for i in $(ls *R1*Pair*.fastq | sort -t'_' -k3); do a=$((a + 1)); b=$(echo $i | cut -d'.' -f1 | cut -d'_' -f5); c=$(ls *_R2*$b*.fastq); echo $i" "$c ; /home/imateus/software/bowtie2-2.3.2/bowtie2 -x /home/imateus/Documents/V_cholerae/test1/ref_VCA0107/reference_template/All_chromosomes-Sa5Y-VCA0107-frt-kan-frt_db  -1 $i -2 $c -S $i'_VCA0107.sam' --no-unal ;done


- Compare mapping quality of different thresholds

     - Make alignment with different mapping quality
  
             #no filter
            for i in $(ls *.sam); do echo $i; /home/imateus/software/samtools-1.4/samtools view -bS $i | samtools sort -o $(echo $i | cut -d'.' -f1)_sorted.bam;done
            #Q10
            for i in $(ls *.sam); do echo $i; /home/imateus/software/samtools-1.4/samtools view -bSq 10 $i | samtools sort -o $(echo $i | cut -d'.' -f1)_sortedQ10.bam;done
            #Q20
            for i in $(ls *.sam); do echo $i; /home/imateus/software/samtools-1.4/samtools view -bSq 20 $i | samtools sort -o $(echo $i | cut -d'.' -f1)_sortedQ20.bam;done
            #Q30
            for i in $(ls *.sam); do echo $i; /home/imateus/software/samtools-1.4/samtools view -bSq 30 $i | samtools sort -o $(echo $i | cut -d'.' -f1)_sortedQ30.bam;done


     - Prepare files for further analysis
            
            # add RG in tag
            for i in $(ls *.bam); do echo $i; ~/software/bamaddrg/bamaddrg -b $i -r "gr."$(echo $i | cut -d"_" -f3) > $(echo $i | cut -d"." -f1)"RG.bam"; done
            # create index of bam file
            for i in $(ls *RG.bam); do echo $i; /home/imateus/software/samtools-1.4/samtools index $i $(echo $i | cut -d"." -f1)".bai" ; done

     - Coverage calculation per site
      
            for i in $(ls *RG.bam); do echo $i; /home/imateus/software/bedtools2/bin/bedtools genomecov -ibam $i -d -split > $(echo $i | cut -d'_' -f1,6)_detail.txt ; done

     - Compare Coverage per site for different analysis in R.
     
- Final analysis with good parameters

     - Mapping with Q30 convert SAM to BAM

            a=0;for folder in $(ls -d */); do cd $folder; a=0; for i in $(ls *.sam); do echo $i;bsub -q normal -L /bin/bash -J sam$i -u ivan.mateusgonzalez@epfl.ch  -N  " module add UHTS/Analysis/samtools/1.3; samtools view -bSq 30 $i | samtools sort -o $(echo $i | cut -d'.' -f1| cut -d'_' -f1)_sortedQ30.bam ";  done ; cd ..; done

     - Add header bam

            for i in $(ls *.bam); do echo $i; ~/software/bamaddrg/bamaddrg -b $i -r "gr."$(echo $i | cut -d"_" -f3) > $(echo $i | cut -d"." -f1)"RG.bam"; done

     - Make index

            for i in $(ls *RG.bam); do echo $i; /home/imateus/software/samtools-1.4/samtools index $i $(echo $i | cut -d"." -f1)".bai" ; done

     - Coverage_calculation

            for i in $(ls *RG.bam); do echo $i; /home/imateus/software/bedtools2/bin/bedtools genomecov -ibam $i -d -split > $(echo $i | cut -d'_' -f1,6)_detail.txt ; done

     - Zip new coverage file (big files)
            
            for i in $(ls *detail.txt); do echo $i ; gzip $i ;done

     - SNP calling 

            a=0; for i in $(ls *_Trimm_Pair_sortedQ30RG.bam); do echo $i;a=$((a + 1)); freebayes -f /home/imateus/Documents/V_cholerae/test1/ref_VCA0107/reference_template/All_chromosomes-Sa5Y-VCA0107-frt-kan-frt.fasta --min-coverage 10 -F 0.3 -p 1 -K -u -v $(echo $i | cut -d'.' -f1)'.vcf' -b $i ; done
