#!/bin/bash
############################################################################################
#                                                                                          
#  Title:    Automated cfTCR Data Processing Pipeline                                       
#  Authors:  Dong Xu, Xinyue Liang                                                         
#  Date:     2025-11-07                                                                    
#  Version:  2024-05-29                                                                    
#                                                                                          
#  Description:                                                                             
#    This pipeline automates the processing of cfTCR sequencing data, including quality     
#    control, adapter trimming, UMI extraction, read merging, genome mapping (hg38 & IMGT), 
#    extraction of poorly mapped reads, non-specific gene filtering, and downstream TCR/BCR 
#    reconstruction using TRUST4 and MiXCR.                                                 
#                                                                                          
#  Updates:                                                                                 
#    - v1  (2024-05-07):  Initial pipeline build for automated cfTCR processing             
#    - v2  (2024-05-09):  Improved MiXCR handling                                           
#    - v3  (2024-05-13):  Added merging and bug fixes                                       
#    - v4  (2024-05-18):  Adjusted ratio calculation logic                                  
#    - v5  (2024-05-20):  Added FASTQ merging strategy                                      
#    - v6  (2024-05-21):  Removed UMI preprocessing                                         
#    - v7  (2024-05-23):  Kept only TRUST4 and MiXCR steps                                  
#    - v8  (2024-05-25):  Paused MiXCR + VDJtools                                           
#    - v9  (2024-05-28):  Fixed mixcr deduplication logic                                   
#    - v10 (2024-05-29):  Consolidated pipeline and removed redundant genes                 
#    - v11 (2024-06-03):  Added environment setup                                           
#    - v12 (2024-06-05):  Added cdhit clustering & consensus selection                      
#    - v13 (2024-06-11):  Integrated visualization-ready outputs                            
#    - v14 (2024-08-19):  Retained poorly mapped reads only                                 
#                                                                                          
#  Notes:                                                                                  
#    - Ensure all paths are correct and accessible                                          
#    - Confirm reference files and indexes exist                                            
#    - Adjust CPU and memory allocation per system                                          
#    - Compatible with SLURM job scheduling system                                          
#                                                                                          
############################################################################################

# -------------------------------
# I. Environment Initialization
# -------------------------------
unset __conda_setup
conda activate base

# Update PATH for dependencies
export PATH=/user/miniconda3/bin:$PATH
export PATH=/user/miniconda3/Trust4:$PATH
export PATH=/user/miniconda3/CDhit/cd-hit:$PATH
export PATH=/user/miniconda3/trimmomatic:$PATH
export PATH=/user/miniconda3/picard:$PATH
export PATH=/user/miniconda3/fgbio:$PATH
export PATH=/user/miniconda3/mixcr-3.0.13:$PATH
export PATH=/user/miniconda3/vdjtools:$PATH

# Resource configuration (auto-detect from SLURM if available)
t_cpu=${SLURM_CPUS_PER_TASK:-101}
M_zhao=$((${SLURM_MEM_PER_NODE:-160000}))
echo "CPU threads: $t_cpu"
echo "Memory (MB): $M_zhao"

start_time=$(date +%s)
echo "Script started at: $(date)"
echo "----- Starting cfTCR Processing Pipeline -----"

# -------------------------------
# II. Define Paths and References
# -------------------------------
trust4="/user/miniconda3/Trust4/TRUST4-master/run-trust4"
trimmomatic="java -Xmx30G -jar /user/miniconda3/trimmomatic/trimmomatic.jar"
picard="java -Xmx30G -jar /user/miniconda3/picard/picard.jar.1"
fgbio="java -Xmx30G -jar /user/miniconda3/fgbio/fgbio-2.2.1.jar"
mixcr="java -Xmx30G -jar /user/miniconda3/mixcr-3.0.13/mixcr.jar"
vdjtools="java -Xmx30G -jar /user/miniconda3/vdjtools/vdjtools-1.2.1.jar"
cdhit="/user/miniconda3/CDhit/cdhit/cd-hit-est"

# Project directories
datadir="/file_path/Rawdata/DATA"
codedir="/prepare_file_path"
workdir="/file_path/Samples/DATA"

# Reference files
hg38_ref="/prepare_file_path/download/hg38.fa"
imgt_ref="/user/miniconda3/Trust4/TRUST4-master/human_IMGT+C.fa"
adapter_path="/user/miniconda3/trimmomatic/adapters/TruSeq3-PE.fa"

mkdir -p $workdir/
cd $datadir || { echo "ERROR: Cannot access $datadir"; exit 1; }

# -------------------------------
# III. Main Loop — Sample Processing
# -------------------------------
for sample_dir in ${datadir}/*/; do    
    sample_id=$(basename $sample_dir)
    echo ">>> Processing sample: ${sample_id}"
    cd $workdir || { echo "ERROR: Cannot access $workdir"; exit 1; }
    mkdir -p "$workdir/$sample_id/outfile" "$workdir/$sample_id/output"
    outdir="$workdir/$sample_id/output"
    outtmp="$workdir/$sample_id/outfile"

    if [ ! -f "$outtmp/${sample_id}_R1.fastq" ]; then
        # ===============================================================
        # Step 1: FASTQC Quality Control
        # ===============================================================
        echo "[1] Running FastQC..."
        fastqc -o "$outtmp" "${sample_dir}${sample_id}_raw_1.fq.gz" "${sample_dir}${sample_id}_raw_2.fq.gz"
        echo "[1] FastQC completed."
        
        # ===============================================================
        # Step 2: Adapter Trimming with Trimmomatic
        # ===============================================================
        echo "[2] Trimming adapters and low-quality reads..."
        $trimmomatic PE -phred33 -threads "$t_cpu" \
            "${sample_dir}${sample_id}_raw_1.fq.gz" "${sample_dir}${sample_id}_raw_2.fq.gz" \
            "$outtmp/${sample_id}_raw_1_val_1.fq.gz" "$outtmp/${sample_id}_raw_1_unpaired.fq.gz" \
            "$outtmp/${sample_id}_raw_2_val_2.fq.gz" "$outtmp/${sample_id}_raw_2_unpaired.fq.gz" \
            ILLUMINACLIP:${adapter_path}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        echo "[2] Trimming complete."
        fastqc -o "$outtmp" "$outtmp/${sample_id}_raw_1_val_1.fq.gz" "$outtmp/${sample_id}_raw_2_val_2.fq.gz"
        
        # ===============================================================
        # Step 3: Convert FASTQ to uBAM and Extract UMIs
        # ===============================================================
        echo "[3] Converting FASTQ to uBAM and extracting UMIs..."
        $picard FastqToSam \
            F1="$outtmp/${sample_id}_raw_1_val_1.fq.gz" \
            F2="$outtmp/${sample_id}_raw_2_val_2.fq.gz" \
            PL=illumina LB=Human RG=Human \
            O="$outtmp/${sample_id}.ubam" \
            SM=$sample_id TMP_DIR="$outtmp/tmp_dir"

        $fgbio --tmp-dir="$outtmp/tmp_dir" ExtractUmisFromBam \
            -i "$outtmp/${sample_id}.ubam" \
            -o "$outtmp/${sample_id}.umi.ubam" \
            -r 5M2S+T 5M2S+T -s RX -t ZA ZB
        echo "[3] UMI extraction completed."

        # ===============================================================
        # Step 4: Convert SAM + Add UMI tags + Convert back to FASTQ
        # ===============================================================
        echo "[4] Adding UMI tags and converting SAM to FASTQ..."
        samtools view -h "$outtmp/${sample_id}.umi.ubam" > "$outtmp/${sample_id}.umi.sam"

        # Append UMI tags to read names
        awk 'BEGIN{FS=OFS="\t"} {
            for(i=1;i<=NF;i++) if($i ~ /^RX:Z:/){rx=substr($i,6); $1=$1"_"rx; break}
            print
        }' "$outtmp/${sample_id}.umi.sam" > "$outtmp/${sample_id}.umis.sam"

        $picard SamToFastq \
            INPUT="$outtmp/${sample_id}.umis.sam" \
            FASTQ="$outtmp/${sample_id}_R1.fastq" \
            SECOND_END_FASTQ="$outtmp/${sample_id}_R2.fastq" \
            INCLUDE_NON_PF_READS=true TMP_DIR="$outtmp/tmp_dir"
        echo "[4] UMI tagging and FASTQ conversion complete."
    else
        echo "[✓] Steps 1–4 already processed. Skipping..."
    fi

    # ===============================================================
    # Step 5: Merge paired-end reads using FLASH
    # ===============================================================
    echo "[5] Merging paired reads..."
    if [ ! -f "$outtmp/${sample_id}_merged_reads.extendedFrags.fastq" ]; then
        flash -t "$t_cpu" -m 3 -M 150 \
            -o "$sample_id/outfile/${sample_id}_merged_reads" \
            "$outtmp/${sample_id}_R1.fastq" "$outtmp/${sample_id}_R2.fastq"
        outputfile="$outtmp/${sample_id}_merged_reads.extendedFrags.fastq"
        echo "[5] Merge complete."
    else
        echo "[✓] Merged FASTQ exists. Skipping merge."
    fi
    outputfile="$outtmp/${sample_id}_merged_reads.extendedFrags.fastq" # Data from a single experiment

    # # *************************************************************************************↑Merging two-end processing for downstream analysis
    # # *************************************************************************************↓Combine data from two experiments
    # input_file1="$workdir/${sample_ids}/outfile/${sample_ids}_merged_reads.extendedFrags.fastq"
    # input_file2="$workdir/${sample_ids}-2/outfile/${sample_ids}-2_merged_reads.extendedFrags.fastq"
    # if [ -f "$input_file1" ] && [ -f "$input_file2" ]; then    
    #     # sample_id=$(basename $sample_dir)
    #     outputfile="$workdir/${sample_ids}/outfile/${sample_ids}_merged_trimmed_reads.fastq"
    #     # Merge two FASTQ files
    #     cat "$input_file1" "$input_file2" > "$outputfile"
    #     echo "The files exist, merge 1+2 successfully!"
    # else
    #     echo "The files do not exist complete, skip the processing merge 1+2."
    # fi
    # outputfile="$workdir/${sample_id}/outfile/${sample_id}_merged_trimmed_reads.fastq" # Combine data from two experiments
    # *************************************************************************************↓Process Analysis Entry Point

    # ===============================================================
    # Step 6: Mapping to hg38 & IMGT, Extracting Unmapped Reads
    # ===============================================================
    echo "[6] Mapping reads to hg38 reference..."
    # (BWA indexing, mapping, MAPQ analysis, median-based thresholding)
    # Followed by extraction of poorly mapped/unmapped reads,
    # intersection with IMGT database, and calculation of overlap metrics.
    #
    # -- Details Omitted for Brevity --
    #
    # Outputs: *_T0_stats.txt, *_IMGT_aligned_results.bam, *_results_summary.txt
    echo "6.1/ Mapping to hg38 for ${sample_id}..."
    # =============== 检查索引文件是否存在 ===============
    if [ ! -f "${hg38_ref}.bwt" ] || [ ! -f "${hg38_ref}.pac" ] || [ ! -f "${hg38_ref}.ann" ] || [ ! -f "${hg38_ref}.amb" ] || [ ! -f "${hg38_ref}.sa" ]; then
        echo "[INFO] Index file not found. Start indexing..."
        bwa index "${hg38_ref}"
    else
        echo "[INFO] Index file already exists. Skipping indexing."
    fi

    # =============== Step ①: 初次比对（T=0） ===============
    T0_sam="$outdir/${sample_id}_aligned_T0.sam"
    T0_primary_sam="$outdir/${sample_id}_T0_primary.sam"
    MAPQ_txt="$outdir/${sample_id}_MAPQ_T0.txt"
    T0_stats="$outdir/${sample_id}_T0_stats.txt"

    if [ ! -f "$T0_stats" ]; then
        echo "[RUNNING] Initial mapping with T=0 ..."
        bwa mem -T 0 -M -t "$t_cpu" -p "$hg38_ref" "$outputfile" > "$T0_sam"
        samtools view -F 0x100 -h "$T0_sam" > "$T0_primary_sam"

        echo "[INFO] Extracting MAPQ distribution..."
        grep "AS:i:" "$T0_primary_sam" | \
            awk '{for(i=1;i<=NF;i++) if($i ~ /^AS:i:/) print $i}' | \
            cut -d':' -f3 | sort | uniq -c | sort -bnr -k2 > "$MAPQ_txt"

        echo "[INFO] Calculating statistics (median, quartile)..."
        awk '{
            for (i = 1; i <= $1; i++) a[++n] = $2;
            sum += $1 * $2
        }
        END {
            ave = sum / n;
            if (n % 2) median = a[(n + 1) / 2]; else median = (a[n/2] + a[n/2 + 1]) / 2;
            print "Median:", median > "'$T0_stats'"
            print "Average:", ave >> "'$T0_stats'"
        }' "$MAPQ_txt"
        # =============== Step ②: 提取 median 并进行重新比对 ===============
        median=$(grep 'Median:' "$T0_stats" | awk '{print $2}')
        if [ -z "$median" ]; then
            echo "[ERROR] Median not found. Please re-check MAPQ stats."
            exit 1
        fi
        echo "[INFO] Using median cutoff T=$median"
    else
        echo "[ERROR] No T0_stats file found. Exiting..."
        exit 1
    fi

    aligned_T_median="$outdir/${sample_id}_aligned_T_median.sam"
    primary_T_median="$outdir/${sample_id}_aligned_T_median_primary.sam"
    unmap_bam="$outdir/${sample_id}_unmap_TCRloci.bam"

    if [ ! -f "$unmap_bam" ]; then
        bwa mem -t "$t_cpu" -M -T "$median" "$hg38_ref" "$outputfile" > "$aligned_T_median"
        samtools view -F 0x100 -h "$aligned_T_median" > "$primary_T_median"
        samtools view -bS -f 4 "$primary_T_median" > "$unmap_bam"
        samtools view "$unmap_bam" | cut -f1 > "$outdir/${sample_id}_unmatched_reads.txt"

        echo "[INFO] Extracting unmatched reads..."
        seqkit grep -f "$outdir/${sample_id}_unmatched_reads.txt" "$outputfile" > "$outdir/${sample_id}_unmatched.fastq"
        total_unmatched_count=$(wc -l < "$outdir/${sample_id}_unmatched.fastq" | awk '{print $1/4}')
        echo "Total unmatched reads: $total_unmatched_count"
    fi

    # =============== Step ③: Mapping to IMGT ===============
    echo "6.2/ Mapping to IMGT for ${sample_id}..."
    if [ ! -f "${imgt_ref}_cleaned.fa.bwt" ]; then
        echo "[INFO] IMGT index not found. Creating..."
        bwa index "${imgt_ref}_cleaned.fa"
    else
        echo "[INFO] IMGT index exists. Skipping."
    fi

    imgt_sam="$outdir/${sample_id}_IMGT_aligned.sam"
    imgt_bam="$outdir/${sample_id}_IMGT_aligned_results.bam"

    if [ ! -f "$imgt_bam" ]; then
        bwa mem -t "$t_cpu" -M "${imgt_ref}_cleaned.fa" "$outputfile" > "$imgt_sam"
        sed -i '/^@SQ.*SN:IGH/d' "$imgt_sam"
        samtools view -buS "$imgt_sam" | samtools sort -o "$outdir/${sample_id}_aligned_to_imgt_primary_only_sorted.bam"
        samtools view -b -F 4 "$outdir/${sample_id}_aligned_to_imgt_primary_only_sorted.bam" > "$imgt_bam"
    fi

    # =============== Step ④: 汇总结果统计 ===============
    echo "6.3/ Calculating results summary for ${sample_id}..."
    total_IMGT_aligned_count=$(samtools view -c "$imgt_bam")
    samtools view "$imgt_bam" | cut -f1 > "$outdir/${sample_id}_imgt_reads.txt"

    awk 'NR==FNR{a[$0];next} $0 in a' \
        "$outdir/${sample_id}_imgt_reads.txt" \
        "$outdir/${sample_id}_unmatched_reads.txt" \
        > "$outdir/${sample_id}_paired_awk.txt"

    overlap_count_awk=$(wc -l < "$outdir/${sample_id}_paired_awk.txt")
    total_reads_count=$(cat "$outputfile" | wc -l | awk '{print $1/4}')
    overlap_ratio=$(awk -v o="$overlap_count_awk" -v t="$total_unmatched_count" 'BEGIN{printf "%.4f",o/t}')
    unmatched_hg38_ratio=$(awk -v u="$total_unmatched_count" -v t="$total_reads_count" 'BEGIN{printf "%.4f",u/t}')

    echo "Overlap ratio: $overlap_ratio"
    echo "Unmatched ratio: $unmatched_hg38_ratio"

    # =============== Step ⑤: GTF注释分析 ===============
    bedtools bamtobed -i "$outdir/${sample_id}_aligned_T_median_primary.bam" > "$outdir/${sample_id}_aligned_T_median_primary.bed"
    bedtools intersect -a "$outdir/${sample_id}_aligned_T_median_primary.bed" \
                    -b "/prepare_file_path/download/gencode.v45.annotation.gtf" \
                    -wa -wb > "$outdir/${sample_id}_annotated_reads.bed"

    awk '$9=="gene"{print $1"\t"$2"\t"$3"\t"$4"\t"$19"\t"$20}' "$outdir/${sample_id}_annotated_reads.bed" | \
        awk '{$6=gensub(/[";]/,"","g",$6);print}' > "$outdir/${sample_id}_awkchr1.bed"

    awk '$6 ~ /^TR[ABDG][VDJ]/ {print}' $outdir/${sample_id}_awkchr1.bed > $outdir/${sample_id}_tcr.bed
    awk '$6 !~ /^TR[ABDG][VDJ]/ {print}' "$outdir/${sample_id}_awkchr1.bed" > "$outdir/${sample_id}_non_TCRgene.bed"
    awk '{count[$6]++;chr[$6]=$1} END{for(g in count) print chr[g],g,count[g]}' "$outdir/${sample_id}_non_TCRgene.bed" \
        | sort -k1,1V > "$outdir/${sample_id}_genomicgene.txt"

    total_reservedGENE_except_TCR=$(awk '{sum+=$3}END{print sum}' "$outdir/${sample_id}_genomicgene.txt")
    reservedGENE_ratio=$(awk -v r="$total_reservedGENE_except_TCR" -v t="$total_reads_count" 'BEGIN{printf "%.4f",r/t}')

    samtools view $outdir/${sample_id}_IMGT_aligned_results.bam 
        | awk '{print $1}' > $outdir/${sample_id}_temp_sequence_names_from_bam.txt
    awk '{print $4}' $outdir/${sample_id}_non_TCRgene.bed > $outdir/${sample_id}_temp_sequence_names_from_bed.txt 
    imgt_reservedGENE_except_TCR_overlap=$(grep -Fwf $outdir/${sample_id}_temp_sequence_names_from_bed.txt $outdir/${sample_id}_temp_sequence_names_from_bam.txt | wc -l)
    overlap_IMGT_ratio=$(awk -v overlap="$imgt_reservedGENE_except_TCR_overlap" -v IMGT="$total_IMGT_aligned_count" 'BEGIN {printf "%.4f", overlap / IMGT}')

    # =============== Step ⑥: 结果输出 ===============
    {
        echo "Sample: $sample_id"
        echo "Total reads count: $total_reads_count"
        echo "Total hg38 unmatched count: $total_unmatched_count"
        echo "Total IMGT aligned count: $total_IMGT_aligned_count"
        echo "Overlap ratio(overlap_count_awk:total_unmatched_hg38_count): $overlap_ratio"
        echo "Unmatched ratio: $unmatched_hg38_ratio"
        echo "Non-specific capture rate: $reservedGENE_ratio"
        echo "Non-specific capture rate(total_reservedGENE_except_TCR:total_reads_count): $reservedGENE_ratio"
        echo "Overlap IMGT(imgt&reservedGENE_except_TCR_overlap:total_IMGT_aligned_count): $overlap_IMGT_ratio"
    } > "$outdir/${sample_id}_results_summary.txt"
    echo "[SUCCESS] All steps for ${sample_id} completed!"

    # =====================================================================
    # Step 7: FASTQ cleanup - remove reserved genes and deduplicate reads
    # Description:
    #   1. Ensure each record starts with '@'
    #   2. Deduplicate reads based on full sequence + UMI tag
    #   3. Output deduplicated FASTQ and a record of duplicates
    # ======================================================================

    # Step ①: Ensure each record starts with '@'
    # Convert original combined.fastq → combined_fixed.fastq
    echo "Step 7.1: Fixing record headers for ${sample_id}..."
    awk 'BEGIN {FS="\n"; RS="@"} NR>1 {print "@"$0}' \
        "${outdir}/${sample_id}_combined.fastq" > \
        "${outdir}/${sample_id}_combined_fixed.fastq"

    # Step ②: Deduplicate based on sequence + UMI
    echo "Step 7.2: Removing duplicates for ${sample_id}..."

    input_fastq="${outdir}/${sample_id}_combined_fixed.fastq"
    output_fastq="${outdir}/${sample_id}_combined_deduped.fastq"
    output_dup="${outdir}/${sample_id}_dup.txt"

    # Ensure output directory exists
    mkdir -p "$(dirname "$output_fastq")"

    # Run AWK to deduplicate
    awk '
        BEGIN {
            FS = "\n";
            RS = "@";  # Use "@" as record separator.
        }
        NR > 1 {  # Skip the first empty record
            id = $1;
            seq = $2;
            plus = $3;
            qual = $4;

            # Extract the last part of the identifier (e.g., UMI tag)
            split(id, parts, "_");
            tag = parts[length(parts)];

            # Create a unique key using sequence + tag
            key = seq "_" tag;

            if (!seen[key]++) {
                # First occurrence of this key → write to deduped FASTQ
                print "@" id "\n" seq "\n" plus "\n" qual > "'$output_fastq'"
            }
        }
        END {
            # Write duplicate records summary
            for (key in seen) {
                if (seen[key] > 1) {
                    print key, seen[key] > "'$output_dup'"
                }
            }
        }
    ' "$input_fastq"

    echo "Step 7 completed for ${sample_id}."
    echo "  - Deduplicated FASTQ: ${output_fastq}"
    echo "  - Duplicate records:  ${output_dup}"

    # Rename deduplicated FASTQ for downstream use
    mv "${output_fastq}" "${outdir}/${sample_id}_combined.fq"

    echo "7/ Processing overall FASTQ files for ${sample_id} (ATCG_ONLY) successful!"

    # ===============================================================
    # Step 8: TRUST4 TCR/BCR Reconstruction
    # ===============================================================
    if [ ! -f "$outdir/${sample_id}_trust4/TRUST_${sample_id}_combined_report.tsv" ]; then
        echo "[8] Running TRUST4..."
        $trust4 \
            -u "$outdir/${sample_id}_combined.fq" \
            -f "/user/miniconda3/Trust4/TRUST4-master/hg38_bcrtcr.fa" \
            --ref "${imgt_ref}_cleaned.fa" \
            --od "$outdir/${sample_id}_trust4" \
            -t "$t_cpu"
        echo "[8] TRUST4 processing complete."
    else
        echo "[✓] TRUST4 output exists. Skipping."
    fi

    # ===============================================================
    # Step 9: MiXCR Clonotype Analysis(Optional)
    # ===============================================================
    if [ ! -f "$outdir/${sample_id}_mixcr_results/${sample_id}.clonotypes.ALL.txt" ]; then
        echo "[9] Running MiXCR..."
        mkdir -p "$outdir/${sample_id}_mixcr_results"
        $mixcr analyze shotgun \
            --species hsa \
            --threads "$t_cpu" \
            --force-overwrite \
            --starting-material dna \
            "$outdir/${sample_id}_combined.fq" \
            "$outdir/${sample_id}_mixcr_results/${sample_id}"
        echo "[9] MiXCR completed successfully."
    else
        echo "[✓] MiXCR results exist. Skipping."
    fi

    echo ">>> Completed processing of $sample_id"
    echo "----------------------------------------------"
done

# ===============================================================
# IV. Completion and Runtime Summary
# ===============================================================
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Pipeline completed successfully at: $(date)"
echo "Total runtime: ${elapsed_time} seconds"
echo "All samples processed!"
############################################################################################
