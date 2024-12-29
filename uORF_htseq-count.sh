

# Create output directories for counts
mkdir -p counts/TEup counts/TEnc counts/TEdown

# For TEnc
for bam in TEnc/LZT10*_uniq_sort_TEnc.bam; do
    sample=$(basename $bam _uniq_sort_TEnc.bam)
    htseq-count -f bam -r name -s yes \
        -t sequence_feature \
        -i ID \
        -m union --nonunique none \
        $bam /path/to/uorf.gff > counts/TEnc/${sample}_uorf_counts.txt
done

# For TEup
for bam in TEup/LZT10*_uniq_sort_TEup.bam; do
    sample=$(basename $bam _uniq_sort_TEup.bam)
    htseq-count -f bam -r name -s yes \
        -t sequence_feature \
        -i ID \
        -m union --nonunique none \
        $bam /path/to/uorf.gff > counts/TEup/${sample}_uorf_counts.txt
done

# For TEdown
for bam in TEdown/LZT10*_uniq_sort_TEdown.bam; do
    sample=$(basename $bam _uniq_sort_TEdown.bam)
    htseq-count -f bam -r name -s yes \
        -t sequence_feature \
        -i ID \
        -m union --nonunique none \
        $bam /path/to/uorf.gff > counts/TEdown/${sample}_uorf_counts.txt
done
