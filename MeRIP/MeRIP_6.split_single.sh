set -ue

# sh   bam_split_single_end.sh      target.bam   ../split_bam/
# Get the bam file from the command line
BAM=$1
TARGET_D=$2

FILE=$(basename $BAM)  # remove path
NAME=${FILE%.*}        # remove ".bam"
BAMF=${TARGET_D}/${NAME}_fwd.bam
BAMR=${TARGET_D}/${NAME}_rev.bam

# Forward strand.
# 0x10 - SEQ being reverse complemented
samtools view -bh -F 16 $BAM > $BAMF
samtools index $BAMF

# Reverse strand
samtools view -bh -f 16 $BAM > $BAMR
samtools index $BAMR





