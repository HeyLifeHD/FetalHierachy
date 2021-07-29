###run Mbias and methylation calling
conda activate MethylDackel

#specify genome reference
genome='/icgc/dkfzlsdf/analysis/C010/genomes/Mmusculus/GRCm38/seq/GRCm38mm10_PhiX_Lambda.fa'
#specify input and output
outdir='/omics/groups/OE0219/internal/FetalHierachy'
input='/icgc/dkfzlsdf/project/OE0219/mouse_hematopoiesis_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/*fetal*/*/paired/merged-alignment/.merging_0/*_merged.mdup.bam'

##Check for mbias
mkdir -p ${outdir}/WGBS/QC/methylDackel/mbias
for file in $input; do
    fname=`basename $file`
    MethylDackel mbias -@ 6 /icgc/dkfzlsdf/analysis/C010/genomes/Mmusculus/GRCm38/seq/GRCm38mm10_PhiX_Lambda.fa $file \
    ${outdir}/WGBS/QC/methylDackel/mbias/${fname}
    echo $fname
done 

##Extract methylation values and adjust for mbias
mkdir ${outdir}/WGBS/QC/methylDackel/methylationCalls/
for file in $input; do
    fname=`basename $file`
    dirname=`echo $fname | awk -F\. '{print $1}'`
    mkdir  ${outdir}/WGBS/QC/methylDackel/methylationCalls/${dirname}
    MethylDackel extract -@ 5  --OT 7,0,6,140 --OB 0,145,12,150 \
    $genome $file \
    -o  ${outdir}/WGBS/QC/methylDackel/methylationCalls//${dirname}/${fname}
    echo $dirname
    echo $fname
done 

##add chr to seqnames
for file in ${outdir}/WGBS/QC/methylDackel//methylationCalls/*/*.bedGraph; do
awk '{print "chr" $0}' $file > $file.chr.bedgraph
echo $file 
done 