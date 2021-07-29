source activate MethylDackel

#Check for mbias
#/icgc/dkfzlsdf/analysis/C010/genomes/Mmusculus/mm10/seq/mm10.fa
mkdir -p /omics/groups/OE0219/internal/FetalHierachy/WGBS/QC/methylDackel/mbias
for file in ls /icgc/dkfzlsdf/project/OE0219/mouse_hematopoiesis_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/*fetal*/*/paired/merged-alignment/.merging_0/*_merged.mdup.bam; do
    fname=`basename $file`
    MethylDackel mbias -@ 6 /icgc/dkfzlsdf/analysis/C010/genomes/Mmusculus/GRCm38/seq/GRCm38mm10_PhiX_Lambda.fa $file \
    /omics/groups/OE0219/internal/FetalHierachy/WGBS/QC/methylDackel/mbias/${fname}
    echo $fname
done 

#Extract methylation values with adjusting for mbias
mkdir /omics/groups/OE0219/internal/FetalHierachy/WGBS/QC/methylDackel/methylationCalls/
for file in ls /icgc/dkfzlsdf/project/OE0219/mouse_hematopoiesis_PBAT/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/*fetal*/*/paired/merged-alignment/.merging_0/*_merged.mdup.bam; do
    fname=`basename $file`
    dirname=`echo $fname | awk -F\. '{print $1}'`
    mkdir  /omics/groups/OE0219/internal/FetalHierachy/WGBS/QC/methylDackel/methylationCalls/${dirname}
    MethylDackel extract -@ 5 --OT 6,0,6,140 --OB 2,145,12,150  \
    /icgc/dkfzlsdf/analysis/C010/genomes/Mmusculus/GRCm38/seq/GRCm38mm10_PhiX_Lambda.fa \
    $file \
    -o  /omics/groups/OE0219/internal/FetalHierachy/WGBS/QC/methylDackel/methylationCalls//${dirname}/${fname}
    echo $dirname


    echo $fname
done 

#add chr to seqnames
for file in /omics/groups/OE0219/internal/FetalHierachy/WGBS/QC/methylDackel//methylationCalls/*/*.bedGraph; do
awk '{print "chr" $0}' $file > $file.chr.bedgraph
echo $file 
done 