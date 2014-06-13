#for i in *R1*gz;do java -jar /usr/local/bin/trimmomatic-0.30.jar PE $i ${i/R1/R2} s1_pe s1_se s2_pe s2_se ILLUMINACLIP:/usr/local/share/adapters/TruSeq3-PE.fa:2:30:10; interleave-reads.py s1_pe s2_pe > ${i/R1/interleaved}; rm s*pe; cat s*se > ${i/R1/singletons};rm s*se;done

#for i in *interleaved*gz;do mv $i ${i/.gz/};done
#for i in *singleton*gz;do mv $i ${i/.gz/};done

#chmod u-w *fastq

#for i in *interleaved* *singleton*;do fastq_quality_filter -Q33 -q 30 -p 50 -i $i -o ${i/fastq/qc.fastq};done

#for i in *interleaved*qc.fastq; do extract-paired-reads.py $i;done

#for i in A_MOcculta_*.pe; do cat ${i/interleaved_001.qc.fastq.pe/singletons_001.fastq} >> ${i/interleaved_001.qc.fastq.pe/interleaved_001.qc.fastq.se};done

#rm *qc.fastq

#for i in *interleaved*.pe; do normalize-by-median.py -p -C 20 -k 20 -N 4 -x 4e9 --savetable $i.kh $i;done

#for i in *interleaved*.pe; do normalize-by-median.py -C 20 -k 20 -N 4 -x 4e9 --loadtable $i.kh --savetable $i.kh ${i/.pe/.se};done

#for i in *interleaved*.pe; do filter-abund.py $i.kh -V $i.keep;done
#for i in *interleaved*.pe; do filter-abund.py $i.kh -V ${i/pe/se.keep};done

#rm *kh

#for i in *pe*abundfilt; do extract-paired-reads.py $i;done

#for i in *se.keep.abundfilt; do cat $i >> ${i/se.keep.abundfilt/pe.keep.abundfilt.se};done

#for i in *.abundfilt.pe
#do
#    normalize-by-median.py -p -C 5 -k 20 -N 4 -x 4e9 --savetable $i.kh $i
#    normalize-by-median.py -C 5 -k 20 -N 4 -x 4e9 --loadtable $i.kh ${i/abundfilt.pe/abundfilt.se}
#done

#rm *kh

trim_singletons=*abundfilt.se.keep.trim
velveth trim_output 29,69,10 -fastq -shortPaired A_MOcculta_807_GCCAAT_L008_interleaved_001.qc.fastq.pe.keep.abundfilt.pe.keep.trim -shortPaired2 A_MOcculta_683_GCCAAT_L007_interleaved_001.qc.fastq.pe.keep.abundfilt.pe.keep.trim -shortPaired3 A_MOcculta_342_CGATGT_L006_interleaved_001.qc.fastq.pe.keep.abundfilt.pe.keep.trim -short $trim_singletons
for((k=29;k<=69;k=k+10));do velvetg trim_output_"$k" -cov_cutoff auto -exp_cov auto -ins_length 810 -ins_length2 680 -ins_length3 340 -min_contig_lgth 200;done

#for i in *.keep;do 