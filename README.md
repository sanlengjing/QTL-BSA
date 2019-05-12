# QTL-BSA
Bulked Sequence Analysis(BSA) Using QTLv1.0.pl

QTLv1.0.pl, an anlysis pipeline for Bulked Sequence Analysis based on highthroughput whole genome sequencing technology to rapidly locate the quantitative trait loci from two populations, developed in October,2013. If you want to know more about the principle of the method, please refer to the paper: "Takagi Hiroki, Abe Akira, Yoshida Kentaro et,al.QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations(2013),Plant.J." It contains four step: sequence alignment, sequence realignment, snp calling, and snp dentity graph plotting. 

About installation
If you need to anlysis from the step of the sequence alignment and snp calling, please make sure your linux operating system have installed the following software: bwa,samtools, Coval,R and Perl module(Bio::seq module). Otherwise, make sure you have installed R.
Note: you should modify the path to your software in the QTLv1.0.pl

About parameter：

-h      get help;
-f      Reference sequence, fasta file
-C      CPU number, default:    2       
-FL     First left reads (paired reads of First bulked file)
-FR     First right reads (paired reads of First bulked file)
-SL     Second left reads (paired reads of second bulked file)
-SR     Second right reads (paired reads of second bulked file)
-O      Output dir (default present directory)
-q      mapping quality less than is filtered(default:1)
-D      coverage depth more than is filtered in the SNP calling(the three fold of the average coverage depth is recommended, default:100)
-freq   the minimum frequency of a non-reference allele in the SNP calling 
              (a fraction of the three-fold of the average coverage depth is recommended, default:0.01)
-m      the minimum frequency of a non-reference allele in the QTL location (default:0.3)
-n      the reads depth less than is filtered in the QTL location
              (the average coverage depth is recommended, default:5)
-fv     the first bulked snp variants
-sv     the second bulked snp variants
-inc    the increment(default unit:Mb, default value:0.01Mb)
-wsize  the sliding window size(default unit:Mb,default value:2Mb)
-fp     the prefix of first bulked file output
-sp     the prefix of second bulked file output
Output: prefix1.sam prefix1.sorted.bam prefix1.snp.txt prefix2.sam prefix2.sorted.bam prefix2.snp.txt ...

About Use
Quick analysis from the start,please use the command:
perl QTLv1.0.pl -f /path_to_ReferSeq/chr.fa -c 4 -FL H_1.fq -FR H_2.fq -SL L_1.fq -SR L_2.fq -D 30 -freq 0.03 -m 0.3 -n 7 -fp H -sp L
Note: when you use Coval for SNP calling, there are many parameters can be setting to reduce the false positives, please refer to the Coval manua
l and add it in the according command.

If you have obtained SNP variant information from other SNP calling software such as samtools,gatk,.. you can skiped the sequence alignment and S
NP calling step using -fv and -sv. the input file shoud be the following format:
RName   POS     RefBase Variant Variant_Support ReadsCoverage Variant_frequency
Chr1    12038   A       G       2       25      0.08
Chr1    233485  G       T       7       11      0.636363636363636
Chr1    401692  A       T       1       15      0.0666666666666667
Chr1    404436  A       G       4       28      0.142857142857143
RName: the reference sequence name
Pos: the variant occured in the position of reference sequence(1-based leftmost).
RefBase: the reference base
Variant: the variant.
Variant_Support: number of the reads that support the variant in the position.
ReadsCoverage: number of the reads that coveraged in the position.
Variant_frequency: the frequency of the variant(the Variant_Support/ReadsCoverage);
Use the command: 
perl GetQTLv1.pl -f /path_to_ReferSeq/ref.fa -fv VariantFile1.txt -sv VariantFile2.txt -D 30 -freq 0.03 -m 0.3 -n 7 -fp H -sp L

It also can analysis for single bulked genome sequencing data or from single SNP variant information.
for the former,use the command:
perl QTLv1.0.pl -f /path_to_ReferSeq/ref.fa -c 4 -FL H_1.fq -FR H_2.fq -D 30 -freq 0.03 -m 0.3 -n 7 -fp H 
perl QTLv1.0.pl -f /path_to_ReferSeq/ref.fa -c 4 -SL L_1.fq -SR L_2.fq -D 30 -freq 0.03 -m 0.3 -n 7 -sp L

for the latter,use the command: 
perl QTLv1.0.pl -f /path_to_ReferSeq/ref.fa -fv VariantFile1.txt -freq 0.03 -m 0.3 -n 7 -fp H 
perl QTLv1.0.pl -f /path_to_ReferSeq/ref.fa -sv VariantFile2.txt -freq 0.03 -m 0.3 -n 7 -sp L
Note: the parameter -fp or -sp must be assigned according to your input data.


If you have any problem during the using, please contact: 290360262@qq.com
