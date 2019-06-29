#!/usr/bin/perl -w
#Modified:30,Nov
#Add:one graph for each chromosome and correct the output order of each chromosome
use strict;
use warnings;
use Cwd;
use Carp;
use Data::Dumper;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

my $usage = <<_USAGE_;

###############################################################################################################################################################
#
# Usage: perl QTLv1.1.pl [options]
#       -f      Reference sequence, fasta file
#       -xmax  The maximum length of chromosome in reference sequence, set the upper limit of x-coordinate for plotting(default unit:Mb,default:50)
#       -C      CPU number, default: 2  
#       -FL     First left reads (paired reads of First bulked file)
#       -FR     First right reads (paired reads of First bulked file)
#       -SL     Second left reads (paired reads of second bulked file)
#       -SR     Second right reads (paired reads of second bulked file)
#       -O      Output dir (default present directory)
#       -q      mapping quality less than is filtered(default:1)
#       -D      coverage depth more than is filtered in the SNP calling(the three fold of the average coverage depth is recommended, default:100)
#       -freq   the minimum frequency of a non-reference allele in the SNP calling
#                                       (a fraction of the three-fold of the average coverage depth is recommended, default:0.01)
#       -m      the minimum Frequency of a non-reference allele in the QTL location (default:0.3) 
#       -rm     the reads coverage depth more than is filtered in the QTL location(default:100)
#       -rn     the reads coverage depth less than is filtered in the QTL location
#                                       (the average coverage depth is recommended, default:5)
#       -fv     the first bulked snp variants
#       -sv     the second bulked snp variants
#       -inc    the increment (default unit:Mb,default value: 0.01Mb)
#       -wsize  the sliding window size (default unit:Mb,default value: 2Mb)
#       -fp     the prefix of first bulked file output
#       -sp     the prefix of second bulked file output
#       -cov    the average sequencing depth(if using any other software for mapping and SNP calling, make sure this option has been assigned. de
fault:0)
# Output: prefix1.sam prefix1.sorted.bam prefix1.snp.txt prefix2.sam prefix2.sorted.bam prefix2.snp.txt ...
#
# Examples: perl QTLv1.1.pl -f /path_to_ReferSeq/chr.fa -C 4 -FL H_1.fq -freqR H_2.fq -SL L_1.fq -SR L_2.fq -D 30 -freq 0.03 -m 0.3 -rn 7 -fp H -
sp L
#           perl QTLv1.1.pl -f /path_to_ReferSeq/chr.fa -fv 18-weisui -sv 18-yisui -D 60 -m 0.3 -rn 7 -fp wei -sp yi -cov 20
#
#################################################################################################################################################
##############

_USAGE_
        ;

my $help_flag;
my $reference;
my $cpu;
my $f1_lef;
my $f1_rig;
my $f2_lef;
my $f2_rig;
my $workdir="";
my $output_directory="";
my $startdir="";
my $mq;
my $depth;
my $freq;
my $readsmax;
my $readsmin;
my $prefix1;
my $prefix2;
my $minf;
my %SeqIndex;
my $bulkedsnp1;
my $bulkedsnp2;
my $increment;
my $wsize;
my $cov;
my $xmax;

unless (@ARGV){
        die $usage;
        }

&GetOptions(    'h' => \$help_flag,
                                                'f=s' => \$reference,
                                                'C=i'   => \$cpu,
                                                'FL=s'  => \$f1_lef,
                                                'FR=s'  => \$f1_rig,
                                                'SL=s'  => \$f2_lef,
                                                'SR=s'  => \$f2_rig,
                                                'O=s'   => \$output_directory,
                                                'q=i'   => \$mq,
                                                'D=i'   => \$depth,
                                                'freq=f'=> \$freq,
                                                'm=f'   => \$minf,
                                                'rm=i'  => \$readsmax,
                                                'rn=i'  => \$readsmin,
                                                'fv=s'  => \$bulkedsnp1,
                                                'sv=s'  => \$bulkedsnp2,
                                                'inc=f' => \$increment,
                                                'wsize=f'=>\$wsize,
                                                'fp=s'  => \$prefix1,
                                                'sp=s'  => \$prefix2,
                                                'cov=i' => \$cov,
                                                'xmax=i'=> \$xmax,
                                        );

if ($help_flag){ die $usage; }

sub process_cmd{
      my ($cmd) = @_;
      print STDERR "CMD:$cmd\n";
      my $ret = system($cmd);
      if($ret){
          confess "Error, cmd: $cmd died with ret $ret\n";
        }
        return ;
}       

# parameter: $reference,$cpu,$f1_lef,$f1_rig,$workdir,$mq,$freq,$depth,$readsmin,$readsmax,$prefix1
sub alignment{                     
                my($RF,$c,$lef_fq,$rig_fq,$outdir,$maq,$fre,$d,$prefix)=@_;
                my $cmd = "bwa index -a is $RF";
                &process_cmd($cmd) unless (-e "$reference.bwt");
                        $cmd = "bwa aln -t $c $RF $lef_fq > $outdir/$prefix-lef.sai";
                &process_cmd($cmd) unless (-e "$outdir/$prefix-lef.sai");
                        $cmd = "bwa aln -t $c $RF $rig_fq > $outdir/$prefix-rig.sai";
                &process_cmd($cmd) unless (-e "$outdir/$prefix-rig.sai");
                        $cmd = "bwa sampe $RF $outdir/$prefix-lef.sai $outdir/$prefix-rig.sai $lef_fq $rig_fq > $outdir/$prefix.sam";
                &process_cmd($cmd) unless (-e "$outdir/$prefix.sam");
                        $cmd = "samtools view -bS $outdir/$prefix.sam > $outdir/$prefix.bam";
                &process_cmd($cmd) unless (-e "$outdir/$prefix.bam");
                        $cmd = "samtools sort $outdir/$prefix.bam $outdir/$prefix.sorted";
                &process_cmd($cmd) unless (-e "$outdir/$prefix.sorted.bam");
                        $cmd = "samtools depth $outdir/$prefix.sorted.bam > $prefix.depth";
                &process_cmd($cmd) unless (-e "$outdir/$prefix.depth");
                        $cmd = "coval refine -r $RF -p $prefix.refined -mq $maq $outdir/$prefix.sorted.bam -s";
                &process_cmd($cmd) unless (-e "$prefix.refined");
                        $cmd = "coval call-sam $prefix.refined.sam -m $d -f $fre -r $RF -p $prefix";
                &process_cmd($cmd) unless (-e "$prefix-snp-C.txt");             
        }

#parameter:$prefix1.snp.txt,$readsmin,$readsmax,$minf,$prefix2.snp.txt,$out1,$out2,$c
sub FilteredSNP{
                my ($fh,$rmin,$rmax,$min,$rh,$out1,$out2,$c)=@_;
                my $cmd = "bwa index -a is $RF";
                &process_cmd($cmd) unless (-e "$reference.bwt");
                        $cmd = "bwa aln -t $c $RF $lef_fq > $outdir/$prefix-lef.sai";
                &process_cmd($cmd) unless (-e "$outdir/$prefix-lef.sai");
                        $cmd = "bwa aln -t $c $RF $rig_fq > $outdir/$prefix-rig.sai";
                &process_cmd($cmd) unless (-e "$outdir/$prefix-rig.sai");
                        $cmd = "bwa sampe $RF $outdir/$prefix-lef.sai $outdir/$prefix-rig.sai $lef_fq $rig_fq > $outdir/$prefix.sam";
                &process_cmd($cmd) unless (-e "$outdir/$prefix.sam");
                        $cmd = "samtools view -bS $outdir/$prefix.sam > $outdir/$prefix.bam";
                &process_cmd($cmd) unless (-e "$outdir/$prefix.bam");
                        $cmd = "samtools sort $outdir/$prefix.bam $outdir/$prefix.sorted";
                &process_cmd($cmd) unless (-e "$outdir/$prefix.sorted.bam");
                        $cmd = "samtools depth $outdir/$prefix.sorted.bam > $prefix.depth";
                &process_cmd($cmd) unless (-e "$outdir/$prefix.depth");
                        $cmd = "coval refine -r $RF -p $prefix.refined -mq $maq $outdir/$prefix.sorted.bam -s";
                &process_cmd($cmd) unless (-e "$prefix.refined");
                        $cmd = "coval call-sam $prefix.refined.sam -m $d -f $fre -r $RF -p $prefix";
                &process_cmd($cmd) unless (-e "$prefix-snp-C.txt");             
        }

#parameter:$prefix1.snp.txt,$readsmin,$readsmax,$minf,$prefix2.snp.txt,$out1,$out2,$c
sub FilteredSNP{
                my ($fh,$rmin,$rmax,$min,$rh,$out1,$out2,$c)=@_;
                my %hash=();
                open (FH, "<$fh") or die ("Cannot open the file:$fh!\n");
                open (RH, "<$rh") or die ("Cannot open the file:$rh!\n");               
                while(<RH>){
                                        my ($id,$loc,$cov,$rate)=(split)[0,1,5,6];
                                        $loc = $loc/1000000;
                                        $hash{$id}{$loc}="$rate\t"."$cov";
                                }
                                                                
                while(<FH>){
                                        chomp;
                                        my ($chr,$pos,$reads,$rate)=(split)[0,1,5,6];
                                        if($reads >= $rmin && $reads <= $rmax && $rate >= $min){
                                                $pos = $pos/1000000;
                                                print $out1 "$chr\t$pos\t$rate\t$reads\n";
                                                if(defined $hash{$chr}{$pos}){
                                                        print $out2 "$chr\t$pos\t$hash{$chr}{$pos}\n";
                                                }else{
                                                        print $out2 "$chr\t$pos\t0\t$c\n";
                                                        }
                                                }                                       
                        }
        close FH;
        close RH;
        }

sub FilteredSingle{
        my ($fh,$rmin,$rmax,$min,$out) = @_;
                my %hash = ();
                open (FH,"<$fh") or die ("Cannot open the file:$fh!\n");
                while(<FH>){
                my ($id,$loc,$reads,$rate)=(split)[0,1,5,6];
                if($reads >= $rmin && $reads <= $rmax && $rate >= $min){
                                 $loc = $loc/1000000;
                                 print $out "$id\t$loc\t$rate\t$reads\n";
                }
            }
}

#parameter: bulked-SNP1,bulked-SNP1.sorted
sub SortedSNP{
                my ($fh,$out)=@_;
                my %hash;
                my $chrsignal=0;
                open (FH, "<$fh") or die ("Cannot open the file:$fh!\n");
                while(<FH>){
                                chomp;
                                my ($chr,$pos,$rate,$cov)=(split)[0,1,2,3];
                                if($chr =~ /^C/){
                                         $chrsignal=1;
                                }
                                $chr =~ s/[Cc]hr(\d+)/$1/;
                                $hash{$chr}{$pos}="$rate\t"."$cov";
                        }
                        foreach my $key (sort {$a <=> $b} keys %hash){
                                foreach my $p (sort {$a <=> $b} keys %{$hash{$key}}){
                                  if($chrsignal == 1 ){
                                        print $out "Chr$key\t$p\t$hash{$key}{$p}\n";
                                        }else{
                                        print $out "chr$key\t$p\t$hash{$key}{$p}\n";
                                        #print $out "$key\t$p\t$hash{$key}{$p}\n";
                                        }
                                        
                                        }
                                }
        close FH;       
        }
        
#parameter:bulked-SNP1.sorted directory(1/2)
sub SplitSNP{
                my ($fh,$dir) = @_;
                my $cmd = "cut -f 1 $fh | uniq > tit";
                &process_cmd($cmd);
                open (RH,"<tit") or die "Cannot open the file:$!\n";
                while(<RH>){
                                chomp;
                        my $cmd = "mkdir $dir";
                        &process_cmd($cmd) unless (-d "$dir");
                         $cmd = "grep '\\b$_\\b' $fh > $dir/$_";
                        &process_cmd($cmd);
                        }               
                close RH;
        }
        
#parameter: %SeqIndex Chr1... $dir
sub SnpDentity{
        my ($PH,$increment,$wsize,%hash) = @_; 
        my %num = ();
        opendir(DH, $PH) or die "Cannot open the directory:$!\n";
        my @dir = readdir DH;
        print "reading dir...\n";
        foreach my $file (@dir){
                        next if($file =~ /\./);
                        my %stat=();
                        open FH, "<$PH/$file" or die "Cannot open the file:$!\n";
                        open my $out, ">$PH/$file.out" or die;
                        foreach (keys %hash){
                                next if($_ ne "$file");
                                        for(my $i=0;$i<=$hash{$_};$i+=$increment){
                                        $stat{$_}{$i}=0;
                                        }
                        }   
                while(<FH>){
                chomp;
                        next if(/INDEL/);
                my ($chr,$pos,$rate)=(split)[0,1,2];
                        foreach my $key (keys %{$stat{$chr}}){
                                if($pos >= $key && $pos < $key+$wsize){
                                                                $stat{$chr}{$key}+=$rate;
                                                                        $num{$chr}{$key}++;
                                        }
                        }
                }
                print "computing...\n";
                foreach my $p (keys %stat){
                        foreach my $k (sort {$a<=>$b} keys %{$stat{$p}}){
                                if(defined $num{$p}{$k} && $num{$p}{$k} >= 10 ){
                                                $stat{$p}{$k}=$stat{$p}{$k}/$num{$p}{$k};
                                
                                        print $out "$p\t$k\t$stat{$p}{$k}\n";
                        }
                }
                }

        }
}

sub substract{
                my ($fh1, $fh2)=@_;
                open (FH1, "<$fh1") or die ("Cannot open the file:$fh1!\n");
                open (FH2, "<$fh2") or die ("Cannot open the file:$fh2!\n");
                open my $out, ">", "Delta-SNP" or die "Cannot ioen the file:$!\n";
                my %hash;
                while(<FH1>){
                                chomp;
                                my ($chr, $pos, $rate1)=(split)[0,1,2];
                                $hash{$chr}{$pos}=$rate1;
                        }
                        while(<FH2>){
                                chomp;
                                my ($id, $loc, $rate2)=(split)[0,1,2];
                                my $sub = $hash{$id}{$loc} - $rate2;
                                        print $out "$id\t$loc\t$sub\n";
                                }
        }
        
main:{
                $startdir = cwd();
                if($output_directory =~ /^\//){
                                $workdir = $output_directory;
                        }elsif($output_directory ne ""){
                                        $workdir = "$startdir/$output_directory";
                                }else{
                                                $workdir = "$startdir";
                                        }
                if(! defined $cpu){
                                $cpu = 2; 
                                }
                if(! defined $freq){
                                $freq = 0.01; 
                                }
                if(! defined $depth){
                                $depth = 100; 
                                }
                if(! defined $increment){
                        $increment = 0.01;
                                }
                if(! defined $wsize){
                        $wsize = 2;
                }
                if(! defined $mq){
                                $mq = 1; 
                        }
                if(! defined $readsmin){
                                $readsmin = 5; 
                        }
                if(! defined $readsmax){
                                $readsmax = 100;
                }
                if(! defined $cov){
                                $cov = 0;
                }
                if(! defined $minf){
                                $minf = 0.3;
                }
                if(! defined $xmax){
                                $xmax = 50;
                }
        die "Error:The output directory isn't exists!" unless (-e "$workdir");


#bwa alignment
if(defined $f1_lef && defined $f2_lef ){
                &alignment($reference,$cpu,$f1_lef,$f1_rig,$workdir,$mq,$freq,$depth,$prefix1) unless (-e "$prefix1-snp-C.txt");
                &alignment($reference,$cpu,$f2_lef,$f2_rig,$workdir,$mq,$freq,$depth,$prefix2) unless (-e "$prefix2-snp-C.txt");
                my $SNPFile1="$prefix1-snp-C.txt";
                my $SNPFile2="$prefix2-snp-C.txt";

                # snp index statistics
                open my $OUT1, ">>", "bulked-SNP1" or die "Cannot open the file:$!\n";
                open my $OUT2, ">>", "bulked-SNP2" or die "Cannot open the file:$!\n";
                &FilteredSNP($SNPFile1,$readsmin,$readsmax,$minf,$SNPFile2,$OUT1,$OUT2,$cov);
                &FilteredSNP($SNPFile2,$readsmin,$readsmax,$minf,$SNPFile1,$OUT2,$OUT1,$cov);
                close $OUT1;
                close $OUT2;
                print "*************** SNP index statistics have done! ************************\n";
}elsif((defined $f1_lef && !defined $f2_lef) || (defined $f2_lef && !defined $f1_lef) || (defined $bulkedsnp1 && !defined $bulkedsnp2) || (!defin
ed $bulkedsnp1 && defined $bulkedsnp2)){
                if(defined $f1_lef || defined $bulkedsnp1){
                        open my $OUT1,">","bulked-SNP1" or die "Cannot open the file:$!\n";
                        if(defined $f1_lef){
                                &alignment($reference,$cpu,$f1_lef,$f1_rig,$workdir,$mq,$freq,$depth,$prefix1) unless (-e "$prefix1-snp-C.txt");
                                my $SNPFile1 = "$prefix1-snp-C.txt";
                                &FilteredSingle($SNPFile1,$readsmin,$readsmax,$minf,$OUT1);
                        }else{
                                &FilteredSingle($bulkedsnp1,$readsmin,$readsmax,$minf,$OUT1);
                        }
                        close $OUT1;
                        print "*************** SNP index statistics have done! ************************\n";
                }else{
                        open my $OUT2,">","bulked-SNP2" or die "Cannot open the file:$!\n";
                        if(defined $f2_lef){
                                &alignment($reference,$cpu,$f2_lef,$f2_rig,$workdir,$mq,$freq,$depth,$prefix2) unless (-e "$prefix2-snp-C.txt");
                                my $SNPFile2 = "$prefix2-snp-C.txt";    
                                &FilteredSingle($SNPFile2,$readsmin,$readsmax,$minf,$OUT2);
                        }else{
                                &FilteredSingle($bulkedsnp2,$readsmin,$readsmax,$minf,$OUT2);
                                }
                                close $OUT2;
                print "*************** SNP index statistics have done! ************************\n";
                }
}else{
                open my $OUT1,">>","bulked-SNP1" or die "Cannot open the file:$!\n";
                open my $OUT2,">>","bulked-SNP2" or die "Cannot open the file:$!\n";
            &FilteredSNP($bulkedsnp1,$readsmin,$readsmax,$minf,$bulkedsnp2,$OUT1,$OUT2,$cov);
                &FilteredSNP($bulkedsnp2,$readsmin,$readsmax,$minf,$bulkedsnp1,$OUT2,$OUT1,$cov);
                close $OUT1;
                close $OUT2;
                print "*************** SNP index statistics have done! ************************\n";
}
my $file3;
my $file4;
my $dir1;
my $dir2;
my $dir3;
# sorted bulked SNP
if(-e "bulked-SNP1"){
            my $file1 = "bulked-SNP1";
                open my $OUT3, ">", "bulked-SNP1.sorted" or die "Cannot open the file:$!\n";
                &SortedSNP($file1,$OUT3)unless(-s "bulked-SNP1.sorted");
                close $OUT3;
print "*************** sorted the bulked SNP have done! ************************\n";
                $file3 = "bulked-SNP1.sorted";
                my $cmd = "mkdir $prefix1";
                &process_cmd($cmd) unless (-d "$startdir/$prefix1");
                 $dir1 = "$startdir/$prefix1";
                &SplitSNP($file3,$dir1) unless (-e "$dir1/Chr1");
                print "*************** split the bulked SNP have done! *************************\n";
}
if(-e "bulked-SNP2"){
                my $file2 = "bulked-SNP2";
                open my $OUT4, ">", "bulked-SNP2.sorted" or die "Cannot open the file:$!\n";
                &SortedSNP($file2,$OUT4) unless (-s "bulked-SNP2.sorted");
                close $OUT4;
print "*************** sorted the bulked SNP have done! ************************\n";
                $file4 = "bulked-SNP2.sorted";
                my $cmd = "mkdir $prefix2";
                &process_cmd($cmd) unless (-d "$startdir/$prefix2");
                $dir2 = "$startdir/$prefix2";
                &SplitSNP($file4,$dir2) unless (-e "$dir2/Chr1");
                print "*************** split the bulked SNP have done! *************************\n";
}

# Get Delta(SNP-index)
if(-e "bulked-SNP1" && -e "bulked-SNP2"){
                &substract($file3,$file4) unless (-e "Delta-SNP" && -s "Delta-SNP");
                my $file5 = "Delta-SNP";
                my $cmd = "mkdir $prefix1-$prefix2";
                &process_cmd($cmd) unless (-d "$startdir/$prefix1-$prefix2");
                $dir3 = "$startdir/$prefix1-$prefix2";
                &SplitSNP($file5,$dir3) unless (-e "$startdir/$prefix1-$prefix2/Chr1");
                print "*************** split the bulked SNP have done! *************************\n";
}

#get snp density
my $seqio = Bio::SeqIO->new(-file => $reference, -format => 'fasta');
while(my $seq = $seqio->next_seq()){
        my $id = $seq->display_id;
        my $len = $seq->length;
        $SeqIndex{$id}=$len/1000000;
        }
                
&SnpDentity($dir1,$increment,$wsize,%SeqIndex) unless (!defined $prefix1 || -e "$dir1/Chr1.out");   
&SnpDentity($dir2,$increment,$wsize,%SeqIndex) unless (!defined $prefix2 || -e "$dir2/Chr2.out");
&SnpDentity($dir3,$increment,$wsize,%SeqIndex) unless (!defined $prefix1 || !defined $prefix2 || -e "$dir3/Chr1.out");

#Get coverage depth and confidence interval
#parameter: tit,$prefix.depth
sub SplitCoverage{
        my ($index,$file) = @_;
        my $d;
        if($file=~/(.+)\.depth/){$d = $1;}
                open (RH,"<$index") or die "Cannot open the file:$!\n";
                while(<RH>){
                                chomp;
                         my $cmd = "grep '\\b$_\\b' $file > $d/$_.depth";
                        &process_cmd($cmd);
                        }               
                close RH;
}
#Get the depth of SNP locus
#parameter:$dir
sub AverageCoverage{
                my $fh = shift;
                my %hash;
                opendir(IH, $fh) or die "Cannot open the directory:$!\n";
                my @dir = readdir IH;
                foreach my $file (@dir){
                next if($file =~/\./);
                open IN, "<","$fh/$file" or die;
                open PH, "<","$fh/$file.depth" or die;
                open my $OUT,">","$fh/$file.snpdepth" or die;
                while(<IN>){
                                chomp;
                                my ($chr,$pos,$rate)=(split)[0,1,2];
                                $hash{$chr}{$pos}=$rate;
                }
                while(<PH>){
                                my ($chr,$pos,$depth)=(split)[0,1,2];
                                $pos=$pos/1000000;
                                if(defined $hash{$chr}{$pos}){
                                                        print $OUT "$chr\t$pos\t$hash{$chr}{$pos}\t$depth\n";
                                }
                }
        }
        
close IN;
close PH;
close IH;
}

sub TwoAverage{
        my ($fh1,$fh2,$fh3)=@_;
  opendir (TH, $fh3) or die "Cannot open the directory:$!\n";
  my @dir = readdir TH;
        foreach my $file (@dir){
                next if($file=~/\./);
                open IN1,"<$fh1/$file.snpdepth" or die;
                open IN2,"<$fh2/$file.snpdepth" or die;
                open OUT, ">$fh3/$file.snpdepth" or die;
                my %hash;
                my $aver;
                while(<IN1>){
        chomp;
        my @arr = split/\s+/,$_;
        $hash{$arr[0]}{$arr[1]}=$arr[3];
                                }
                        while(<IN2>){
        chomp;
        my ($chr,$pos,$depth)=(split)[0,1,3];
        if(defined $hash{$chr}{$pos} ){
                $aver=($hash{$chr}{$pos}+$depth)/2;
              print OUT "$chr\t$pos\t$aver\n";
        }

        }
close IN1;
close IN2;
close OUT;
}
closedir TH;
}

sub VcfTwoAverage{
        my ($fh1,$fh2,$fh3)=@_;
  opendir (TH, $fh3) or die "Cannot open the directory:$!\n";
  my @dir = readdir TH;
        foreach my $file (@dir){
                next if($file=~/\./);
                open IN1,"<$fh1/$file" or die;
                open IN2,"<$fh2/$file" or die;
                open OUT, ">$fh3/$file.snpdepth" or die;
                my %hash;
                my $aver;
                while(<IN1>){
        chomp;
        my @arr = split/\s+/,$_;
        $hash{$arr[0]}{$arr[1]}=$arr[3];
                                }
                        while(<IN2>){
        chomp;
        my ($chr,$pos,$depth)=(split)[0,1,3];
        if(defined $hash{$chr}{$pos} ){
                $aver=($hash{$chr}{$pos}+$depth)/2;
              print OUT "$chr\t$pos\t$aver\n";
        }

        }
close IN1;
close IN2;
close OUT;
}
closedir TH;
}

#parameter:$dir,$file,$out2,$%SeqIndex
sub GetConfidenceInterval{
    my ($dir,$fh,$increment,$wsize,%hash)=@_;
                my %stat=();
                my %num=();
                my $key;
                if($fh =~ /(.+)\.snpdepth/){
                        $key = $1;
                }
                foreach (keys %hash){
                     next if($_ ne $key);
                        for(my $i=0;$i<=$hash{$_};$i+=$increment){
                                $stat{$_}{$i}=0;
                        }
                }
                open SP,"<","$dir/$fh" or die;
                while (<SP>) {
                                chomp;
                                my ($chr,$pos,$depth)=(split)[0,1,2];
                                foreach my $key (keys %{$stat{$chr}}){
                                                if($pos>=$key && $pos < $key+$wsize){
                                                        $stat{$chr}{$key}+=$depth;
                                                        $num{$chr}{$key}++;
                                                }
                                }
                }
                open OUT1, ">$dir/$fh.avercov" or die;  
                open OUT2, ">$dir/$fh.c01" or die;
                open OUT3, ">$dir/$fh.c005" or die;
                open OUT4, ">$dir/$fh.c001" or die;
                open IN, "<$startdir/ConfiInterval" or die;
                my @arr;
                my %hash01;
                my %hash001;
                my %hash005;
                while(<IN>){
                                chomp;
                                next if (/read_depth/);
                                @arr = split/\s+/,$_;
                                if($arr[3] ne '-'){
                                        $hash01{$arr[3]}=$arr[0];
                                }
                                if($arr[2] ne '-'){
                                        $hash005{$arr[2]}=$arr[0];
                                }
                                if($arr[1] ne '-'){
                                $hash001{$arr[1]}=$arr[0];
                                }
                }
                foreach my $p (keys %stat){
                        foreach my $k (sort {$a<=>$b} keys %{$stat{$p}}){
                                        if(defined $num{$p}{$k}){
                                                        $stat{$p}{$k}=int($stat{$p}{$k}/$num{$p}{$k}+0.5);
                                                        print OUT1 "$p\t$k\t$stat{$p}{$k}\n";
                                        }
                                        if($stat{$p}{$k} !=0 && defined $hash01{$stat{$p}{$k}} ){
                                                        print OUT2 "$p\t$k\t$hash01{$stat{$p}{$k}}\t-$hash01{$stat{$p}{$k}}\n";
                                        }
                                        if($stat{$p}{$k} !=0 && defined $hash005{$stat{$p}{$k}} ){
                                                        print OUT3 "$p\t$k\t$hash005{$stat{$p}{$k}}\t-$hash005{$stat{$p}{$k}}\n";
                                        }
                                        if($stat{$p}{$k} !=0 && defined $hash001{$stat{$p}{$k}} ){
                                                        print OUT4 "$p\t$k\t$hash001{$stat{$p}{$k}}\t-$hash001{$stat{$p}{$k}}\n";
                                        }
                        }
                }
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close IN;
close SP;
}
my $index="$startdir/tit";
my $fdepth1;
my $fdepth2;
if(defined $prefix1){
         $fdepth1 = "$prefix1.depth";
}
if(defined $prefix2){
         $fdepth2 = "$prefix2.depth";
}

#compute single bulked confidence Interval

#Compute paired-bulked confidence Interval
if(defined $f1_lef && defined $f2_lef){
        &SplitCoverage($index,$fdepth1);
        &SplitCoverage($index,$fdepth2);
        &AverageCoverage($dir1);
        &AverageCoverage($dir2);
        &TwoAverage($dir1,$dir2,$dir3);          # Get $dir3 *.snpdepth
        opendir (FH, $dir3 ) or die "Cannot open the directory:$!\n";
        my @dir = readdir FH;
        foreach my $file (@dir){
                        if($file =~/\.snpdepth\b/){
                        #open my $out1, ">", "$dir3/$file.avercov" or die "Cannot open the file:$!\n";
                        &GetConfidenceInterval($dir3,$file,$increment,$wsize,%SeqIndex);  # output the confidence Interval of 0.1,0.05,0.01     
                        
        }
}
closedir FH;
}

if(defined $bulkedsnp1 && defined $bulkedsnp2){
        &VcfTwoAverage($dir1,$dir2,$dir3);
        opendir (FH, $dir3 ) or die "Cannot open the directory:$!\n";
        my @dir = readdir FH;
        foreach my $file (@dir){
                        if($file =~/\.snpdepth\b/){
                        #open my $out1, ">", "$dir3/$file.avercov" or die "Cannot open the file:$!\n";
                        &GetConfidenceInterval($dir3,$file,$increment,$wsize,%SeqIndex);  # output the confidence Interval of 0.1,0.05,0.01     
                        
        }
}
closedir FH;
}

sub GetArr{
        my $dir = shift ;
        my @arr1=();
        my @arr2=();
        my @arr3=();
        my @arr4=();
        my @arr5=();
#       opendir (RH,$dir) or die "Cannot open the file:$!\n";
#       my @dir = readdir RH;
#        print "@dir\n";
        `ls $dir > dir`;
        open IN,"<tit" or die;
my @dir=();
       while(<IN>){
          chomp;
        push @dir,$_;
}
        for my $file (@dir){
                next if($file =~ /\./);
                        push @arr1, "$dir/$file";
                        push @arr2, "$dir/$file.out";
                        push @arr3, "$dir/$file.snpdepth.c01";
                        push @arr4, "$dir/$file.snpdepth.c005";
                        push @arr5, "$dir/$file.snpdepth.c001";
        }
        return (\@arr1,\@arr2,\@arr3,\@arr4,\@arr5);    
}
my @SNPFile1;
my @SNPFile2;
my @SNPFile3;
my @SNPDEN1;
my @SNPDEN2;
my @SNPDEN3;
my @CONF01;
my @CONF001;
my @CONF005;
if(defined $prefix1 ){
                my ($ref1,$ref2,$ref3,$ref4,$ref5)=GetArr($dir1);
                @SNPFile1 = @$ref1;
                @SNPDEN1 = @$ref2;
#       if (!prefix2){
#               @CONF101 = @ref3;
#               @CONF1005 = @ref4;
#               @CONF1001 = @ref5;
#       }
}
if(defined $prefix2 ){
                my ($ref1,$ref2)=GetArr($dir2);
                @SNPFile2 = @$ref1;
                @SNPDEN2 = @$ref2;
#       if(!prefix1){
#               @CONF201 = @ref3;
#               @CONF2005 = @ref4;
#               @CONF2001 = @ref5;
#  }
}
if(defined $dir3){
        my($ref1,$ref2,$ref3,$ref4,$ref5)=GetArr($dir3);
                @SNPFile3 = @$ref1;
                @SNPDEN3 = @$ref2;
                @CONF01 = @$ref3;
                @CONF005 = @$ref4;
                @CONF001 = @$ref5;
}
#plot R.1 graph  
my $R1_script;
my $R2_script;
my $R3_script;
if(defined $prefix1){
                $R1_script = "plot_SNP1.R";     
                open (my $RH, ">$R1_script") or die "Error, cannot write to $R1_script";
                print $RH "files = c(\"".join("\",\"",@SNPFile1)."\");\n";    ##?
                print $RH "dent = c(\"".join("\",\"",@SNPDEN1)."\");\n";
#               if(!prefix2){
#               print $RH "pv = c(\"".join("\",\"",@CONF101)."\");\n";
#               print $RH "pv005 = c(\"".join("\",\"",@CONF1005)."\");\n";
#               print $RH "pv001 = c(\"".join("\",\"",@CONF1001)."\");\n";
#               }
                print $RH "pdf(file=\"$dir1/$prefix1-bulked.pdf\",width=7,height=10)\n";
                print $RH "par(mfrow=c(6,2),bty='l',mar=c(4,4,3,2)+0.1,oma=c(0.1,0.1,0,0),mgp=c(2.5,1,0))\n";
                print $RH "#  png(file=\"my_plot.png\");\n";
                print $RH "for (i in 1:length(files))   {\n";
                print $RH "         sx=strsplit(files[i],split=\"/\");\n";
                print $RH "         x=c(sx[[1]][sapply(sx,length)]);\n";
                print $RH "                     data = read.table(files[i]);\n";
#               print $RH "                     ymin = min(data[,3]); ymax = max(data[,3]);\n";
#               print $RH "                     xmax = max(data[,2]);\n";
                print $RH "pdf(file=paste(\"$dir1/\", as.character(x),\".pdf\", sep=''));\n";
                print $RH "                     plot(data[,2],data[,3],pch=16,col='darkgreen',main=x,xlab='chromosome position(Mb)',ylab='SNP-ind
ex',font.axis=2,font.lab=2,xlim=c(0,$xmax),ylim=c(0,1),type='p',las=1,cex=0.5);\n";
                print $RH "                     aver = read.table(dent[i]);\n";
                print $RH "                     lines(aver[,2],aver[,3],col='red',lwd=2);\n";
                print $RH "                     abline(h=0.5,lty=2,col='black',lwd=1);\n";
                print $RH "dev.off()\n";
#               print $RH "                     xmax = max(data[,2]);\n";
                print $RH "                     plot(data[,2],data[,3],pch=16,col='darkgreen',main=x,xlab='chromosome position(Mb)',ylab='SNP-ind
ex',font.axis=2,font.lab=2,xlim=c(0,$xmax),ylim=c(0,1),type='p',las=1,cex=0.5);\n";
                print $RH "                     aver = read.table(dent[i]);\n";
                print $RH "                     lines(aver[,2],aver[,3],col='red',lwd=2);\n";
                print $RH "                     abline(h=0.5,lty=2,col='black',lwd=1);\n";
                print $RH "}\n";
                print $RH "dev.off()\n";
                
                close $RH;
}
                
#plot R.2 graph
if(defined $prefix2){
                $R2_script = "plot_SNP2.R";     
                open (my $SH, ">$R2_script") or die "Error, cannot write to $R2_script";
                print $SH "files = c(\"".join("\",\"",@SNPFile2)."\");\n";   
                print $SH "dent = c(\"".join("\",\"",@SNPDEN2)."\");\n";
#               if(!prefix2){
#               print $SH "pv = c(\"".join("\",\"",@CONF201)."\");\n";
#               print $SH "pv005 = c(\"".join("\",\"",@CONF2005)."\");\n";
#               print $SH "pv001 = c(\"".join("\",\"",@CONF2001)."\");\n";
#               }
                print $SH "pdf(file=\"$dir2/$prefix2-bulked.pdf\",width=7,height=10)\n";
                print $SH "par(mfrow=c(6,2),bty='l',mar=c(4,4,3,2)+0.1,oma=c(0.1,0.1,0,0),mgp=c(2.5,1,0))\n";
                print $SH "#  png(file=\"my_plot.png\");\n";
                print $SH "for (i in 1:length(files))   {\n";
                print $SH "         sx=strsplit(files[i],split=\"/\");\n";
                print $SH "         x=c(sx[[1]][sapply(sx,length)]);\n";
                print $SH "                     data = read.table(files[i]);\n";
#               print $SH "                     ymin = min(data[,2]); ymax = max(data[,2]);\n";
#               print $SH "                     xmin = min(data[,1]); xmax = max(data[,1]);\n";
                print $SH "pdf(file=paste(\"$dir2/\", as.character(x),\".pdf\", sep=''));\n";
                print $SH "                     plot(data[,2],data[,3],pch=16,col='darkgoldenrod1',main=x,xlab='chromosome position(Mb)',ylab='SN
P-index',font.axis=2,font.lab=2,xlim=c(0,$xmax),ylim=c(0,1),type='p',las=1,cex=0.5);\n";
                print $SH "                     aver = read.table(dent[i]);\n";
                print $SH "                     lines(aver[,2],aver[,3],col='red',lwd=2);\n";
                print $SH "                     abline(h=0.5,lty=2,col='black',lwd=1);\n";
                print $SH "dev.off()\n";
                print $SH "                     plot(data[,2],data[,3],pch=16,col='darkgoldenrod1',main=x,xlab='chromosome position(Mb)',ylab='SN
P-index',font.axis=2,font.lab=2,xlim=c(0,$xmax),ylim=c(0,1),type='p',las=1,cex=0.5);\n";
                print $SH "                     aver = read.table(dent[i]);\n";
                print $SH "                     lines(aver[,2],aver[,3],col='red',lwd=2);\n";
                print $SH "                     abline(h=0.5,lty=2,col='black',lwd=1);\n";
                print $SH "}\n";
                print $SH "dev.off()\n";
                
                close $SH;
}

#plot R.3 graph
if(defined $dir3){
                $R3_script = "plot_SNP3.R";     
                open (my $CH, ">$R3_script") or die "Error, cannot write to $R3_script";
                print $CH "files = c(\"".join("\",\"",@SNPFile3)."\");\n";    
                print $CH "dent = c(\"".join("\",\"",@SNPDEN3)."\");\n";
                print $CH "pv = c(\"".join("\",\"",@CONF01)."\");\n";
                print $CH "pv005 = c(\"".join("\",\"",@CONF005)."\");\n";
                print $CH "pv001 = c(\"".join("\",\"",@CONF001)."\");\n";
                print $CH "pdf(file=\"$dir3/$prefix1-$prefix2.pdf\",width=7,height=10)\n";
                print $CH "par(mfrow=c(6,2),bty='l',mar=c(4,4,3,2)+0.1,oma=c(0.1,0.1,0,0),mgp=c(2.5,1,0))\n";
                print $CH "#  png(file=\"my_plot.png\");\n";
                print $CH "for (i in 1:length(files))   {\n";
                print $CH "         sx=strsplit(files[i],split=\"/\");\n";
                print $CH "         x=c(sx[[1]][sapply(sx,length)]);\n";
                print $CH "                     data = read.table(files[i]);\n";
#               print $CH "                     ymin = min(data[,2]); ymax = max(data[,2]);\n";
#               print $CH "                     xmin = min(data[,1]); xmax = max(data[,1]);\n";
                print $CH "pdf(file=paste(\"$dir3/\", as.character(x),\".pdf\", sep=''));\n";
                print $CH "                     plot(data[,2],data[,3],pch=16,col='darkblue',main=x,xlab='chromosome position(Mb)',ylab=expressio
n(Delta(SNP-index)),font.axis=2,font.lab=2,xlim=c(0,$xmax),ylim=c(-1,1),type='p',las=1,cex=0.5);\n";
                print $CH "                     aver = read.table(dent[i]);\n";
                print $CH "                     lines(aver[,2],aver[,3],col='red',lwd=2);\n";
                print $CH "                     abline(h=0.5,lty=2,col='black',lwd=1);\n";
                print $CH "                     pv1=read.table(pv[i]);\n";
                print $CH "                     lines(pv1[,2],pv1[,3],col='darkgoldenrod4',lwd=1);\n";
                print $CH "                     lines(pv1[,2],pv1[,4],col='darkgoldenrod4',lwd=1);\n";
                print $CH "                     pv2=read.table(pv005[i]);\n";
                print $CH "                     lines(pv2[,2],pv2[,3],col='darkolivegreen',lwd=1);\n";
                print $CH "                     lines(pv2[,2],pv2[,4],col='darkolivegreen',lwd=1);\n";
                print $CH "                     pv3=read.table(pv001[i]);\n";
                print $CH "                     lines(pv3[,2],pv3[,3],col='pink',lwd=1);\n";
                print $CH "                     lines(pv3[,2],pv3[,4],col='pink',lwd=1);\n";
                print $CH "dev.off()\n";
                print $CH "                     plot(data[,2],data[,3],pch=16,col='darkblue',main=x,xlab='chromosome position(Mb)',ylab=expressio
n(Delta(SNP-index)),font.axis=2,font.lab=2,xlim=c(0,$xmax),ylim=c(-1,1),type='p',las=1,cex=0.5);\n";
                print $CH "                     aver = read.table(dent[i]);\n";
                print $CH "                     lines(aver[,2],aver[,3],col='red',lwd=2);\n";
                print $CH "                     abline(h=0.5,lty=2,col='black',lwd=1);\n";
                print $CH "                     pv1=read.table(pv[i]);\n";
                print $CH "                     lines(pv1[,2],pv1[,3],col='darkgoldenrod4',lwd=1);\n";
                print $CH "                     lines(pv1[,2],pv1[,4],col='darkgoldenrod4',lwd=1);\n";
                print $CH "                     pv2=read.table(pv005[i]);\n";
                print $CH "                     lines(pv2[,2],pv2[,3],col='darkolivegreen',lwd=1);\n";
                print $CH "                     lines(pv2[,2],pv2[,4],col='darkolivegreen',lwd=1);\n";
                print $CH "                     pv3=read.table(pv001[i]);\n";
                print $CH "                     lines(pv3[,2],pv3[,3],col='pink',lwd=1);\n";
                print $CH "                     lines(pv3[,2],pv3[,4],col='pink',lwd=1);\n";
                print $CH "}\n";
                print $CH "dev.off()\n";
                
                close $CH;
}

                &process_cmd("R --vanilla -q < $R1_script") unless (!defined $prefix1 || -e "$prefix1-bulked.pdf");
                &process_cmd("R --vanilla -q < $R2_script") unless (!defined $prefix2 || -e "$prefix2-bulked.pdf");
                &process_cmd("R --vanilla -q < $R3_script") unless (!defined $dir3 || -e "$prefix1-$prefix2.pdf");
                
                exit(0);
                
        }

