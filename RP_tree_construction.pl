#! usr/bin/perl
#written by xiaoyuzhu 2021-10-27

##step1: copy ribosomal protein sequences to a new folder from the orthofinder folder based on eggnog-mapper result
##step2: count the number of genomes each genes exist in and select the proper genes for tree buliding
##step3: retain only one sequence of multi-copy genes
##step4: conduct multi sequences alignment for each ribosomal protein (MAFFT)
##step5: trim aligned sequences (TrimAl)
##step6: Concatenate all protein sequences for each genome
##step7: construct maximum likelihood (ML) phylogenetic tree using IQtree or Fasttree 

#Info of 16 bacterial ribosomal proteins
# COG_class	Description
# COG00090	Ribosomal protein L2
# COG00087	Ribosomal protein L3
# COG00088	Ribosomal protein L4
# COG00094	Ribosomal protein L5
# COG00097	Ribosomal protein L6P/L9E
# COG00093	Ribosomal protein L14
# COG00200	Ribosomal protein L15
# COG00197	Ribosomal protein L16/L10E
# COG00256	Ribosomal protein L18
# COG00091	Ribosomal protein L22
# COG00198	Ribosomal protein L24
# COG00092	Ribosomal protein S3
# COG00096	Ribosomal protein S8
# COG00051	Ribosomal protein S10
# COG00186	Ribosomal protein S17
# COG00185	Ribosomal protein S19


#Info of 56 archaeal ribosomal proteins
# arCOG_class   COG	 	Description
# arCOG04289	COG0081	Ribosomal protein L1
# arCOG04288	COG0244	Ribosomal protein L10
# arCOG04372	COG0080	Ribosomal protein L11
# arCOG04287	COG2058	Ribosomal protein L12E/L44/L45/RPP1/RPP2
# arCOG04242	COG0102	Ribosomal protein L13
# arCOG04095	COG0093	Ribosomal protein L14
# arCOG04167	COG0093	Ribosomal protein L14
# arCOG00779	COG0200	Ribosomal protein L15
# arCOG04209	COG1632	Ribosomal protein L15E
# arCOG04113	COG0197	Ribosomal protein L16/L10AE
# arCOG04088	COG0256	Ribosomal protein L18
# arCOG00780	COG1727	Ribosomal protein L18E
# arCOG04089	COG2147	Ribosomal protein L19E
# arCOG04067	COG0090	Ribosomal protein L2
# arCOG04175	COG2157	Ribosomal protein L20A (L18A)
# arCOG04129	COG2139	Ribosomal protein L21E
# arCOG04098	COG0091	Ribosomal protein L22
# arCOG04072	COG0089	Ribosomal protein L23
# arCOG04094	COG0198	Ribosomal protein L24
# arCOG01950	COG2075	Ribosomal protein L24E
# arCOG00785	COG0255	Ribosomal protein L29
# arCOG04070	COG0087	Ribosomal protein L3
# arCOG04086	COG1841	Ribosomal protein L30/L7E
# arCOG01752	COG1841	Ribosomal protein L30/L7E
# arCOG04473	COG2097	Ribosomal protein L31E
# arCOG00781	COG1717	Ribosomal protein L32E
# arCOG04208	COG1997	Ribosomal protein L37AE/L43A
# arCOG04071	COG0088	Ribosomal protein L4
# arCOG04109	COG1631	Ribosomal protein L44E
# arCOG04092	COG0094	Ribosomal protein L5
# arCOG04090	COG0097	Ribosomal protein L6P/L9E
# arCOG01751	COG1358	Ribosomal protein L7Ae or related RNA K-turn-binding protein
# arCOG01758	COG0051	Ribosomal protein S10
# arCOG04240	COG0100	Ribosomal protein S11
# arCOG04255	COG0048	Ribosomal protein S12
# arCOG01722	COG0099	Ribosomal protein S13
# arCOG04185	COG0184	Ribosomal protein S15P/S13E
# arCOG04096	COG0186	Ribosomal protein S17
# arCOG01885	COG1383	Ribosomal protein S17E
# arCOG04099	COG0185	Ribosomal protein S19
# arCOG01344	COG2238	Ribosomal protein S19E (S16A)
# arCOG04245	COG0052	Ribosomal protein S2
# arCOG04182	COG2004	Ribosomal protein S24E
# arCOG04183	COG1998	Ribosomal protein S27AE
# arCOG04108	COG2051	Ribosomal protein S27E
# arCOG04314	COG2053	Ribosomal protein S28E/S33
# arCOG04097	COG0092	Ribosomal protein S3
# arCOG04186	COG1890	Ribosomal protein S3AE
# arCOG04239	COG0522	Ribosomal protein S4 or related protein
# arCOG04093	COG1471	Ribosomal protein S4E
# arCOG04087	COG0098	Ribosomal protein S5
# arCOG01946	COG2125	Ribosomal protein S6E (S10)
# arCOG04254	COG0049	Ribosomal protein S7
# arCOG04091	COG0096	Ribosomal protein S8
# arCOG04154	COG2007	Ribosomal protein S8E
# arCOG04243	COG0103	Ribosomal protein S9

use POSIX qw(strftime); #time
use Getopt::Long;
use Term::ANSIColor qw(:constants); #color
use strict;

my ($proteins, $outdir, $tree, $annotation, $min, $lineage, $help);

GetOptions
(
	"p|proteins=s" => \$proteins,   
	"a|annotation=s" => \$annotation,        # string
	"o|outdir=s" => \$outdir,
	"t|tree=s" => \$tree,
	"m|min=s" => \$min,
	"s|lineage=s" => \$lineage,
	"h|help" => \$help                 # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Author: xiaoyuzhu 2021-10-27
Github: 

Function: 
This script is aimed to construct ribosomal protein (RP) tree based on 16 bacterial RPs [1] or 56 archaeal RPs [2]. 
[1] Hug, L.A., et al., A new view of the tree of life. Nat Microbiol, 2016. 1: p. 16048.
[2] Seitz, K.W., Dombrowski, N., Eme, L. et al. Asgard archaea capable of anaerobic hydrocarbon cycling. Nat Commun 10, 1822 (2019).

Dependency: MAFTT, TrimAl, Fastree, IQtree and seqtk. Please ensure they have been installed ahead.

Example: perl $0 -p Orthogroup_Sequences -o RP_tree -t fasttree -m 0.4 -a emapper.annotations
Options:
	-p -proteins <dir>           required, the path to the Orthogroup_Sequences folder generated by orthofinder
	-a -annotation <file>		 required, the path to the result of eggnog annotation (e.g. eggnog.annotations)
	-o -outdir <dir>             required, the output directory
	-t -tree <tree menthod>      required, two options: "iqtree" or "fasttree"
	-m -min <num>                minmum percentage of genome containing each RP, default: 0.4. Range:[0-1]
	-s -lineage <lineage>	     required, two options: "arch" or "bact"
	-h -help                     print help information

Requirement of input files:
1. For protein sequences, the sequence names must be composed of genome name and gene id with "__" as separator:
   >MAG6666__gene00001
   >MAG6666__gene00002
   To meet this requirement, you had better add the genome names in front of each gene before conduct the orthofinder pipline.
   You can simply conduct the following command in your directory of faa files:
   "for filename in *.faa; do base=\$(basename \$filename .faa); echo \$base; sed -i "s/^>/>\${base}__/g" \$filename; done"
2. The query name in eggnog.annotations file (first line) must be consistent with the name of each Orthogroup Sequences files, like OG0000001.

Note: 
When you apply iqtree to build tree, it may inform you:
"OMP: Info #270: omp_set_nested routine deprecated, please use omp_set_max_active_levels instead."
This message can be ignored.

Output files: 
01.RP_tree(.nwk/.treefile)         tree file 
02.Concatenated.faa                concatenated RP sequences for constructing tree
03.RP_number_in_each_genome        distribution of RP in each genome
xxxoriginal.faa                    all RPs from each genome
xxx(single_copy/mafft/trimal).faa  intermediate files of tree building

USAGE
############################## usage end #######################################

if ($help || !(defined $proteins) || !(defined $outdir) || !(defined $tree)){
	print BRIGHT_YELLOW, $usage."\n", RESET;
	exit;
}
if (!(defined $min)){
	$min = 0.4;
}

my @RP_genes;
if ($lineage eq "arch") {
	@RP_genes = qw/arCOG04086 arCOG04255 arCOG04254 arCOG01758 arCOG04245 arCOG04372 arCOG04289 arCOG04070
 arCOG04071 arCOG04072 arCOG04067 arCOG04098 arCOG04097 arCOG04095 arCOG04167 arCOG04092 arCOG04091 arCOG04090 
 arCOG04087 arCOG01722 arCOG04240 arCOG04242 arCOG04243 arCOG04185 arCOG04099 arCOG04096 arCOG04113 arCOG04094 
 arCOG00779 arCOG04288 arCOG00785 arCOG04088 arCOG04239 arCOG01751 arCOG01885 arCOG04093 arCOG04109 arCOG04209 
 arCOG00781 arCOG00780 arCOG01752 arCOG04186 arCOG04208 arCOG04183 arCOG04182 arCOG04154 arCOG04108 arCOG04314 
 arCOG04287 arCOG01950 arCOG04473 arCOG01946 arCOG04129 arCOG04089 arCOG04175 arCOG01344/;
}elsif ($lineage eq "bact") {
	@RP_genes = qw/COG0090 COG0087 COG0088 COG0094 COG0097 COG0093 COG0200 COG0197 COG0256 COG0091 
	COG0198 COG0092 COG0096 COG0051 COG0186 COG0185/;
}else {
	print "Parameter -s error: Please specify the lineage, \"bact\" or \"arch\"";
	exit;
}

my $start_time=time;
print strftime("\nStart time is %Y-%m-%d %H:%M:%S\n", localtime(time));

#grasp genes from genomes
system "mkdir $outdir";
print "starting to identify RPs and cp them to the $outdir ...\n";
foreach my $genes (@RP_genes) {
	system "for i in `grep \"$genes\" $annotation | cut -f 1`; do cat $proteins/\$i.fa >> $outdir/$genes\_original.faa;done";
}

my @dir_files = <$outdir/*>;
if ( @dir_files ) {
    print scalar @dir_files." kind of RPs have been found in all genomes.\n\n";    
}else{
    print "Didn't find any RPs. Please check whether the query names in eggnog.annotations file (first line) are consistent with the names of each Orthogroup Sequences files.\n";
    exit;
}

#generate the genome array
chomp (my $genomes= `cat $proteins/*fa | grep \">\" | sed 's/>//g'| sed 's/__.*//g' | sort | uniq`);
my @genomes = split/\n/, $genomes;

#sum the number of ribosomal proteins in each genome
print "starting to select proper genes for tree building ...\n";
open INFO, '>>', "$outdir/03.RP_number_in_each_genome.txt" || die "03.RP_number_in_each_genome.txt\n";
my @IDs;
push @IDs, "RP_Genes"."\t";
foreach (@genomes) {
	push @IDs, $_."\t";
}
print INFO (@IDs, "\n");
print "##################################################################################################################\n";
print "The following genes presenting in at least ". $min*100 ."% of genomes are selected to construct ML tree:\n";
my @selected_genes;
foreach my $genes (@RP_genes) {
	my @content;
	my $genome_count = 0;
	push @content, $genes . "\t";
	foreach my $genomes (@genomes){
		chomp(my $num= `grep -s \"$genomes\" $outdir/$genes\_original.faa | wc -l`);
		push @content, $num . "\t";
		if ($num > 0) {
			$genome_count ++;
		}
	}
	my $percentage = $genome_count / scalar @genomes;
	if ($percentage >= $min ) {
		push @selected_genes, $genes;
		print "$genes\n";
	}else{
		system "rm -rf $outdir/$genes\_original.faa";
	}
	print INFO (@content, "\n");
}
print "\n".scalar @selected_genes." of total ". @RP_genes ." RP genes were selected to construct tree\n";
print "##################################################################################################################\n\n";

#identify multi-copy genes and retain only one sequence
print "starting to prepare the single-copy genes files ...\n\n";
foreach my $genes (@selected_genes) {
	my @arr;
	open OUTPUT, '>>', "$outdir/$genes\_temp.faa" || die "Cannot creat the file: $outdir/$genes\_temp.faa\n";
	open RP, '<', "$outdir/$genes\_original.faa" || die "Cannot open the file: $genes\_original.faa\n";
	$/ = ">";
	while (<RP>) {
		chomp;
		my $gene_list = (split/__/,(split/\n/,$_)[0])[0];  ##基因组名按两个下划线分隔
		if ((grep {$_ eq $gene_list} @arr) || !(defined $gene_list)){
			next;
		}else{
			print OUTPUT ">" . $_;
			push @arr, $gene_list;
		}	
	}
}

#add empty sequences names to single-copy files
system "touch COGxxxx_temp.faa";  #不知道为什么最后一个temp.faa总是无法处理，所以创了一个临时的空文件，等循环跑完再删除
foreach my $genes (@selected_genes) {
	my %hash;
	open TEMP, '<', "$outdir/$genes\_temp.faa";
	$/ = ">";
	while (<TEMP>) {
		chomp;	
		my $gene_list = (split/__/,(split/\n/,$_)[0])[0]; 
		$hash{$gene_list} = $_;
	}
	open OUT, '>>', "$outdir/$genes\_single_copy.faa";
	foreach my $genome (@genomes) {
		if (exists $hash{$genome}) {      
			print OUT ">" . $hash{$genome};	
		}else{
			print OUT ">"."$genome\__ribosomal_protein_$genes\n";
		}
	}
}
system "rm -rf $outdir/*_temp.faa $outdir/COGxxxx*";


##align and trim
print "starting to align sequences ...\n\n";
foreach my $genes (@selected_genes) {
!system "mafft --quiet --maxiterate 1000 --localpair  $outdir/$genes\_single_copy.faa > $outdir/$genes\_mafft.faa" || die "mafft error";
}
print "starting to trim sequences ...\n\n";
foreach my $genes (@selected_genes) {
!system "trimal -in $outdir/$genes\_mafft.faa -out $outdir/$genes\_temp1.faa -automated1 -keepseqs" || die "trimal error"; #-automated1 is suitable for ML tree
system "seqtk seq -l 0 $outdir/$genes\_temp1.faa >> $outdir/$genes\_trimal.faa && rm -rf $outdir/*_temp1.faa";
}


##Concatenate sequences
print "starting to concatenate sequences ...\n\n";
my %newhash;
foreach my $genes (@selected_genes) {
	open FINAL, '<', "$outdir/$genes\_trimal.faa" || die "Cannot open the file: $outdir/$genes\_trimal.faa\n";
	$/ = ">";
	while (<FINAL>) {
		chomp;	
		my $gene_list = (split/__/,(split/\n/,$_)[0])[0]; 
		my $sequence = (split/\n/,$_)[1];
		$newhash{$gene_list}.= $sequence;
	}
}
open CONCAT, '>>', "$outdir/02.Concatenated.faa" || die "Cannot creat the file: 02.Concatenated.faa\n";
foreach (@genomes) {
	print CONCAT ">".$_."\n";
	print CONCAT $newhash{$_}."\n";

}


##Construct trees using Fastree or IQtree
if ($tree eq "iqtree") {
	print "starting to constucting RP tree with iqtree ...\n";
	system "cd $outdir && iqtree -s 02.Concatenated.faa -alrt 1000 -bb 1000 --quiet --prefix 01.RP_iqtree";
}elsif ($tree eq "fasttree") {
	print "\nstarting to constucting RP tree with fasttree ...\nFasttree log:\n";
	system "FastTreeMP -gamma $outdir/02.Concatenated.faa > $outdir/01.RP_fasttree.nwk";
}


##Record the program running time!
my $duration_time=time-$start_time;
print strftime("\n\nEnd time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Running time: $duration_time seconds\.\n";
print "Congratulations!!!\n\n"
