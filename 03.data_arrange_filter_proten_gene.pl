#!/data2/wangb/tools/perl
#!/usr/bin/perl

####################################################
# $ARGV[0] - 20160128-LUSC-RNAseqGene.txt
# $ARGV[1] - 20160128-LUSC-Clinical.txt
# $ARGV[2] - filter.list
# example
# perl 03.data_arrange_filter_proten_gene.pl 20160128-LUSC-RNAseqGene.txt 20160128-LUSC-Clinical.txt filter.list
####################################################

# open raw file
open IN, "$ARGV[0]";
open IN2, "$ARGV[1]";
open IN3, "$ARGV[2]";

# create output file
open OUT, ">01.inputdata.txt";
open OUT2, ">01.coldata.txt";
print OUT2 "sample_id\tcondition\n";
open OUT3, ">01.clinical.txt";

# load filter.list file 
while ($line = <IN3>) {
	chomp($line);
	$protein_gene{$line} = 1;
}
close(IN3);

$line = <IN>;
chomp($line);
$line =~  tr/a-z-/A-Z_/;
@RNASeqSample = split("\t",$line);
shift @RNASeqSample;
for($i=0;$i<@RNASeqSample;$i++){
	$RNASeqSample[$i] =~ s/(......_....)_..._..._...._../$1/g;
}
close(IN);

$line = <IN2>;
chomp($line);
$line =~  tr/a-z-/A-Z_/;
$clinical_head = $line ;
@ClinicalSample = split("\t",$line);
shift @ClinicalSample;
@mergeSample = (@ClinicalSample, @RNASeqSample);
%unionSample = {};
foreach $sample_id(@mergeSample){
	if(exists($sampleHash{$sample_id})){
		$unionSample{$sample_id} = 1;
	}else{
		$sampleHash{$sample_id} = 1;
	}
}



############################################################
###  prepare
############################################################
open IN, "$ARGV[0]";
$line = <IN>;
chomp($line);
@arr = split("\t",$line);
for($i=0;$i<@arr;$i++){
	$arr[$i] =~ tr/a-z-/A-Z_/;
	if($i == 0){
		push @list, $i;
		push @head, '';
		next;
	}
	$arr[$i] =~ /(.*?_.*?_.*?)_..._.*?_.*?_.*?/;
	if(!exists($hash{$1})){
		if($arr[$i] =~ /(.*?_.*?_.*?)_..B_.*?_.*?_.*?/){
			$unionSample{$1} = 0;
			next;
		}
		else{
			if(exists($unionSample{$1})){
				if($unionSample{$1} == 1){
					push @list, $i;
					$hash{$1} = 1;
					push @head, $arr[$i];
					push @head2, $1;			
					next;
				}
			}
		}
	}
}

############################################################
###  01.colData.txt
############################################################
$sample_num = 0;
for($i=1;$i<@head;$i++)
{
	$head[$i] =~ /(.*?_.*?_.*?)_(..)._.*?_.*?_.*/;
	if($2 < 10){
		$sample_num ++;
		print OUT2 $1."\tT\n";
	}
	elsif($2 >= 10){
		$sample_num ++;
		print OUT2 $1."\tN\n";
	}
	else{
		print $1."\t MOVE OUT \n";
	}
}
print "01.colData.txt \t $sample_num samples\n";



############################################################
###  01.inputdata.txt
############################################################
$sample_num = @head2;
print "01.inputdata.txt \t $sample_num samples\n";

$head_line = join("\t",@head2);
print OUT "\t$head_line\n";
$line = <IN>;
while($line = <IN>){
	chomp($line);
	@arr = split("\t",$line);
	@tmp = split("\\|",$arr[0]);
	if(exists($protein_gene{$tmp[0]}))
	{
		for($i=0;$i<@list;$i++){
			push @row, $arr[$list[$i]];
		}
		print OUT join("\t",@row)."\n";
		@row =();
	}
}


############################################################
###  01.clinical.txt
############################################################
@list = ();
push @list, 0;
for($i=0;$i<@ClinicalSample;$i++){
	if(exists($unionSample{$ClinicalSample[$i]})){
		if($unionSample{$ClinicalSample[$i]} == 1){
			push @list, $i+1;
		}
	}
}
$sample_num = @list - 1;
print "01.clinical.txt \t $sample_num samples\n";

@arr =split("\t",$clinical_head);
for($i=0;$i<@list;$i++){
	push @row, $arr[$list[$i]];
}
$row[0] =~s/\s/_/g;
print OUT3 join("\t",@row)."\n";
@row =();

while($line = <IN2>)
{
	chomp($line);
	@arr = split("\t",$line);
	for($i=0;$i<@list;$i++){
		push @row, $arr[$list[$i]];
	}
	print OUT3 join("\t",@row)."\n";
	@row =();
}