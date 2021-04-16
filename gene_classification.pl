#!/data2/wangb/tools/perl
# usage
# perl gene_classification.pl module_color_gene.list;


`rm -r Gene_Moudle; mkdir Gene_Moudle`;
open IN,$ARGV[0]; # module_color_gene.list;
$line = <IN>;
while ($line = <IN>) {
	chomp($line);
	@arr = split("\t",$line);
	$arr[2] =~ s/[^a-zA-Z0-9]//g;
	$group = "$arr[1]-$arr[2]";
	$arr[0] =~ s/\|.*//g;
	push @{$hash{$group}}, $arr[0];
}
foreach $key (keys %hash)  {
	$num = @{$hash{$key}};
	open OUT,">>Gene_Moudle/$key-$num.list";
	foreach $item(@{$hash{$key}}){
		print OUT "$item\n";
	}
	close(OUT);
}

