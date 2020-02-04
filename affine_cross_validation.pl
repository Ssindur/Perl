$training_file=shift;
open(IN,$training_file);
@aa=<IN>;
chomp(@aa);

$data_length=scalar @aa;
#print($data_length);
@ge=(-2,-1,-.5,-.1);
@go=(-16,-14,-12,-10,-8,-6,-4,-2);
chomp(@ge);
chomp(@go);
$length_ge=scalar @ge;
$length_go=scalar @go;
$best_accuracy=-1;
%error={};
$error{$go[0]}{$ge[0]}=0;
$s=0;
$n=0;

for ($i=0; $i< $data_length; $i++)
{	
chomp($aa[$i]);
#opening the amino acids in the file
open(my $handle,$aa[$i]);
chomp(my @lines=<$handle>);
	$t_seq1=$lines[1];
	$t_seq2=$lines[3];
	close $handle;


for($j=0; $j < $length_go; $j++)
  {
	for ($k=0; $k < $length_ge; $k++)
		{	
		$gap=$go[$j];
		$gap_extension=$ge[$k];
		print("The gap_extension value = ".$ge[$k]." The gap open value is = ".$go[$j]."\n");
		#system("perl affine_modified.pl $aa[$i] $gap $gap_extension");
#with blossom matrix
#		(qx(perl affine_modified.pl $aa[$i] $gap $gap_extension  > output.txt));
#withut bl0ssom matrix		
(qx(perl affine_modified_asgn2.pl $aa[$i] $gap $gap_extension  > output1.txt));
		
#opening of output text file into an array 		
		open(FILE,"<","output.txt");
		@output_data=<FILE>;
		chomp(@output_data);
		close(FILE);

$len=scalar @output_data;

#printing the NW_affine alignment output

#for ($a=0; $a<$len; $a++)
#	{ 		
#	 print($output_data[$a]."\n");		
#	}


#for calculating same seq from prefab file
$str=$aa[$i];
$index=index($str,'.unaligned');
$true_alignment_path=substr($str,0,$index);

#print($true_alignment_path."\n");
#with blossom aligned matrix
#$accuracy_nw_affine=qx(perl alignment_accuracy.pl  $true_alignment_path output.txt);
#without blossom aligned matrix
$accuracy_nw_affine=qx(perl alignment_accuracy.pl  $true_alignment_path output1.txt );

print("The Accuracy of sequence ".$true_alignment_path." is = ".$accuracy_nw_affine."\n");

$error{$go[$j]}{$ge[$k]} = $error{$go[$j]}{$ge[$k]} + $accuracy_nw_affine;
print("the error for each line ".$error{$go[$j]}{$ge[$k]}."\n");
#Best alignment accuracy calculation

}

}
}
$max=-1;
for($i=0; $i < $length_go; $i++)
{
	for ($j=0; $j < $length_ge; $j++)
	{
		$gap=$go[$i];
		$gap_ext=$ge[$j];
#	print(" the error before averaging = ".$error{$gap}{$gap_ext}."\n");
	$error{$gap}{$gap_ext}= $error{$gap}{$gap_ext}/$data_length;
#	print("For gap open = ".$go[$i]." and for gap extension = ".$ge[$j]." the error is = ".$error{$gap}{$gap_ext}."\n");
	if($error{$gap}{$gap_ext} > $max)
	{
	 	$max=$error{$gap}{$gap_ext};
		$best_gap=$gap;
		$best_gap_ext=$gap_ext;

	
	}	

}

} 
#	print(" the alignement accuracy without  blossom matrix "."\n");
	print(" the alignement accuracy with blossom matrix "."\n");
	print(" the BEST GAP = ".$best_gap."\n");
	print(" the BEST GAP EXTENSION  = ".$best_gap_ext."\n");
	print(" the ERROR = ".$error{$best_gap}{$best_gap_ext}."\n");
	print(" the lowest ERROR = ".(1-$max)."\n");

