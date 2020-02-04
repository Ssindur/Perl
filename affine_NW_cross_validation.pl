$training_file=shift;
open(IN,$training_file);
@aa=<IN>;
chomp(@aa);

$data_length=scalar @aa;
print($data_length);

@go=(-16,-14,-12,-10,-8,-6,-4,-2);

chomp(@go);
$length_go=scalar @go;

%error={};
$error{$go[0]}=0;
for ($i=0; $i< $data_length; $i++)
{	
chomp($aa[$i]);

#opening the amino acids in the file
open(my $handle,$aa[$i]);
chomp(my @lines=<$handle>);
	$t_seq1=$lines[1];
	$t_seq2=$lines[3];
	close $handle;
#print("the sequence input from Affine_NW_cross_validation"."\n");
#print("The sequence 1 ".$t_seq1."\n");
#print("The sequence 2 ".$t_seq2."\n");

for($j=0; $j < $length_go; $j++)
  {
		$gap=$go[$j];
		print("The gap open value is = ".$go[$j]."\n");
#With blossom matrix	
	(qx(perl needleman_modified.pl $aa[$i] $gap > output_nw1.txt));
#Without blossom matrix
		(qx(perl needleman_modified_asgn2.pl $aa[$i] $gap > output_nw2.txt));
		
#opening of output text file into an array 		
		open(FILE,"<","output_nw.txt");
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


#with blossom matrix alignment
$accuracy_nw_affine1=qx(perl alignment_accuracy.pl $true_alignment_path output_nw1.txt);

#without blossom matrix alignment
$accuracy_nw_affine2=qx(perl alignment_accuracy.pl $true_alignment_path output_nw2.txt);

print("The Accuracy of sequence with blossom ".$true_alignment_path." is = ".$accuracy_nw_affine1."\n");
print("The Accuracy of sequence without blossom ".$true_alignment_path." is = ".$accuracy_nw_affine2."\n");

$error{$gap}= $error{$gap} + $accuracy_nw_affine1;
$error2{$gap}= $error2{$gap} + $accuracy_nw_affine2;
print("The cumulative  error with blossom for ".$true_alignment_path." is = ".$error{$gap}." and the gap is =  ".$gap."\n");
print("The cumulative  error without blossom for ".$true_alignment_path." is = ".$error2{$gap}." and the gap is =  ".$gap."\n");

}
	

}

$best_accuracy=-1;

#Best alignment accuracy calculation
#Average gap theory

for($i=0;$i< scalar(@go); $i++)
{	
	$g=$go[$i];
#	print("For gap = ".$g." Before Averaging, The error is = ".$error{$g}."\n");
	$error{$g} =$error{$g}/$data_length ;
	print("For gap = ".$g."The error is  =".$error{$g}."\n");

	if ($error{$g} > $best_accuracy)
	{
	$best_accuracy=$error{$g};
	$best_gap=$g;
	

	}
	

} 

print("The best alignment is as follows for NW with blossom matrix"."\n");

print("The gap open :- ".$best_gap."\n");
print("The best accuracy is :- ".$best_accuracy."\n");
print("The lowest error : ".(1-$best_accuracy)."\n");
$best_accuracy2=-1;

#Best alignment accuracy calculation
#Average gap theory

for($i=0;$i< scalar(@go); $i++)
{	
	$g=$go[$i];
#	print("For gap = ".$g." Before Averaging, The error is = ".$error{$g}."\n");
	$error2{$g} =$error2{$g}/$data_length ;
	print("For gap = ".$g."The error is  =".$error2{$g}."\n");

	if ($error2{$g} > $best_accuracy2)
	{
	$best_accuracy2=$error2{$g};
	$best_gap2=$g;
	

	}
	

} 

print("The best alignment is as follows for NW without  blossom matrix"."\n");

print("The gap open :- ".$best_gap2."\n");
print("The best accuracy is :- ".$best_accuracy2."\n");
print("The lowest error : ".(1-$best_accuracy2)."\n");

