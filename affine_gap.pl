$file=shift;
open(FILE,$file);
@data = <FILE>;

chomp(@data);

$seq1 = $data[1];
$seq2 = $data[3];
$aligned_seq1="";
$aligned_seq2="";
chomp($seq1);
chomp($seq2);
@score_matrix;#V matrix
$score_matrix[0][0] = 0;
@M;@E;@F;
$M[0][0] = 0;
$E[0][0] = $F[0][0] = -Inf;
$match = 5;
$mismatch = -4;
$gap = -10;
$gap_extension = -.1;
@T;
$T[0][0]=0;

##INITIALIZATION
for (my $i=1; $i <= length($seq1); $i++)
{
$score_matrix[$i][0]  = $gap + ($i)*$gap_extension;
$E[$i][0] = -Inf;
$F[$i][0] = $gap + ($i)*$gap_extension;
$M[$i][0] = -Inf;
$T[$i][0] = "U";
}

for (my $i=1; $i <= length($seq2); $i++)
{
$score_matrix[0][$i] =  $gap + ($i)*$gap_extension;
$E[0][$i] = $gap + ($i)*$gap_extension;
$F[0][$i] = -Inf;
$M[0][$i] = -Inf;
$T[0][$i] = "L";
}

##Recurrence
for (my $i=1; $i<=length($seq1); $i++)
{
for (my $j=1; $j<=length($seq2); $j++)
{
 	$letter_seq1 = substr($seq1,$i-1,1);
	$letter_seq2 = substr($seq2,$j-1,1);
#	print("Seq 1".$letter_seq1."\n");
#	print("Seq 2".$letter_seq2."\n");

#calculatin the M matrix
	if($letter_seq1 eq $letter_seq2)
	{
	   $M[$i][$j] = $score_matrix[$i-1][$j-1] + $match;
	    
	}
	else			
	{
	
	  $M[$i][$j] = $score_matrix[$i-1][$j-1] + $mismatch;
		
	}

#calculating the F matrix value 																																									
	if($M[$i-1][$j]+$gap >= $F[$i-1][$j]+$gap_extension)
	{
		$F[$i][$j] = $M[$i-1][$j] + $gap;
	}
	else
	{
		$F[$i][$j] = $F[$i-1][$j] + $gap_extension;
	}
#calculating the E matrix
	if($M[$i][$j-1]+$gap >= $E[$i][$j-1]+ $gap_extension)
	{
		$E[$i][$j] = $M[$i][$j-1] +$gap;
	}
	else
	{
 
		$E[$i][$j] = $E[$i][$j-1] + $gap_extension;
	}
		
#print("M = ".$M[$i][$j]."\n");
#print("E = ".$E[$i][$j]."\n");
#print("F = ".$F[$i][$j]."\n");

##Maximum value selection
if($M[$i][$j] >= $E[$i][$j] && $M[$i][$j] >= $F[$i][$j])
	{
		$score_matrix[$i][$j] = $M[$i][$j];
		$T[$i][$j] = "D";
	}
elsif($F[$i][$j] >= $M[$i][$j] && $F[$i][$j] >= $E[$i][$j])
	{
		$score_matrix[$i][$j] = $F[$i][$j];
		$T[$i][$j] = "U";
	}
else 
	{
		$score_matrix[$i][$j] = $E[$i][$j];
		$T[$i][$j] = "L";
	}

}
}

print("The sequence 1 id :- ".$seq1."\n");
print("The sequence 2 id :- ".$seq2."\n");

for ($i=0; $i<=length($seq1); $i++)
{
for ($j=0; $j<=length($seq2); $j++)
{ 
print($score_matrix[$i][$j]."  ");
}
print("\n");
}

print("The Traceback matirx"."\n");
for ($i=0; $i<=length($seq1); $i++)
{
for ($j=0; $j<=length($seq2); $j++)
{
print($T[$i][$j]);
}
print("\n");
}
##Traceback
$aligned_seq1="";
$aligned_seq2="";
$i=length($seq1);
$j=length($seq2);
while ($i != 0 || $j != 0)
{
 if ($T[$i][$j] eq "L")
	{
	 $aligned_seq1="-".$aligned_seq1;
	 $aligned_seq2=substr($seq2,$j-1,1).$aligned_seq2;
	 $j =$j-1;
	}
elsif ($T[$i][$j] eq "U")
 {
   $aligned_seq2="-".$aligned_seq2;
   $aligned_seq1=substr($seq1,$i-1,1).$aligned_seq1;
   $i=$i-1;
  }
else
 {
  $aligned_seq1=substr($seq1,$i-1,1).$aligned_seq1;
  $aligned_seq2=substr($seq2,$j-1,1).$aligned_seq2;
  $i=$i-1;
  $j=$j-1;
}
}
print("The aligned sequence 1 is ".$aligned_seq1."\n");
print("The aligned sequence 2 is ".$aligned_seq2."\n");

