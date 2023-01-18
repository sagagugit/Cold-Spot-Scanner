#!/usr/bin/perl

while ($ARGV[$a])           {
	open(IN, $ARGV[$a]);
	$pdb[$a]=$ARGV[$a];

	$nombre=rindex($ARGV[$a],p);
	$numero=$nombre-1;
	$preout=substr($pdb[$a],0,$numero);
	$pdbfile=$preout."_FH.pdb";

	`perl remediator.pl $pdb[$a] > a`;
	`reduce a > b`;
	`perl remediator.pl b -oldout > $pdbfile`;
	$a++;
}

`rm a b -rf`
 


