#! /usr/bin/perl
# Simon K. Schaefer
# Version 1.2

use strict;
use warnings;

#global variabels for handling input and output files
my $input_file;
my $output_file;
my $header = ">Generic Header Placeholder";


#global variables for calculation and cleanup process
my @file;
my @fasta;
my $dir;


#Primer counter
my $IgG1_2 = 0;
my $IgG3 = 0;
my $IgM = 0;
my $IgA = 0;

#calculation variables
my $RID;
my $FID;
my $RID_FID;
my $UMI;
my $Iso;
my $primary_key;
my $all;
my $i_counter;
my $corrected;
my $loopkey;

my $one;
my $two;
my $three;	


my %FID_hash;
my %RID_hash;
my %FID_RID_hash;
my %FID_hash_correted;
my %RID_hash_correted;
my %FID_RID_hash_correted;
my %sorting_hash;
my @sorted_hash;

#input via commandline
$input_file = $ARGV[0];

#open file and test for filetype
if ($input_file =~ m/.*\.fasta$/i){
open(READ,"$input_file") || die "error:cant open file";
	@file = <READ>;
	close(READ);
}else{
	die "error:wrong file type";
}


#create output files from input file
$output_file = substr($input_file,0,int((length($input_file)))-6)."_umified_norm.fasta";

#create output folder
$dir = substr($input_file,0,int((length($input_file)))-6)."_umified";
mkdir $dir;

#open output file
open (FILE1, "> $dir/log_$output_file") || die "problem opening $output_file\n";

#iterates over all lines in the file and creates fasta format with header and content
foreach (@file){
##############
#all
##############
if($_ !~ m/^>.*/i){
$all++;
}

$Iso = "NA";

#Primer for IgA:
if ($_ =~ m/^[ACT]{5}ACA[ACT]{5}ACA[ACT]{4}[ACGT]ATTGAGCTCGTGGGAGTGTCAGTG/i){
	$IgA++;
	$Iso = "IgA";

}
#Primer for IgG1/2:
if ($_ =~ m/^[ACT]{5}ACA[ACT]{5}ACA[ACT]{4}[ACGT]ATTCCTTTGACAAGGCATCC/i){
	$IgG1_2++;
	$Iso = "IgG1_2";
}
#Primer for IgG3:
if ($_ =~ m/^[ACGT]{1}[ACT]{3}CA[ACT]{3}[ACGT]{1}TTCGGGGAAGTAGTCCTTGAC/i){
	$IgG3++;
	$Iso = "IgG3";
}
#Primer for IgM:

if ($_ =~ m/^[ACT]{5}ACA[ACT]{5}ACA[ACT]{4}[ACGT]ATTCCATGGCCACCAGATTCTT/i){
	$IgM++;
	$Iso = "IgM";
}

#all FID
if ($_ =~ m/^[ACT]{5}ACA[ACT]{5}ACA[ACT]{4}[ACGT]ATT/i){
	$FID++;
}
#all RID
if ($_ =~ m/(CTGC[ACGT][AGT]{3}GT[AGT]{4}GT[AGT]{4})$/i){
	$RID++;
}

##############
#FID && RID
##############
#all RID + FID
if ($_ =~ m/^([ACT]{5}ACA[ACT]{5}ACA[ACT]{4}[ACGT]ATT)(.*)(CTGC[ACGT][AGT]{3}GT[AGT]{4}GT[AGT]{4})$/i){
	$RID_FID++;
	$one = $1;
	$two = $2;
	$three = $3;	
	
	if (length($2) > 49){
		push @{$FID_hash{$Iso."\|".$one}}, "$two";
		push @{$RID_hash{$Iso."\|".$three}}, "$two";
		push @{$FID_RID_hash{$Iso."\|".$one."\|".$three}}, "$two";
	}
}


}


#open output file
open (FILE2, "> $dir/$output_file") || die "problem opening $output_file\n";
#FID clean up
foreach (keys %FID_hash){
	#print "$_\n";
	$UMI = $_;
	foreach (@{$FID_hash{$_}}){
	#print "$_\n";
	$sorting_hash{$_}++;		
	}
	@sorted_hash = sort { $sorting_hash{$b} <=> $sorting_hash{$a} } keys(%sorting_hash);
	$corrected = $sorted_hash[0];	#

	foreach (keys @{$FID_hash{$_}}){	#
	$i_counter++;	#
	 #
	@sorted_hash = ();
	%sorting_hash = ();
	}
	$FID_hash_correted{$UMI."-$i_counter"} = $corrected;
	$i_counter = 0;		#
}

foreach (sort keys %FID_hash_correted){
$primary_key++;
print FILE2 ">$primary_key\|$_\n";
print FILE2 $FID_hash_correted{$_}."\n";
}
$primary_key = 0;
#close output file
close (FILE2);

#open output file
#open (FILE2, "> $dir/RID_$output_file") || die "problem opening $output_file\n";
#RID clean up
#foreach (keys %RID_hash){
#	#print "$_\n";
#	$UMI = $_;
#	foreach (@{$RID_hash{$_}}){
#	#print "$_\n";
#	$sorting_hash{$_}++;		
#	}
#	@sorted_hash = sort { $sorting_hash{$b} <=> $sorting_hash{$a} } keys(%sorting_hash);
#	$corrected = $sorted_hash[0];	#
#
#	foreach (keys @{$RID_hash{$_}}){	#
#	$i_counter++;	#
#	@sorted_hash = ();
#	%sorting_hash = ();
#	}			#
#	$RID_hash_correted{$UMI."-$i_counter"} = $corrected;
#	$i_counter = 0;		#
#}

#foreach (sort keys %RID_hash_correted){
#$primary_key++;
#print FILE2 ">$primary_key\|$_\n";
#print FILE2 $RID_hash_correted{$_}."\n";
#}
#$primary_key = 0;
#close output file
#close (FILE2);

#open output file
#open (FILE2, "> $dir/FID_RID_$output_file") || die "problem opening $output_file\n";
#FID/RID clean up
#foreach (keys %FID_RID_hash){
#	#print "$_\n";
#	$UMI = $_;
#	foreach (@{$FID_RID_hash{$_}}){
#	#print "$_\n";
#	$sorting_hash{$_}++;		
#	}
#	@sorted_hash = sort { $sorting_hash{$b} <=> $sorting_hash{$a} } keys(%sorting_hash);
#	$corrected = $sorted_hash[0];	#
#
#	foreach (keys @{$FID_RID_hash{$_}}){	#
#	$i_counter++;	#
#	@sorted_hash = ();
#	%sorting_hash = ();
#	}			
#	$FID_RID_hash_correted{$UMI."-$i_counter"} = $corrected; #
#	$i_counter = 0;		#
#}

#foreach (keys %FID_RID_hash_correted){
#$primary_key++;
#print FILE2 ">$primary_key\|$_\n";
#print FILE2 $FID_RID_hash_correted{$_}."\n";
#}
#$primary_key = 0;
#close output file
#close (FILE2);

#output to log file

#Iso Types
print FILE1 "IgG1/IgG2:	$IgG1_2\n";
print FILE1 "IgG3:	$IgG3\n";
print FILE1 "IgM:	$IgM\n";
print FILE1 "IgA:	$IgA\n";


#FID/RID
print FILE1 "all:	$all\n";
print FILE1 "FID:	$FID\n";
print FILE1 "RID:	$RID\n";
print FILE1 "RID\+FID:	$RID_FID\n";

#close output files
close (FILE1);


#isotype specific output
open (FILE1, "> $dir/IgG12_$output_file") || die "problem opening $output_file\n";
foreach (sort keys %FID_hash_correted){
$primary_key++;
$loopkey=$_;
if ($loopkey =~ m/^(IgG1_2)\|.*?$/i){
	if ($1 eq "IgG1_2"){
	print FILE1 ">$primary_key\|$loopkey\n";
	print FILE1 $FID_hash_correted{$loopkey}."\n";
	}
}
}
$primary_key = 0;
#close output file
close (FILE1);

open (FILE1, "> $dir/IgG3_$output_file") || die "problem opening $output_file\n";
foreach (sort keys %FID_hash_correted){
$primary_key++;
$loopkey=$_;
if ($loopkey =~ m/^(IgG3)\|.*?$/i){
	if ($1 eq "IgG3"){
	print FILE1 ">$primary_key\|$loopkey\n";
	print FILE1 $FID_hash_correted{$loopkey}."\n";
	}
}
}
$primary_key = 0;
#close output file
close (FILE1);

open (FILE1, "> $dir/IgM_$output_file") || die "problem opening $output_file\n";
foreach (sort keys %FID_hash_correted){
$primary_key++;
$loopkey=$_;
if ($loopkey =~ m/^(IgM)\|.*?$/i){
	if ($1 eq "IgM"){
	print FILE1 ">$primary_key\|$loopkey\n";
	print FILE1 $FID_hash_correted{$loopkey}."\n";
	}
}
}
$primary_key = 0;
#close output file
close (FILE1);

open (FILE1, "> $dir/IgA_$output_file") || die "problem opening $output_file\n";
foreach (sort keys %FID_hash_correted){
$primary_key++;
$loopkey=$_;
if ($loopkey =~ m/^(IgA)\|.*?$/i){
	if ($1 eq "IgA"){
	print FILE1 ">$primary_key\|$loopkey\n";
	print FILE1 $FID_hash_correted{$loopkey}."\n";
	}
}
}
$primary_key = 0;
#close output file
close (FILE1);

open (FILE1, "> $dir/NA_$output_file") || die "problem opening $output_file\n";
foreach (sort keys %FID_hash_correted){
$primary_key++;
$loopkey=$_;
if ($loopkey =~ m/^(NA)\|.*?$/i){
	if ($1 eq "NA"){
	print FILE1 ">$primary_key\|$loopkey\n";
	print FILE1 $FID_hash_correted{$loopkey}."\n";
	}
}
}
$primary_key = 0;
#close output file
close (FILE1);

