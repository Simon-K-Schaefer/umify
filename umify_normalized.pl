#! /usr/bin/perl
# Simon K. Schaefer
# Version 1.1

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


#Primer
my $IgG = 0;
my $IgM = 0;
my $IgA = 0;

my $RID;
my $FID;
my $RID_FID;
my $UMI;
my $Iso;
my $primary_key;
my $all;
my $i_counter;
my $corrected;

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
$output_file = substr($input_file,0,int((length($input_file)))-12)."_umified_norm.fasta";

#create output folder
$dir = substr($input_file,0,int((length($input_file)))-12)."_umified";
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
#NHHHCAHHHNTTGGGGCTGGTCGGGGATGC
if ($_ =~ m/^[ACGT]{1}[ACT]{3}CA[ACT]{3}[ACGT]{1}TTGGGGCTGGTCGGGGATGC/i){
	$IgA++;
	$Iso = "IgA";
#Primer for IgG:
}
#NHHHCAHHHNTTTTCGGGGAAGTAGTCCTTGAC
if ($_ =~ m/^[ACGT]{1}[ACT]{3}CA[ACT]{3}[ACGT]{1}TTCGGGGAAGTAGTCCTTGAC/i){
	$IgG++;
	$Iso = "IgG";
}
#Primer for IgM:
#NHHHCAHHHNGGTTGGGGCGGATGCACTCC
if ($_ =~ m/^[ACGT]{1}[ACT]{3}CA[ACT]{3}[ACGT]{1}GGTTGGGGCGGATGCACTCC/i){
	$IgM++;
	$Iso = "IgM";
}

#all FID
if ($_ =~ m/^[ACGT]{1}[ACT]{3}CA[ACT]{3}[ACGT]{1}/i){
	$FID++;
}
#all RID
if ($_ =~ m/[ACGT]{4}$/i){
	$RID++;
}

##############
#FID && RID
##############
#all RID + FID
if ($_ =~ m/^([ACGT]{1}[ACT]{3}CA[ACT]{3}[ACGT]{1})(.*)([AGCT]{4})$/i){
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
open (FILE2, "> $dir/FID_$output_file") || die "problem opening $output_file\n";
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
open (FILE2, "> $dir/RID_$output_file") || die "problem opening $output_file\n";
#RID clean up
foreach (keys %RID_hash){
	#print "$_\n";
	$UMI = $_;
	foreach (@{$RID_hash{$_}}){
	#print "$_\n";
	$sorting_hash{$_}++;		
	}
	@sorted_hash = sort { $sorting_hash{$b} <=> $sorting_hash{$a} } keys(%sorting_hash);
	$corrected = $sorted_hash[0];	#

	foreach (keys @{$RID_hash{$_}}){	#
	$i_counter++;	#
	@sorted_hash = ();
	%sorting_hash = ();
	}			#
	$RID_hash_correted{$UMI."-$i_counter"} = $corrected;
	$i_counter = 0;		#
}

foreach (sort keys %RID_hash_correted){
$primary_key++;
print FILE2 ">$primary_key\|$_\n";
print FILE2 $RID_hash_correted{$_}."\n";
}
$primary_key = 0;
#close output file
close (FILE2);

#open output file
open (FILE2, "> $dir/FID_RID_$output_file") || die "problem opening $output_file\n";
#FID/RID clean up
foreach (keys %FID_RID_hash){
	#print "$_\n";
	$UMI = $_;
	foreach (@{$FID_RID_hash{$_}}){
	#print "$_\n";
	$sorting_hash{$_}++;		
	}
	@sorted_hash = sort { $sorting_hash{$b} <=> $sorting_hash{$a} } keys(%sorting_hash);
	$corrected = $sorted_hash[0];	#

	foreach (keys @{$FID_RID_hash{$_}}){	#
	$i_counter++;	#
	@sorted_hash = ();
	%sorting_hash = ();
	}			
	$FID_RID_hash_correted{$UMI."-$i_counter"} = $corrected; #
	$i_counter = 0;		#
}

foreach (keys %FID_RID_hash_correted){
$primary_key++;
print FILE2 ">$primary_key\|$_\n";
print FILE2 $FID_RID_hash_correted{$_}."\n";
}
$primary_key = 0;
#close output file
close (FILE2);


#Iso Types
print FILE1 "IgG:	$IgG\n";
print FILE1 "IgM:	$IgM\n";
print FILE1 "IgA:	$IgA\n";


#FID/RID
print FILE1 "all:	$all\n";
print FILE1 "FID:	$FID\n";
print FILE1 "RID:	$RID\n";
print FILE1 "RID\+FID:	$RID_FID\n";

#close output files fore cleaned PDB (FILE1) and Log for changes (FILE3)
close (FILE1);

