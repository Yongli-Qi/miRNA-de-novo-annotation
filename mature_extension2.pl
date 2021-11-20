#!/bin/perl -w
use strict;

my $fileIn   = $ARGV[0];
my $fileOut  = $ARGV[1];


####################
### fasta              
open SRC, "<$fileIn";
open TGT, ">$fileOut";

while (my $line = <SRC>){
    chomp $line;
    if($line=~/^\>/){  
    	  chomp $line;
    	  my $id_full = substr $line,1; 	
    	  my @temp = split(/\s+/,$id_full);
    	  my $id = $temp[0];
    	  $line=<SRC>;
    	  chomp $line;    		
    	  my $seq = 'AAAAAAAAAAAAAAAAAAAAAAAAA'.$line.'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA';    
    	  
          $seq =~ s/U/T/g;
    	  print TGT '>',$id,"\n";
    	  	  	  print TGT $seq,"\n";
    	  	  
    } ### 
 }
 
 
