#!/bin/perl -w
use strict;
my $filein  = $ARGV[0];
my $fileout = $ARGV[1];


open(SRC,"$filein")||die "can not open file:";

my %clu2hit = ();

while (my $line=<SRC>) {
    chomp $line;
    my @temp=split(/\t/,$line);
    my $start = $temp[2];
    
    my $dis = abs($start-26);
    push @temp, $dis;
    my $clu = $temp[1];
    
    if($dis <=1000 ){
    	if(exists $clu2hit{$clu}){
    		push @{$clu2hit{$clu}},[@temp];
    	}
    	else{
    		my @ttt = ();
    		$clu2hit{$clu} = \@ttt ;
    		push @{$clu2hit{$clu}},[@temp];
    	}
    }
}

### output file
open(TGT,">$fileout")||die "can not open file:";


for my $clu (sort keys %clu2hit){
	my @test =@{$clu2hit{$clu}};
	
	@test = sort {$a->[-1] <=> $b->[-1]} @test ;
	@test = sort {$a->[-2] <=> $b->[-2]} @test ;
	

	my @line = @{$test[0]};

	my $id = $line[0];
	my $fam = '';
	if($id =~ m/miR(\d+)/){
		$fam = $1;
	}
	
	push @line, $fam;
	print TGT join("\t",@line),"\n";
	
}

