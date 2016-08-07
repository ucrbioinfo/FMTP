#! /usr/bin/perl -w

use strict;

my $usage = "$0 original_ctg_clone_table GLPK_solution_file output_new_ctg_clone_table\n";


unless (@ARGV){
    die $usage;
}

unless (open(CTG_CLONE, $ARGV[0])){
    die "cannot open file\n";
}

unless (open(GLPK_SOL, $ARGV[1])){
    die "cannot open file\n";
}

unless (open(OUTPUT_FILE, ">>$ARGV[2]")){
    die "cannot open output_file\n";
}

my %clone_to_ctg = (); #key:clonename value:ctg_id.. stores the original clone_ctg associatioin..
my @solution = <GLPK_SOL>;
my %output_ctg_to_clone = (); #key ctg_id, value: selected clones in this contig.. (array)

#parse the original ctg_clone table and generate clone_to_ctg map...

foreach my $ln (<CTG_CLONE>){
    chomp($ln);

    my @columns = split("\t", $ln);
    my @columns_clones = split(" ", $columns[2]);
    foreach my $clone (@columns_clones){
	chomp($clone);
	$clone_to_ctg{$clone} = $columns[0];
    } 
}


#parse the GLPK solution file..

foreach my $ln (@solution){
    chomp($ln);
    if ($ln =~ /x\[/){
	my @columns = split(" ", $ln);
	#if this clone is in solution.. 
	if ($columns[3]){
	    my $clonename = substr($columns[1], 2, length($columns[1])-3);
	    #print $clonename."\n";
	    my $ctg_id = $clone_to_ctg{$clonename};
	    if (exists ${output_ctg_to_clone{$ctg_id}}){
		my @t = @{${output_ctg_to_clone{$ctg_id}}};
		push(@t, $clonename);
		${output_ctg_to_clone{$ctg_id}} = \@t;
	    }else{
		my @t = ();
		push(@t, $clonename);
		${output_ctg_to_clone{$ctg_id}} = \@t;
	    }
	}
    }
}
	    
#output_the_results..

foreach my $k (sort {$a <=> $b} keys %output_ctg_to_clone){
    my @t = @{${output_ctg_to_clone{$k}}};
    print OUTPUT_FILE "$k\t".@t."\t@t\n";
}
