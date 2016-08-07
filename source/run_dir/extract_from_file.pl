#! /usr/bin/perl -w

use strict;

my $DESCRIPTION = "Given a file in tabular format and another file with list of keys, extract the rows whose ith column is in keys file...\n";
my $USAGE = "$0 input_file  input_file (contains only list_of_keys) column_number(which column in the input file contains the keys..) output_file\n";

sub trim($);

unless(@ARGV){
    print "$USAGE\n";
    exit;
}

#my $chr_all_file;
#open input and output files...

unless (open(INPUT_FILE, $ARGV[0])){
    print "Cannot open file\n";
    exit;
}

unless (open(INPUT_FILE_KEYS, $ARGV[1])){
    print "Cannot open file\n";
    exit;
}

my $column_no = $ARGV[2]; #key value in input file is at this column...

unless (open(OUTPUT_FILE, ">$ARGV[3]")){
    print "Cannot open file\n";
    exit;
}

my @input_file = <INPUT_FILE>;
my @requested_keys = <INPUT_FILE_KEYS>;

my $key = "";



my %is_requested = ();

#undef %is_requested;

for (@requested_keys) {$_ = trim($_); $is_requested{$_} = 1}

foreach my $ln (@input_file){
    chomp($ln);
    
    my @columns = split("\t", $ln);
    $key = $columns[$column_no];
    $key = trim($key);

    if ($is_requested{$key}){
	print OUTPUT_FILE $ln."\n";
    }
}
################################################################################
# Subroutines
################################################################################



# get_next_record
#
#   - given GenBank record, get annotation and DNA
 
sub get_next_record {
 
    my($fh) = @_;
 
    my($offset);
    my($record) = '';
    my($save_input_separator) = $/;
 
    $/ = "\n";
 
    $record = <$fh>;
 
    $/ = $save_input_separator;
 
    return $record;
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
