#!/usr/bin/env perl

use local::lib;
use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use Getopt::Long;
use Config::Abstract::Ini;

my $paramfile;
my $samplefile;

my $versionString = "1.0";

my $miniManual = qq"
Usage: perl analyzeRNA.pl [-h] -p <parameters> -d <sampleTable>

Options:
   [-v|--version]                                 : Print out just the version.
   [-h|--help]                                    : Print this miniManual.
   [-d|--sampleTable]   <tab-delimited file>        : Sample table with all samples 
   [-p|--parameters]      <parameters file>            : Parameters file.
	
Version: $versionString	
Release date: 2020
Copyright: Alex M. Mawla. 2014-Present. All Rights Reserved. 
For questions and comments, please contact Alex.
";

if ((@ARGV == 0) || ($ARGV[0] eq "-h")||($ARGV[0] eq "--help")) {
    print "$miniManual\n";
    exit();
}

my $version=0;

GetOptions( 
	"help|h",
	"version|v" => \$version,
	"sampleTable|d:s" => \$samplefile, 
	"parameters|c:s" => \$paramfile,
) or die "Unknown option\n";

if ($version) {
    print "Version: $versionString\n";
    exit();
}	





open(PARAMETERS, "$paramfile") or die "Cannot open parameters file";
my %params;

while(<PARAMETERS>){
    $_ =~ s/\n//g;
    if($_ =~ /^#/){next;}
    if($_ =~ /^PARAM/){
	my @a = split/\s+/, $_;
	$params{$a[1]} = $a[2];
    }
}close(PARAMETERS);



open(SAMPLES, "$samplefile") or die "Cannot open sample file";
my @samples;
while(<SAMPLES>){
    $_ =~ s/\n//g;
    if($_ =~ /^#/){next;}
#if($_ =~ /^SAMPLE/){
    push(@samples, $_);
#    }
}close(SAMPLES);



my ($sec, $min, $hour, $mday,$mon,$year) = (localtime(time))[0..5];
$mon +=1;
$year += 1900;
print "$mon/$mday/$year $hour:$min:$sec\n";



foreach my $sample(@samples){
	 system("./pipeline.sh $params{cores} $params{mem} $params{part} $params{starCores} $params{bamCores} $params{genome} $params{length} $params{type} $params{rawDir} $params{outDir} $sample $params{adapter} $params{fastqcAdapters} $params{fastqcContaminants} &");
}




