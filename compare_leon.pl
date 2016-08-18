#!/usr/bin/perl -w
# Author : Charles VAN GOETHEM

# Perl general libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper;
use Cwd;
use Time::HiRes;

# Perl libs for given/when (smartmatch is experimental)
use v5.14;
no warnings 'experimental::smartmatch';

use File::stat;

# Perl for signature
use Term::ANSIColor;


##########################################################################################
##########################################################################################


my $path_leon = "/Users/adminbioinfo/Documents/Leon/leon/leon";
my $pwd = getcwd();


##########################################################################################
##########################################################################################

# Mandatory arguments
my @files;
my @directories;

# Optional arguments
my $output = $pwd."/";
my $makeGraphs = '';

# General arguments
my $man 		= 0;
my $help 		= 0;
my $verbosity	= 0;


## Parse options and print usage if there is a syntax error, 
## or if usage was explicitly requested.

GetOptions(
	'f|file=s'	 		=> \@files,
	'd|directory=s'		=> \@directories,
	'o|output=s'		=> \$output,
	'g|makeGraphs'		=> \$makeGraphs,
	'v|verbosity=i'		=> \$verbosity,
	'help|?' 			=> \$help,
	'm|man' 			=> \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

@files = split(/,/,join(',',@files));
@directories = split(/,/,join(',',@directories));

# Test arguments
pod2usage("Error : need file(s) and/or directory(ies).") unless (@files or @directories);

# Check the files
foreach my $file (@files) {
	pod2usage("Error : $file is not a plain text or not exist.") unless (-f $file);
	
	# Check if file is a fastq or fastq.gz
	#	my ($name, $dir, $ext) = fileparse($file,qw(.fastq .fastq.gz));
	my ($name, $dir, $ext) = fileparse($file,qw(.fastq));
	pod2usage("Error : $file is not a fastq or a fastq.gz. (only fastq until now)") unless ($ext);
}

# Check directories and get fastq and fastq.gz
foreach my $directory (@directories) {
	pod2usage("Error : $directory is not a directory or not exist.") unless (-d $directory);
	
	# Complete the list of files with the fastq and fastq.gz in the repertory
	#push(@files, glob "$directory/*.{fastq,fastq.gz}");
	push(@files, glob "$directory/*.fastq");
}

# print Dumper @files;
# exit();

my $filename = $output."/report.tab";
open(OUT , ">$filename") or die "Could not open file '$filename' $!";

######## Header file (on the same line)
#	File_ID	Size_fastQ
#	Time_comp_leon_lossy	Size_comp_leon_lossy
#	Time_comp_leon_lossless	Size_comp_leon_lossless
#	Time_comp_gzip1	Size_gzip1
#	Time_comp_gzip6	Size_gzip6
#	Time_comp_gzip9	Size_gzip9
#	Time_uncomp_leon_lossy	Time_uncomp_leon_lossless
#	Time_uncomp_gzip1	Time_uncomp_gzip6	Time_uncomp_gzip9
print OUT "File_ID\tSize_fastQ\t";
print OUT "Time_comp_leon_lossy\tSize_comp_leon_lossy\t";
print OUT "Time_comp_leon_lossless\tSize_comp_leon_lossless\t";
print OUT "Time_comp_gzip1\tSize_gzip1\t";
print OUT "Time_comp_gzip6\tSize_gzip6\t";
print OUT "Time_comp_gzip9\tSize_gzip9\t";
print OUT "Time_uncomp_leon_lossy\tTime_uncomp_leon_lossless\t";
print OUT "Time_uncomp_gzip1\tTime_uncomp_gzip6\tTime_uncomp_gzip9\n";


##########################################################################################
##########################################################################################


my $i=1;
my $line ="";
foreach my $file (@files) {
	my ($name, $dir, $ext) = fileparse($file,qw(.fastq .fastq.gz));
	
	# uncompress if necessary 
	given($ext) {
		when(".fastq.gz") { # Uncompress then launch the bench
			print "file : \"".$dir.$name.$ext."\" is a fastQ.gz\n";
		}
		#when(".fastq") { # Launch the fastq bench function
		#	print "file : \"".$dir.$name.$ext."\" is a fastQ\n";
		#}
	}
	
	my $size_origin = stat($file)->size;
	$line = "File_".$i."\t".$size_origin."\t";
	
	##################################
	######### Compressions
	
	##################################
	##### LEON
	## lossy (best)
	## lossless (depreciated)
	
	my $file_leon_lossy = $output.$name.".lossy-leon.fastq";
	my $file_leon_lossless = $output.$name.".lossless-leon.fastq";
	
	system("ln -s ".$dir.$name.".fastq ".$file_leon_lossy);
	system("ln -s ".$dir.$name.".fastq ".$file_leon_lossless);
	
	my $time_comp_leon_lossy = timer_cmd_bash($path_leon." -file ".$file_leon_lossy." -c");
	my $time_comp_leon_lossless = timer_cmd_bash($path_leon." -file ".$file_leon_lossless." -c -lossless");
	
	my $size_lossy = size_for_leon_files($file_leon_lossy);
	my $size_lossless = size_for_leon_files($file_leon_lossless);
	
	$line .= $time_comp_leon_lossy."\t".$size_lossy."\t";
	$line .= $time_comp_leon_lossless."\t".$size_lossless."\t";
	
	##################################
	##### Gzip
	## level 1 (fastest)
	## level 6 (default)
	## level 9 (best)
	
	my $file_gzip1 = $output.$name."_level-1.fastq.gz";
	my $file_gzip6 = $output.$name."_level-6.fastq.gz";
	my $file_gzip9 = $output.$name."_level-9.fastq.gz";
	
	my $time_comp_gzip1 = timer_cmd_bash("gzip -ck1 ".$file." > ".$file_gzip1);
	my $time_comp_gzip6 = timer_cmd_bash("gzip -ck6 ".$file." > ".$file_gzip6);
	my $time_comp_gzip9 = timer_cmd_bash("gzip -ck9 ".$file." > ".$file_gzip9);
	
	my $size_gzip1 = stat($file_gzip1)->size;
	my $size_gzip6 = stat($file_gzip6)->size;
	my $size_gzip9 = stat($file_gzip9)->size;
	
	$line .= $time_comp_gzip1."\t".$size_gzip1."\t";
	$line .= $time_comp_gzip6."\t".$size_gzip6."\t";
	$line .= $time_comp_gzip9."\t".$size_gzip9."\t";
	
	##################################
	######### Uncompresses
	
	##################################
	##### LEON
	
	my $time_uncomp_leon_lossy = timer_cmd_bash($path_leon." -file ".$file_leon_lossy.".leon -d");
	my $time_uncomp_leon_lossless = timer_cmd_bash($path_leon." -file ".$file_leon_lossless.".leon -d");
	
	$line .= $time_uncomp_leon_lossy."\t".$time_uncomp_leon_lossless."\t";
	
	##################################
	##### Gzip
	
	my $time_uncomp_gzip1 = timer_cmd_bash("gunzip ".$file_gzip1);
	my $time_uncomp_gzip6 = timer_cmd_bash("gunzip ".$file_gzip6);
	my $time_uncomp_gzip9 = timer_cmd_bash("gunzip ".$file_gzip9);
	
	$line .= $time_uncomp_gzip1."\t".$time_uncomp_gzip6."\t".$time_uncomp_gzip9."\n";
	
	##################################
	######### Write .tab file
	
	print OUT $line;
	
	system("rm $output/$name*");
	
	$i++;
}

### LAUNCH make_Awesome_graphs.R
if($makeGraphs) {
	system("Rscript --vanilla make_Awesome_graphs.R ".$output."report.tab ".$output);
}

##########################################################################################
##########################################################################################


sub timer_cmd_bash {
	my $cmd = shift(@_);
	
	my $timer_S;
	my $timer_E;
	my $diff;
	
	$timer_S = time();
	system($cmd);
	#print $cmd."\n";
	$timer_E = time();
	
	$diff = $timer_E - $timer_S;
	
	return $diff;
}

sub size_for_leon_files {
	my $file = shift(@_);

	my $stat_leon = stat($file.".leon");
	my $stat_qual = stat($file.".qual");
	
	my $size_total = ($stat_leon->size) + ($stat_qual->size);
	
	return $size_total;
}


##########################################################################################
##########################################################################################

__END__

=pod

=encoding UTF-8

=head1 NAME

compare_leon.pl - Compare Leon compression algorithm with gzip (generally used) for fastQ.

=head1 VERSION

version 0.01

=head1 SYNOPSIS

compare_leon.pl  -f file.fastq (-f file2.fastq,file3.fastq) -d directory/with/some/fastq

=head1 DESCRIPTION

This script compress and uncompress automatically some fastQ file.

=head1 OPTIONS

=head2 General

	-h,--help		Print this help
	-m,--man		Open man page
	-v,--verbosity		Level of verbosity

=head2 Mandatory arguments

	-f,--file=file.fastq			Specify the fastq file you want use (possible multiple file) 
	-d,--directory=path/to/directory	Specify a directory contains some fastq you want use (possible multiple directory)

=head2 Optional arguments

	-o,--output=repertory			You can specify the output repertory (default Current)
	-g,--makeGraphs					Generate automatically awesome graphs with R and ggplot2

=head1 AUTHORS

=over 4

=item -
Charles VAN GOETHEM

		      _ __,~~~/
		,~~`( )_( )-\|
		    |/|  `--.
		    ! !  !

=back

=cut
