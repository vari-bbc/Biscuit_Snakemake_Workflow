#!/usr/bin/perl -w
use strict;
use List::Util qw( min max );
use Statistics::Basic qw(:all);
=head1
Author: Ian Beddows
Date: 3/9/18 

Short description:	
		
		Takes a directory (-indir) and loads all *.bed.gz files (these files are the output from biscuit pileup).
		These files are first loaded and the average observed methylation foreach window is calculated (bookended windows)
		Window averages are calculated only for positions where the depth is greater than or equal to -covFilter (set to 0 to include all positions)
		Finally, the window averages are output to a file (-out) in tab-separed file format (with header):
		chromosome	windowStart	windowEnd	sample1_average_methylation sample2_avg_meth ....
		
		
		
=cut
use Getopt::Long;

my $usage = <<EOF;
OPTIONS:
-indir (directory with bed.gz files)
-windowSize (e.g. 10, 1000, 10000000)
-out
-covFilter (positions >=covFilter are considered)
OPTIONAL:
-h|help = print usage
-b = uppercase arguments
EOF
my($help,$bold,$indir,$outfile,$windowSize,$covFilter);
#=======================================================================
GetOptions(
	'outfile=s' => \$outfile,	# string
	'indir=s' => \$indir,	# string
	#~ 'windowSize=s' => \$windowSize,	# string
    'windowSize=i' => \$windowSize,	# integer
    'covFilter=i' => \$covFilter,	# integer
	#~ '=f' => \$,	# floating point number
    #~ 'bold' => \$bold,	# flag
	'h|help' => \$help	# flag for help
);

if (defined($help)){die print "HELP:\n",$usage;}
if (!defined($outfile)){die print "HELP:\n",$usage;}
if (!defined($indir)){die print "HELP:\n",$usage;}
if (!defined($covFilter)){die print "HELP:\n",$usage;}
if (!defined($windowSize)){die print "HELP:\n",$usage;}

open(my $out,'>',$outfile)|| die

my %data = ();
my %chr = ();
my %samples = ();

#~ my @files = `ls $indir|grep .bed.gz\$ |head -2`;
my @files = `ls $indir|grep .bed.gz\$`;
chomp(@files);
for(my $i=0;$i<=$#files;$i++){
	my $f = $files[$i];
	my $sample = $f;
	$sample=~s/\.bed\.gz//;
	$samples{$sample}=();
	print "Opening $f (sample $sample)\n";
	#~ open(my $in, "gunzip -c $indir$f |head -20000 |") || die "can't open pipe";
	open(my $in, "gunzip -c $indir$f  |") || die "can't open pipe";
	while(<$in>){
		chomp;
		my($chr,$start,$stop,$meth,$cov) = split("\t",$_);
		#~ print("$sample\t$chr\t$meth\n");
		if($cov>=$covFilter){
			${data}{$chr}{$start}{$sample} = $meth;
			$chr{$chr}=();
		} # else act like there is no observation there
	}
	#~ if($i==2){last;} # testing
}

my $nSamples = scalar keys %samples;
print STDOUT "\nData for ",$nSamples," samples loaded.\n";

# do the sliding window
my %window_avg=();
my %outdata = ();


# print header to outfile
print $out "chr\twindowPos";
foreach my $s (sort keys %samples){
	print $out "\t$s";
}
print $out "\n";

foreach my $chr (sort keys %chr){
	
	#
	#	This is a sliding window approach with -window bp
	#
	print STDOUT "Getting window averages for $chr.\n";
	#my $m = max(keys $data{$chr});
	#~ #print "$chr\t$m\n";
	#~ for(my $i=0;$i<=$m;$i+=$windowSize){
		#~ # foreach of the samples, get the average within the window
		#~ for(my $j=$i;$j<($i+$windowSize);$j++){
			#~ if(exists $data{$chr}{$j}){
				#~ #print "$chr\twindow # $i\tposition $j\n";
				#~ #print "\t\t",scalar keys %{$data{$chr}{$j}}," sample(s)\n";
				#~ #print "\t\t",join("\t",keys %{$data{$chr}{$j}},"\n");
				#~ foreach my $sample (keys %{$data{$chr}{$j}}){
					#~ push(@{$window_avg{$i}{$sample}},$data{$chr}{$j}{$sample});
				#~ }
			#~ }
			
		#~ }
	#~ }
	
	# this is an approach for a -window size of # CpGs (e.g. 5 CpG windows)
	my $nCpG=0;
	my $windowN=0;
	foreach my $pos (sort {$a<=>$b} keys %{$data{$chr}}){
		$nCpG++;
		# foreach of the samples, get the average within the window
		foreach my $sample (keys %{$data{$chr}{$pos}}){
			push(@{$window_avg{$windowN}{$sample}},$data{$chr}{$pos}{$sample});
		}
		if($nCpG==$windowSize){
			# go to next windowN
			$windowN++;
			$nCpG=0;
		}
	}
	

	# now we have the averages, so report those
	# but only if $numberNA flag is equal to number of samples (no missing data)
	foreach my $window (sort {$a<=>$b} keys %window_avg){
		my $nNA=0;
		my $outString = "$chr\t$window";
		#~ print STDOUT "$chr\t$window\n";
		foreach my $sample (sort keys %samples){
			if(exists $window_avg{$window}{$sample}){
				my $nObs = scalar @{$window_avg{$window}{$sample}};
				my $mean_meth = mean(@{$window_avg{$window}{$sample}});
				#~ print STDOUT "\t$sample\t$nObs\t$mean_meth\n";
				$outString .= "\t$mean_meth";
			}else{
				$nNA++;
				$outString .= "\tNA"
			}
		}
		
		if($nNA!=$nSamples){
			print $out "$outString\n";
		}else{
			print STDOUT "Skipping window $window because $nNA NAs\n";
			# do nothing
		} 
	}
	#~ last; # testing
}









print "\nFinished $0\n";
#=======================================================================
#( Subroutines                  )
# ------------------------------------ 
#  o
#   o   \_\_    _/_/
#    o      \__/
#           (oo)\_______
#           (__)\       )\/\
#               ||----w |
#               ||     ||
