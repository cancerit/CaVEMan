#! /usr/bin/perl

############## LICENSE ##############
# Copyright (c) 2014-2015 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of CaVEMan.
#
# CaVEMan is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#    1. The usage of a range of years within a copyright statement contained within
#    this distribution should be interpreted as being equivalent to a list of years
#    including the first and last year specified and all consecutive years between
#    them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
#    2009, 2011-2012’ should be interpreted as being identical to a statement that
#    reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
#    statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
#    identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
#    2009, 2010, 2011, 2012’."
#
############## LICENSE ##############

use strict;
use warnings FATAL=>'all';

use autodie qw(:all);
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Carp;

use Data::Dumper;

{
  	my $opts = setup();
  	my $baselookup = undef;
  	my $IN;
  	my $OUT;
  	open($OUT,'>',$opts->{'o'}) or croak("Error trying to open output VCF file '".$opts->{'o'}."': $!");
			open($IN,'<',$opts->{'f'}) or croak("Error opening input VCF file '".$opts->{'f'}."': $!");
				while(<$IN>){
					my $line = $_;
					if($line =~ m/^\s*#/){
						print $OUT $line;
					}else{
						chomp($line);
						my ($chr,$pos,$tmpa,$ref_base,$tmpb,$tmpc,$tmpd,$dptheq,$format,@samples) = split /\t/, $line;
						my ($gtf,@ftbases) = split /:/, $format;
						my ($d,$dpth) = split /=/, $dptheq ;
						if(!defined($baselookup)){
							%$baselookup = ( 	"A" => 0,
																"C" => 0,
																"G" => 0,
																"T" => 0);
							for(my $i=0;$i<scalar(@ftbases);$i++){
								my $bs = substr($ftbases[$i],0,1);
								foreach my $ky(keys %$baselookup){
									if($bs eq $ky){
										$baselookup->{$bs} = $i;
									}
								}
							}
						}
						my $check = 0;
						my $total = 0;
						for my $samp(@samples){
							next if($samp eq "-");
							my ($gt,@base_counts) = split /:/, $samp;
							$total += $base_counts[$baselookup->{$ref_base}];
						}
						if($dpth != $total){
							print $OUT $line,"\n";
						}
					}
				}
			close($IN);
  	close($OUT);
}

sub setup {
  my %opts;
  GetOptions( 	'h|help' => \$opts{'h'},
								'f|file=s' => \$opts{'f'},
								'o|output=s' => \$opts{'o'},
  ) or pod2usage(2);

  pod2usage(0) if(defined $opts{'h'});
	pod2usage(1) unless(defined($opts{'f'}) && defined($opts{'o'}));
	pod2usage(1) unless(-e $opts{'f'} && -f $opts{'f'});
	return \%opts;
}

__END__

=head1 NAME

removeVCFPAnelRefLines.pl - Remove the lines where samples only show reference alleles from VCF unmatched normal panel files.

=head1 SYNOPSIS

removeVCFPAnelRefLines.pl [options]

  	Required parameters:
    -output       -o   File to write results to.
    -file         -f   input VCF file

	  Other:
    -help     -h   Brief help message.

=cut
