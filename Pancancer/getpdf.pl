#!/usr/bin/perl -w
=cut
to download the pdf format picture from UCSC for the desired region
we should set up the configurations in the webbrowser at first
the $hgsid should be the same as the local webbrowser 
Can only be run in desktop since you need input hgsid manually from internet explorer. 
This hgsid is located in the end of the UCSC browser after "hgsid="
input format： 
APC   chr2   2003485487  2003485687  cg564845541
awk '{print $2,$3-5000,$4+5000,$5,$6}' OFS="\t" TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.FinallMarker_E_6.bed  | grep chr | sort -u > TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.FinallMarker_E_6.bed.getpdf.bed

=cut

use strict;

use LWP 5.64;
use HTTP::Cookies;
my $cookie_jar = HTTP::Cookies->new( file => "cookies.lwp", );

my $browser = LWP::UserAgent->new;
# $browser->proxy( 'http', 'http://192.168.217.197:808' );
$browser->cookie_jar($cookie_jar);
 
my $hgsid = "718562257_uriICVZWLNlHhQwwa5vOAhCv6m7q";  # change it each time when you open your chorsome

my $DIR = '/gpfs/home/guosa/hpc/methylation/Pancancer';       
chdir $DIR or die "Cannot change to the $DIR:$!";
my $file = "getpdf.V2.hg19.bed";    #the coord_file to be download sequence
my $address;

#set the font size
#my $address = "http://genome.ucsc.edu/cgi-bin/hgTracks?textSize=medium&hgsid=$hgsid";#text size:tiny,small,medium,large,huge
#print "$address\n";
#my $tempory = $browser->get($address);

open FH, $file or die "Cannot open the file:$!";    #NM_001029935	chr13	23790020	23790578 
my @Coordination = <FH>;
my $length = scalar(@Coordination);
for my $Coord (@Coordination) {
	next if $Coord =~ /^\s|Name/;
	my ( $ID, $chr, $start, $end, $sample1 ) = split /\t/, $Coord;
	# mkdir $sample1 if !( -e $sample1 );
    #other features on the ucsc site
	my $UCSC_address = "&refGene=pack&refGene.priority=0.012";
	my @ucscFeature = qw (cons46wayViewphyloP.phyloP.vis=full cutters=dense cons46wayViewalign.align.vis=full xenoRefGene=pack cpgIslandExt=pack rmsk=dense cons46way=full );
	for my $ucscFeature (@ucscFeature) {
		( my $temp, undef ) = split /=/, $ucscFeature;
		$UCSC_address .= "&$ucscFeature&$temp.priority=0.02";
	}
	$address = "http://genome.ucsc.edu/cgi-bin/hgTracks?pix=1500&textSize=medium&hgsid=$hgsid$UCSC_address";
    print "change display mode\n$address\n";
	my $tempory = $browser->get($address);

	my $addition = 2000; 
	my $coor = $chr . ":" . $start . "-" . $end;
	getPDF( $coor, "_0K", $ID,$sample1);
	my ($coor2K) =  $chr . ":" . ( $start - $addition ) . "-" . ( $end + $addition );
	getPDF( $coor2K, "_2K",$ID, $sample1);
	print "\n";
}

sub getPDF {
	my $add = shift;
	my $extension = shift;
	my $ID = shift;
	my $sample1 = shift;
	$address ="http://genome.ucsc.edu/cgi-bin/hgTracks?position=$add&hgsid=$hgsid";
	print "The Address is:\n\t", $address, "\n";
	my $content = $browser->get($address);

	$content = $browser->get("http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=$hgsid&hgt.psOutput=on");
	my @content = split /\n/, $content->content;

	foreach my $temp (@content) {
		chomp $temp;
		next if ( $temp =~ /^\<\// );
		if ( $temp =~ /(hgt_genome.+?pdf)/ ) {
			my $pdf_add = "http://genome.ucsc.edu/trash/hgt/" . "$1";
			print "The location of the pdf is:\n\t$pdf_add\n";
			my $pdf = $browser->get( $pdf_add, ':content_file' => "$ID$extension.pdf" );#
		}
	}

}

sub getPDFWithBed {
	my $add = shift;
	my $extension = shift;
	my $ID = shift;
	my $sample1 = shift;
	$address ="http://genome.ucsc.edu/cgi-bin/hgTracks?position=$add&hgsid=$hgsid";
	print "The Address is:\n\t", $address, "\n";
	my $content = $browser->get($address);

	$content = $browser->get("http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=$hgsid&hgt.psOutput=on");
	my @content = split /\n/, $content->content;

	foreach my $temp (@content) {
		chomp $temp;
		next if ( $temp =~ /^\<\// );
		if ( $temp =~ /(hgt_genome.+?pdf)/ ) {
			my $pdf_add = "http://genome.ucsc.edu/trash/hgt/" . "$1";
			print "The location of the pdf is:\n\t$pdf_add\n";
			my $pdf = $browser->get( $pdf_add, ':content_file' => "$ID$extension.pdf" );#
		}
	}

}

sub getPDFWithPCRPRIMER {
	my $add = shift;
	my $extension = shift;
	my $ID = shift;
	my $sample1 = shift;
	$address ="http://genome.ucsc.edu/cgi-bin/hgTracks?position=$add&hgsid=$hgsid";
	print "The Address is:\n\t", $address, "\n";
	my $content = $browser->get($address);

	$content = $browser->get("http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=$hgsid&hgt.psOutput=on");
	my @content = split /\n/, $content->content;

	foreach my $temp (@content) {
		chomp $temp;
		next if ( $temp =~ /^\<\// );
		if ( $temp =~ /(hgt_genome.+?pdf)/ ) {
			my $pdf_add = "http://genome.ucsc.edu/trash/hgt/" . "$1";
			print "The location of the pdf is:\n\t$pdf_add\n";
			my $pdf = $browser->get( $pdf_add, ':content_file' => "$ID$extension.pdf" );#
		}
	}

}


