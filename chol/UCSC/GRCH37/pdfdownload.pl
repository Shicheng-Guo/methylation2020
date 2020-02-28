use strict;
use Cwd;
use LWP 5.64;
use HTTP::Cookies;
my $cookie_jar = HTTP::Cookies->new(
   file     => "cookies.lwp",
);
my $browser = LWP::UserAgent->new;
$browser->cookie_jar( $cookie_jar );
my $hgsid = "686065221_uoaMviuXpomSlykn0VM0jEEDRIA2";
my $DIR =getcwd;
chdir $DIR or die "Cannot change to the $DIR:$!";
open FH, "tobedownload.txt" or die "Cannot open the file:$!";
my @Coordination = <FH>;
my $length = scalar(@Coordination);
for my $Coord (@Coordination) {
	next if $Coord =~ /^\s|Name/;
	chop $Coord;
	my ($ID,$chr,$start,$end)= split /[\s+:-]/, $Coord;
	my $addition = 200;
	my($coor) = $chr.":".($start-$addition)."-".($end+$addition);
    my $filename=$chr.".".($start-$addition)."-".($end+$addition);
	my $address = "http://genome.ucsc.edu/cgi-bin/hgTracks?position=$coor&hgsid=$hgsid";
	print "The Address is:\n\t",$address,"\n";
	my $content =  $browser->get($address);
	$content = $browser->get("http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=$hgsid&hgt.psOutput=on");
	my @content = split /\n/, $content->content;
	foreach my $temp (@content){
		chomp $temp;
		next if ($temp =~ /^\<\//);
		if ($temp =~ /(hgt_genome.+?pdf)/){
			my $pdf_add = "http://genome.ucsc.edu/trash/hgt/"."$1";
			print "The location of the pdf is:\n\t$pdf_add\n";
			if ($addition == 200){
				my $pdf = $browser->get($pdf_add,':content_file'=>"$ID.$filename.2k.pdf");
			}
			else{
				my $pdf = $browser->get($pdf_add,':content_file'=>"$ID.$filename.pdf");
			}	
		}
	}
}

