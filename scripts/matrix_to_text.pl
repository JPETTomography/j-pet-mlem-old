#!/usr/bin/perl -w 

=pod
/// uint32_t magic = 'PETt'
/// uint32_t n_pixels_2 
/// uint32_t n_emissions // per pixel
/// while (!eof)
///   uint16_t lor_a, lor_b // pair
///   uint32_t pixel_pair_count
///   for(.. count ..)
///     uint16_t t_pixel_x, t_pixel_y
///     uint32_t pixel_hits
=cut


use strict;


while(my $file=shift @ARGV) {
    print STDERR $file,"\n";
    open(FILE,"<$file") or warn "cannot open file `$file'  for reading.";
    my $buffer;
    read FILE,$buffer,4;
    my $magic=unpack "a4", $buffer;

    print "magic=$magic\n";

    my ($pixels_2,$n_emissions,$n_detectors);
    if($magic eq "PETP"  ) {
	print "PETP\n";
	read FILE,$buffer,12;
	($pixels_2,$n_emissions,$n_detectors)=unpack "III", $buffer;
    } elsif($magic eq "PETs" ) {
	print "PETs\n";
	read FILE,$buffer,8;
	($pixels_2,$n_emissions)=unpack "II", $buffer;	
	$n_detectors="unknown";
    } 


    print "#$magic $pixels_2 $n_emissions $n_detectors\n";
    while(!(eof FILE)) {
	read FILE,$buffer,8;
	my ($lor_i,$lor_j,$count)=unpack "SSI",$buffer;
	print "($lor_i ,  $lor_j) $count\n";  
	for(my $i=0;$i<$count;++$i) {
	    read FILE, $buffer,8;
	    my ($pix_x,$pix_y,$pix_count)=unpack "SSI",$buffer;
	    print "$pix_x $pix_y $pix_count\n";
	}
    }
}
