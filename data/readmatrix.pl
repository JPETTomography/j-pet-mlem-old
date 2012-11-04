#!/usr/bin/perl -w 

=pod
/// uint32_t magic = 'PETt'
/// uint32_t n_pixels_2 
/// while (!eof)
///   uint16_t lor_a, lor_b // pair
///   uint32_t pixel_pair_count
///   for(.. count ..)
///     uint16_t t_pixel_x, t_pixel_y
///     uint32_t pixel_hits
=cut


use strict;

use PDL;

use PDL::Graphics::PGPLOT::Window;
use PDL::Graphics::LUT;

my $win=pgwin(Device=>'/xs',Aspect=>1,WidowWidth=>10);


while(my $file=shift @ARGV) {
    print STDERR $file,"\n";
    open(FILE,"<$file") or warn "cannot open file `$file'  for reading.";
    my $buffer;
    read FILE,$buffer,8;
    my ($magic,$pixels_2)=unpack "a4I", $buffer;

    my $mat=zeroes($pixels_2,$pixels_2);
    print "#$magic $pixels_2\n";
    while(!(eof FILE)) {
	read FILE,$buffer,8;
	my ($lor_i,$lor_j,$count)=unpack "SSI",$buffer;

	for(my $i=0;$i<$count;++$i) {
	    read FILE, $buffer,8;
	    my ($pix_x,$pix_y,$pix_count)=unpack "SSI",$buffer;
	 #   print "$pix_x $pix_y $pix_count\n";
	    my $val=at $mat, ($pix_x,$pix_y);
	    set $mat, ($pix_x,$pix_y),$val+$pix_count;
	    $val=at $mat, ($pix_y,$pix_x);
	    set $mat, ($pix_y,$pix_x),$val+$pix_count;
	}
    }
    
    $win->ctab(lut_data('idl5'));
    $win->imag($mat,{DrawWedge=>1});

}
