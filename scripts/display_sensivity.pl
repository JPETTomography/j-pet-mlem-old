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

use PDL;

use PDL::Graphics::PGPLOT::Window;
use PDL::Graphics::LUT;

my $radius=0.450;
my $width =0.006;
my $height=0.020;


my $win=pgwin(Device=>'/cps',Aspect=>1,WindowWidth=>10,Recording=>1);


my $n_detectors=shift @ARGV;

while(my $file=shift @ARGV) {
    print STDERR $file,"\n";
    open(FILE,"<$file") or warn "cannot open file `$file'  for reading.";
    my $buffer;
    read FILE,$buffer,12;
    my ($magic,$pixels_2,$n_emissions)=unpack "a4II", $buffer;

    my $mat=zeroes($pixels_2,$pixels_2);
    print "#$magic $pixels_2 $n_emissions\n";
    while(!(eof FILE)) {
	read FILE,$buffer,8;
	my ($lor_i,$lor_j,$count)=unpack "SSI",$buffer;

	for(my $i=0;$i<$count;++$i) {
	    read FILE, $buffer,8;
	    my ($pix_x,$pix_y,$pix_count)=unpack "SSI",$buffer;
	    my $val=at $mat, ($pix_x,$pix_y);
	    set $mat, ($pix_x,$pix_y),$val+$pix_count;
	    $val=at $mat, ($pix_y,$pix_x);
	    set $mat, ($pix_y,$pix_x),$val+$pix_count;
	}
    }

    my $pixel_size=$radius/($pixels_2*sqrt(2.0));
    my $r_in=$radius/$pixel_size;
    my $r_out=$r_in+$height/$pixel_size;
    my $r_mid=($r_in+$r_out)/2;
    my $r_fov=$r_in/sqrt(2.0);

    my $pw=$width/$pixel_size;
    my $ph=$height/$pixel_size;


    $win->rectangle($r_out/2,$r_out/2,$r_out,$r_out,{FILLTYPE=>'OUTLINE'});
    $win->hold();
    $win->circle(0,0,$r_out,{FILLTYPE=>'SOLID',COLOR=>'Black'});

    $win->ctab(lut_data('idl5'));
    $win->imag($mat);
    $win->draw_wedge();
    my $pi=3.1415926;
    for(my $angle=0;$angle<=$pi/2.0+0.0001; $angle+=2*$pi/$n_detectors) {
	$win->rectangle($r_mid*cos($angle),$r_mid*sin($angle),
			$ph,$pw,$angle,{FILLTYPE=>'SOLID',COLOR=>'RED'});
    }

    $win->circle(0,0,$r_out,{FILLTYPE=>'OUTLINE'});
    $win->circle(0,0,$r_in,{FILLTYPE=>'OUTLINE'});
    $win->circle(0,0,$r_fov,{FILLTYPE=>'OUTLINE'});


    $win->release();

}
