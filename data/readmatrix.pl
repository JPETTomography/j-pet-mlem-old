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


while(my $file=shift @ARGV) {
    print STDERR $file,"\n";
    open(FILE,"<$file") or warn "cannot open file `$file'  for reading.";
    my $buffer;
    read FILE,$buffer,8;
    my ($magic,$pixels_2)=unpack "a4I", $buffer;
    print "#$magic $pixels_2\n";
    while(!(eof FILE)) {
	read FILE,$buffer,8;
	my ($lor_i,$lor_j,$count)=unpack "SSI",$buffer;
	print "$lor_i $lor_j $count\n";
	read FILE, $buffer,8*$count;
    }

}
