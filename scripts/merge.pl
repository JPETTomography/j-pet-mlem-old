#!/usr/bin/perl -w 

use strict;


my @files = @ARGV;
my @file_data=();
my $max_id=0;
foreach my $file (@files) {
	open(FILE,"<$file") or die "Cannot open file `$file'";
	my %data = ();
	while(<FILE>) {
		my ($id, $time) = split(' ',$_,3);
	 	#print "$file $id $time<-\n";
		$data{$id} = $time;
			
		$max_id=$id if($max_id<$id);
		
	}
	push @file_data, \%data;
	close(FILE);
}

#print "$max_id\n";

for( my $id = 0;$id<=$max_id; $id++) {
	my $sid = sprintf("%05d",$id);
	
	my $line = "$sid ";
	my $valid = 1;
	foreach my $fd (@file_data) {
		if(exists($fd->{$sid}) ) {
			$line.="$fd->{$sid} ";
		} else {
			$valid=0;
		}
		
	}
	print "$line\n" if($valid);
}