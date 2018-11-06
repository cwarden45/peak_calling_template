use warnings;
use strict;
use Bio::SeqIO;

my $bed_file = "../Results/i-cisTarget_motif/input_files/[compID]_[direction].bed";
my $ref = "/path/to/ref.fa";
my $peak_seq = "[compID]_[direction].fa";

#if peaks must be same size, select set distance from middle of peak (unless set to "-1", then keep full sequence)
my $set_length = -1;

my $seqio_obj = Bio::SeqIO->new(-file => $ref, 
								-format => "fasta" );
my %seq_hash;

while (my $seq_obj = $seqio_obj->next_seq){
	my $chr_name = $seq_obj->id;
	my $chr_seq = $seq_obj->seq;
	$seq_hash{$chr_name}=$chr_seq;
}#end while (my $seq_obj = $seqio_obj->next_seq)

open(OUT, "> $peak_seq") || die("Could not open  $peak_seq!");
								
open(IN, $bed_file)|| die("Cannot open $bed_file!\n");

my $line_count = 0;

while(<IN>){
	my $line = $_;
	chomp $line;
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	
	$line_count++;
	
	if($line_count > 1){
		my @line_info = split("\t",$line);
		
		my $peak_chr = $line_info[0];
		my $peak_start = $line_info[1];
		my $peak_stop = $line_info[2];
		my $peak_name = "$peak_chr:$peak_start-$peak_stop";
		
		my $chr_seq = $seq_hash{$peak_chr};
		
		my $peak_seq = "";
		
		if($set_length == -1){
			$peak_seq = uc(substr($chr_seq, $peak_start-1,$peak_stop-$peak_start+1));
		}else{
			#print "Centered, Fixed-Length Sequence\n";
			my $peak_center = int(($peak_stop+$peak_start)/2);
			my $flank_distance = int($set_length/2);
			$peak_start=$peak_center-$flank_distance;
			$peak_stop=$peak_center+$flank_distance;
			$peak_seq = uc(substr($chr_seq, $peak_start-1,$peak_stop-$peak_start+1));
			#print "$peak_center:$flank_distance:$peak_start-$peak_stop\n";
			#print ">$peak_name\n$peak_seq\n";
		}#end else
		
		print OUT ">$peak_name\n$peak_seq\n";
	}#end if($line_count > 1)
}#end while(<IN>)

close(IN);

close(OUT);						

exit;