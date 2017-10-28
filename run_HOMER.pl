use warnings;
use strict;

my $genome = "hg19";
my $bed = "/path/to/merged_peaks_2kb_[comp ID]_[direction].bed";
my $output_folder = "/path/to/HOMER_motifs/merged_peaks_2kb_[comp ID]_[direction]";
my $extra_param = "-size given";

my $command = "findMotifsGenome.pl $bed $genome $output_folder $extra_param";
system($command);

#While this does seem to be the better way of annotating specific peaks, this only work for de-novo motifs (without corrected p-values)
my $de_novo_folder = "$output_folder/homerResults";

opendir DH, $de_novo_folder or die "Failed to open $de_novo_folder: $!";
my @files = readdir(DH);
foreach my $file (@files){
	if($file =~ /motif\d+.motif/){
		my ($motif_num) = ($file =~ /motif(\d+).motif/);
		print "$file|$motif_num\n";
		my $example_motif = "$de_novo_folder/$file";
		my $motif_hits = "$de_novo_folder/motif$motif_num.txt";
		$command = "annotatePeaks.pl $bed $genome -m $example_motif > $motif_hits";
		system($command);
	}#end if(-d ("$inputfolder/$folder"))
}#end foreach my $folder (@folders)

closedir(DH);

exit;
