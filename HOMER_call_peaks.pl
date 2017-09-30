use warnings;
use strict;

my %finished_samples = ();
my $model_type = "histone";
my $alignment_folder = "../Alignment_Folder";
my $bam_suffix = "\_filtered.bam";
my $tag_folder = "Tag_Directories";
my $peak_folder = "Peaks";
my $HOMER = "/path/to/HOMER/bin/";

opendir DH, $alignment_folder or die "Failed to open $alignment_folder: $!";
my @files= readdir(DH);

foreach my $file (@files){
	if($file =~ /$bam_suffix$/){
		my ($sampleID) = ($file =~ /(.*)$bam_suffix$/);
		print "$sampleID\n";
		my $bam = "$alignment_folder/$sampleID\_filtered.bam";

		my $sample_tags = "$tag_folder/$sampleID";

		my $command = "makeTagDirectory $sample_tags $bam -format sam ";
		system($command);

		my $raw_peaks = "$peak_folder/$sampleID.peaks";
		$command = "findPeaks $sample_tags -style $model_type -o $raw_peaks";
		system($command);

		my $peak_bed = "$peak_folder/$sampleID.bed";

		open(IN,$raw_peaks)||die("Cannot open $raw_peaks\n");
		open(OUT,"> $peak_bed")||die("Cannot open $peak_bed\n");

		while(<IN>){
			my $line = $_;
			chomp $line;
			$line =~ s/\r//g;
			$line =~ s/\n//g;
			
			if(!($line =~ /^#/)){
				my @line_info = split("\t",$line);
				print OUT "$line_info[1]\t$line_info[2]\t$line_info[3]\n";
			}#end if(!($line =~ /^#/))
		}#end while(<IN>)

		close(OUT);
		close(IN);
	}#end if($file =~ /$bam_suffix$/)
}#end foreach my $file (@files)
closedir(DH);

exit;