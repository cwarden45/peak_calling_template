use warnings;
use strict;

#have to first run extract_peak_sequences.pl

my $peakFA = "[compID]_[direction].fa";
my $dbRoot = "/path/to/meme_4.12.0/motif_databases";
my @motifDBs = ("$dbRoot/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme",
				"$dbRoot/EUKARYOTE/jolma2013.meme",
				"$dbRoot/EUKARYOTE/macisaac_theme.v1.meme",
				"$dbRoot/JASPAR/JASPAR_CORE_2016_vertebrates.meme");
my $outputfolder = "[compID]_[direction]";

my $command = "meme-chip -o $outputfolder -db ".join(" -db ",@motifDBs)." $peakFA";
system($command);

exit;