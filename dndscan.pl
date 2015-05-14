#!/usr/bin/perl

use strict;

my $input = "";
my $window = 3;
my $step   = 3;

my ($help, $input, $window, $step);

usage() if( @ARGV == 0);

GetOptions(
	'help|h'     => \$help,
	'input|i=s'  => \$input,
	'window|w=i' => \$window,
	'step|s=i'   => \$step
) or die usage();

die "ERROR: The window length should be divisible by 3\n" if ( $window % 3 != 0 );
die "ERROR: The step length should be divisible by 3\n" if ( $step % 3 != 0 );
die "ERROR: No input file\n" if($input == "");


my $LOG10 = log(10);

my ( $length, @list_seq ) = readFasta($input);

my %codons = (
	"TCA" => "S",    #Serine
	"TCC" => "S",    #Serine
	"TCG" => "S",    #Serine
	"TCT" => "S",    #Serine
	"TTC" => "F",    #Phenylalanine
	"TTT" => "F",    #Phenylalanine
	"TTA" => "L",    #Leucine
	"TTG" => "L",    #Leucine
	"TAC" => "Y",    #Tyrosine
	"TAT" => "Y",    #Tyrosine
	"TAA" => "*",    #Stop
	"TAG" => "*",    #Stop
	"TGC" => "C",    #Cysteine
	"TGT" => "C",    #Cysteine
	"TGA" => "*",    #Stop
	"TGG" => "W",    #Tryptophan
	"CTA" => "L",    #Leucine
	"CTC" => "L",    #Leucine
	"CTG" => "L",    #Leucine
	"CTT" => "L",    #Leucine
	"CCA" => "P",    #Proline
	"CCC" => "P",    #Proline
	"CCG" => "P",    #Proline
	"CCT" => "P",    #Proline
	"CAC" => "H",    #Histidine
	"CAT" => "H",    #Histidine
	"CAA" => "Q",    #Glutamine
	"CAG" => "Q",    #Glutamine
	"CGA" => "R",    #Arginine
	"CGC" => "R",    #Arginine
	"CGG" => "R",    #Arginine
	"CGT" => "R",    #Arginine
	"ATA" => "I",    #Isoleucine
	"ATC" => "I",    #Isoleucine
	"ATT" => "I",    #Isoleucine
	"ATG" => "M",    #Methionine
	"ACA" => "T",    #Threonine
	"ACC" => "T",    #Threonine
	"ACG" => "T",    #Threonine
	"ACT" => "T",    #Threonine
	"AAC" => "N",    #Asparagine
	"AAT" => "N",    #Asparagine
	"AAA" => "K",    #Lysine
	"AAG" => "K",    #Lysine
	"AGC" => "S",    #Serine
	"AGT" => "S",    #Serine
	"AGA" => "R",    #Arginine
	"AGG" => "R",    #Arginine
	"GTA" => "V",    #Valine
	"GTC" => "V",    #Valine
	"GTG" => "V",    #Valine
	"GTT" => "V",    #Valine
	"GCA" => "A",    #Alanine
	"GCC" => "A",    #Alanine
	"GCG" => "A",    #Alanine
	"GCT" => "A",    #Alanine
	"GAC" => "D",    #Aspartic Acid
	"GAT" => "D",    #Aspartic Acid
	"GAA" => "E",    #Glutamic Acid
	"GAG" => "E",    #Glutamic Acid
	"GGA" => "G",    #Glycine
	"GGC" => "G",    #Glycine
	"GGG" => "G",    #Glycine
	"GGT" => "G"     #Glycine
);

my @bases = ( 'A', 'C', 'G', 'T' );

my $nb_seq   = scalar(@list_seq);
my $doublets = $nb_seq * ( $nb_seq - 1 ) * 0.5;

print "pos,%id,%Nid,S,NS,S-rf1,NS-rf1,S/NS-rf1,TS,TV,TS/TV,"
  . "A-pos1,C-pos1,G-pos1,T-pos1,A-pos2,C-pos2,G-pos2,T-pos2,A-pos3,C-pos3,G-pos3,T-pos3\n";

my ( @NS, @S );

# process the alignment step by step
for ( my $i = 0 ; $i < $length ; $i += $step ) {
	my $id = 0;
	my ( $NS, $S, $TS, $TV ) = (0) x 4;
	my ( $NSn, $Sn, $SnNSn, $TSn, $TVn ) = (0) x 5;
	my @freq = ();
	
	# process each window
	for ( my $w = 0 ; $w < $window ; $w += 3 ) {
		
		# process each doublet
		for ( my $j = 0 ; $j < $nb_seq ; $j++ ) {
			my $codon1 = substr( $list_seq[$j]->{SEQ}, $i + $w, 3 );
			freq( $codon1, \@freq );
			for ( my $k = $j + 1 ; $k < $nb_seq ; $k++ ) {
				my $codon2 = substr( $list_seq[$k]->{SEQ}, $i + $w, 3 );
				baseStats( $codon1, $codon2, \$id, \$TV, \$TS );
				next if ( $codon1 eq $codon2 );
				permute( $codon1, $codon2, \$S, \$NS );
				permute( $codon2, $codon1, \$S, \$NS );
			}
		}
	}
	my $idp  = ( $id / ( $doublets * $window ) ) * 100;
	my $nidp = 100 - $idp;

	# Normaliza NS and S scores
	if ( $idp != 100 ) {
		( $NSn, $Sn )  = normalize( $NS, $S,  $nidp );
		( $TSn, $TVn ) = normalize( $TS, $TV, $nidp );
	}
	
	# Keep track of the NS and S scores for the Zipf test
	push( @NS, $NSn );
	push( @S,  $Sn );

	# format results
	$idp  = sprintf( "%.2f", $idp );
	$nidp = sprintf( "%.2f", $nidp );

	eval { $SnNSn = sprintf( "%.2f", $Sn / $NSn ); };
	$SnNSn = "NaN" if $@;

	my $TSnTVn = "NaN";
	eval { $TSnTVn = sprintf( "%.2f", $TSn / $TVn ); };

	my $Sn  = sprintf( "%.2f", $Sn );
	my $NSn = sprintf( "%.2f", $NSn );
	my $TSn = sprintf( "%.2f", $TSn );
	my $TVn = sprintf( "%.2f", $TVn );

	# print results
	print eval( $i + 1 ) . ",$idp,$nidp,$S,$NS,$Sn,$NSn,$SnNSn,$TSn,$TVn,$TSnTVn";

	for ( my $i = 0 ; $i < 3 ; $i++ ) {
		map {
			print( defined $freq[$i]->{$_} ? sprintf( ",%.2f",( $freq[$i]->{$_} * 100 ) / ( $nb_seq * ( $window / 3 ) ) ): ",0"	);
		} @bases;
	}
	print ",\n";
}

zipf( \@S, \@NS );

#####################################
# SUBROUTINES
####################################

# Compute base frequencies for 1 codon
sub freq($\@) {
	my ( $seqa, $freq ) = @_;
	my @seq = split( //, $seqa );
	for ( my $i = 0 ; $i < 3 ; $i++ ) {
		$$freq[$i]{ $seq[$i] }++;
	}
}

# Do permuations of 2 codons
sub permute($$\$\$) {
	my ( $seqa, $seqb, $S, $NS ) = @_;
	for ( my $pos = 0 ; $pos < 3 ; $pos++ ) {
		my @seq1 = split( //, $seqa );
		my @seq2 = split( //, $seqb );
		next if ( $seq1[$pos] eq $seq2[$pos] );
		my $aa = $codons{ join( "", @seq2 ) };
		$seq2[$pos] = $seq1[$pos];
		my $newaa = $codons{ join( "", @seq2 ) };
		if ( $newaa eq $aa ) { $$S++; }
		else { $$NS++; }
	}
}

# Compute base stats: TS, TV and homology between 2 codons
sub baseStats($$\$) {
	my ( $seqa, $seqb, $id, $TV, $TS ) = @_;
	my @seq1 = split( //, $seqa );
	my @seq2 = split( //, $seqb );
	for ( my $i = 0 ; $i < 3 ; $i++ ) {
		if($seq1[$i] =~ /[actg]/i and $seq2[$i] =~ /[actg]/i){
			if ( $seq1[$i] ne $seq2[$i] ) {
				if ( isTransition( $seq1[$i], $seq2[$i] ) ) { $$TS++; }
				else { $$TV++; }
			}
			else {
				$$id++;
			}
		}
	}
}

# Test if it is a transition
sub isTransition($$) {
	my ( $a, $b ) = @_;
	if (   ( $a eq 'A' and $b eq 'G' )
		or ( $a eq 'G' and $b eq 'A' )
		or ( $a eq 'C' and $b eq 'T' )
		or ( $a eq 'T' and $b eq 'C' ) )
	{
		return 1;
	}
	return 0;
}

# Normalize S and NS scores
sub normalize($$$) {
	my ( $a, $b, $tot ) = @_;
	my $an = ( $tot * $a ) / ( $a + $b );
	my $bn = $tot - $an;
	return ( $an, $bn );
}

# Compute the NS, S, log NS and log S bins for the Zipf test
sub zipf(\@\@) {
	my ( $Sn, $NSn ) = @_;
	my $max = @$Sn[0];
	map { $max = $_ if ( $_ > $max ) } @$Sn;
	map { $max = $_ if ( $_ > $max ) } @$NSn;
	my $nBins      = int( $max + 0.5 ) + 1;
	my $numEntries = scalar(@$Sn);
	my $binMin     = 0;

	my @NSbinCount = (0) x ( $nBins + 1 );
	my @SbinCount  = (0) x ( $nBins + 1 );

	for ( my $i = 0 ; $i < $numEntries ; $i++ ) {
		my $NSidx = int( 1.0 * ( @$NSn[$i] ) * $nBins / $max );
		my $Sidx  = int( 1.0 * ( @$Sn[$i] ) * $nBins / $max );
		$NSbinCount[$NSidx]++ unless $NSidx < 0;
		$SbinCount[$Sidx]++   unless $Sidx < 0;
	}

	# print the data
	print "\nClass,S-rf1,logS-rf1,NS-rf1,logNS-rf1\n";
	for my $i ( 0 .. $nBins ) {
		my $logS  = 0;
		my $logNS = 0;
		if ( $SbinCount[$i] == 1 ) {$logS = 0.3;}
		else {
			eval {
				$logS = log( $SbinCount[$i] ) / $LOG10;
				$logS = sprintf( "%.2f", $logS );
			};
		}

		if ( $NSbinCount[$i] == 1 ) {$logNS = 0.3;}
		else {
			eval {
				$logNS = log( $NSbinCount[$i] ) / $LOG10;
				$logNS = sprintf( "%.2f", $logNS );
			};
		}
		print "$i,$SbinCount[$i],$logS,$NSbinCount[$i],$logNS\n";
	}
}

# Read FASTA file
sub readFasta($) {
	my ($input) = @_;
	my $rec     = {};
	my $length  = 0;
	open( IN, $input ) or die "Cannot open $input $!\n";
	while (<IN>) {
		s/\r?\n$//;
		next if (/^\s*$/);
		if (/^>/) {
			if ( $rec->{SEQ} ) {
				die "Sequences of different length\n" if ( length( $rec->{SEQ} ) != $length and $length != 0 );
				$length = length( $rec->{SEQ} );
				$rec->{SEQ} =~ tr/U/T/;
				push( @list_seq, $rec );
				$rec = {};
			}
			$rec->{NAME} = substr( $_, 1 );
		}
		else {
			$rec->{SEQ} .= uc($_);
		}
	}
	die "Sequences of different length\n" if ( length( $rec->{SEQ} ) != $length );
	push( @list_seq, $rec );
	return ( $length, @list_seq );
}

sub usage{
	print "\nperl dndscan.pl - Perl implemention of the DnDscan program described in: Gibbs MJ, Wayper P, Fourment ML, Wood JT, Ohshima K, Armstrong JS, Gibbs AJ. The variable codons of H3 influenza A virus haemagglutinin genes. Arch Virol. 2007 \n\n";
	print "Usage: perl dndscan.pl [-ws] [file]\n";
	print "Example: perl dndscan.pl -w 3 -s 3 file.fa\n\n";
	exit 0;
}
