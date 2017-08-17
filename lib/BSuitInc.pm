package main;
use strict;
use warnings;
use Digest::SHA;
use IPC::Open2;
use Galaxy::Data;

sub getFilesHash(@) {
	my $fileStr = join( ',', @_ );
	return Digest::SHA::sha1_base64($fileStr);
}

sub getRef2char($$) {
	my ( $HostRefName, $VirusRefName ) = @_;
	my $HostChar  = substr $HostRefName,  0, 1;
	my $VirusChar = substr $VirusRefName, 0, 1;
	my $tlen      = length($VirusRefName) - 1;
	my $i         = 0;
	while ( ( $VirusChar eq $HostChar ) and ( ++$i <= $tlen ) ) {
		my $tmpChar = substr $VirusRefName, $i, 1;
		$VirusChar = $tmpChar if $tmpChar =~ /\w/;
		warn "$i - $tmpChar $VirusChar -\n";
	}
	$VirusChar = 'V' if $VirusChar eq $HostChar;
	return $HostChar . $VirusChar;
}

sub sortChrPos($$) {
	my ( $ChrA, $PosA ) = split /\t/, $_[0];
	my ( $ChrB, $PosB ) = split /\t/, $_[1];
	if ( $ChrA eq $ChrB ) {
		return $PosA <=> $PosB;
	}
	return 1  if exists $main::VirusChrIDs{$ChrA};
	return -1 if exists $main::VirusChrIDs{$ChrB};
	$ChrA cmp $ChrB
	  || $PosA <=> $PosB;
}

sub cigar2poses($) {
	my ($cigar) = @_;
	my @cigar = $cigar =~ /(\d+)(\w)/g;
	my ( $reflen, $maxM, $readlen ) = ( 0, 0, 0 );
	while (@cigar) {
		my ( $len, $op ) = splice( @cigar, 0, 2 );
		if ( $op eq 'M' ) {
			$reflen  += $len;
			$readlen += $len;
			$maxM = $len if $maxM < $len;
			next;
		}
		$reflen  += $len if $op eq 'D';
		$readlen += $len if $op eq 'I';
	}
	return ( $reflen, $readlen );
}

sub mergeIn($$$$$$) {
	my ( $isHost, $rChrRange, $aRead, $rStore, $cid, $i ) = @_;
	my ( $reflen, $readlen ) = cigar2poses( $aRead->[5] );
	my $thisehPos = $aRead->[3] + $reflen;
	my $ret;
	if ( $i > 1 or $isHost == 0 or scalar( keys %{$rChrRange} ) == 0 ) {
		$ret = 0;
		if ( keys %{$rChrRange} and exists $rChrRange->{ $aRead->[2] } ) {
			if ( $aRead->[3] <= $rChrRange->{ $aRead->[2] }->[0] ) {
				$rChrRange->{ $aRead->[2] }->[0] = $aRead->[3];
			}
			if ( $thisehPos >= $rChrRange->{ $aRead->[2] }->[1] ) {
				$rChrRange->{ $aRead->[2] }->[1] = $thisehPos;
			}
		}
		else {
			$rChrRange->{ $aRead->[2] } = [ $aRead->[3], $thisehPos, 0 ];
		}
	}
	elsif ( exists $rChrRange->{ $aRead->[2] } ) {
		if ( $aRead->[3] <= $rChrRange->{ $aRead->[2] }->[1] ) {
			$rChrRange->{ $aRead->[2] }->[1] = $thisehPos if $thisehPos >= $rChrRange->{ $aRead->[2] }->[1];
			$ret = 0;
		}
		else {
			$ret = 1;
		}
	}
	else {
		$ret = 1;
	}
	unless ($ret) {
		my $tid = join( "\n", $cid, $i );
		push @{$rStore}, $tid;
		++$rChrRange->{ $aRead->[2] }->[2];
	}
	return $ret;
}

sub formatChrRange($) {
	my ($rChrRange) = @_;
	my @ret;
	for my $c ( sort keys %{$rChrRange} ) {
		push @ret, "${c}:" . $rChrRange->{$c}->[0] . '-' . $rChrRange->{$c}->[1] . ':' . $rChrRange->{$c}->[2];
	}
	return join( ',', @ret );
}

sub revcom($) {
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	my $rev = reverse $str;
	$rev =~ tr/[](){}<>/][)(}{></;
	return $rev;
}

sub guessMethyl($) {
	my ($seq) = @_;
	my %BaseCnt = (
		A => 0,
		G => 0,
		C => 0,
		T => 0,
		N => 0,
	);
	for ( split //, $seq ) {
		++$BaseCnt{ uc $_ };
	}
	my $seqlen = length($seq) - $BaseCnt{'N'};
	return 'N' if $seqlen == 0;
	my @Cnts = sort { $BaseCnt{ uc $b } <=> $BaseCnt{ uc $a } } keys %BaseCnt;
	if ( $BaseCnt{'C'} <= $seqlen * $main::methly3BaseErrRate and $BaseCnt{'T'} > 0 ) {
		return '1CT';
	}
	elsif ( $BaseCnt{'G'} <= $seqlen * $main::methly3BaseErrRate and $BaseCnt{'A'} > 0 ) {
		return '2GA';
	}
	else {
		return '0Raw';
	}
}

my ( $match, $methylmatch, $mismatch, $INDEL ) = ( 2, 1, -3, -4 );
our ( %MatrixO, %MatrixFct, %MatrixRga );
my @Bases = sort qw(A C G T);
for my $i ( 0 .. $#Bases ) {
	for my $j ( $i .. $#Bases ) {
		my $str = $Bases[$i] . $Bases[$j];
		if ( $i == $j ) {
			$MatrixO{$str} = $MatrixFct{$str} = $MatrixRga{$str} = $match;
		}
		else {
			$MatrixO{$str} = $mismatch;
			if ( $str =~ /^[CT]+$/ ) {
				$MatrixFct{$str} = $methylmatch;
			}
			else {
				$MatrixFct{$str} = $mismatch;
			}
			if ( $str =~ /^[GA]+$/ ) {
				$MatrixRga{$str} = $methylmatch;
			}
			else {
				$MatrixRga{$str} = $mismatch;
			}
		}
	}
}

sub doAlign($$$) {
	my ( $AssemHRef, $retHostARef, $retVirusARef ) = @_;
	my $retHost  = $$retHostARef[0];
	my $retVirus = $$retVirusARef[0];
	my @froDat   = sort keys %{$AssemHRef};
	my $fro0     = shift @froDat;
	my $result   = [ $AssemHRef->{$fro0}->[0]->[1] ];
	for my $fro (@froDat) {
		$result = mergeAln( $result->[0], $AssemHRef->{$fro}->[0]->[1] );
	}
	my $HostResult = mergeAln( $retHost->[1], $result->[0] );
	my @RefInfo = split /[:-]/, $retHost->[0];
	my ( $RefCut, $Left, $Insert, @Rcuts, @Vcuts ) = (0);
	if ( $HostResult->[1] =~ /^(D*)([MmR]+)(I+)/ ) {
		( undef, $Left, $Insert ) = ( $1, $2, $3 );
		push @Rcuts, ( $RefInfo[1] + length($Left) + 1 );
	}
	if ( $HostResult->[1] =~ /^(I*)([MmR]+)(D+)/ ) {
		( $Left, $Insert, undef ) = ( $1, $2, $3 );
		push @Rcuts, ( $RefInfo[1] + length($Left) + length($Insert) + 1 );
	}
	$RefCut = join( ';', @Rcuts );
	my $VirusResult = mergeAln( $retVirus->[1], $result->[0] );
	my @VirusInfo = split /[:-]/, $retVirus->[0];
	my ( $VirCut, $VirLeft, $VirInsert ) = ( 0, 0, 0 );
	if ( $VirusResult->[1] =~ /^(D*)([MmR]+)(I+)/ ) {
		( undef, $VirLeft, $VirInsert ) = ( $1, $2, $3 );
		push @Vcuts, ( $VirusInfo[1] + length($VirLeft) + 1 );
	}
	if ( $VirusResult->[1] =~ /^(I*)([MmR]+)(D+)/ ) {
		( $Left, $Insert, undef ) = ( $1, $2, $3 );
		push @Vcuts, ( $VirusInfo[1] + length($Left) + length($Insert) + 1 );
	}
	$VirCut = join( ';', @Vcuts );
	return [ $RefInfo[0], $RefCut, $VirusInfo[0], $VirCut, length($VirInsert) ];
}

sub dynAln($$$) {
	my ( $ref, $query, $MatrixR ) = @_;
	my @a = ( '', split //, $ref );
	my @b = ( '', split //, $query );
	my ( @matrix, @path );
	my ( $i, $j, $sc );
	$matrix[$_][0] = 0 for ( 0 .. $#a );
	$matrix[0][$_] = 0 for ( 0 .. $#b );
	$path[$#a][$#b] = 0;

	for ( $i = 1 ; $i <= $#a ; $i++ ) {

		for ( $j = 1 ; $j <= $#b ; $j++ ) {
			my $str = join( '', sort ( $a[$i], $b[$j] ) );
			$matrix[$i][$j] = $MatrixR->{$str} + $matrix[ $i - 1 ][ $j - 1 ];
			$path[$i][$j]   = 1;
			$sc             = $matrix[ $i - 1 ][$j] + $INDEL;
			if ( $sc > $matrix[$i][$j] ) {
				$matrix[$i][$j] = $sc;
				$path[$i][$j]   = 2;
			}
			$sc = $matrix[$i][ $j - 1 ] + $INDEL;
			if ( $sc > $matrix[$i][$j] ) {
				$matrix[$i][$j] = $sc;
				$path[$i][$j]   = 3;
			}
			if ( $matrix[$i][$j] < 0 ) {
				$matrix[$i][$j] = 0;
				$path[$i][$j]   = 0;
			}
		}
	}
	my ( $count, @ax, @ay, @as ) = (0);
	--$i;
	--$j;
	my ( $len, $mv, $mi, $mj ) = (0);
	while ( $i >= 0 and $j >= 0 ) {
		( $mv, $mi, $mj ) = @{ &getmax( \@matrix, $i, $j ) };
		while ( $mi >= 0 and $mj >= 0 ) {
			$as[$count] = $matrix[$mi][$mj] unless $as[$count];
			++$len;
			my $path = $path[$mi][$mj];
			if ( !defined $path ) {
				--$len;
				if ( $len > 0 ) {
					@{ $ax[$count] } = reverse @{ $ax[$count] };
					@{ $ay[$count] } = reverse @{ $ay[$count] };
					++$count;
					$len = 0;
				}
				--$mi;
				--$mj;
				$i = $mi;
				$j = $mj;
				( $mv, $mi, $mj ) = @{ &getmax( \@matrix, $i, $j ) };
			}
			elsif ( $path == 1 ) {
				push @{ $ax[$count] }, $a[$mi];
				push @{ $ay[$count] }, $b[$mj];
				--$mi;
				--$mj;
			}
			elsif ( $path == 2 ) {
				push @{ $ax[$count] }, $a[$mi];
				push @{ $ay[$count] }, "-";
				--$mi;
			}
			elsif ( $path == 3 ) {
				push @{ $ax[$count] }, "-";
				push @{ $ay[$count] }, $b[$mj];
				--$mj;
			}
			else {
				die;
			}
		}
	}
	for ( $i = 1 ; $i <= $count ; $i++ ) {
		print "No. $i:\n";
		print 'X:[', @{ $ax[ $i - 1 ] }, "]\n", 'Y:[', @{ $ay[ $i - 1 ] }, "]\n", 'Score:', $as[ $i - 1 ], "\n";
	}
}

sub getmax($$$) {
	my ( $arref, $ii, $ij ) = @_;
	my ( $mv, $mi, $mj, $i, $j ) = ( -1, -1, -1 );
	for ( $i = 1 ; $i <= $ii ; $i++ ) {
		for ( $j = 1 ; $j <= $ij ; $j++ ) {
			if ( $mv == $$arref[$i][$j] ) {
				if ( ( $mi + $mj ) < ( $i + $j ) ) {
					$mi = $i;
					$mj = $j;
				}
			}
			elsif ( $mv < $$arref[$i][$j] ) {
				$mv = $$arref[$i][$j];
				$mi = $i;
				$mj = $j;
			}
		}
	}
	return [ $mv, $mi, $mj ];
}

sub getMaxM(@) {
	my %Dat;
	for my $c (@_) {
		my @cigar = $c =~ /(\d+)(\w)/g;
		my $maxM = 0;
		while (@cigar) {
			my ( $len, $op ) = splice( @cigar, 0, 2 );
			if ( $op eq 'M' ) {
				$maxM = $len if $maxM < $len;
			}
		}
		push @{ $Dat{$maxM} }, $c;
	}
	my @keys = sort { $b <=> $a } keys %Dat;
	return $Dat{ $keys[0] };
}

sub doAln($$$) {
	my ( $refile, $query, $dir ) = @_;
	my $str    = "$RealBin/bin/water $refile <(echo $query) -scircular2 -gapopen 10 -gapextend 0.5 stdout -aformat3 sam";
	my $revstr = "$RealBin/bin/water $refile <(echo $query) -sreverse1 -scircular2 -gapopen 10 -gapextend 0.5 stdout -aformat3 sam";
	my %Result;
	my @Ret = `bash -c \'$str\' 2>/dev/null`;
	my @tmp = split /\t/, $Ret[2];
	$tmp[0] = '+';
	push @{ $Result{ $tmp[5] } }, [@tmp];
	@Ret    = `bash -c \'$revstr\' 2>/dev/null`;
	@tmp    = split /\t/, $Ret[2];
	$tmp[0] = '-';
	push @{ $Result{ $tmp[5] } }, [@tmp];
	my $t = getMaxM( keys %Result );
	my $k = $t->[0];
	my $ret;

	if ( $dir == 1 ) {
		if ( ( ( $Result{$k}->[0][0] eq '+' and $k =~ /^(\d+)M/ ) or ( $Result{$k}->[0][0] eq '-' and $k =~ /(\d+)M$/ ) ) and $1 >= $main::minVirMapLen ) {
			my $p = $Result{$k}->[0][3];
			$ret = [ $Result{$k}->[0][0], $p, $1 ];
		}
	}
	else {
		if ( ( ( $Result{$k}->[0][0] eq '+' and $k =~ /(\d+)M$/ ) or ( $Result{$k}->[0][0] eq '-' and $k =~ /^(\d+)M/ ) ) and $1 >= $main::minVirMapLen ) {
			my $p = $Result{$k}->[0][3];
			if ( $Result{$k}->[0][0] eq '+' ) {
				my $rlen = bam_cigar2rlen($k);
				$p += $rlen;
			}
			elsif ( $Result{$k}->[0][0] eq '-' ) {
				my $rlen = bam_cigar2rlen($k);
				$p -= $rlen;
			}
			$ret = [ $Result{$k}->[0][0], $p, $1 ];
		}
	}
	return $ret;
}

sub mergeAln($$) {
	my ( $ref, $query ) = @_;
	my $pid = open2( \*READER, \*WRITER, "$RealBin/bin/alnmethly" );
	WRITER->autoflush();
	my @Dat = ( [ $query, undef ], [ revcom($query), undef ] );
	$_->[1] = guessMethyl( $_->[0] ) for @Dat;
	@Dat = sort { $a->[1] cmp $b->[1] } @Dat;
	unshift @Dat, [ $ref, guessMethyl($ref) ];
	for (@Dat) {

		if ( $_->[1] eq '1CT' ) {
			$_->[0] =~ s/[CT]/Y/ig;
		}
		elsif ( $_->[1] eq '2GA' ) {
			$_->[0] =~ s/[GA]/R/ig;
		}
		$_->[0] = 'N' if $_->[1] eq 'N';
	}
	print WRITER join( "\n", map { $_->[0] } @Dat ), "\n";
	my %Result;
	while (<READER>) {
		chomp;
		if (/^Path(\d): ([IDMmR]+),(\d+)$/) {
			$Result{$1} = [ $2, $3 ];
		}
	}
	waitpid( $pid, 0 );
	close READER;
	my @Resu = sort { $Result{$b}->[1] <=> $Result{$a}->[1] } keys %Result;
	my @ResDat   = split //, $Result{ $Resu[0] }->[0];
	my @QuaryDat = split //, $Dat[ $Resu[0] ]->[0];
	my @Refdat   = split //, $ref;
	my @AlnDat;
	my ( $Refp, $Querp ) = ( 0, 0 );

	for my $p ( 0 .. $#ResDat ) {
		if ( $ResDat[$p] =~ /^[DMmR]$/ ) {
			++$Querp if $ResDat[$p] eq 'D';
			my $theBases = $IUB{ $Refdat[ $p - $Refp ] } or die;
			push @$theBases, @$theBases if $#$theBases == 0;
			++$AlnDat[$p]->{$_} for @$theBases;
		}
		if ( $ResDat[$p] =~ /^[MmR]$/ ) {
			my $theBases = $IUB{ $QuaryDat[ $p - $Querp ] };
			push @$theBases, @$theBases if $#$theBases == 0;
			++$AlnDat[$p]->{$_} for @$theBases;
		}
		elsif ( $ResDat[$p] eq 'I' ) {
			++$Refp;
			my $theBases = $IUB{ $QuaryDat[ $p - $Querp ] };
			push @$theBases, @$theBases if $#$theBases == 0;
			++$AlnDat[$p]->{$_} for @$theBases;
		}
	}
	my $reseq;
	for (@AlnDat) {
		my %t = %{$_};
		my @k = sort { $t{$a} <=> $t{$b} } keys %t;
		my $iub;
		if ( @k > 1 and $t{ $k[0] } == $t{ $k[1] } ) {
			$iub = join( '', sort @k[ 0, 1 ] );
		}
		else {
			$iub = $k[0];
		}
		$reseq .= $REV_IUB{$iub};
	}
	return [ $reseq, $Result{ $Resu[0] }->[0] ];
}

sub warnFileExist(@) {
	my %NotFound;
	for (@_) {
		next if /^\s*$/;
		++$NotFound{$_} unless -f $_;
	}
	my @NF = sort keys %NotFound;
	if ( @NF > 0 ) {
		warn "[!!!] File NOT Found:[", join( '],[', @NF ), "]\n";
	}
	return join( ' ', @_ );
}

sub bam_cigar2rlen($) {
	my @cigar = $_[0] =~ /(\d+\w)/g;
	my $pos = 0;
	for (@cigar) {
		if (/(\d+)([MDN=X])/) {
			$pos += $1;
		}
	}
	return $pos;
}

sub bam_cigar2qlen($) {
	my @cigar = $_[0] =~ /(\d+\w)/g;
	my $pos = 0;
	for (@cigar) {
		if (/(\d+)([MIS=X])/) {
			$pos += $1;
		}
	}
	return $pos;
}

sub getSeqCIGAR($$) {
	my ( $minLeft, $i ) = @_;
	my $firstSC   = 0;
	my $deltraPos = $i->[3] - $minLeft;
	my @cigar     = $i->[5] =~ /(\d+[MIDNSHP=XB])/g;
	my $pos       = $i->[3];
	my $offset    = $pos - $minLeft;
	my ( $cursorQ, $cursorR ) = ( 0, 0 );
	my $seqCIGAR = 'B' x $offset;

	for (@cigar) {
		/(\d+)([MIDNSHP=XB])/ or die;
		if ( $2 eq 'S' ) {
			if ( $cursorQ == 0 ) {
				substr $seqCIGAR, -$1, $1, $2 x $1;
				$firstSC = $1;
			}
			else {
				$cursorQ += $1;
				$seqCIGAR .= $2 x $1;
			}
		}
		elsif ( $2 eq 'M' ) {
			$cursorQ += $1;
			$seqCIGAR .= $2 x $1;
			for my $p ( ( $cursorQ - $1 ) .. ( $cursorQ - 1 ) ) {
				my $refpos = $pos + $cursorR + $p;
			}
			$cursorR += $1;
		}
		elsif ( $2 eq 'I' ) {
			$cursorQ += $1;
			$seqCIGAR .= $2 x $1;
		}
		elsif ( $2 eq 'D' ) {
			$cursorR += $1;
		}
		else { die; }
	}
	return ( $firstSC, $seqCIGAR );
}

sub grepmerge($$) {
	my $DEBGUHERE = 0;
	my ( $minLeft, $maxLS, @clipReads, %relPoses, %relPosesFR, $chr ) = ( 5000000000, 0 );
	my %Results;
	for my $i ( @{ $_[0] } ) {
		my @cigar = $i->[5] =~ /(\d+[MIDNSHP=XB])/g;
		my $flag = 0;
		for (@cigar) {
			if (/(\d+)S/) {
				$flag = 1 if $1 >= $main::minSoftClip;
			}
		}
		if ($flag) {
			unless ( defined $chr ) {
				$chr = $i->[2];
			}
			else {
				next if $chr ne $i->[2];
			}
			push @clipReads, $i;
			my $left = $i->[3];
			if ( $cigar[0] =~ /(\d+)S/ ) {
				$maxLS = $1 if $maxLS < $1;
				$left -= $1;
			}
			$minLeft = $left if $minLeft > $left;
		}
	}
	return \%Results unless @clipReads;
	my @ReadsCIGAR;
	print "$minLeft\n" if $DEBGUHERE;
	for my $i (@clipReads) {
		my ( $firstSC, $seqCIGAR ) = getSeqCIGAR( $minLeft, $i );
		my $offset = $i->[3] - $minLeft;
		my @Poses = getDeriv( "M", "B", $seqCIGAR );
		for (@Poses) {
			++$relPosesFR{$_};
			if ( $_ >= 0 ) {
				++$relPoses{ $_ + 0.5 };
			}
			else {
				++$relPoses{ -$_ - 0.5 };
			}
		}
		if ($DEBGUHERE) {
			my @thePoses = sort { $relPoses{$b} <=> $relPoses{$a} } keys %relPoses;
			my @absPoses;
			for ( @Poses, @thePoses ) {
				if ( $_ < 0 ) {
					push @absPoses, -( $minLeft - $_ );
				}
				else {
					push @absPoses, ( $minLeft + $_ );
				}
			}
			print "$i->[5]\t$offset\t@{$i}[0..4]; @Poses, @thePoses -> @absPoses\n$seqCIGAR\n";
			print 'B' x ( $offset - $firstSC ), $i->[9], "\n";
			my $tmpstr = ' ' x length($seqCIGAR);
			for my $p (@Poses) {
				if ( $p >= 0 ) {
					substr $tmpstr, $p, 1, ']';
				}
				else {
					substr $tmpstr, -$p, 1, '[';
				}
			}
			print $tmpstr, "\n";
		}
		push @ReadsCIGAR, $seqCIGAR;
	}
	my @thePosesA = sort { $relPoses{$b} <=> $relPoses{$a} } keys %relPoses;
	my $thePos = int( $thePosesA[0] );
	my @usingPoses;
	push @usingPoses, $thePos      if exists $relPosesFR{$thePos};
	push @usingPoses, -1 - $thePos if exists $relPosesFR{ -1 - $thePos };
	if ( @usingPoses == 1 ) {
		if ( $usingPoses[0] < 0 ) {
			for my $p ( ( $thePos - $main::posAround ) .. $thePos ) {
				if ( exists $relPosesFR{$p} ) {
					push @usingPoses, $p;
				}
			}
			if ( @usingPoses > 2 ) {
				my @t = @usingPoses;
				shift @t;
				@t = sort { $relPosesFR{$b} <=> $relPosesFR{$a} || $b <=> $a } @t;
				@usingPoses = ( $usingPoses[0], $t[0] );
			}
		}
		elsif ( $usingPoses[0] > 0 ) {
			for my $p ( ( $thePos + 1 ) .. ( $thePos + $main::posAround + 1 ) ) {
				if ( exists $relPosesFR{ -$p } ) {
					push @usingPoses, -$p;
				}
			}
			if ( @usingPoses > 2 ) {
				my @t = @usingPoses;
				shift @t;
				@t = sort { $relPosesFR{$b} <=> $relPosesFR{$a} || $a <=> $b } @t;
				@usingPoses = ( $usingPoses[0], $t[0] );
			}
		}
	}
	my %Bases;
	for my $i (@clipReads) {
		my $mtype;
		if ( $_[1] eq 'bwa-meth' ) {
			my ($YC) = grep /^YC:Z:/, @$i;
			if ( $i->[1] & 16 ) {
				if ( $YC eq 'YC:Z:CT' ) {
					$mtype = 'GA';
				}
				elsif ( $YC eq 'YC:Z:GA' ) {
					$mtype = 'CT';
				}
				else { die; }
			}
			else {
				if ( $YC eq 'YC:Z:CT' ) {
					$mtype = 'CT';
				}
				elsif ( $YC eq 'YC:Z:GA' ) {
					$mtype = 'GA';
				}
				else { die; }
			}
		}
		elsif ( $_[1] eq 'BSseeker2' ) {
			my ($XO) = grep /^XO:Z:/, @$i;
			if ( $XO eq 'XO:Z:+FR' ) {
				$mtype = 'CT';
			}
			elsif ( $XO eq 'XO:Z:-FR' ) {
				$mtype = 'GA';
			}
			elsif ( $XO eq 'XO:Z:+RF' ) {
				$mtype = 'GA';
			}
			elsif ( $XO eq 'XO:Z:-RF' ) {
				$mtype = 'CT';
			}
			else { die; }
		}
		elsif ( $_[1] eq 'bwa' ) {
			$mtype = '-';
		}
		my ( $firstSC, $seqCIGAR ) = getSeqCIGAR( $minLeft, $i );
		my $offset = $i->[3] - $minLeft - $firstSC;
		for my $p (@usingPoses) {
			my ( $tlen, $tmp, $vseq, $vqual, $t ) = (0);
			my $readlen = length $i->[9];
			if ( $_[1] eq 'BSseeker2' and $i->[10] eq '*' ) {
				$i->[10] = '@' x $readlen;
			}
			if ( $p > 0 and $p < length($seqCIGAR) ) {
				$tmp = substr( $seqCIGAR, $p + 1 ) or next;

				$tmp =~ /^(S+)/;
				if ( defined $1 ) {
					$tlen = length $1;
					$t    = $p + 1 - $offset;
					my $tl = $readlen - $t;
					if ( abs($t) >= $readlen ) {
						print STDERR ']';
						next;
					}
					$vseq  = substr $i->[9],  $t, $tlen;
					$vqual = substr $i->[10], $t, $tlen;
				}
			}
			elsif ( $p < 0 and ( -$p - 1 ) < length($seqCIGAR) ) {
				$tmp = substr( $seqCIGAR, 0, -$p - 1 ) or next;
				$tmp =~ /(S+)$/;
				if ( defined $1 ) {
					$tlen = length $1;
					$t    = -$p - 1 - $tlen - $offset;
					if ( abs($t) >= $readlen ) {
						print STDERR '[';
						next;
					}
					$vseq  = substr $i->[9],  $t, $tlen;
					$vqual = substr $i->[10], $t, $tlen;
				}
			}
			if ( $vseq and length($vseq) >= $main::minVirLen ) {
				push @{ $Bases{$p} }, [ $vseq, $vqual, $mtype ];
			}
		}
	}
	ddx \%relPoses, \%relPosesFR, \@usingPoses, \%Bases if $DEBGUHERE;
	for my $p (@usingPoses) {
		my ( @theReads, $absPos );
		if ( $p < 0 ) {
			$absPos = $p - $minLeft;
			for ( @{ $Bases{$p} } ) {
				my $x = reverse $_->[0];
				my $y = reverse $_->[1];
				push @theReads, [ $x, $y, $_->[2] ];
			}
		}
		else {
			$absPos = $p + $minLeft;
			if ( $_[1] eq 'BSseeker2' ) {
				next unless exists $Bases{$p};
			}
			@theReads = @{ $Bases{$p} };
		}
		my $mergstr = mergeStr( \@theReads );
		$mergstr = reverse $mergstr if $p < 0;
		my $depth = @theReads;
		$Results{$absPos} = [ $depth, $mergstr ];
	}
	print "@usingPoses\t", '-' x 25, "\n" if $DEBGUHERE;
	return ( \%Results );
}

sub mergeStr($) {
	my @Strs = @{ $_[0] };
	my ( $maxLen, $merged ) = (0);
	for (@Strs) {
		my $l = length $_->[0];
		$maxLen = $l if $maxLen < $l;
	}
	for my $p ( 0 .. ( $maxLen - 1 ) ) {
		my ( %Col, %ColBp ) = ();
		for (@Strs) {
			my $str  = $_->[0];
			my $qual = $_->[1];
			if ( $p < length($str) ) {
				my $c = substr $str,  $p, 1;
				my $q = substr $qual, $p, 1;
				if ( $_->[2] eq 'CT' and $c eq 'T' ) {
					$c = 'Y';
				}
				elsif ( $_->[2] eq 'GA' and $c eq 'A' ) {
					$c = 'R';
				}
				++$ColBp{$c};
				if ( $c eq 'Y' ) {
					$Col{$_} += $main::Qual2LgP{$q}->[4] for qw(G A);
					$Col{$_} += $main::Qual2LgP{$q}->[3] for qw(C T);
				}
				elsif ( $c eq 'R' ) {
					$Col{$_} += $main::Qual2LgP{$q}->[4] for qw(C T);
					$Col{$_} += $main::Qual2LgP{$q}->[3] for qw(G A);
				}
				else {
					$Col{$_} += $main::Qual2LgP{$q}->[2] for qw(A T C G);
					$Col{$c} += $main::Qual2LgP{$q}->[1];
				}
			}
		}
		my ( $res, @Bps );
		if ( keys(%ColBp) == 1 ) {
			$res = ( keys(%ColBp) )[0];
		}
		else {
			@Bps = sort { $Col{$a} <=> $Col{$b} } keys %Col;
			if ( $Col{ $Bps[1] } - $Col{ $Bps[0] } < 0.02 ) {
				my @t = sort { $a cmp $b } @Bps[ 0, 1 ];
				$res = $REV_IUB{ $t[0] . $t[1] };
			}
			else {
				$res = $Bps[0];
			}
		}
		$merged .= $res;
	}
	return $merged;
}

sub getDeriv($$$) {
	my ( $interest, $bypass, $str ) = @_;
	my $len      = length $str;
	my $Inserted = 0;
	my $Deleted  = 0;
	my @interestedPoses;
	for my $i ( 0 .. $len - 2 ) {
		my $thischar = substr $str, $i, 1;
		my $nextchar = substr $str, $i + 1, 1;
		if ( $thischar eq $nextchar or $thischar =~ /[$bypass]/ ) {
			next;
		}
		elsif ( $thischar eq $interest ) {
			push @interestedPoses, $i;
		}
		elsif ( $nextchar eq $interest ) {
			push @interestedPoses, -( $i + 1 );
		}
		else {
			next;
		}
	}
	return @interestedPoses;
}

sub mergehash($$) {
	my ( $sink, $in ) = @_;
	return unless ref($in);
	return unless ref($in) eq "HASH";
	for my $k ( keys %{$in} ) {
		if ( exists $sink->{$k} ) {
			$sink->{$k} += $in->{$k};
		}
		else {
			$sink->{$k} = $in->{$k};
		}
	}
}

sub sortWsum {
	if    ( $a eq '=Sum=' ) { return -1; }
	elsif ( $b eq '=Sum=' ) { return 1; }
	else                    { return $a cmp $b; }
}

sub do_patch {
	my $Refilename = warnFileExist( $main::RefConfig->{$main::RefFilesSHA}->{'Refilename'} );
	my ( %tID, %tFH );
	for ( @{ $main::Config->{$main::FileData}->{'='} } ) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
	}
	for my $k ( keys %tID ) {
		my $outf  = "$main::RootPath/${main::ProjectID}_analyse/.$k.analyse";
		my $outf2 = "$main::RootPath/${main::ProjectID}_analyse/$k.analyse";
		open IN,  '<', $outf  or die;
		open OUT, '>', $outf2 or die;
		print OUT join( "\t", qw[ID Chr BreakPoint1 BreakPoint2 Virus Strand Start End] ), "\n";
		my %Results;
		while (<IN>) {
			chomp;
			my @dat = split /\t/;
			$Results{ $dat[1] }{ $dat[2] } = \@dat;
		}
		for my $chr ( sort keys %Results ) {
			my @Poses = sort { $a <=> $b } keys %{ $Results{$chr} };
			if ( @Poses == 1 ) {
				my @tmp   = @{ $Results{$chr}{ $Poses[0] } };
				my $last1 = pop @tmp;
				my $last2 = pop @tmp;
				my $last3 = pop @tmp;
				my @Virus;
				my @tVr = split /,/, $last1;
				for ( @tVr, $last2, $last3 ) {
					push @Virus, $_ if $_ != -1;
				}
				print OUT join( "\t", @tmp, $Virus[0], $Virus[-1] ), "\n";
				next;
			}

			my @TTT;
			for my $i ( 0 .. $#Poses ) {
				if ( ( $i == 0 ) or ( $Poses[$i] - $Poses[ $i - 1 ] <= 20 ) ) {
					push @TTT, $Results{$chr}{ $Poses[$i] };
				}
				else {
					if (@TTT) {
						my ( @Virus, @Hum );
						for my $tt (@TTT) {
							push @Hum,   $tt->[2];
							push @Hum,   $tt->[3] if $tt->[3] != -1;
							push @Virus, $tt->[6] if $tt->[6] != -1;
							push @Virus, $tt->[7] if $tt->[7] != -1;
							my @tVr = split /,/, $tt->[8];
							for (@tVr) {
								push @Virus, $_ if $_ != -1;
							}
						}
						@Hum   = sort { $a <=> $b } @Hum;
						@Virus = sort { $a <=> $b } @Virus;
						push @Hum, -1 if scalar @Hum == 1;
						if ( scalar @Virus == 1 ) {
							push @Virus, -1;
						}
						else {
							my $vlen = $Virus[1] - $Virus[0];
							if ( $vlen < $main::MinVirusLength ) {
								next;
							}
						}
						print OUT join( "\t", $Results{$chr}{ $Hum[0] }->[0], $chr, $Hum[0], $Hum[-1], $Results{$chr}{ $Hum[0] }->[4], $Results{$chr}{ $Hum[0] }->[5], $Virus[0], $Virus[-1] ), "\n";
					}
					@TTT = ( $Results{$chr}{ $Poses[$i] } );
				}
			}
			if (@TTT) {
				my @tmp   = @{ $TTT[0] };
				my $last1 = pop @tmp;
				my $last2 = pop @tmp;
				my $last3 = pop @tmp;
				my @Virus;
				my @tVr = split /,/, $last1;
				for ( @tVr, $last2, $last3 ) {
					push @Virus, $_ if $_ != -1;
				}
				print OUT join( "\t", @tmp, $Virus[0], $Virus[-1] ), "\n";
			}
		}
		close IN;
		close OUT;
	}
}

1;
