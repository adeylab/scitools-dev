package sci_commands::atac_count;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_count");

sub atac_count {

#defaults
$format = "S";

@ARGV = @_;
getopts("s:b:O:BC:XF:uI", \%opt);

$die2 = "
scitools atac-count [options] [bam file] [peaks bed file]
   or    count(s)

Creates a counts or binary matrix file and a values file
with the fraction of reads on target for each cell. In addition if matrix defined instead of bed file
can add cells to counts matrix if cells are new (at moment merge bam if cell names do not overlap).
Using a matrix as input is not currently compatible with the sparse matrix option.

Options:
   -O   [STR]   Output prefix (default is peaks prefix)
                (adds .counts.matrix)
   -b   [STR]   Bedtools call (def = $bedtools)
   -s   [STR]   Samtools call (def = $samtools)
   -B           Make it a binary matrix instead of counts matrix
                (file name will end in .binary.matrix)
   -C   [STR]   Complexity file (speeds up on-target calc)
   -F   [S/D]   Format: D = dense, S = sparse (def = $format)
   -u           Do not gzip output (defauly = yes)
   -X           Remove temp files
   -I           Index mode (experimental, sparse matrix only)

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.bed//};
if (defined $opt{'F'}) {$format = $opt{'F'}};
if ($format !~ /(S|D)/i) {die "ERRROR: Formats must be S or D\n$die2"};

#if ends in .matrix then uses peaks from matrix to do counts or binary
if ($ARGV[1] =~ /\.matrix$/) { # THIS BLOCK IS DEPROCATED

	if ($format =~ /S/i) {
		die "ERROR: Sparse matrix functionality not supported if matrix is input.\n";
	}

	print "USING MATRIX FILE, ADDING CELLS: \n";
	#read in matrix
	read_matrix($ARGV[1]);
	$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.binary.matrix//; $opt{'O'} =~ s/\.counts.matrix//; $opt{'O'} =~ s/\.matrix//;
	# create temporary bed file from matrix pre defined peaks
	$tembedfilename=$ARGV[0];
	$tembedfilename=~ s/\.bam//;
	$tembedfilename=~ s/.*\///;
	open OUT, ">$tembedfilename.temp.bed";
	foreach $peak (@MATRIX_ROWNAMES) {
		print OUT join("\t",split("_",$peak))."\n";
	}
	close(OUT);

	#create options string
	$common_opts = "";
	if (defined $opt{'b'}) {$common_opts .= "-b $opt{'b'} "};
	if (defined $opt{'s'}) {$common_opts .= "-s $opt{'s'} "};
	if (defined $opt{'B'}) {$common_opts .= "-B $opt{'B'} "};
	if (defined $opt{'C'}) {$common_opts .= "-C $opt{'C'} "};
	if (defined $opt{'u'}) {$common_opts .= "-u "};
	$common_opts =~ s/\s$//;
	#call same script but on temp bed
	system("scitools atac-count $common_opts $ARGV[0] $tembedfilename.temp.bed");
	#join matrixes and remove temp matrix and bed file if X is not defined
	#DEV: merging matrix files if overlap in cell names see notes in merge-matrix
	if (defined $opt{'B'}) {
		system("scitools matrix-merge -O $opt{'O'}_$tembedfilename $ARGV[1] $tembedfilename.temp.binary.matrix");
		if (!defined $opt{'X'}) {system("rm -f $tembedfilename.temp.bed $tembedfilename.binary.temp.matrix")};
	} else {
		system("scitools matrix-merge -O $opt{'O'}_$tembedfilename $ARGV[1] $tembedfilename.temp.counts.matrix");
		if (!defined $opt{'X'}) {system("rm -f $tembedfilename.temp.bed $tembedfilename.counts.temp.matrix")};
	}
	
} else {

	if (defined $opt{'C'}) {
		read_complexity($opt{'C'});
	} else {
		open IN, "$samtools view $ARGV[0] |";
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			if ($P[2] !~ /(M|Y|L|K|G|Un|Random|Alt)/i) {
				($cellID,$null) = split(/:/, $P[0]);
				$CELLID_uniq_reads{$cellID}++;
			}
		} close IN;
	}
	
	if (defined $opt{'I'}) { # EXPERIMENTAL MODE! SPARSE MATRIX ONLY
	
		if (!defined $opt{'u'}) {
			if (defined $opt{'B'}) {
				open CELLS, "| $gzip > $opt{'O'}.binary.sparseMatrix.cols.gz";
				open SITES, "| $gzip > $opt{'O'}.binary.sparseMatrix.rows.gz";
				open VALS, "| $gzip > $opt{'O'}.binary.sparseMatrix.values.gz";
			} else {
				open CELLS, "| $gzip > $opt{'O'}.counts.sparseMatrix.cols.gz";
				open SITES, "| $gzip > $opt{'O'}.counts.sparseMatrix.rows.gz";
				open VALS, "| $gzip > $opt{'O'}.counts.sparseMatrix.values.gz";
			}
		} else {
			if (defined $opt{'B'}) {
				open CELLS, ">$opt{'O'}.binary.sparseMatrix.cols";
				open SITES, ">$opt{'O'}.binary.sparseMatrix.rows";
				open VALS, ">$opt{'O'}.binary.sparseMatrix.values";
			} else {
				open CELLS, ">$opt{'O'}.counts.sparseMatrix.cols";
				open SITES, ">$opt{'O'}.counts.sparseMatrix.rows";
				open VALS, ">$opt{'O'}.counts.sparseMatrix.values";
			}
		}
		
		$bai = $ARGV[0]; $bai =~ s/\.bam$/\.bai/;
		if (-e "$ARGV[0].bai" || -e $bai) {} else {
			system("$samtools index $ARGV[0]");
		}
		
		%CELLID_onTarget = ();
		
		$cellNum = 0; $siteID = 0;
		open IN, "$ARGV[1]";
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			$siteID++; $siteName = "$P[0]_$P[1]_$P[2]";
			print SITES "$siteName\n";
			
			%CELL_VALUES = ();
			open INT, "$samtools view $ARGV[0] $P[0]:$P[1]-$P[2] |";
			while ($int_l = <INT>) {
				chomp $int_l;
				@INT_P = split(/\t/, $int_l);
				$cellID = $INT_P[0]; $cellID =~ s/:.+$//;
				if (!defined $CELLID_number{$cellID}) {
					$cellNum++;
					$CELLID_number{$cellID} = $cellNum;
					print CELLS "$cellID\n";
				}
				$CELL_VALUES{$cellID}++;
				$CELLID_onTarget{$cellID}++;
			} close INT;
			
			foreach $cellID (keys %CELL_VALUES) {
				print VALS "$siteID\t$CELLID_number{$cellID}\t$CELL_VALUES{$cellID}\n";
			}
		} close IN;
		close CELLS; close VALS; close SITES;
		
	} else {

		if ($format =~ /D/i) {
		
			open IN, "$bedtools intersect -abam $ARGV[0] -b $ARGV[1] -bed -wa -wb |";
			while ($l = <IN>) {
				chomp $l;
				@P = split(/\t/, $l);
				$cellID = $P[3]; $cellID =~ s/:.+$//;
				$siteID = $P[12]."_".$P[13]."_".$P[14];
				if (!defined $SITE_CELL_STATUS{$siteID}{$cellID}) {
					$SITE_CELL_STATUS{$siteID}{$cellID} = 1;
				} else {
					if (defined $opt{'B'}) {
						$SITE_CELL_STATUS{$siteID}{$cellID} = 1;
					} else {
						$SITE_CELL_STATUS{$siteID}{$cellID}++;
					}
				}
				$CELLID_onTarget{$cellID}++;
				if (!defined $CELL_IDS{$cellID}) {
					push @CELL_ID_LIST, $cellID;
					$CELL_IDS{$cellID} = 1;
				}
			} close IN;

			if (!defined $opt{'u'}) {
				if (defined $opt{'B'}) {
					open OUT, "| $gzip > $opt{'O'}.binary.matrix.gz";
				} else {
					open OUT, "| $gzip > $opt{'O'}.counts.matrix.gz";
				}
			} else {
				if (defined $opt{'B'}) {
					open OUT, ">$opt{'O'}.binary.matrix";
				} else {
					open OUT, ">$opt{'O'}.counts.matrix";
				}
			}
			
			print OUT "$CELL_ID_LIST[0]";
			for ($i = 1; $i < @CELL_ID_LIST; $i++) {
				print OUT "\t$CELL_ID_LIST[$i]";
			} print OUT "\n";
			
			if (!defined $opt{'n'}) {
				foreach $siteID (sort keys %SITE_CELL_STATUS) {
					print OUT "$siteID";
					for ($i = 0; $i < @CELL_ID_LIST; $i++) {
						if ($SITE_CELL_STATUS{$siteID}{$CELL_ID_LIST[$i]}>0.5) {
							if (defined $opt{'B'}) {
								print OUT "\t1";
							} else {
								print OUT "\t$SITE_CELL_STATUS{$siteID}{$CELL_ID_LIST[$i]}";
							}
						} else {
							print OUT "\t0";
						}
					}
					print OUT "\n";
				}
			} else {
				foreach $siteID (keys %SITE_CELL_STATUS) {
					print OUT "$siteID";
					for ($i = 0; $i < @CELL_ID_LIST; $i++) {
						if ($SITE_CELL_STATUS{$siteID}{$CELL_ID_LIST[$i]}>0.5) {
							if (defined $opt{'B'}) {
								print OUT "\t1";
							} else {
								print OUT "\t$SITE_CELL_STATUS{$siteID}{$CELL_ID_LIST[$i]}";
							}
						} else {
							print OUT "\t0";
						}
					}
					print OUT "\n";
				}
			}
			
			close OUT;
			
		} elsif ($format =~ /S/i) {
		
			if (!defined $opt{'u'}) {
				if (defined $opt{'B'}) {
					open CELLS, "| $gzip > $opt{'O'}.binary.sparseMatrix.cols.gz";
					open SITES, "| $gzip > $opt{'O'}.binary.sparseMatrix.rows.gz";
				} else {
					open CELLS, "| $gzip > $opt{'O'}.counts.sparseMatrix.cols.gz";
					open SITES, "| $gzip > $opt{'O'}.counts.sparseMatrix.rows.gz";
				}
			} else {
				if (defined $opt{'B'}) {
					open CELLS, ">$opt{'O'}.binary.sparseMatrix.cols";
					open SITES, ">$opt{'O'}.binary.sparseMatrix.rows";
				} else {
					open CELLS, ">$opt{'O'}.counts.sparseMatrix.cols";
					open SITES, ">$opt{'O'}.counts.sparseMatrix.rows";
				}
			}
			
			$cellNum = 1; $siteNum = 1;
			open IN, "$bedtools intersect -abam $ARGV[0] -b $ARGV[1] -bed -wa -wb |";
			while ($l = <IN>) {
				chomp $l;
				@P = split(/\t/, $l);
				$cellID = $P[3]; $cellID =~ s/:.+$//;
				$siteID = $P[12]."_".$P[13]."_".$P[14];
				
				if (!defined $SITEID_number{$siteID}) {
					$SITEID_number{$siteID} = $siteNum;
					print SITES "$siteID\n";
					$siteNum++;
				}
				
				if (!defined $CELLID_number{$cellID}) {
					$CELLID_number{$cellID} = $siteNum;
					print CELLS "$cellID\n";
					$cellNum++;
				}
				
				$SITENUM_CELLNUM_count{$SITEID_number{$siteID}}{$CELLID_number{$cellID}}++;
				$CELLID_onTarget{$cellID}++;
				
			} close IN; close CELLS; close SITES;

			if (!defined $opt{'u'}) {
				if (defined $opt{'B'}) {
					open VALS, "| $gzip > $opt{'O'}.binary.sparseMatrix.values.gz";
				} else {
					open VALS, "| $gzip > $opt{'O'}.counts.sparseMatrix.values.gz";
				}
			} else {
				if (defined $opt{'B'}) {
					open VALS, ">$opt{'O'}.binary.sparseMatrix.values";
				} else {
					open VALS, ">$opt{'O'}.counts.sparseMatrix.values";
				}
			}

			foreach $siteNum (keys %SITENUM_CELLNUM_count) {
				foreach $cellNum (keys %{$SITENUM_CELLNUM_count{$siteNum}}) {
					if (defined $opt{'B'}) {
						print VALS "$siteNum\t$cellNum\t1\n";
					} else {
						print VALS "$siteNum\t$cellNum\t$SITENUM_CELLNUM_count{$siteNum}{$cellNum}\n";
					}
				}
			}

			close VALS;
			
		}
	}

	open OUT, ">$opt{'O'}.fracOnTarget.values";
	for ($i = 0; $i < @CELL_ID_LIST; $i++) {
		if ($CELLID_uniq_reads{$cellID}>0) {
			$cellID = $CELL_ID_LIST[$i];
			$frac = sprintf("%.3f", $CELLID_onTarget{$cellID}/$CELLID_uniq_reads{$cellID});
			print OUT "$cellID\t$frac\n";
		} 
	} close OUT;
}

}
1;
