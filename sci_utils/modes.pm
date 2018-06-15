package modes;

sub new {

	$class = shift;
	$alias = shift;
	$file = shift;
	
	$self = {};
	bless($self,$class);
	
	# alias to search for and the input file
	$self->{alias} = $alias;
	$self->{file} = $file;
	
	# set all read inclusions as false - switching to true is global (ie all modalities)
	$self->{read1} = "false";
	$self->{read2} = "false";
	$self->{index1} = "false";
	$self->{index2} = "false";
	
	# set multimodal to false - true will be global
	$self->{multimodal} = "false";
	
	open FILE, "$file";
	while ($l = <FILE>) {
		if ($l !~ /^#/) {
			chomp $l;
			if ($l =~ /^>/) {
			
				if ($include>0) {
					last;
				}
			
				$l =~ s/^>//;
				@aliases = split(/,/, $l);
				$name = $aliases[0];
				$include = 0;
				
				foreach $listed_alias (@aliases) {
					if (uc($listed_alias) eq uc($alias)) {
						$include = 1;
						$self->{name} = $name;
						@indexes_present = (); @output_reads = (); @modalities = ();
						$modality_ID = 0;
						$modalities[$modality_ID] = "modality".$modality_ID;
						$self->{$modalities[$modality_ID]} = {};
						$self->{$modalities[$modality_ID]}->{ins} = {};
						$self->{$modalities[$modality_ID]}->{outs} = {};
					}
				}
				
			} elsif ($include>0) {
				
				if ($l =~ /\&/) {
					# specify multimodal
					$self->{multimodal} = "true";
					
					# close previous modality
					$self->{$modalities[$modality_ID]}->{indexes} = [@indexes_present];
					$self->{$modalities[$modality_ID]}->{outputs} = [@output_reads];
					
					# setup next modality
					$modality_ID++;
					@indexes_present = (); @output_reads = ();
					$modalities[$modality_ID] = "modality".$modality_ID;
					$self->{$modalities[$modality_ID]} = {};
					$self->{$modalities[$modality_ID]}->{ins} = {};
					$self->{$modalities[$modality_ID]}->{outs} = {};
					
				} else {
					
					($read_name,$read_parts) = split(/=/, $l);
					
					if ($read_name =~ /name/) {
						$modalities[$modality_ID] = $read_parts;
					} else {
					
						if ($read_name !~ /(read1|read2|index1|index2)/) {
							die "ERROR - new mode - Read name: $read_name in $file is not a viable read name. Options are read1, read2, index1, or index2.\n";
						}
						
						$self->{$read_name} = "true";
						
						@parts_list = split(/,/, $read_parts);
						$offset = 0;
						for ($i = 0; $i < @parts_list; $i++) {
							if ($parts_list[$i] =~ /read/) { # output read
								($out_name,$length) = split(/:/, $parts_list[$i]);
								$out_name =~ s/^read_//;
								push @output_reads, $out_name;
								$self->{$modalities[$modality_ID]}->{outs}->{$out_name} = {
									"read"   =>   $read_name,
									"offset" =>   $offset,
									"length" =>   $length
								};
								if ($length !~ /(end|all)/) {
									$offset += $length;
								}
							} else { # index read or umi
								($index_name,$length) = split(/:/, $parts_list[$i]);
								push @indexes_present, $index_name;
								$self->{$modalities[$modality_ID]}->{ins}->{$index_name} = {
									"read"   =>   $read_name,
									"offset" =>   $offset,
									"length" =>   $length
								};
								$offset += $length;
							}
						}
					}
				}
			}
		}
	} close FILE;
	
	# close last
	$self->{$modalities[$modality_ID]}->{indexes} = [@indexes_present];
	$self->{$modalities[$modality_ID]}->{outputs} = [@output_reads];
	$self->{modalities} = [@modalities];
	
	return $self;

}

sub modalities {
	$self = shift;
	return @{$self->{modalities}};
}

sub check_multimodal {
	$self = shift;
	return $self->{multimodal};
}

sub indexes {
	$self = shift;
	$modality = shift;
	return @{$self->{$modality}->{indexes}};
}

sub index_length {
	$self = shift;
	$modality = shift;
	$index_name = shift;
	return $self->{$modality}->{ins}->{$index_name}->{length};
}

sub outputs {
	$self = shift;
	$modality = shift;
	return @{$self->{$modality}->{outputs}};
}

sub pull_index {

	$self = shift;
	$modality = shift;
	$index = shift;
	$read_set = shift; # read_set is an object specific to fastq_dump.pm
	
	$index_read = $self->{$modality}->{ins}->{$index}->{read};
	$index_offset = $self->{$modality}->{ins}->{$index}->{offset};
	$index_length = $self->{$modality}->{ins}->{$index}->{length};
	
	$read_seq = $read_set->{$index_read}->{seq};
	
	$index_seq = substr($read_seq,$index_offset,$index_length);
	
	return $index_seq;
	
}

sub pull_read {
	
	$self = shift;
	$modality = shift;
	$read = shift;
	$read_set = shift;
	
	$out_read = $self->{$modality}->{outs}->{$read}->{read};
	$out_offset = $self->{$modality}->{outs}->{$read}->{offset};
	$out_length = $self->{$modality}->{outs}->{$read}->{length};
	
	$read_seq = $read_set->{$out_read}->{seq};
	
	if ($out_length =~ /all/) {
		$out_seq = $read_seq;
	} elsif ($out_length =~ /end/) {
		$out_seq = substr($read_seq,$out_offset);
	} else {
		$out_seq = substr($read_seq,$out_offset,$out_length);
	}
	
	return $out_seq;
	
}

sub read_check {
	$self = shift;
	$read = shift;
	return $self->{$read};
}

sub alias {
	$self = shift;
	return $self->{alias};
}

sub name {
	$self = shift;
	return $self->{name};
}

sub file {
	$self = shift;
	return $self->{file};
}

1;