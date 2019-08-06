#!/usr/local/bin/perl -w


### for ch6
#my $awk_std = "'".'$3~/Chr07/'."'";
### for CHS & GRFs
##my $target_chr = "'".'@SQ\|Chr04\|Chr06\|Chr07'."'";
#my $awk_std = "'".'$3~/Chr04/&&$4>613972200&&$4<61441827||$3~/Chr06/&&$4>55387501&&$4<55421315||$3~/Chr07/&&$4>60501046&&$4<60544815'."'";

### three step for aligning reads download to GenBank
#1). download then split the reads into small chunk
#2). After the first step finished, run each chunk individually
#3). After all chunk finished, merge bams into a single bam file


use strict;
my $bwa = 'bwa';
my $btrim  = '/home/lixr/NGS_bam/btrim64';
my $PE_pl = '/home/lixr/NGS_bam/paired_end_trim.pl';
#my $filter_pl = '/home/lixr/PollenSeq/scripts/filter_unique_reads_bwa.pl';
#my $genome_ref_4_bwa = '/work/jmyu/lixr/ZmNAM_bams/Refs/B73_v4';
#my $genome_ref_4_samtools = '/work/agron_1/Sb/Refs/Sb_V2_all.fa';
my ($pre_dir, $sample, $genome_ref_4_bwa, $k, $i) = @ARGV;
#my $pre_dir = '/work/jmyu/lixr/SbNC_bams/';

my $qsub_dir = '/home/lixr/NGS_bam/temp/';
mkdir $qsub_dir unless -e $qsub_dir;


#my $awk_std = "'".'$3~/Chr07/&&$4>56430000&&$4<56475000'."'";
my $awk_std = "'".'$3~/Chr04/&&$4>61315396&&$4<63318779||$3~/Chr02/&&$4>6884092&&$4<8889052||$3~/Chr01/&&$4>11271177&&$4<13277291||$3~/Chr01/&&$4>6505747&&$4<8509193'."'";

my $sam_header_file = '/work/jmyu/gb_refs/SbV3_samheader';

		my $sample_dir = $pre_dir.$sample.'/';
		mkdir $sample_dir unless -e $sample_dir;
		my $merged_bam = $pre_dir.$sample.'_gb.bam';
		unless (-e $merged_bam) {

#			next if $k > 1;
			my $strain = $sample.'_'.$k;
#
			my $left_chunk_prex = $sample_dir.$strain.'_1_';
			my $right_chunk_prex = $sample_dir.$strain.'_2_';
			my $run_quality_file = $sample_dir.$strain.'_fastq_sample';
			my $quality_code = Quality_Code($run_quality_file);

#				next unless $i == 0;
#				print $i." ";
			my $suffix = sprintf "%03d", $i;
			my $left_chunk_file = $left_chunk_prex.$suffix.'.gz';
			my $right_chunk_file = $right_chunk_prex.$suffix.'.gz';

			my $trim_left_chunk = $left_chunk_file.'.fq_trm';
			my $trim_left_chunk_s = $trim_left_chunk.'_sum';
			my $trim_right_chunk = $right_chunk_file.'.fq_trm';
			my $trim_right_chunk_s = $trim_right_chunk.'_sum';

			my $trim_PE_1_chunk = $trim_left_chunk.'.pe';
			my $trim_PE_1_chunk_s = $sample_dir.$strain.'_1_'.$suffix.'.pe_sum';

			my $trim_PE_2_chunk = $trim_right_chunk.'.pe';
			my $trim_PE_2_chunk_s = $sample_dir.$strain.'_2_'.$suffix.'.pe_sum';

			 ## -e $trim_PE_1_chunk && -e $trim_PE_2_chunk

			my $trim_SE_1_chunk = $trim_left_chunk.'.se';
			my $trim_SE_2_chunk = $trim_right_chunk.'.se';
			my $trim_SE_chunk = $sample_dir.$strain.'_0_'.$suffix.'.se';

			my $PE_target_bam = $sample_dir.$strain.'_pe_Target.bam_'.$suffix;
			my $SE_target_bam = $sample_dir.$strain.'_se_Target.bam_'.$suffix;

			my $trim_PE_chunk_sam = $sample_dir.$strain.'_pe_all.sam_'.$suffix;
			my $trim_SE_chunk_sam = $sample_dir.$strain.'_se_all.sam_'.$suffix;

			my $trim_PE_chunk_uni_sam = $sample_dir.$strain.'_uni_pe.sam_'.$suffix;
			my $trim_SE_chunk_uni_sam = $sample_dir.$strain.'_uni_se.sam_'.$suffix;

			my $trim_PE_chunk_uni_bam = $sample_dir.$strain.'_pe_uni.bam_'.$suffix;
			my $trim_SE_chunk_uni_bam = $sample_dir.$strain.'_se_uni.bam_'.$suffix;
			my $chunk_uni_bam = $sample_dir.$strain.'_uni.bam_'.$suffix;

			my $chunk_uni_srt_bam = $sample_dir.$strain.'_Target_'.$suffix.'.bam';
#			unlink $left_chunk_file if -e  $chunk_uni_srt_bam;
#			unlink $right_chunk_file if -e  $chunk_uni_srt_bam;
#			next if -e $chunk_uni_srt_bam;
			unlink $trim_PE_chunk_sam if -e $trim_PE_chunk_sam;
			unlink $trim_SE_chunk_sam if -e $trim_SE_chunk_sam;

			unless (-e $trim_left_chunk_s) {
				if ($quality_code eq 'I') {
					system("$btrim -i -q -t $left_chunk_file -Z -l 40 -o $trim_left_chunk -s $trim_left_chunk_s");
					system("rm $left_chunk_file");
					system("$btrim -i -q -t $right_chunk_file -Z -l 40 -o $trim_right_chunk -s $trim_right_chunk_s");
					system("rm $right_chunk_file");
					}
					elsif ($quality_code eq 'S') {
						system("$btrim -q -t $left_chunk_file -Z -l 40 -o $trim_left_chunk -s $trim_left_chunk_s");
						system("rm $left_chunk_file");
						system("$btrim -q -t $right_chunk_file -Z -l 40 -o $trim_right_chunk -s $trim_right_chunk_s");
						system("rm $right_chunk_file");
					}
			}
			system("$PE_pl $trim_left_chunk_s $trim_right_chunk_s $trim_left_chunk $trim_right_chunk");
			system("rm $trim_left_chunk_s $trim_right_chunk_s $trim_left_chunk $trim_right_chunk");

			system("cat $trim_SE_1_chunk $trim_SE_2_chunk > $trim_SE_chunk");
			system("rm $trim_SE_1_chunk $trim_SE_2_chunk");

#			system("bwa mem -t 2 -v 0 $genome_ref_4_bwa $trim_PE_1_chunk $trim_PE_2_chunk | samtools view - -q 10 -Sb | samtools sort - -o $PE_target_bam");
#			system("rm $trim_PE_1_chunk $trim_PE_2_chunk");
#
#			system("bwa mem -t 2 -v 0 $genome_ref_4_bwa $trim_SE_chunk | samtools view - -q 10 -Sb | samtools sort - -o $SE_target_bam");
#			system("rm $trim_SE_chunk");
			### for targeted candidate genes
			system("bwa mem -t 2 -v 0 $genome_ref_4_bwa $trim_PE_1_chunk $trim_PE_2_chunk | awk $awk_std >> $trim_PE_chunk_sam");
			system("cat $sam_header_file $trim_PE_chunk_sam | samtools view - -q 10 -Sb | samtools sort - -o $PE_target_bam")	;
			system("rm $trim_PE_1_chunk $trim_PE_2_chunk $trim_PE_chunk_sam");

			system("bwa mem -t 2 -v 0 $genome_ref_4_bwa $trim_SE_chunk | awk $awk_std >> $trim_SE_chunk_sam");
			system("cat $sam_header_file $trim_SE_chunk_sam | samtools view - -q 10 -Sb | samtools sort - -o $SE_target_bam")	;
			system("rm $trim_SE_chunk $trim_SE_chunk_sam");

			system("samtools merge $chunk_uni_bam $PE_target_bam $SE_target_bam");

			system("rm $PE_target_bam $SE_target_bam");

			system("samtools sort $chunk_uni_bam -o $chunk_uni_srt_bam");
			system("rm $chunk_uni_bam");
		}
#			my $target_bams = $sample_dir.'*Target*bam';
#			system("samtools merge $merged_bam $target_bams -f");
#			system("samtools index $merged_bam");

###########################################
#sub SRA_list_DIR {
#	my ($bwa_group, $f) = @_;
#	my (%hash, @array);
#	open (F, $f) || die;
#	while (<F>) {
#		chomp;
#		my @t = split /\t/;
#		next if $t[0] eq 'bwa_Group';
##		  next unless $#t == 1;
#		next unless $t[0] == $bwa_group;
#		push @array, $t[1];
#		  @{ $hash{$t[1]} } = @t[3..$#t];
#		}
#	close F;
#	return (\@array, \%hash);
#	}

sub Quality_Code {
	my ($f) = @_;
	my ($j, $flag, %hash, $code) ;
	open (F, $f) || die;
	while (<F>) {
		chomp;
		my $line = $_;
		$flag = 0 if $line =~ /^\+/;
		$flag ++;
		$j ++;
		if ($flag) {
			my @t = split //, $line;
			foreach my $t(@t) {
				my $n = ord($t);
				  $hash{$n} ++;
				}
			}
		last if $j > 1000
		}
	close F;
	if (exists $hash{90}) { $code = 'I'} elsif (exists $hash{48}) {$code = 'S'};
	return $code;
	}


#sub SRR_list_by_SRX {
##	my ($bwa_group) = @_;
##	my (%hash, @array);
#	open (F, '/home/lixr/Sb/NC_SRAs') || die;
#	open (O, '>/home/lixr/Sb/NC_SRAs_wSRRs') || die;
#	my $ftp_site = 'ftp-trace.ncbi.nih.gov';
#	my $ftp = Net::FTP->new($ftp_site);
#     $ftp->login();
#	while (<F>) {
#		chomp;
#		my $line = $_;
#		my @t = split /\t/, $line;
#		next if $t[0] eq 'bwa_Group';
###		  next unless $#t == 1;
##		next unless $t[-1] == $bwa_group;
##		push @array, $t[1];
#	my $srx = $t[2];
#		my ($pre_srs, $x) = $srx =~ /(SRX\d\d\d)(\d+)/;
#		my $remote_ftp_dir = '/sra/sra-instant/reads/ByExp/sra/SRX/'.$pre_srs.'/'.$srx;
#		  $ftp->cwd($remote_ftp_dir);
#		my @remote_files = $ftp->ls();
#		print O $line;
#		print O "\t".$_ foreach (@remote_files);
#		print O "\n";
#		}
#	close F;
##	return (\@array, \%hash);
#	   $ftp->quit;
#
#	}
#
#sub SRR_list_by_SRS {
##	my ($bwa_group) = @_;
##	my (%hash, @array);
#	open (F, '/home/lixr/Sb/SbNAM_founders_1') || die;
#	open (O, '>/home/lixr/Sb/SbNAM_founders_wSRRs') || die;
#	my $ftp_site = 'ftp-trace.ncbi.nih.gov';
#	my $ftp = Net::FTP->new($ftp_site);
#     $ftp->login();
#	while (<F>) {
#		chomp;
#		my $line = $_;
#		my @t = split /\t/, $line;
#		next if $t[0] eq 'bwa_Group';
###		  next unless $#t == 1;
##		next unless $t[-1] == $bwa_group;
##		push @array, $t[1];
#	my $srx = $t[1];
#		my ($pre_srs, $x) = $srx =~ /(SRS\d\d\d)(\d+)/;
#		my $remote_ftp_dir = '/sra/sra-instant/reads/ByStudy/sra/SRS/'.$pre_srs.'/'.$srx;
#		  $ftp->cwd($remote_ftp_dir);
#		my @remote_files = $ftp->ls();
#		print O $line;
#		print O "\t".$_ foreach (@remote_files);
#		print O "\n";
#		}
#	close F;
##	return (\@array, \%hash);
#	   $ftp->quit;
#
#	}	
