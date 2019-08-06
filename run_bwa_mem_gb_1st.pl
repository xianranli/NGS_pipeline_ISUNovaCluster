#!/usr/local/bin/perl -w


### for ch6
#my $awk_std = "'".'$3~/Chr07/'."'";
### for CHS & GRFs
##my $target_chr = "'".'@SQ\|Chr04\|Chr06\|Chr07'."'";
#my $awk_std = "'".'$3~/Chr04/&&$4>613972200&&$4<61441827||$3~/Chr06/&&$4>55387501&&$4<55421315||$3~/Chr07/&&$4>60501046&&$4<60544815'."'";
## make sure load "module load sra-toolkit/2.9.6-ub7kz5h"

### three step for aligning reads download to GenBank
#1). download then split the reads into small chunk
#2). After the first step finished, run each chunk individually
#3). After all chunk finished, merge bams into a single bam file

## SC1103.sra 5.8G => 4.4G * 2 => 25
## SC35.sra 7.9G => 6G * 2
## SC971.sra 8.5G => 6.4G * 2
## SC283.sra 5.7G => 4.4G * 2
use strict;
my $bwa = 'bwa';
#my $btrim  = '/home/lixr/NGS_bam/btrim64';
#my $PE_pl = '/home/lixr/NGS_bam/paired_end_trim.pl';
my $pl_2nd = '/home/lixr/NGS_bam/run_bwa_mem_gb_2nd_target.pl';
#my $genome_ref_4_bwa = '/work/jmyu/lixr/SbNAM_bams/Refs/Sb_V3_all_bwa';
my $genome_ref_4_bwa = '/work/jmyu/gb_refs/Sb_V3_all_bwa';


my $pre_dir = '/work/jmyu/lixr/Tan12Sh1/';
#my $pigz = '/work/jmyu/softwares/bin/pigz/';

my $qsub_dir = '/home/lixr/NGS_bam/temp/';
mkdir $qsub_dir unless -e $qsub_dir;
my $split_lines = 8000000; #8000000; ###8_000_000; ## total of BW350-1_1 is 353_607_940. 45 files
my $pigz_cpus = 8;
my $split_line_pre = 'split -l '.$split_lines.' -d -a 3 --filter='."'".'/work/jmyu/softwares/bin/pigz -p '.$pigz_cpus. ' > $FILE.gz'."' ";
#my $awk_std = "'".'$3~/Chr07/&&$4>56430000&&$4<56475000'."'";
#my $awk_std = "'".'$3~/Chr04/&&$4>613972200&&$4<61441827||$3~/Chr06/&&$4>55387501&&$4<55421315||$3~/Chr07/&&$4>60501046&&$4<60544815'."'";

my $srr_list_file = $pre_dir.'SbGRNC_wSRRs'; ##SbNAM_founders_wSRRs';
#my $bwa_group = 3;
for (my $bwa_group = 1; $bwa_group <= 10; $bwa_group ++) { ## 9, 2, 7, 3-4, 1
#	next unless $bwa_group == 1;
	my ($sra_arrayref, $sra_hashref) = SRA_list_DIR($bwa_group, $srr_list_file);

	my @srs = @$sra_arrayref; #keys %$sra_hashref;
#	print $_." A\n" foreach (@srs);
#	my $sam_header_file = '/home/lixr/Sb/Refs/Sb_SAM_header';
#	my $local_SRA_chunk_dir = '/home/lixr/Sb/SRAs/';
	for (my $i = 0; $i <= $#srs; $i ++) { ## 0 - 3; 4 - 6;
		my $sample = $srs[$i];
#		next unless $sample eq 'PE08503';
#		next unless $sample eq 'PE08171' || $sample eq 'PE08482' || $sample eq 'PE08503';

		my $sample_dir = $pre_dir.$sample.'/';
		mkdir $sample_dir unless -e $sample;
#		my $merged_bam = $pre_dir.$sample.'.bam';
#		next if -e $merged_bam;
#		print $sample."\n";
		my @remote_files = @{ $$sra_hashref{$sample}};

#	my @remote_files = qw/SRR2759161/;
		my $k;
		foreach my $run (@remote_files) {

			$k ++;
#			next unless $run eq 'SRR4028763';
			my $squene_000_file    = $qsub_dir.$sample.'_'.$k.'_000';
			my $squene_1st_file_A = $qsub_dir.$sample.'_'.$k.'_1st_A';
			my $squene_1st_file_B = $qsub_dir.$sample.'_'.$k.'_1st_B';

			my $squene_2nd_file = $qsub_dir.$sample.'_'.$k.'_2nd';
			my $squene_3rd_file = $qsub_dir.$sample.'_'.$k.'_3rd';
#			next if $k > 1;
			my ($pre_srs, $x) = $run =~ /(SRR\d\d\d)\d+(\d)$/;
			my $strain = $sample.'_'.$k;
			my $run_dir = $sample_dir.$run.'/';
			mkdir $run_dir unless -e $run_dir;
#			my $ori_run_sra_file = $run_dir.$run.'.sra';
			my $ori_run_sra_file = '/work/jmyu/ncbi/sra/'.$run.'.sra';
			my $ln_run_sra_file = $sample_dir.$strain.'.sra';
#
			my $left_file = $sample_dir.$strain.'_1.fastq';
			my $right_file = $sample_dir.$strain.'_2.fastq';
#			my $left_gz_file = $left_file.'.gz';
#			my $right_gz_file = $right_file.'.gz';
			my $left_chunk_prex = $sample_dir.$strain.'_1_';
			my $right_chunk_prex = $sample_dir.$strain.'_2_';
			my $left_000_gz_file = $sample_dir.$strain.'_1_000.gz';
			my $run_quality_file = $sample_dir.$strain.'_fastq_sample';
#			my $run_quality_file = $local_SRA_chunk_dir.$sample.'/'.$strain.'_fastq_sample';

			my $merged_gb_bam = $pre_dir.$strain.'_Tan12Sh1.bam';
			next if -e $merged_gb_bam;
			print $merged_gb_bam."\n";
			my $target_bams = $sample_dir.$strain.'*Target*bam';
			my $remain_gzs = $sample_dir.$strain.'*gz';
#			my $remote_ftp_file = $remote_ftp_pre.$run.'/'.$run.'.sra';
			my $chunk_num = 35;

			open (SQ_0, '>'.$squene_000_file) || die;
			print SQ_0 '#!/bin/bash'."\n".'#SBATCH --time=6:00:00'."\n".'#SBATCH --nodes=1'."\n".'#SBATCH --ntasks-per-node=6'."\n".'#SBATCH --mem=4G'."\n";

			if ($k > 0) {
				print SQ_0 "module load sra-toolkit/2.9.6-ub7kz5h\n";
				my $size_info = `vdb-dump --info $run | grep SEQ`;
				my @t1 = split /\s+/, $size_info;
				my $size = $t1[-1];
				   $size =~ s/\,//g;
				  $chunk_num = int(4 * $size / $split_lines);
				  $chunk_num = $chunk_num > 120 ? 120 : $chunk_num;
##				print $size." $chunk_num \n";
				print SQ_0 "prefetch $run -X 100G \n";
				print SQ_0 "ln -s  $ori_run_sra_file $ln_run_sra_file\n"  unless -e $ln_run_sra_file;
				print SQ_0 "fasterq-dump  $ln_run_sra_file -o $strain -O $sample_dir -t /work/jmyu/temp/ \n" unless -e $left_file && -e $right_file;
				print SQ_0 "rm $ln_run_sra_file\n";
				print SQ_0 "rm $ori_run_sra_file -f\n";
				print SQ_0 "echo nothing\n";
			}
			close SQ_0;
			chdir $qsub_dir;
			my $jobnum_info_0 = `sbatch $squene_000_file`;
			chop $jobnum_info_0; my @t0 = split /\s+/, $jobnum_info_0;
			my $jobnum_0 = $t0[-1];

			open (SQ_1st_A, '>'.$squene_1st_file_A) || die;
			print SQ_1st_A '#!/bin/bash'."\n".'#SBATCH --time=6:00:00'."\n".'#SBATCH --nodes=1'."\n".'#SBATCH --ntasks-per-node='."$pigz_cpus\n".'#SBATCH --mem=4G'."\n";
			print SQ_1st_A $split_line_pre."$left_file $left_chunk_prex\n";
			print SQ_1st_A "zless $left_000_gz_file | head -n 4000 > $run_quality_file\n" ;
			print SQ_1st_A "rm $left_file\n";
			print SQ_1st_A "echo nothing\n";
			close SQ_1st_A;
			my $jobnum_info_1stA = `sbatch --dependency=afterok:$jobnum_0 $squene_1st_file_A`;
			chop $jobnum_info_1stA; my @t1A = split /\s+/, $jobnum_info_1stA;
			my $jobnum_1stA = $t1A[-1];

			open (SQ_1st_B, '>'.$squene_1st_file_B) || die;
			print SQ_1st_B '#!/bin/bash'."\n".'#SBATCH --time=6:00:00'."\n".'#SBATCH --nodes=1'."\n".'#SBATCH --ntasks-per-node=1'."$pigz_cpus\n".'#SBATCH --mem=4G'."\n";
			print SQ_1st_B $split_line_pre."$right_file $right_chunk_prex\n";
			print SQ_1st_B "rm $right_file\n";
			print SQ_1st_B "echo nothing\n";
			close SQ_1st_B;
			my $jobnum_info_1stB = `sbatch --dependency=afterok:$jobnum_0 $squene_1st_file_B`;
			chop $jobnum_info_1stB; my @t1B = split /\s+/, $jobnum_info_1stB;
			my $jobnum_1stB = $t1B[-1];

			open (SQ_2nd, '>'.$squene_2nd_file) || die;
			print SQ_2nd '#!/bin/bash'."\n".'#SBATCH --time=0:40:00'."\n".'#SBATCH --nodes=1'."\n".'#SBATCH --mem=4G'."\n".'#SBATCH --ntasks-per-node=2'."\n";
			print SQ_2nd '#SBATCH --array=0-'.$chunk_num."\n"; ## 00%20
			print SQ_2nd "module load samtools/1.9-k6deoga\n";
			print SQ_2nd "module load bwa/0.7.17-zhcbtza\n";
			print SQ_2nd "perl $pl_2nd $pre_dir $sample $genome_ref_4_bwa $k ".'$SLURM_ARRAY_TASK_ID'."\n";
			close SQ_2nd;
#			my $jobnum_info_2nd = `sbatch  $squene_2nd_file`;
			my $jobnum_info_2nd = `sbatch --dependency=afterok:$jobnum_1stA:$jobnum_1stB $squene_2nd_file`;
			chop $jobnum_info_2nd; my @t2 = split /\s+/, $jobnum_info_2nd;
			my $jobnum_2nd = $t2[-1];

			open (SQ_3rd, '>'.$squene_3rd_file) || die;
			print SQ_3rd '#!/bin/bash'."\n".'#SBATCH --time=0:40:00'."\n".'#SBATCH --nodes=1'."\n".'#SBATCH --ntasks-per-node=1'."\n".'#SBATCH --mem=4G'."\n".'#SBATCH --ntasks-per-node=2'."\n";
			print SQ_3rd "module load samtools/1.9-k6deoga\n";
			print SQ_3rd 'samtools merge '.$merged_gb_bam.' '.$target_bams."\n";
			print SQ_3rd 'samtools index '.$merged_gb_bam."\n";
			print SQ_3rd 'rm '.$target_bams.' '.$remain_gzs."\n";
			close SQ_3rd;
			my $last_chunk = $jobnum_2nd.'_'.$chunk_num;
#			print "sbatch --dependency=afterok:$jobnum_2nd $squene_3rd_file\n";
			system("sbatch --dependency=afterok:$jobnum_2nd $squene_3rd_file");
			sleep(60)
			}
		}
	}

###########################################
sub SRA_list_DIR {
	my ($bwa_group, $f) = @_;
	my (%hash, @array);
	open (F, $f) || die;
	while (<F>) {
		chomp;
		my @t = split /\t/;
		next if $t[0] eq 'bwa_Group';
#		  next unless $#t == 1;
		next unless $t[0] == $bwa_group;
		push @array, $t[1];
		  @{ $hash{$t[1]} } = @t[3..$#t];
		}
	close F;
	return (\@array, \%hash);
	}

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
