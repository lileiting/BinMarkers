#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use PDL::LiteF;
use PDL::Stats::Distr;

#--Preset--#
my $author = "Leiting Li";
my $version = "20150514";

sub message {local $\ = "\n"; print STDERR @_}
sub usage{
    print <<EOF;

perl binmarkers.pl [OPTIONS] <MARKER_MATRIX>

  MARKER_MATRIX  Scafold/contig and position information
                 as marker name, '-' or '_' as seperator,

  OPTIONS
  
  -t,--threshold NUM 
        Maximum difference allowed within a block
        Default: 5

  -0,--letter_for_0_0 CODE
        Specify code for 0/0 SNP
        Default: a
        
  -1,--letter_for_0_1 CODE
        Specify code for 0/1 SNP
        Default: h
        
  -2,--letter_for_1_1 CODE
        Specify code for 1/1 SNP
        Default: b

  -m,--missing CODE
        Specify code for missing
        Default: -
  
  --error_rate_for_0_0 DECIMAL
        Specify the error rate for 0/0 SNP
        Default: 0.04

  --error_rate_for_0_1 DECIMAL
        Specify the error rate for 0/1 SNP
        Default: 0.03

  --error_rate_for_1_1 DECIMAL
        Specify the error rate for 0/1 SNP
        Default: 0.01
     
  -h,--help      
        print this usage message
  
EOF
    exit;
}

my $infile;
my $threshold = 5;
my $help;

my $letter_for_0_0 = 'a'; # a_letter
my $letter_for_0_1 = 'h'; # h_letter
my $letter_for_1_1 = 'b'; # b_letter
my $letter_for_missing = '-';

my $error_rate_for_0_0 = 0.04; # a_error
my $error_rate_for_0_1 = 0.03; # h_error
my $error_rate_for_1_1 = 0.01; # b_error

GetOptions(
  "t|threshold=i" => \$threshold,
  "h|help" => \$help,
  "0|letter_for_0_0=s" => \$letter_for_0_0,
  "1|letter_for_0_1=s" => \$letter_for_0_1,
  "2|letter_for_1_1=s" => \$letter_for_1_1,
  "m|missing=s" => $letter_for_missing,
  "error_rate_for_0_0=i" => \$error_rate_for_0_0,
  "error_rate_for_0_1=i" => \$error_rate_for_0_1,
  "error_rate_for_1_1=i" => \$error_rate_for_1_1
) or die("Error in command line arguments");

usage if $help;

$infile = shift @ARGV;
usage unless $infile;
die "$infile: NOT EXIST!" unless -e $infile;

my $a_letter = $letter_for_0_0;
my $h_letter = $letter_for_0_1;
my $b_letter = $letter_for_1_1;
my $a_error = $error_rate_for_0_0;
my $h_error = $error_rate_for_0_1;
my $b_error = $error_rate_for_1_1;

# Confirm error rate range is in [0,1]
my @error_rates = ($a_error, $h_error, $b_error);
map{die unless $_ <= 1 and $_ >= 0}@error_rates;

my @heterozygous_genotype = ($letter_for_0_1);
my @homologous_genotype = ($letter_for_0_0, $letter_for_1_1);
my @missing_genotype = ($letter_for_missing);
my @valid_genotype = (@homologous_genotype, @heterozygous_genotype);
my @all_genotypes = (@valid_genotype, @missing_genotype);
my %is_homologous = map{$_ => 1}@homologous_genotype;
my %is_heterozygous = map{$_ => 1}@heterozygous_genotype;
my %is_valid = map{$_ => 1}@valid_genotype;
my %is_missing = map{$_ => 1}@missing_genotype;
my %is_genotype = map{$_ => 1}@all_genotypes;

my %matrix;
my @marker_index;
my $num_of_markers = 0;
my %bin_marker;
my $num_of_prob_select = 0;
my $num_of_aaabbb_type = 0;
my $num_of_rand_select = 0;

my $log_file = "$infile.bin_markers_$threshold.log";
my $bin_markers_file = "$infile.bin_markers_$threshold.map";
my $manual_checking_file = "$infile.bin_markers_$threshold.check.txt";
my $prob_file = "$infile.bin_markers_$threshold.prob.txt";
open my $prob_fh, ">", $prob_file or die;

my $status_memory = 0;
sub print_status{
    my ($current, $total) = @_;
    if($current / $total * 10 > $status_memory){
        $status_memory++;
        message "${status_memory}0%";
    }
}

sub are_diff_genotypes{
    my($first, $second) = @_;
    return 0 if $is_missing{$first};
    return 0 if $is_missing{$second};
    return 0 if $first eq $second;
    return 1;
}

my %bin_marker_threshold;
my $total_diff_cal_times;
sub difference{
    my ($block_id, $first, $second) = @_;
    $total_diff_cal_times++;
    my @first = @{$matrix{$first}->{array}};
    my @second = @{$matrix{$second}->{array}};
    for my $index (1..$#first){
         my $a = $first[$index];
         my $b = $second[$index];
         if(are_diff_genotypes($a, $b)){
            $bin_marker_threshold{$block_id}->{$index}++;
         }
    }
    my $difference = keys %{$bin_marker_threshold{$block_id}};
    #message "Block $block_id: Difference for $a and $b is $difference";
    return $difference;
}

sub with_title{
    my $array = shift;
    my $first_genotype = $array->[1];
    return $is_genotype{$first_genotype} ? 0 : 1;
}

sub checking_genotypes{
    my $array = shift;
    my @genotypes = @{$array}[1..$#{$array}];
    for (@genotypes){
        die "CAUTION: $_ is not an valid genotype!" 
            unless $is_genotype{$_};
    }
}

sub count_missing{
    my $array = shift;
    my @genotypes = @{$array}[1..$#{$array}];
    my $missing = 0;
    for(@genotypes){
        $missing++ if $is_missing{$_}
    }
    return $missing;
}

open my $log_fh, "> $log_file" or die;
sub print_log{
    print $log_fh @_,"\n";
}

my $title;
sub load_marker_matrix{
    my $file = shift;
    my $line_count;
    open my $fh, "<", $file or die "$file: $!";
    while(my $marker = <$fh>){
        $line_count++;
        chomp $marker;
        my $ref = [split /\s+/, $marker];
        if($line_count > 1 and with_title $ref){
            die "CAUTION: Undefined genotype codes for marker: $marker\n",
                "Current defined genotypes: @all_genotypes\n";
        }
        if($line_count == 1 and with_title $ref){
            $title = "$marker\n";
            next;
        }
        my ($scaffold, $position) = split /[\-_]/, $ref->[0];
        die unless ($position =~ /^(\d+).*/);
        $position = $1;
        checking_genotypes $ref;
        my $index = ++$num_of_markers;
        push @marker_index, $index;
        $matrix{$index}->{array} = $ref;
        $matrix{$index}->{missing} = count_missing $ref;
        $matrix{$index}->{scaffold} = $scaffold;
        $matrix{$index}->{position} = $position;
    }
}

sub random_select{
    my $n = scalar(@_);
    die unless $n > 0;
    my $random_index = int(rand($n));
    return $_[$random_index];
}

sub convert_h_to_a_or_b{
    my $ref = shift;
    my $first = random_select($letter_for_0_0, $letter_for_1_1);
    my $second = $first eq $letter_for_0_0 ? $letter_for_1_1 : $letter_for_0_0;
    my @genotypes;
    my $i = 0;
    for(@$ref){
        my $choice = do{
            if(   $_ eq $letter_for_0_1 and $i % 2 == 0){
                $first;
            }elsif($_ eq $letter_for_0_1 and $i % 2 == 1){
                $second;
            }else{$_}
        };
        $i++;
        push @genotypes, $choice;
    }
    return @genotypes;
}

sub count_valid_genotypes{
    my $n = 0;
    map{$n++ if $is_valid{$_}}@_;
    return $n;
}

sub count_b{
    my $n = 0;
    map{$n++ if $_ eq $b_letter}@_;
    return $n;
}

sub no_het_genotype{
    map{return 0 if $is_heterozygous{$_}}@_;
    return 1;
}

sub select_genotype{
    my ($genotypes_array_ref, $valid_genotypes_array_ref, $countif_hash_ref) = @_;
    my @genotypes = grep{$is_valid{$_}}@$genotypes_array_ref;
    my $genotypes_str = join('', @genotypes);
    my @valid_genotypes = @{$valid_genotypes_array_ref};
    my %countif = %{$countif_hash_ref};
    my ($first, $second, $third) = 
        sort{$countif{$b} <=> $countif{$a}} @valid_genotypes;
    
    my $equal_case_choose = 0;
    if(@valid_genotypes == 2 and no_het_genotype(@valid_genotypes) and  
       $genotypes_str =~ 
       /^(($a_letter|$b_letter)\2*)(($a_letter|$b_letter)\4*)$/
    ){
        $num_of_aaabbb_type++;
        my($first_part,$first_letter, $second_part, $second_letter) = ($1,$2,$3,$4);
        $equal_case_choose = do{
            if(length($first_part) == length($second_part)){
                $num_of_rand_select++;
                random_select($first_letter, $second_letter);
            }elsif(length($first_part) > length($second_part)){
                $first_letter
            }else{
                $second_letter
            }
        };
        #die "$genotypes_str => $equal_case_choose";
    } 
    
    my @converted_genotypes = convert_h_to_a_or_b($genotypes_array_ref);
    my $num_of_valid_genotypes = count_valid_genotypes(@converted_genotypes);
    my $num_of_b = count_b(@converted_genotypes);
    my $a_ex_prob = pmf_binomial($num_of_b, $num_of_valid_genotypes, $a_error);
    my $h_ex_prob = pmf_binomial($num_of_b, $num_of_valid_genotypes, 0.5 + $h_error/ 2);
    my $b_ex_prob = pmf_binomial($num_of_b, $num_of_valid_genotypes, 1 - $b_error);
#    my %prob = ($a_letter => $a_ex_prob, 
#                $h_letter => $h_ex_prob, 
#                $b_letter => $b_ex_prob);
#    my $best_prob_genotype = (sort{$prob{$b} <=> $prob{$a}}(keys %prob))[0];
    my $best_prob_genotype = $h_ex_prob >= $a_ex_prob ? $h_letter : $a_letter;
    my $tmp_prob = $h_ex_prob >= $a_ex_prob ? $h_ex_prob : $a_ex_prob;
    $best_prob_genotype = $tmp_prob >= $b_ex_prob ? $best_prob_genotype : $b_letter;

    print $prob_fh "$genotypes_str => ", 
           join('',@converted_genotypes),
           " => $equal_case_choose| $best_prob_genotype\n",
           "Valid genotypes: ", scalar(@valid_genotypes),
           " N(b): ", $num_of_b, 
           " N: ", $num_of_valid_genotypes, 
           " P(a): ", $a_ex_prob, 
           " P(h): ", $h_ex_prob,
           " P(b): ", $b_ex_prob, "\n";
    if($equal_case_choose){return $equal_case_choose}
    else{
        $num_of_prob_select++;
        return $best_prob_genotype;
    }
}

sub judge_genotype{
    my @genotypes = @{$_[0]};
    my %countif;
    map{$countif{$_}++}@genotypes;
    my @all_genotypes = keys %countif;
    my @valid_genotypes = grep{$is_valid{$_}}@all_genotypes;
    my @missing_genotypes = grep{$is_missing{$_}}@all_genotypes;
    my $result;
    if(@all_genotypes == 1){
        return $all_genotypes[0];
    }elsif(@valid_genotypes == 1){
        return $valid_genotypes[0];
    }elsif(@valid_genotypes == 2 or @valid_genotypes == 3){
        return select_genotype(\@genotypes, \@valid_genotypes, \%countif);
    }else{die "More than three genotypes? @genotypes"}
}

sub consensus_marker{
    my $block_id = shift;
    my @marker_indexes = @_;
    my @consensus_marker;

    my $max_array_index = $#{$matrix{$marker_indexes[0]}->{array}};
    for(my $i = 1; $i <= $max_array_index; $i++){
        my @genotypes = map{$matrix{$_}->{array}->[$i]}@marker_indexes;
        push @consensus_marker, judge_genotype(\@genotypes);
    }
    return @consensus_marker;
}

my %block_contents;
sub cluster_markers{
    #my @sorted_marker_index = sort
    #    {$matrix{$a}->{scaffold} cmp $matrix{$b}->{scaffold} or
    #     $matrix{$a}->{position} <=> $matrix{$b}->{position}
    #    }keys %matrix;
    my @sorted_marker_index = (1..$num_of_markers);
    
    my $block_id = 0;
    for (my $i = 0; $i <= $#sorted_marker_index; $i++){
        my $marker_index = $sorted_marker_index[$i];
        print_status($i, $num_of_markers);
        my $start_scaffold = $matrix{$marker_index}->{scaffold};
        my $end_scaffold;
        my $start_position = $matrix{$marker_index}->{position};
        my $end_position = $start_position;

        #-- Start dynanmic searching 
        $block_id++;
        my @block;
        push @block, $marker_index;
        LABEL: for(my $j = $i + 1; $j <= $#sorted_marker_index; $j++){
            my $new_marker_index = $sorted_marker_index[$j];
            $end_scaffold = $matrix{$new_marker_index}->{scaffold};
            
            # Stop searching with these conditions
            last LABEL if $end_scaffold ne $start_scaffold;
            for my $index_in_block (@block){
                my $d = difference($block_id, $index_in_block, $new_marker_index);
                last LABEL if $d > $threshold;
            }
            
            # Add a marker to a block
            push @block, $new_marker_index;
            $i++;
            $end_position = $matrix{$new_marker_index}->{position};
        }
        #-- End dynamic searching
        
        map{$matrix{$_}->{cluster} = $block_id}@block;
        $block_contents{$block_id} = [@block];
        my $block_size = scalar(@block);
        my $marker_name = $block_size == 1 ? 
                "${start_scaffold}_$start_position" : 
                "${start_scaffold}_$start_position-$end_position($block_size)";
        $bin_marker{$block_id} = [$marker_name, consensus_marker($block_id, @block)];

        print_log("$marker_name: ", join(", ", map{$matrix{$_}->{array}->[0]}@block));
    }
    return $block_id;
}

# use $out_file and %bin_marker
sub print_bin_markers{
    my $out_file = shift;
    open my $out_fh, ">", $out_file or die;
    print $out_fh $title if $title;
    for(sort{$a <=> $b}(keys %bin_marker)){
        print $out_fh join("\t", @{$bin_marker{$_}}), "\n";
    }
    close $out_fh;
}

# use %matrix and 
sub print_marker_matrix{
    my $out_file = shift;
    open my $out_fh, ">", $out_file or die;
    print $out_fh $title if $title;
    for my $block_id (sort {$a <=> $b} keys %block_contents){
        for my $index (@{$block_contents{$block_id}}){
            my $bin_marker_name = $bin_marker{$block_id}->[0];
            print $out_fh "$bin_marker_name|",join("\t", @{$matrix{$index}->{array}}), "\n";
        }
        print $out_fh join("\t", @{$bin_marker{$block_id}}), "\n";
    }
    close $out_fh;
}

sub hr{message '-' x 60}

sub main{
    message "Loading marker matrix from $infile ...";
    load_marker_matrix($infile);
    message "Done! $num_of_markers markers";

    message "Start clustering ...";
    my $num_of_bin_markers = cluster_markers;
    message "Done! $num_of_bin_markers clusters!";

    message "Print original markers append bin marker name";
    print_marker_matrix($manual_checking_file);
    message "Done!";

    message "Print bin markers ...";
    print_bin_markers($bin_markers_file);
    message "Done!";

    hr;
    message "Some log information is in file: $log_file";
    message "Marker matrix for manual checking is in file: $manual_checking_file";
    message "Final bin markers for further analysis if in file: $bin_markers_file";
    hr;
    message "Total number of markers: $num_of_markers";
    message "Bin markers: $num_of_bin_markers";
    message "$num_of_prob_select genotypes were determined based on probability";
    message "$num_of_aaabbb_type genotypes were aaabbb style";
    message "In which, $num_of_rand_select genotypes were random determined!";
    message "Marker comparison times: $total_diff_cal_times"; 
    hr;
}

main();
