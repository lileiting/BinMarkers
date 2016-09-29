#!/usr/bin/env perl

use warnings;
use strict;
use FindBin;
use Getopt::Long;
use PDL::LiteF;
use PDL::Stats::Distr;

############################################################
# Part 1. Usage
############################################################

sub message {
    my @messages = @_;
    local $\ = "\n";
    print STDERR @messages;
    return 1;
}

sub usage {
    print <<"EOF";

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

  Source: https://github.com/lileiting/BinMarkers

EOF
    exit;
}

############################################################
# Part 2. Read commands
############################################################

sub read_commands {

    my $infile;
    my $threshold = 5;
    my $help;

    my $letter_for_0_0     = 'a';    # a_letter
    my $letter_for_0_1     = 'h';    # h_letter
    my $letter_for_1_1     = 'b';    # b_letter
    my $letter_for_missing = '-';

    my $error_rate_for_0_0 = 0.04;    # a_error
    my $error_rate_for_0_1 = 0.03;    # h_error
    my $error_rate_for_1_1 = 0.01;    # b_error

    GetOptions(
        "t|threshold=i"        => \$threshold,
        "h|help"               => \$help,
        "0|letter_for_0_0=s"   => \$letter_for_0_0,
        "1|letter_for_0_1=s"   => \$letter_for_0_1,
        "2|letter_for_1_1=s"   => \$letter_for_1_1,
        "m|missing=s"          => \$letter_for_missing,
        "error_rate_for_0_0=f" => \$error_rate_for_0_0,
        "error_rate_for_0_1=f" => \$error_rate_for_0_1,
        "error_rate_for_1_1=f" => \$error_rate_for_1_1
    ) or die("Error in command line arguments");

    usage if $help;

    $infile = shift @ARGV;
    usage unless $infile;
    die "$infile: NOT EXIST!" unless -e $infile;

    # Confirm error rate range is in [0,1]
    map { die unless $_ <= 1 and $_ >= 0 }
      ( $error_rate_for_0_0, $error_rate_for_0_1, $error_rate_for_1_1 );

    my %para = (
        infile    => $infile,
        threshold => $threshold,
        a_letter  => $letter_for_0_0,
        h_letter  => $letter_for_0_1,
        b_letter  => $letter_for_1_1,
        missing   => $letter_for_missing,
        a_error   => $error_rate_for_0_0,
        h_error   => $error_rate_for_0_1,
        b_error   => $error_rate_for_1_1,
        is_valid  => {
            $letter_for_0_0 => 1,
            $letter_for_0_1 => 1,
            $letter_for_1_1 => 1
        },
        is_genotype => {
            $letter_for_0_0 => 1,
            $letter_for_0_1 => 1,
            $letter_for_1_1 => 1,
            $letter_for_missing => 1
        }
    );
    return \%para;
}

############################################################
# Part 3. Load marker matrix from file
############################################################

sub load_marker_matrix {
    my $para          = shift;
    my @all_genotypes = (
        $para->{a_letter}, $para->{h_letter},
        $para->{b_letter}, $para->{missing}
    );

    my $file = $para->{infile};
    my @matrix;
    $matrix[0] = "title\tindividuals\n";    # First element is the title

    message "Loading marker matrix from $file ...";
    my $line_count     = 0;
    my $num_of_markers = 0;
    open my $fh, "<", $file or die "$file: $!";
    while ( my $marker = <$fh> ) {
        $line_count++;
        warn "Loaded $line_count lines ...\n" if $line_count % 1000 == 0;
        chomp $marker;
        my $ref = [ split /\s+/, $marker ];

        # Check if title was present
        if ( $line_count == 1 and exists $para->{is_genotype}->{$ref->[1]} ) {
            $matrix[0] = "$marker\n";
            next;
        }

        my ( $scaffold, $position ) = split /[\-_]/, $ref->[0];
        die unless $position =~ /^(\d+).*/;
        $position = $1;
        for (@{$ref}[1..$#{$ref}]) {
            die "CAUTION: $_ is not an defined genotype!"
              unless exists $para->{is_genotype}->{$_};
        }

        $num_of_markers++;
        my $index = $num_of_markers;

        #push @marker_index, $index;
        $matrix[$index]->{array}    = $ref;
        $matrix[$index]->{scaffold} = $scaffold;
        $matrix[$index]->{position} = $position;
    }
    close $fh;

    message "Done! $num_of_markers markers";
    $para->{stats}->{markers} = $num_of_markers;
    return \@matrix;
}

############################################################
# Part 4. Clustering
############################################################

sub convert_h_to_a_or_b {
    my ( $ref, $para ) = @_;

    my $a_letter = $para->{a_letter};
    my $h_letter = $para->{h_letter};
    my $b_letter = $para->{b_letter};

    my $first = ( $a_letter, $b_letter )[ int( rand( 2 ) ) ];
    my $second = $first eq $a_letter ? $b_letter : $a_letter;
    my @genotypes;
    my $i = 0;
    for (@$ref) {
        my $choice = do {
            if ( $_ eq $h_letter and $i % 2 == 0 ) {
                $first;
            }
            elsif ( $_ eq $h_letter and $i % 2 == 1 ) {
                $second;
            }
            else { $_ }
        };
        $i++;
        push @genotypes, $choice;
    }
    return @genotypes;
}

sub count_valid_genotypes {
    my @genotypes = @_;
    my $para = pop @genotypes;

    my $n = 0;
    map { $n++ if exists $para->{is_valid}->{$_} } @genotypes;
    return $n;
}

sub count_b {
    my @genotypes = @_;
    my $para      = pop @genotypes;

    my $b_letter = $para->{b_letter};
    my $n        = 0;
    map { $n++ if $_ eq $b_letter } @genotypes;
    return $n;
}

sub no_het_genotype {
    my @genotypes = @_;
    my $para = pop;

    map { return 0 if is_heterozygous( $_, $para ) } @genotypes;
    return 1;
}

sub prob_select_genotype {
    my ( $genotypes_array_ref, $valid_genotypes_array_ref, $countif_hash_ref,
        $para )
      = @_;

    my $a_letter = $para->{a_letter};
    my $h_letter = $para->{h_letter};
    my $b_letter = $para->{b_letter};
    my $a_error  = $para->{a_error};
    my $h_error  = $para->{h_error};
    my $b_error  = $para->{b_error};

    my @genotypes = grep { exists $para->{is_valid}->{$_} } @$genotypes_array_ref;
    my $genotypes_str   = join( '', @genotypes );
    my @valid_genotypes = @{$valid_genotypes_array_ref};
    my %countif         = %{$countif_hash_ref};
    my ( $first, $second, $third ) =
      sort { $countif{$b} <=> $countif{$a} } @valid_genotypes;

    my $equal_case_choose = '';
    if (    @valid_genotypes == 2
        and no_het_genotype( @valid_genotypes, $para )
        and $genotypes_str =~
        /^(($a_letter|$b_letter)\2*)(($a_letter|$b_letter)\4*)$/ )
    {
        $para->{stats}->{aaabbb_type}++;
        my ( $first_part, $first_letter, $second_part, $second_letter ) =
          ( $1, $2, $3, $4 );
        $equal_case_choose = do {
            if ( length($first_part) == length($second_part) ) {
                $para->{stats}->{rand_select}++;
                ( $first_letter, $second_letter ) [int( rand( 2 ) )];
            }
            elsif ( length($first_part) > length($second_part) ) {
                $first_letter;
            }
            else {
                $second_letter;
            }
        };

        #die "$genotypes_str => $equal_case_choose";
    }

    my @converted_genotypes =
      convert_h_to_a_or_b( $genotypes_array_ref, $para );
    my $num_of_valid_genotypes =
      count_valid_genotypes( @converted_genotypes, $para );
    my $num_of_b = count_b( @converted_genotypes, $para );

    die "ERROR: Negative values were passed to  pmf_binomial!"
      if $num_of_b < 0 or $num_of_valid_genotypes < 0;
    my $a_ex_prob =
      pmf_binomial( $num_of_b, $num_of_valid_genotypes, $a_error );
    my $h_ex_prob =
      pmf_binomial( $num_of_b, $num_of_valid_genotypes, 0.5 + $h_error / 2 );
    my $b_ex_prob =
      pmf_binomial( $num_of_b, $num_of_valid_genotypes, 1 - $b_error );

   #    my %prob = ($a_letter => $a_ex_prob,
   #                $h_letter => $h_ex_prob,
   #                $b_letter => $b_ex_prob);
   #    my $best_prob_genotype = (sort{$prob{$b} <=> $prob{$a}}(keys %prob))[0];
    my $best_prob_genotype = $h_ex_prob >= $a_ex_prob ? $h_letter  : $a_letter;
    my $tmp_prob           = $h_ex_prob >= $a_ex_prob ? $h_ex_prob : $a_ex_prob;
    $best_prob_genotype =
      $tmp_prob >= $b_ex_prob ? $best_prob_genotype : $b_letter;

    unless ( $para->{prob_fh} ) {
        open my $prob_fh, ">", prob_file_name($para);
        $para->{prob_fh} = $prob_fh;
    }
    print { $para->{prob_fh} } "$genotypes_str => ",
      join( '', @converted_genotypes ),
      " => $equal_case_choose| $best_prob_genotype\n",
      "Valid genotypes: ", scalar(@valid_genotypes),
      " N(b): ",           $num_of_b,
      " N: ",              $num_of_valid_genotypes,
      " P(a): ",           $a_ex_prob,
      " P(h): ",           $h_ex_prob,
      " P(b): ",           $b_ex_prob, "\n";
    if ($equal_case_choose) { return $equal_case_choose }
    else {
        $para->{stats}->{prob_select}++;
        return $best_prob_genotype;
    }
}

sub judge_genotype {
    my ( $genotypes, $para ) = @_;

    my @genotypes = @$genotypes;
    my %countif;
    map { $countif{$_}++ } @genotypes;
    my @all_genotypes = keys %countif;
    my @valid_genotypes = grep { exists $para->{is_valid}->{$_} } @all_genotypes;
    my $result;
    if ( @all_genotypes == 1 ) {
        return $all_genotypes[0];
    }
    elsif ( @valid_genotypes == 1 ) {
        return $valid_genotypes[0];
    }
    elsif ( @valid_genotypes == 2 or @valid_genotypes == 3 ) {
        return prob_select_genotype( \@genotypes, \@valid_genotypes, \%countif,
            $para );
    }
    else { die "More than three genotypes? @genotypes" }
}

sub create_consensus_marker {
    my ( $bin_contents, $matrix, $para ) = @_;
    my @consensus_marker;

    my $max_array_index = $#{ $matrix->[ $bin_contents->[0] ]->{array} };
    for ( my $i = 1 ; $i <= $max_array_index ; $i++ ) {
        my @genotypes = map { $matrix->[$_]->{array}->[$i] } @$bin_contents;
        push @consensus_marker, judge_genotype( \@genotypes, $para );
    }
    return @consensus_marker;
}

sub exist_different_genotypes {
    my @genotypes       = @_;
    my $para            = pop @genotypes;
    my @column          = @genotypes;
    my %genotypes       = map { $_ => 1 } @column;
    my @valid_genotypes = grep { exists $para->{is_valid}->{$_} } ( keys %genotypes );
    return scalar(@valid_genotypes) > 1 ? 1 : 0;
}

sub not_meet_the_threshold {
    my ( $candidate_marker_index, $bin_contents, $matrix, $para ) = @_;
    my $max_array_index =
      scalar( @{ $matrix->[$candidate_marker_index]->{array} } ) - 1;
    my $difference = 0;
    my $threshold  = $para->{threshold};
    for ( my $i = 1 ; $i <= $max_array_index ; $i++ ) {
        my @column = map { $matrix->[$_]->{array}->[$i] }
          ( $candidate_marker_index, @$bin_contents );
        $difference++ if exist_different_genotypes( @column, $para );
    }
    return $difference > $threshold ? 1 : 0;
}

sub print_status {
    my ( $current, $max, $para ) = @_;
    $para->{clustering_status} = 0 unless defined $para->{clustering_status};
    if ( $current / $max * 10 > $para->{clustering_status} ) {
        $para->{clustering_status}++;
        message $para->{clustering_status} . "0%";
    }
    return 1;
}

sub cluster_markers {

    # This is the main subroutine for clustering
    # Input: markers matrix[array], stats[hash], para[hash]
    # Output: bin_marker[hash]
    #my %block_contents;
    my ( $matrix, $para ) = @_;

    message "Start clustering ...";
    my $bin_markers = [];
    $bin_markers->[0] = $matrix->[0];    # First element is title

    my $log_file = log_file_name($para);
    open my $log_fh, ">", $log_file or die;

    #my @sorted_marker_index = sort
    #    {$matrix{$a}->{scaffold} cmp $matrix{$b}->{scaffold} or
    #     $matrix{$a}->{position} <=> $matrix{$b}->{position}
    #    }keys %matrix;
    my $num_of_markers      = $para->{stats}->{markers};
    my @sorted_marker_index = ( 1 .. $num_of_markers );

    my $bin_id        = 0;
    my $status_memory = 0;
    for ( my $i = 0 ; $i <= $#sorted_marker_index ; $i++ ) {
        my $first_marker_index = $sorted_marker_index[$i];

        #message "i = $first_marker_index";
        print_status( $i, $num_of_markers, $para );
        my $start_scaffold = $matrix->[$first_marker_index]->{scaffold};
        my $end_scaffold   = $start_scaffold;
        my $start_position = $matrix->[$first_marker_index]->{position};
        my $end_position   = $start_position;

        #-- Start dynanmic searching
        $bin_id++;

        #message "Bin ID: $bin_id";
        my @bin_contents = ($first_marker_index);
        for ( my $j = $i + 1 ; $j <= $#sorted_marker_index ; $j++ ) {
            my $second_marker_index = $sorted_marker_index[$j];

            #message "j = $second_marker_index";
            $end_scaffold = $matrix->[$second_marker_index]->{scaffold};

            # Stop searching with these conditions
            last if $end_scaffold ne $start_scaffold;
            last
              if not_meet_the_threshold( $second_marker_index, \@bin_contents,
                $matrix, $para );

            # Add a marker to a block
            push @bin_contents, $second_marker_index;
            $i++;
            $end_position = $matrix->[$second_marker_index]->{position};
        }

        #-- End dynamic searching

        # Assign bin ID to each marker
        map { $matrix->[$_]->{cluster} = $bin_id } @bin_contents;
        $bin_markers->[$bin_id]->{contents} = [@bin_contents];

        my $bin_size = scalar(@bin_contents);
        my $marker_name =
          $bin_size == 1
          ? "${start_scaffold}_$start_position"
          : "${start_scaffold}_$start_position-$end_position($bin_size)";
        $bin_markers->[$bin_id]->{consensus} = [
            $marker_name,
            create_consensus_marker( \@bin_contents, $matrix, $para )
        ];

        print $log_fh "$marker_name: ",
          join( ", ", map { $matrix->[$_]->{array}->[0] } @bin_contents ),
          "\n";
    }
    $para->{stats}->{bin_markers} = $bin_id;
    message "Done! $bin_id clusters!";
    return $bin_markers;
}

###########################################################
# Part 5. Print results
###########################################################

sub print_bin_markers {
    my ( $matrix, $bin_markers, $para ) = @_;

    my $bin_markers_file = bin_markers_file_name($para);
    open my $out_fh, ">", $bin_markers_file or die;
    print $out_fh $matrix->[0];
    for my $index ( 1 .. $#{$bin_markers} ) {
        print $out_fh join( "\t", @{ $bin_markers->[$index]->{consensus} } ),
          "\n";
    }
    close $out_fh;
    return 1;
}

sub print_marker_matrix {
    my ( $matrix, $bin_markers, $para ) = @_;

    my $out_file = check_file_name($para);
    open my $out_fh, ">", $out_file or die;
    print $out_fh $matrix->[0];    # title
    for my $bin_id ( 1 .. $#{$bin_markers} ) {
        for my $index ( @{ $bin_markers->[$bin_id]->{contents} } ) {
            my $bin_marker_name = $bin_markers->[$bin_id]->{consensus}->[0];
            print $out_fh "$bin_marker_name|",
              join( "\t", @{ $matrix->[$index]->{array} } ),
              "\n";
        }
        print $out_fh join( "\t", @{ $bin_markers->[$bin_id]->{consensus} } ),
          "\n";
    }
    close $out_fh;
    return 1;
}

sub get_file_name {
    my ( $suffix, $para ) = @_;
    return join( '', $para->{infile}, '.t', $para->{threshold}, $suffix );
}

sub log_file_name {
    my $para = shift;
    return get_file_name( ".log", $para );
}

sub check_file_name {
    my $para = shift;
    return get_file_name( ".check.txt", $para );
}

sub bin_markers_file_name {
    my $para = shift;
    return get_file_name( ".map", $para );
}

sub prob_file_name {
    my $para = shift;
    return get_file_name( ".prob.txt", $para );
}

sub hr { message '-' x 60; return 1; }

sub print_stats {
    my $para                 = shift;
    my $log_file             = log_file_name($para);
    my $manual_checking_file = check_file_name($para);
    my $bin_markers_file     = bin_markers_file_name($para);
    my $prob_file            = prob_file_name($para);
    my $num_of_markers       = $para->{stats}->{markers};
    my $num_of_bin_markers   = $para->{stats}->{bin_markers};
    my $num_of_prob_select   = $para->{stats}->{prob_select} // 0;
    my $num_of_aaabbb_type   = $para->{stats}->{aaabbb_type} // 0;
    my $num_of_rand_select   = $para->{stats}->{rand_select} // 0;

    hr;
    message "Some log information is in file: $log_file";
    message
      "Marker matrix for manual checking is in file: $manual_checking_file";
    message
      "Final bin markers for further analysis if in file: $bin_markers_file";
    message "Probability selection inforation in file $prob_file";
    hr;
    message "Total number of markers: $num_of_markers";
    message "Bin markers: $num_of_bin_markers";
    message
      "$num_of_prob_select genotypes were determined based on probability";
    message "$num_of_aaabbb_type genotypes were aaabbb style";
    message "In which, $num_of_rand_select genotypes were random determined!";
    hr;
    return 1;
}

sub main {
    my $para        = read_commands;
    my $matrix      = load_marker_matrix($para);
    my $bin_markers = cluster_markers( $matrix, $para );

    message "Print original markers append bin marker name";
    print_marker_matrix( $matrix, $bin_markers, $para );
    message "Done!";

    message "Print bin markers ...";
    print_bin_markers( $matrix, $bin_markers, $para );
    message "Done!";

    print_stats($para);

    return 1;
}

main() unless caller;

__END__
