#!/usr/bin/perl

use warnings;
use strict;

#--Preset--#
my $author = "Leiting Li";
my $version = "20150511";

sub message {local $\ = "\n"; print STDERR @_}
sub usage{
    message "\nperl $0 <Marker_matrix> [NUM]";
    message "\n  NUM  allowed number of genotype differences\n";
    exit;
}
usage unless @ARGV;

my ($file,$threshold) = @ARGV;
$threshold = 1 unless $threshold;

my @valid_genotype = qw/ll lm nn np hh hk kk h a b/;
my @missing_genotype = qw/-- .. - ./;
my @all_genotypes = (@valid_genotype, @missing_genotype);
my %is_valid = map{$_, 1}@valid_genotype;
my %is_missing = map{$_, 1}@missing_genotype;
my %is_genotype = map{$_, 1}@all_genotypes;

my %matrix;
my @marker_index;
my $num_of_markers;
my @seeds;
my $count;
my %cal_times;

my $log_file = "$file.bin_markers_$threshold.log";
my $out_file = "$file.bin_markers_$threshold.txt";
my $seed_file = "$file.bin_markers_$threshold.seeds.txt";

my %seeds;
sub new_cluster{
    my $cluster = shift;
    die if $seeds{$cluster};
    $seeds{$cluster}++;
    push @seeds, $cluster;
    return $cluster;
}

my $status_memory = 0;
sub print_status{
    my ($current, $total) = @_;
    if($current / $total * 10 > $status_memory){
        $status_memory++;
        message "${status_memory}0%";
    }
}

sub difference{
    my ($a, $b) = @_;
    $cal_times{$b}++;
    my $difference = 0;
    my @a = @{$matrix{$a}->{array}};
    my @b = @{$matrix{$b}->{array}};
    for my $index (1..$#a){
         my $element_a = $a[$index];
         my $element_b = $b[$index];
         #--#
         #next unless $is_valid{$element_a};
         #next unless $is_valid{$element_b};
         #next unless $element_a eq $element_b;
         #----#
         if($is_missing{$element_a} or 
             $is_missing{$element_b} or 
             $element_a ne $element_b){
            $difference++;
         }
    }
    #p "Difference for $a and $b is $difference";
    return $difference;
}

sub with_title{
    my $array = shift;
    if($is_genotype{$array->[1]}){
        return 0
    }else{
        return 1
    }
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

sub print_markers_to_file{
    my $file = shift;
    my @indexes = @_;
    open my $fh, "> $file" or die;
    my $num_of_ind = scalar(@{$matrix{1}->{array}}) - 1;
    print $fh "marker_index, cluster, marker_name, ",
              join(", ", map{"ind_$_"}(1..$num_of_ind)),
              "\n";
    for my $marker_index(@indexes){
        print $fh join(", ", $marker_index, 
                  $matrix{$marker_index}->{cluster}, 
                  @{$matrix{$marker_index}->{array}}), 
                  "\n";
    }
}

open my $log_fh, "> $log_file" or die;
sub print_log{
    print $log_fh @_,"\n";
}

sub load_marker_matrix{
    my $file = shift;
    open my $fh, $file or die;
    while(my $marker = <$fh>){
        chomp $marker;
        my $ref = [split /\t/, $marker];
        next if with_title $ref;
        checking_genotypes $ref;
        my $index = ++$num_of_markers;
        push @marker_index, $index;
        $matrix{$index}->{array} = $ref;
        $matrix{$index}->{missing} = count_missing $ref;
    }
}

sub cluster_markers{
    my @sorted_marker_index = sort
        {$matrix{$b}->{missing} <=> $matrix{$a}->{missing}}
        keys %matrix;

    for my $marker_index (@sorted_marker_index){
        $count++;
        print_status($count, $num_of_markers);
        my $cluster = 0;
        if (@seeds){
            my @check_seeds = @seeds;
            for my $seed (@check_seeds){
                if(&difference($seed, $marker_index) <= $threshold){
                    # Redundant
                    $cluster = $seed;
                    last;
                }
            }
        }
        if(@seeds == 0 or $cluster == 0){
            $cluster = new_cluster $marker_index;
            $cal_times{$marker_index} = 0 
                unless $cal_times{$marker_index};
        }
        $matrix{$marker_index}->{cluster} = $cluster;
        print_log("$count: Num_of_seeds => ", scalar(@seeds),
                  ", Cal_times => $cal_times{$marker_index}",
                  ", Cluster of $marker_index => ", 
                  $matrix{$marker_index}->{cluster});
    }
}

sub main{
    message "Loading marker matrix ...";
    load_marker_matrix($file);
    message "Done! $num_of_markers markers";

    message "Start clustering ...";
    cluster_markers;
    message qq/Done! /. @seeds." clusters!";

    message qq/Print all markers .../;
    print_markers_to_file($out_file, @marker_index);
    message qq/Done/;

    message qq/Print seed markers .../;
    print_markers_to_file($seed_file, @seeds);
    message qq/Done/;
}

main();
