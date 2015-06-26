BinMarkers
======

Overview
------

BinMarkers is a program for creating bin markers from a large set of SNP markers based on their positions on scaffolds/contigs. The purpose of this program is to reduce the number of markers for constructing genetic maps using [Joinmap](http://www.kyazma.nl/index.php/mc.JoinMap/), which limited to about six thousands markers.

* Author: [Leiting Li](https://github.com/lileiting), [Chenyong Miao](https://github.com/freemao)
* Email: lileiting@gmail.com
* LICENCE: [BSD](http://opensource.org/licenses/bsd-license.php)

Installation
------

    git clone https://github.com/lileiting/BinMarkers.git

Dependencies
------

[Perl](http://www.perl.org) programming language and [pmf_binomial](http://pdl-stats.sourceforge.net/Distr.htm#pmf_binomial) function in Perl module [PDL::Stats](https://metacpan.org/pod/PDL::Stats)

Dependencies of PDL::Stats

- [PDL](https://metacpan.org/pod/PDL), The Perl Data Language
- [GSL](http://www.gnu.org/software/gsl/), GNU Scientific Library

How to install PDL::Stats

    perl -MCPAN -e shell
    install PDL
    install PDL::Stats

or

    cpan PDL
    cpan PDL::Stats

or (if without permission to install)

    sudo cpan PDL
    sudo cpan PDL::Stats

or

    cpanp -i PDL
    cpanp -i PDL::Stats

or

    cpanm PDL
    cpanm PDL::Stats

How to install GSL

Use [brew](http://brew.sh) in Mac OS

    brew install gsl

How to use BinMarkers
------

Input data format

Input data should be a matrix, one marker per line, one individual per column.
The first line is the title containing individual names, the first column is 
marker names.

Marker names should contain two parts. First is the scaffold/contig/chromosome 
name, second is the position, and they were joined by dash or underscore, e.g.
scaffold1\_12345, scaffold1-12345

The default genotypes are "a" for 0/0, "h" for 0/1, "b" for 1/1, "-" for 
missingthese could be customized using -0, -1, -2 and -m options.

Examples

    perl binmarkers.pl markers.map -t 5
    perl binmarkers.pl markers.map --threshold 5
    perl binmarkers.pl markers.map -t 5 -0 a -1 h -2 b -m '-'
    perl binmarkers.pl markers.map --threshold=5 --letter_for_0_0=a --letter_for_0_1=h --letter_for_1_1=b --missing='-'
    perl binmarkers.pl markers.map -t 5 -0 a -1 h -2 b -m '-' --error_rate_for_0_0 0.04 --error_rate_for_0_1 0.03 --error_rate_for_1_1 0.01
