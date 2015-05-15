BinMarkers
======

Dependencies
------
[pmf_binomial](http://pdl-stats.sourceforge.net/Distr.htm#pmf_binomial) function in Perl module [PDL::Stats](https://metacpan.org/pod/PDL::Stats)

Dependencies of PDL::Stats

- [PDL](https://metacpan.org/pod/PDL), The Perl Data Language
- [GSL](http://www.gnu.org/software/gsl/), GNU Scientific Library

How to install PDL::Stats

    perl -MCPAN -e shell
    install PDL
    install PDL::Stats

How to use BinMarkers
------

Examples

    perl binmarkers.pl markers.map -t 5
    perl binmarkers.pl markers.map --threshold 5
    perl binmarkers.pl markers.map -t 5 -0 a -1 h -2 b -m '-'
    perl binmarkers.pl markers.map --threshold=5 --letter_for_0_0=a --letter_for_0_1=h --letter_for_1_1=b --missing='-'
    perl binmarkers.pl markers.map -t 5 -0 a -1 h -2 b -m '-' --error_rate_for_0_0 0.04 --error_rate_for_0_1 0.03 --error_rate_for_1_1 0.01
