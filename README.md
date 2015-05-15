BinMarkers
======

### Dependencies
- [pmf_binomial](http://pdl-stats.sourceforge.net/Distr.htm#pmf_binomial) function in Perl module [PDL::Stats](https://metacpan.org/pod/PDL::Stats)

#### Dependencies of PDL::Stats
- [PDL](https://metacpan.org/pod/PDL), The Perl Data Language
- [GSL](http://www.gnu.org/software/gsl/), GNU Scientific Library

#### How to install PDL::Stats
- perl -MCPAN -e shell
- install PDL
- install PDL::Stats

### How to use BinMarkers
#### Example
`perl markers.map -t 5`
`perl markers.map -t 5 -0 a -1 h -2 b -m '-'`
`perl markers.map -t 5 -0 a -1 h -2 b -m '-' --error_rate_for_0_0 0.04 --error_rate_for_0_1 0.03 --error_rate_for_1_1 0.01`

#### Usage
`perl binmarkers.pl [OPTIONS] <MARKER_MATRIX>

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
  
`
