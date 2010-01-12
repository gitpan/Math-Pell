# 
# This file is part of Math-Pell
# 
# This software is copyright (c) 2010 by Stefan Petrea.
# 
# This is free software; you can redistribute it and/or modify it under
# the same terms as the Perl 5 programming language system itself.
# 
use strict;
use warnings;
package Math::Pell;
our $VERSION = '0.6';
use Moose;
#use Math::Counting;
use Math::BigInt;
use POSIX qw(ceil floor);
use Carp qw/confess/;
use List::AllUtils qw/any/;

=pod

=head1 NAME

Math::Pell - Solution for Pell equations


=head1 VERSION

version 0.6

=head1 DESCRIPTION

This module solves the Pell equation by finding a minimal solution and making a method available for
generating all solutions to the Pell equation:

    x^2 - Dy^2 = 1

Of particular interest are non-trivial solutions to this equation.
The numbers x and y are integers.



=head1 SYNOPSIS

Finding the first 5 solutions of x^2 - 7y^2 = 1

    use Math::Pell;
    my $p = Math::Pell->new({D=>7});# equation is x^2 - Dy^2 = 1 where x,y are integers
    $p->find_minimal_sol;
    printf "(%d,%d)\n",$p->nth_solution($_)
    for 1..5;


    (8,3)
    (127,48)
    (2024,765)
    (32257,12192)
    (514088,194307)

=head1 SAMPLE TEST OUTPUT

    sqrt[13]=  [3,1,1,1,1,6]
    n= 2 p=4 q=1 p/q = 4
    n= 3 p=7 q=2 p/q = 3
    n= 4 p=11 q=3 p/q = 3
    n= 5 p=18 q=5 p/q = 3
    n= 6 p=119 q=33 p/q = 3
    n= 7 p=137 q=38 p/q = 3
    n= 8 p=256 q=71 p/q = 3
    n= 9 p=393 q=109 p/q = 3
    n= 10 p=649 q=180 p/q = 3
    Minimal solution -> (649,180)
    ok 14 - (649,180) is a working solution 1
    ok 15 - (1093435849,303264540) is a working solution 1
    ok 16 - (1419278889601,393637139280) is a working solution 1
    ok 17 - (1842222905266249,510940703520900) is a working solution 1
    ok 18 - (2391203911756701601,663200639532988920) is a working solution 1
    ok 19 - (3103780835237293411849,860833919173116097260) is a working solution 1
    ok 20 - (4028705132934095091878401,1117361763886065161254560) is a working solution 1
    ok 21 - (5229256158767620191964752649,1450334708690193406192321620) is a working solution 1
    ok 22 - (6787570465375238075075157060001,1882533334518107155172472208200) is a working solution 1
    ok 23 - (8810261234800900253827361899128649,2443526817869794397220462733921980) is a working solution 1
    ok 24 - (11435712295201103154229840669911926401,3171695927061658609485005456158521840) is a working solution 1
    ok 25 - (14843545748909797093290079362183781339849,4116858869799215005317139861631027426340) is a working solution 1
    ok 26 - (19266910946372621425987368782273878267197601,5343679641303454015243038055391617440867480) is a working solution 1
    ....






    sqrt[61]=  [7,1,4,3,1,2,2,1,3,4,1,14]
    ...
    Minimal solution -> (1766319049,226153980)
    ok 14 - (1766319049,226153980) is a working solution 1
    ok 15 - (22042834973108102061352541449,2822295814832482312327709940) is a working solution 1
    ok 16 - (77869358613928486808166555366140995201,9970149719303180503641083029374964080) is a working solution 1
    ....


=head1 CF2fraction($max_iter,$conv_callback,@cont_fraction)

=over

=item * @cont_fraction is the continued fraction form of the square root the non-periodic
and periodic parts are provided here.

=item * $conv_callback is a subref which is called whenever a new convergent is found and it's passed as argument
to the callback

=item * $max_iter is the maximum iteration

=back

=head1 nth_solution($i)

if a minimal solution a+b\sqrt{D} is found , then all solutions can be found by expanding (a+b\sqrt{D})^i = A+B\sqrt{D}
so (A,B) is also a solution.

=head1 cfrac(D)

computes the continued fraction representation of a square root D. it can be proved that the continued fraction will be periodic.

=head1 find_minimal_sol()

returns a minimal solution for the Pell equation by finding the first convergent of \sqrt{D} for which it's numerator and denominator are solutions to the Pell equation.

=head1 BIBLIOGRAPHY

    [1] L. Panaitopol A. Gica - O Introducere in aritmetica si Teoria Numerelor

=head1 SEE ALSO

L<http://en.wikipedia.org/wiki/Pell's_equation>

=head1 AUTHOR

Stefan Petrea, C<< <stefan.petrea at gmail.com> >>


=cut



has minimal_sol => (
    isa =>'ArrayRef[Math::BigInt]',
    is  =>'rw',
    lazy    => 1,
    default => sub {[]},
);

has $_   => (
    isa => 'Int',
    is  => 'rw',
    default => 0,
) for qw/D debug/;


sub BUILDARGS {
    my ($self,$param) = @_;

    confess "D is supposed to be square-free, you passed $param->{D}"
    if
        any { $param->{D} % $_ == 0 } 
        map { $_**2 } 2..$param->{D}/2;

    { D=> $param->{D} };
};


# gets the continued fraction of a square root
sub cfrac {
    my ($self,$D) = @_;
    my ($n,@m,@d,@a);
    my $i_rep;#index where a_n starts to repeat
    my $repeating; #shows if the sequence is repeating or not

    ($repeating,$i_rep,$n,$m[0],$d[0],$a[0]           ) =
    (0         ,    -1, 0,    0,    1,floor(sqrt($D)) );
    
    while(1) {
        ++$n;
        $m[$n] = $d[$n-1]*$a[$n-1] - $m[$n-1];
        $d[$n] = ($D-$m[$n]**2)/$d[$n-1];
        $a[$n] = floor(($a[0]+$m[$n])/$d[$n]);

        for (reverse(1..$n-1)) {
            $repeating=
            $m[$_] == $m[$n] &&
            $d[$_] == $d[$n] &&
            $a[$_] == $a[$n];
            if($repeating) {
                $i_rep = $_;
                last;
            }
        };
        last if $repeating;
    };
    pop @a;
    return @a;
}


# computing a convergent with the continued fraction representation
# using http://en.wikipedia.org/wiki/Fundamental_recurrence_formulas
# code->(p,q) is run for each of the convergents
sub CF2fraction {
    my ($self,$max_iter,$code,@b) = @_;
#    print Dumper \@b;
    my @p = ($b[0],$b[1]*$b[0]+1);
    my @q = (1    ,$b[1]        );

    $p[$_] = Math::BigInt->new($p[$_]) for qw/0 1/;
    $q[$_] = Math::BigInt->new($q[$_]) for qw/0 1/;

    return if 
    $code->($p[0],$q[0]) || 
    $code->($p[1],$q[1])  ;

    my $n=2;
    my $jump=0; # jump over the floor part

    while(1) {
        print "n= $n p=$p[$n-1] q=$q[$n-1] p/q = ".$p[$n-1]/$q[$n-1]."\n"
            if $self->debug;
        if($code) {
            last if $code->($p[$n-1],$q[$n-1]); # stop if code returns something defined or 1
        };
        $p[$n] = Math::BigInt->new(0);
        $q[$n] = Math::BigInt->new(0);

        my $j = ($n+$jump)%@b; # need to jump over the first number which is the floor of the square root
                               # and keep track of the jump over iterations
        
        $j=1,$jump++ if $j==0; # $b[0] is used when $j==0 and we want to avoid that because
                               # that's the first number in the continued fraction representation
                               # and it's not part of the periodic part

        $p[$n] = $b[$j]*$p[$n-1] + $p[$n-2];
        $q[$n] = $b[$j]*$q[$n-1] + $q[$n-2];
        last if $n++ > $max_iter;
    };
}

# returns true if parameters are a solution to the eq and false otherwise
sub is_solution {
    my ($self,$a,$b) = @_;
    return $a**2 - ($self->D * ($b**2)) == 1;
}



sub iterate_convergents {
    my ($self,$lim,$code) = @_;
    $self->CF2fraction(
        $lim,
        sub { $code->(@_) },
        $self->cfrac($self->D)
    );
}

sub find_minimal_sol {
    my ($self) = @_;
    my @min_sol;

    #search a solution to the Pell equation through the first 100 convergents of sqrt(D)
    $self->iterate_convergents(
        100,
        sub {
            if( $self->is_solution(@_) ) { # the Pell equation
                @min_sol = @_;
                return 1;# stop if minimal solution is found
            };
            0;
        }
    );
    $self->minimal_sol([@min_sol]);
    return @min_sol;
}


sub nth_solution {
    my ($self,$i) = @_;
    my ($a,$b) = @{$self->minimal_sol};
    return ($a,$b) if $i==1;

    my ($A,$B) =
    map { Math::BigInt->new($_) }
    ($a,$b);

    for(2..$i){
        my $oldA = $A->copy(); #old values for $A
        $A = $a*$A + $self->D * $b*$B;
        $B = $a*$B + $b*$oldA;
    };
    return ($A,$B);
}

1;