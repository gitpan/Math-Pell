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
# Project Euler , problem 57 http://projecteuler.net/index.php?section=problems&id=57

use lib './lib';
use Math::Pell;

my $result;
Math::Pell->new({D=>2})->iterate_convergents(
    1001,
    sub {
#        print join(',',@_)."\n";
        $result += length($_[0]) > length($_[1]);
        0;
    }
);

print "$result\n";