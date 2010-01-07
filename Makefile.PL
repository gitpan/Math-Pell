
use strict;
use warnings;



use ExtUtils::MakeMaker;

WriteMakefile(
  DISTNAME  => 'Math-Pell',
  NAME      => 'Math::Pell',
  AUTHOR    => "Stefan\ Petrea\ \<stefan\.petrea\@gmail\.com\>",
  ABSTRACT  => "Solution\ for\ Pell\ equations",
  VERSION   => '0.5',
  EXE_FILES => [ qw() ],
  (eval { ExtUtils::MakeMaker->VERSION(6.31) } ? (LICENSE => 'perl') : ()),
  PREREQ_PM    => {
    "Test::More" => '0',
    "Math::BigInt" => '0',
    "List::AllUtils" => '0.02',
    "Test::Deep" => '0.106',
    "Moose" => '0.92',
  },
  test => {TESTS => 't/*.t'}
);

