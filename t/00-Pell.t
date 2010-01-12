use Test::More;
use Test::Deep;
use lib './lib';
use Math::Pell;

my $sols = 20;#number of solutions to test;

sub test_min_sol {
    my ($D,$a,$b) = @{$_[0]};
    cmp_deeply(
        [Math::Pell->new({D=>2})->find_minimal_sol()],
        [map { Math::BigInt->new($_) } (3,2)],
        "D=$D has minimal solution ($a,$b)"
    );
};

test_min_sol($_) for (
    [2,3,2],
    [3,2,1],
    [5,9,4],
    [7,8,3],
    [11,10,3],
    [13,649,180],
    [15,4,1],
    [17,33,8],
    #[28,127,24],
    [29,9801,1820],
    [61,1766319049,226153980],
    [97,62809633,6377352],
    [101,201,20]
);

sub D_eq { 
    my $D = $_[0];
    my $n = Math::Pell->new({D=>$D});

    print "sqrt[$D]=  [".join(',',$n->cfrac($D))."]\n";
    $n->find_minimal_sol;
    print "Minimal solution -> (".(join ',',@{$n->minimal_sol}).")\n";
    for my $i(1..$sols) {
        my ($A,$B) = $n->nth_solution($i);
        if($D==7 && $i==2) {
            ok($n->is_solution(127,49)==0,'is_solution works fine');

            my $no_other = 1;
            for my $x     ($n->minimal_sol->[0]+1..-1+$A) {
                for my $y ($n->minimal_sol->[1]+1..-1+$B) {
                    if($n->is_solution($x,$y)==1) {
                        $no_other = 0;
                        last;
                    }
                }
            };
            ok($no_other,'ensure a solution wasn\'t skipped');
        }
        ok( $n->is_solution($A,$B) , "($A,$B) is a working solution $LHS");
    }
};

D_eq $_ for qw/2 3 5 7 11 13 17 19 43 61/;


#my $p = Math::Pell->new({D=>19});
#
#$p->CF2fraction(
#    40,
#    sub { 0 },
#    $p->cfrac($p->D)
#);


done_testing();



# Table with minimal solutions and D from x^2 - Dy^2 = 1 (Pell equation)
# (can be used for testing)
# source: http://mathworld.wolfram.com/PellEquation.html
#|==========================================================|
#|D       x       y       D       x       y
#|===========================================================
#|2       3       2       54      485     66
#|3       2       1       55      89      12
#|5       9       4       56      15      2
#|6       5       2       57      151     20
#|7       8       3       58      19603   2574
#|8       3       1       59      530     69
#|10      19      6       60      31      4
#|11      10      3       61      1766319049      226153980
#|12      7       2       62      63      8
#|13      649     180     63      8       1
#|14      15      4       65      129     16
#|15      4       1       66      65      8
#|17      33      8       67      48842   5967
#|18      17      4       68      33      4
#|19      170     39      69      7775    936
#|20      9       2       70      251     30
#|21      55      12      71      3480    413
#|22      197     42      72      17      2
#|23      24      5       73      2281249 267000
#|24      5       1       74      3699    430
#|26      51      10      75      26      3
#|27      26      5       76      57799   6630
#|28      127     24      77      351     40
#|29      9801    1820    78      53      6
#|30      11      2       79      80      9
#|31      1520    273     80      9       1
#|32      17      3       82      163     18
#|33      23      4       83      82      9
#|34      35      6       84      55      6
#|35      6       1       85      285769  30996
#|37      73      12      86      10405   1122
#|38      37      6       87      28      3
#|39      25      4       88      197     21
#|40      19      3       89      500001  53000
#|41      2049    320     90      19      2
#|42      13      2       91      1574    165
#|43      3482    531     92      1151    120
#|44      199     30      93      12151   1260
#|45      161     24      94      2143295 221064
#|46      24335   3588    95      39      4
#|47      48      7       96      49      5
#|48      7       1       97      62809633        6377352
#|50      99      14      98      99      10
#|51      50      7       99      10      1
#|52      649     90      101     201     20
#|53      66249   9100    102     101     10
#|-------------------------------------------------
