#include <math.h>
#include <stdio.h>
#include <flint/arb.h>
#include <flint/dirichlet.h>

int main() {
    //init the dirichlet group G
    dirichlet_group_t G;

    //init the dirichlet character chi
    dirichlet_char_t chi;

    //precision set up
    long prec = 200;

    //number of terms to compute
    long len = 10000;

    //set up sigma
    arb_t sigma;
    arb_init(sigma);
    arb_set_d(sigma, 1.1);

    //
    fmpz_t prime;
    fmpz_init(prime);

    //
    fmpz_t pow;
    fmpz_init(pow);

    //loop through all q
    for (long q = 2; q < 100000; q++) {

        dirichlet_group_init(G, q);

        dirichlet_char_init(chi, G);

        //loop through all primitive characters
        dirichlet_char_next_primitive(chi, G);
        do {
            // check if it is quadratic
            if (dirichlet_char_is_real(G, chi)) {
                arb_t sum;

                arb_init(sum);

                //calculate the partial sum

                fmpz_set_ui(prime, 2);

                while (fmpz_cmp_si(prime, len) < 0) {
                    long val = dirichlet_chi(G, chi, fmpz_get_si(prime));
                    if (val != -1) {

                        arb_t num;
                        arb_init(num);
                        arb_log_ui(num, fmpz_get_si(prime), prec);

                        arb_t den;
                        arb_init(den);

                        fmpz_init_set(pow, prime);

                        int sign, cur;
                        if (val == 0) {
                            cur = 1;
                            sign = 1;
                        }else {
                            cur = -1;
                            sign = -1;
                        }
                        while (fmpz_cmp_si(pow, len) < 0) {

                            arb_t b;
                            arb_init(b);
                            arb_set_si(b, fmpz_get_si(pow));

                            arb_pow(den, b, sigma, prec);

                            arb_t term;
                            arb_init(term);
                            arb_div(term, num, den, prec);

                            // printf("term:");
                            // arb_print(term);

                            if (cur == 1) {
                                arb_add(sum, sum, term, prec);
                                // printf("sign:+");
                            } else {
                                arb_neg(term, term);
                                arb_add(sum, sum, term, prec);
                                // printf("sign:-");
                            }
                            // printf("\n");
                            fmpz_mul(pow, pow, prime);
                            cur = cur * sign;
                        }
                    }
                    //
                    fmpz_nextprime(prime, prime, 1);
                }

                //print
                printf("q = %ld @", q);
                dirichlet_char_print(G, chi);
                printf("\n");
                arb_print(sum);
                printf("\n");
            }
        } while (dirichlet_char_next_primitive(chi, G) >= 0);
    }

    return 0;
}
