#include <stdio.h>
#include <flint/arb.h>
#include <flint/dirichlet.h>

/* An alternate version of von mangolt function,
 * it returns the exp of von mangolt lambda
*/
long vonMangolt1(long x) {
    // find the first prime divisor
    for (long i = 2; i*i <= x; i++) {
        // check if it is a prime power
        if (x % i == 0) {
            long remainder = x;
            while (remainder > 1) {
                if (remainder % i == 0) {
                    remainder = remainder / i;
                }else {
                    return 1;
                }
            }
            return i;
        }
    }
    return 1;
}

int main() {

    //init the dirichlet group G
    dirichlet_group_t G;

    //init the dirichlet character chi
    dirichlet_char_t chi;

    //number of terms to compute
    long len = 100000;

    //precision set up
    long prec = 100;

    //set up sigma
    arb_t sigma;
    arb_init(sigma);
    arb_set_d(sigma, 1.1);

    //loop through all q
    for (long q = 2; q < 500; q++) {

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
                for (long i = 1; i <= len; i++) {

                    long val = dirichlet_chi(G, chi, i);

                    if (val != - 1) {
                        long n = vonMangolt1(i);
                        if (n != 1) {
                            arb_t num;
                            arb_init(num);
                            arb_log_ui(num, n,  prec);

                            arb_t den;
                            arb_init(den);

                            arb_t b;
                            arb_init(b);
                            arb_set_si(b, i);

                            arb_pow(den, b, sigma, prec);

                            arb_t term;
                            arb_init(term);
                            arb_div(term, num, den, prec);

                            if (val == 0) {
                                arb_add(sum, sum, term, prec);
                            }else {
                                arb_neg(term, term);
                                arb_add(sum, sum, term, prec);
                            }


                        }
                    }

                }
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

