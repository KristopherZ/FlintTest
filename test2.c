#include <math.h>
#include <stdio.h>
#include <flint/arb.h>
#include <flint/dirichlet.h>
#include <time.h>

int main() {

    clock_t start, end;
    start = clock();  // start timer

    //len of prime
    long lenPrime = 1000000;

    FILE *file;
    long primes[lenPrime];
    int count = 0;

    file = fopen("primes.txt", "r");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Read each long from the file
    while (fscanf(file, "%ld", &primes[count]) == 1 && count < lenPrime) {
        count++;
    }

    fclose(file);

    //init the dirichlet group G
    dirichlet_group_t G;

    //init the dirichlet character chi
    dirichlet_char_t chi;

    //precision set up
    long prec = 200;

    double alpha = 1;

    //set up sigma
    arb_t sigma;
    arb_init(sigma);
    arb_set_d(sigma, 1.1);

    //
    fmpz_t pow1;
    fmpz_init(pow1);

    //
    arb_t sum;
    arb_t num;
    arb_t den;
    arb_t b;
    arb_t term;

    //loop through all q
    for (long q = 2; q < 10000; q++) {

        dirichlet_group_init(G, q);

        dirichlet_char_init(chi, G);

        //number of terms to compute
        long len = pow(q,alpha);

        //loop through all primitive characters
        dirichlet_char_next_primitive(chi, G);
        do {
            // check if it is quadratic
            if (dirichlet_char_is_real(G, chi)) {


                arb_init(sum);

                //calculate the partial sum

                long prime = primes[0];
                long n = 0;

                while (prime < len ) {
                    long val = dirichlet_chi(G, chi, prime);
                    if (val != -1) {

                        arb_init(num);
                        arb_log_ui(num, prime, prec);


                        arb_init(den);

                        fmpz_init(pow1);

                        fmpz_set_si(pow1, prime);

                        int sign, cur;
                        if (val == 0) {
                            cur = 1;
                            sign = 1;
                        }else {
                            cur = -1;
                            sign = -1;
                        }
                        while (fmpz_cmp_si(pow1, len) < 0) {


                            arb_init(b);
                            arb_set_si(b, fmpz_get_si(pow1));

                            arb_pow(den, b, sigma, prec);


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
                            fmpz_mul_si(pow1, pow1, prime);
                            cur = cur * sign;
                        }
                    }
                    //
                    prime = primes[++n];
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

    end = clock();  // end timer

    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Elapsed time: %.3f seconds\n", elapsed_time);

    return 0;
}
