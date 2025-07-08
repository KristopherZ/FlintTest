#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/dirichlet.h>
#include <flint/long_extras.h>
#include <time.h>

int main() {
    clock_t start, end;
    start = clock(); // start timer

    //len of prime
    long lenPrime = 50000;

    FILE *file;
    long *primes = (long *) malloc(lenPrime * sizeof(long));
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

    printf("Elapsed time after file reading: %.3f seconds\n", (double) (clock() - start) / CLOCKS_PER_SEC);

    //precision set up
    long prec = 100;
    double alpha = 0.3;

    //set up lambda
    arb_t lambda;
    arb_init(lambda);
    arb_set_d(lambda, 1);

    // setting up variables
    arb_t logq;
    arb_t sigma;
    arb_t sum;
    arb_t logp;
    arb_t p;
    arb_t psigma;
    arb_t zeta_term;
    arb_t one; // equal to 1
    arb_init(one);
    arb_set_ui(one, 1);
    arb_t temp1; //temp variables for calculations
    arb_t temp2;
    arb_t l_term;
    arb_t term;

    FILE *outfile;
    outfile = fopen("output.csv", "w");
    if (outfile == NULL) {
        perror("Error opening file");
        return 1;
    }

    //loop through all q
    long qMax = 1000000;
    for (long q = -qMax; q <= qMax; q++) {
        if (q % 100000 == 0) {
            printf("Elapsed time when q=%d: %.3f seconds\n", q, (double) (clock() - start) / CLOCKS_PER_SEC);
        }

        if (q % 4 == 0 || q % 4 == 1) {
            //number of terms to compute
            long len = pow(abs(q), alpha);

            // sets sigma= 1 + lambda/log(q)
            arb_init(sigma);
            arb_init(logq);
            arb_init(temp1);
            arb_log_ui(logq, abs(q), prec);
            arb_div(temp1, lambda, logq, prec);
            arb_add(sigma, one, temp1, prec);

            //calculate the partial sum
            arb_init(sum);

            // loop over all primes < q^alpha
            long prime = primes[0];
            int primeIndex = 0;
            while (prime < len) {
                int chi = z_kronecker(q, prime); //the value of chi(prime)

                // Calculating log(p)
                arb_init(logp);
                arb_log_ui(logp, prime, prec);

                // Calculating p^sigma
                arb_init(p);
                arb_set_ui(p, prime);
                arb_init(psigma);
                arb_pow(psigma, p, sigma, prec);

                // The infinite sum from the p terms for zeta is 1/(p^sigma -1)
                arb_init(zeta_term);
                arb_init(temp1);
                arb_sub(temp1, psigma, one, prec);
                arb_inv(zeta_term, temp1, prec);

                // The infinite sum from the p terms for L is chi(p)/(p^sigma -chi(p))
                arb_init(l_term);
                if (chi == 1) {
                    arb_set(l_term, zeta_term);
                } else if (chi == -1) {
                    arb_init(temp1);
                    arb_init(temp2);
                    arb_add(temp1, psigma, one, prec);
                    arb_inv(temp2, temp1, prec);
                    arb_neg(l_term, temp2);
                } else {
                }

                // add log(p)*(zeta_term + l_term) to the sum
                arb_init(term);
                arb_init(temp1);
                arb_add(term, zeta_term, l_term, prec);
                arb_mul(temp1, term, logp, prec);
                arb_add(sum, sum, temp1, prec);

                prime = primes[primeIndex++];
            }

            char *output = (char *) malloc(50 * sizeof(char));
            output = arb_get_str(sum, 40, 0);
            fprintf(outfile, "%ld,%s\n",q,output);
        }
    }

    fclose(outfile);

    end = clock(); // end timer

    double elapsed_time = (double) (end - start) / CLOCKS_PER_SEC;
    printf("Elapsed time: %.3f seconds\n", elapsed_time);

    return 0;
}
