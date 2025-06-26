#include <stdio.h>
#include <flint/arb.h>


long vonMangolt1(long x) {
    for (long i = 2; i*i <= x; i++) {
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

    arb_t sigma;
    arb_init(sigma);
    arb_set_d(sigma, 1.1);

    arb_t sum;
    arb_init(sum);

    long len = 10000000;

    long prec = 100;

    for (long i = 1; i <= len; i++) {

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

            arb_add(sum, sum, term, prec);
        }

    }
    arb_print(sum);
    return 0;
}

