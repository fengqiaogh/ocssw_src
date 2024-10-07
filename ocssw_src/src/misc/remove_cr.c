#include <stdio.h>

int main(int argc, char* argv[]) {
    int tmp;
    while ((tmp = getchar()) != EOF) {
        if (tmp != '\r') {
            putchar(tmp);
        }
    }
    return 0;
}

