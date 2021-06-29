#ifndef AFFICHAGES_H
#define AFFICHAGES_H
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "projet.h"


void print_tab_BYTE(BYTE *tab, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        printf("{");
        for (size_t j = 0; j < 8; ++j)
        {
            printf("%d", (tab[i] >> (7 - j)) & 1);
        }
        printf("} ");
    }
    printf("\n\n");
}

void print_tab(BYTE *tab, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
      if(i%8==0) printf("\n");
        printf("%2d ", tab[i]);
    }
    printf("\n\n");
}

void print_tab_int(int *tab, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        printf("%d ", tab[i]);
    }
    printf("\n\n");
}

void print_byte_binary(BYTE x)
{
    for (int i = 7; i >= 0; i--)
    {
        printf("%d", (x >> i) & 1 ? 1 : 0);
    }
}

void print_tab_binary(BYTE *tab, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        print_byte_binary(tab[i]);
        printf(" ");
    }
    printf("\n");
}

void print_tab_hex(BYTE *tab, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        printf("%x", tab[i]);
        printf(" ");
    }
    printf("\n");
}

void print_tab_tab(BYTE tab[16][48])
{
    printf("\n");
    for (int i = 0; i < 16; ++i)
    {
        for (int j = 0; j < 48; ++j)
        {
            printf("%d ", tab[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void print_tab_tab_hex(BYTE (*tab)[48], int rounds)
{
    printf("\n");
    BYTE res[48];
    for (int i = 0; i < rounds; ++i)
    {
        pack8(res, tab[i]);
        print_tab_hex(res, 6);
    }
    printf("\n\n");
}


#endif //AFFICHAGES_H
