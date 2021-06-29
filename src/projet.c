#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "affichages.h"
#include "projet.h"







/**************************************************************************/
/*****************************   EXERCICE 1   *****************************/
/**************************************************************************/




/*
 * setkey code l'étape de diversification des clefs.
 * A partir d'une clef de 64 bits (key),
 * on remplit le tableau des clefs de tour 48 bits (KS).

 * Tout d'abord, on applique la permutation PC1 :
   elle donne en sortie une clef de 56 bits (on perd les bits de parité).
 * Puis, on sépare ce tableau en deux parties de 28 bits chacunes.
   On effectue alors un décalage à gauche d'un ou deux crans (selon shifts) dans
 chaque partie.
 * Enfin, on applique la permutation CP2, qui renvoie une clef de tour de 48
 bits.
 */
void setkey(BYTE *pkey, int rounds, BYTE (*KS)[48])
{
    // printf("\n\n> Set Key\n");
    int i, j, k, t1, t2;
    static BYTE key[64] = {0};
    static BYTE CD[56] = {0};
    unpack8(pkey, key);

    for (i = 0; i < 56; i++)
    {
        // On remplit CD avec la clef dans l'ordre de la permutation PC1
        CD[i] = key[PC1[i] - 1];
    }

    for (i = 0; i < rounds; i++)
    {
        // on applique le décalage de shifts[i] crans
        for (j = 0; j < shifts[i]; j++)
        {
            t1 = CD[0];  // pointe sur C0
            t2 = CD[28]; // pointe sur D0
            for (k = 0; k < 27; k++)
            {
                // on fait "tourner" en faisant i<--i+1
                CD[k] = CD[k + 1];
                CD[k + 28] = CD[k + 29];
            }
            CD[27] = t1;
            CD[55] = t2;
        }

        for (k = 0; k < 48; k++)
        {
            // KS est le tableau des 16 clefs de tours (chacune de 48 bits)
            KS[i][k] = CD[PC2[k] - 1];
        }
    }
}


/*
 * reverse_key_schedule permet de retrouver l'emplacement des 48 bits
   d'une clef de tour Kr dans la clef d'origine.
 * kSortie[i] indique que le ième bit de la clef d'origine
   est le kSortie[i]ème bit de la clef de tour Kr.

 * Le paramètre Kr est la clef de tour r.
 * IDC est le tableau de l'indice des bits de Kr dans la clef originale
 * kSortie donne la clef potentielle, avec des 0 aux différents place des
   bits "perdus" au cours des permutations PC1 et PC2.

 * La fonction reprend les étapes de set_key, à cela près que la clef d'origine
   utilisée est remplie de 1 à 64, et qu'au lieu de remplir KS, on va chercher
   les indices IDC des bits qui ont été "perdus" au cours de la diversification.
 * Ainsi, après l'action de PC1, on aura perdu les bits de parités,
   donc les multiples de 8 n'apparaîtront plus parmi les indices.
 * Après celle de PC2, on aura perdu 8 autres bits.
 * Ensuite, on recopie dans kSortie[IDC[i]] les indices des bits de Kr dans la
 clef originale.
 * IDC faisant 48 bits, lorsqu'on remplit kSortie, il restera des 0.
   Ces "trous" sont en fait les bits qui ont été perdus après PC1 et PC2 :
   On y compte les 8 bits de parités en positions 7,15,..., et d'autres, qui
 sont donc ceux sur lesquels il faudra lancer une attaque brute force.
 */
void reverse_key_schedule(int rounds,
                          BYTE kSortie[64]) // BYTE KR[48] était en deuxième
                                            // paramètre mais ne servait pas
{
    printf("\n\n> Set Key\n");
    int i, j, k, t1, t2;
    int key[64] = {0};
    BYTE IDC[48] = {0};
    for (i = 0; i < 64; i++)
    {
        key[i] = i + 1;
    }
    BYTE CD[56] = {0};

    for (i = 0; i < 56; i++)
    {
        CD[i] = key[PC1[i] - 1];
    }

    for (i = 0; i < rounds; i++)
    {
        for (j = 0; j < shifts[i]; j++)
        {
            t1 = CD[0];
            t2 = CD[28];
            for (k = 0; k < 27; k++)
            {
                CD[k] = CD[k + 1];
                CD[k + 28] = CD[k + 29];
            }
            // ils sont modifiés ici
            CD[27] = t1;
            CD[55] = t2;
        }
    }
    for (k = 0; k < 48; k++)
    {
        // KS est le tableau des 16 clefs de tours (chacune de 48 bits)
        IDC[k] = CD[PC2[k] - 1];
    }
    printf("IDC : ");
    print_tab(IDC, 48);
    for (i = 0; i < 64; i++)
    {
        kSortie[i] = 0;
    }
    for (i = 0; i < 48; i++)
    {
        // ici peut-être à modifier...(?)
        // kSortie[IDC[i]-1]=KR[i];
        kSortie[IDC[i] - 1] = key[i];
    }
    printf("Clef en termes d'indices dans la clef de tour\n");
    print_tab(kSortie, 64);
}


/*
 * indiceBytesToFindKey permet de stocker dans 'indices' les indices des bits
 perdus dans la diversification de la clef, et qui ne sont pas des bits de
 parité.
 * Ce sont les bits sur lesquels il nous faudra faire une attaque brute force.
 * On lance simplement reverse_key_schedule, qui nous donne une clef "à trous"
   et on retient les indices des 8 bits perdus après l'action de PC2
   (on fait attention à ne pas prendre les bits de parité aux emplacements
 7,15,...).
 */
void indiceBytesToFindInKey(BYTE clefATrou[64], int indices[8])
{
    int i = 0, j = 0;
    while (j < 8 && i < 64)
    {
        if (clefATrou[i] == 0 && (i + 1) % 8 != 0)
        {
            indices[j] = i;
            ++j;
        }
        ++i;
    }
    // indices contient donc les indices de la clef à chercher par force brute
}


/*
 * completerAvecBitsDeParite complète les bits en position 7,15,... pour que
 chaque octet ait une parité impaire.
 */
void completerAvecBitsDeParite(BYTE key[64])
{ // À tester
    for (int i = 0; i < 8; ++i)
    {
        BYTE sum = 0;
        for (int j = 0; j < 7; ++j)
        {
            sum ^= key[8 * i + j];
        }
        key[8 * i + 7] = sum ^ 1;
    }
}



/**************************************************************************/
/*****************************   EXERCICE 2   *****************************/
/**************************************************************************/


/*
 En pratique, le chiffrement consiste à effectuer une permutation initiale IP, puisà lancer des chiffrements de Feistel, et à appliquer la permutation inverse de l’initiale IP⁻¹. Conceptuellement, le chiffré c obtenu à partir du message m se construiten calculant : c= (IP◦F(K0,···,Kn)◦IP−1)(m), où F(K0,···,Kn) représente le chiffrement de Fesitel effectué avec, dans l’ordre, les clefs K0 à Kn.
  
Le déchiffrement consiste donc à appliquer la permutation initiale (qui s'annule avec l'inverse), puis à déchiffrer Feistel, c'est-à-dire, à utiliser le même algorithme que celui du chiffrement de Feistel, en prenant les clefs dans l'ordre inverse. Enfin, il suffit d'appliquer la permutation inverse de l'initiale, pour réobtenir le message original !

On a
(IP◦F(Kn,...,K0)◦IP⁻¹)(c) = (IP◦F(Kn,···,K0)◦IP⁻¹)◦(IP◦F(K0,···,Kn)◦IP⁻¹)(m)
                           = (IP◦F(Kn)◦···◦F(K0)◦IP⁻¹)◦(IP◦F(K0)◦···◦F(Kn)◦IP⁻¹)(m)
                           = (IP◦F(Kn)◦···◦F(K0)◦(IP⁻¹◦IP)◦F(K0)◦···◦F(Kn)◦IP⁻¹)(m)
                           = (IP(◦F(Kn)◦(···(◦F(K0)◦F(K0))◦···)◦F(Kn))◦IP⁻¹)(m)


Or F(Ki)◦F(Ki) =Id.
D’où (IP◦F(Kn,...,K0)◦IP⁻¹)(c) = (IP◦IP⁻¹)(m)
                                = Id(m)
                                = m.

De manière plus précise, en notant F(Ki) = fKi, on suppose que l’on effectue k tours pour les opérations de chiffrement et déchiffrement. De la formule de chiffrement Li=Ri−1, Ri=Li−1⊕fKi(Ri−1), on déduit la formule de déchiffrement Ri−1=Li, Li−1=Ri⊕fKi(Li).
Le chiffrement produit en sortie IP⁻¹*Rk*Lk. On aboutit avec le déchiffrement à IP⁻¹*R0*L0 en appliquant la formule de déchiffrement k fois, puis on applique la permutation IP et l’on échange R et L pour obtenir L0R0.

*/

  
/**************************************************************************/
/*****************************   EXERCICE 3   *****************************/
/**************************************************************************/


/*
 * unpack8 permet de transformer un tableau de 8 octets codant 8 bits
   par un tableau de 64 octets codant 1 bit par octet.
 * packed est un tableau de 8 octets
 * binary est un tableau de 64 octets

 * Il s'agit de mettre chaque bit de packed dans un octet de binary.
 */
void unpack8(BYTE *packed, BYTE *binary)
{
    memset(binary, 0, 64 * sizeof(BYTE));
    int index_b = 0;
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            // printf(" %d ; ", ((packed[i]) >> (7 - j)) & 1);
            binary[index_b] = ((packed[i]) >> (7 - j)) & 1;
            ++index_b;
        }
    }
}

/*
 * pack8 permet de transformer un tableau de 64 octets codant 1 bit par octet
   par un tableau de 8 octets codant 8 bits.
 * packed est un tableau de 8 octets
 * binary est un tableau de 64 octets

 * Il s'agit de mettre chaque octet de binary dans un bit de packed.
 */
void pack8(BYTE *packed, BYTE *binary)
{
    memset(packed, 0, 8 * sizeof(BYTE));
    int index_b = 0;
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            packed[i] = (packed[i] << 1) ^ (binary[index_b] & 1);
            ++index_b;
        }
    }
}


/**************************************************************************/
/*****************************   EXERCICE 4   *****************************/
/**************************************************************************/




/*
 * DES est la fonction qui code le chiffrement DES.
 * Elle effectue trois grandes étapes.
 * Premièrement, on effectue une permutation initiale fixe des blocs (IP).
 * Deuxièmement, pendant tounds tours, on découpe le bloc de 64 bits
   en deux de 32 bits, et on échange les blocs selon un schéma de Feistel.
 * Finalement, on effectue une permutation finale (RFP), qui est l'inverse de
 IP.
 */
void DES(BYTE *in, BYTE *out, int rounds, BYTE (*KS)[48])
{
    int k;
    int t;
    static BYTE block[64] = {0};
    static BYTE LR[64] = {0};
    static BYTE preS[48] = {0};
    static BYTE f[64] = {0};

    unpack8(in, block);

    for (int j = 0; j < 64; j++)
    {
        LR[j] = block[IP[j] - 1]; // permutation initiale
    }

    for (int i = 0; i < rounds; i++) // rounds tours de Feistel
    {
        for (int j = 0; j < 48; j++)
        {
            preS[j] =
                LR[E[j] + 31] ^ KS[i][j]; // XOR avec la sous clé   (la fonction
                                          // fKr du sujet, page 4)
        }

        for (int j = 0; j < 8; j++) // traitement par blocs
        {
            k = 6 * j;

            t = preS[k];
            t = (t << 1) | preS[k + 5];
            t = (t << 1) | preS[k + 1];
            t = (t << 1) | preS[k + 2];
            t = (t << 1) | preS[k + 3];
            t = (t << 1) | preS[k + 4];
            t = Sb[j][t]; // Sbox

            k = 4 * j;

            f[k] = (t >> 3) & 1;
            f[k + 1] = (t >> 2) & 1;
            f[k + 2] = (t >> 1) & 1;
            f[k + 3] = t & 1;
        }

        for (int j = 0; j < 32; j++) // croisement schema de Feistel
        {
            t = LR[j + 32];
            LR[j + 32] = LR[j] ^ f[P[j] - 1];
            LR[j] = t;
        }
    }

    BYTE tmp[32] = {0};

    // Exchange L and R: LR -> RL
    memcpy(tmp, LR, 32);
    memcpy(LR, LR + 32, 32);
    memcpy(LR + 32, tmp, 32);

    for (int j = 0; j < 64; j++)
    {
        block[j] = LR[RFP[j] - 1]; // permutation finale
    }

    pack8(out, block);
}


/*
 * reverse_subkeys permet d'inverser l'ordre des sous-clefs
   obtenues lors de la diversification de la clef (avec set_key).
 */
void reverse_subkeys(BYTE (*KS)[48], int s)
{
    BYTE tmp[48] = {0};
    for (int i = 0; i < s / 2; i++)
    {
        memcpy(tmp, KS[i], 48);
        memcpy(KS[i], KS[s - i - 1], 48);
        memcpy(KS[s - i - 1], tmp, 48);
    }
}


/*
 * des_encrypt est la fonction du chiffrement DES.
 * Après la diversification des clefs, on applique la fonction DES.
 */
void des_encrypt(BYTE *plaintext, BYTE *key, BYTE *ciphertext, int rounds)
{
    BYTE(*KS)[48] = malloc(48 * rounds);
    setkey(key, rounds, KS);
    DES(plaintext, ciphertext, rounds, KS);
    free(KS);
}


/*
 * des_decrypt est la fonction du déchiffrement DES.
 * Après la diversification des clefs, on inverse l'ordre des clefs
 * avant d'appliquer la fonction DES.

 * En effet (comme expliqué p.6 du pdf), l'algorithme de déchiffrment du DES
   est le même que celui de chiffrement en inversant l'ordre des sous-clefs.
 */
void des_decrypt(BYTE *ciphertext, BYTE *key, BYTE *plaintext, int rounds)
{
    BYTE(*KS)[48] = malloc(48 * rounds);
    setkey(key, rounds, KS);
    reverse_subkeys(KS, rounds);
    DES(ciphertext, plaintext, rounds, KS);
    free(KS);
}




/**************************************************************************/
/*****************************   EXERCICE 5   *****************************/
/**************************************************************************/





/*
 * brute_force implémente l'attaque brute force qui permet,
   à partir d'une clef de tour Kr, de retrouver la clef initiale.
 * On récupère d'abord la clef à trous kap à partir de Kr,
   à l'aide de la fonction reverse_key_schedule.
 * Puis, grâce à la fonction indiceBytesToFindInKey,
   on récupère les indices des bits sur lesquels effectuer une brute force.
 * Il y a 2^8 = 256 possibilités. Pour chacune d'entre elles,
   on complète les trous des bits de parité, et on teste si la clef
   est la bonne.
 */
void brute_force(BYTE kap[64], BYTE Kr[48], BYTE clair[8], BYTE chiffre[8],
                 int r)
{
    BYTE kconnu[64] = {0};
    int indices[8] = {0};
    reverse_key_schedule(
        r, kconnu); // On a retiré le Kr qui ne servait pas dans les paramètres
    indiceBytesToFindInKey(kconnu, indices);

    for (int i = 0; i < 64; i++)
    {
        if (kconnu[i] != 0)
        {
            kap[i] = Kr[kconnu[i] - 1];
        }
    }
    // il faut maintenant ajouter des 0 ou des 1 dans les positions "à trou" de
    // kap, et tester le DES pour voir si on obtient bien le même résultat.
    for (int i = 0; i < 256; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            kap[indices[j]] = (i >> j) & 1;
        }
        completerAvecBitsDeParite(kap);
        BYTE pkap[8] = {0};
        pack8(pkap, kap);
        BYTE ciphertext[8] = {0};
        des_encrypt(clair, pkap, ciphertext, r);
        int errorcounter = 0;
        for (int l = 0; l < 8; l++)
        {
            if (ciphertext[l] != chiffre[l])
            {
                errorcounter++;
                break;
            }
        }
        if (errorcounter == 0)
        {
            printf("VICTOIRE !!!! La clef est  : ");
            print_tab_hex(pkap, 8);
            return;
        }
    }
}



/**************************************************************************/
/*****************************   EXERCICE 6   *****************************/
/**************************************************************************/




/* Génération de nombre aléatoire
 */
int rand_num(int n)
{
    assert(n < RAND_MAX);
    return rand() % n;
}


int rand_num256()
{
    return rand_num(256);
}


/* Génération de nombre aléatoire sous forme de uint32_t
 */
uint32_t rand_uint32()
{
    uint32_t ret = 0;
    for (int i = 0; i < 4; i++)
    {
        ret ^= ((uint32_t)rand_num256()) << (i * 8);
    }
    return ret;
}


/* On crée 8 couples de clairs qui respectent la différentielle de l'énoncé pour r=3
 */
void gen_plaintexts_r_eq_3(uint64_t pairs[8][2])
{
    uint32_t l_1;
    uint32_t l_2;
    uint32_t r_1;
    for (int i = 0; i < 8; i++)
    {
        l_1 = rand_uint32();
        l_2 = rand_uint32();
        r_1 = rand_uint32();
        pairs[i][0] = ((uint64_t)l_1 << 32) ^ (uint64_t)r_1;
        pairs[i][1] = ((uint64_t)l_2 << 32) ^ (uint64_t)r_1;
    }
}


/* On crée 600 couples de clairs qui respectent la différentielle de l'énoncé pour r=6
 */
void gen_plaintexts_r_eq_6(uint64_t pairs[600][2])
{
    short i = 0;
    uint32_t l_1;
    uint32_t l_2;
    uint32_t r_1;
    uint32_t r_2;
    while (i <= 599)
    {
        if (i <= 299)
        {
            r_1 = rand_uint32();
            r_2 = r_1 ^ 0x04000000;
            l_1 = rand_uint32();
            l_2 = l_1 ^ 0x40080000;
        }
        else
        {
            r_1 = rand_uint32();
            r_2 = r_1 ^ 0x00000400;
            l_1 = rand_uint32();
            l_2 = l_1 ^ 0x00200008;
        }
        pairs[i][0] = ((uint64_t)l_1 << 32) ^ (uint64_t)r_1;
        pairs[i][1] = ((uint64_t)l_2 << 32) ^ (uint64_t)r_2;
        i++;
    }
}


/* Transforme un uint64_t en tableau de 8 BYTE.
 */
void cast_uint64_to_BYTE_tab(uint64_t in, BYTE out[8])
{
    for (int i = 0; i < 8; i++)
    {
        out[7 - i] = (BYTE)(in >> (i * 8));
    }
}



/* Pour que "l'oracle" qui calcule le DES tronqué à r-1 tours ait accès à la clef
 */
BYTE clef_a_retrouver[8] = {
    0x13, 0x34, 0x57, 0x79, 0x9B,
    0xBC, 0xDF, 0xF1}; 





/* On génère les couples de chiffrés correspondants aux couples de clairs stockés dans le tableau
 * pairsPlain (couples générés avec les fonctions gen_plaintexts_r_eq_3 et gen_plaintexts_r_eq_6).
 * out est un tableau de 8/600 paires de chiffrés i.e. 600 paires de tableau de 8 BYTE.
 */
void gen_cipher_pair(int rounds, uint64_t pairsPlain[600][2],
                     BYTE clairs[600][2][8], BYTE out[600][2][8], BYTE key[8])
{ 
    int tailleTab = 8; // NB de couples par défaut pour un 3 tours
    if (rounds == 6)
    { // NB de couples pour un 6 tours
        tailleTab = 600;
    }
    for (int i = 0; i < tailleTab; i++)
    {
        for (int j = 0; j < 2; j++)
        {
	  // Nous sommes dans pairsPlain[i][j] qui contient un uint64_t
            cast_uint64_to_BYTE_tab(pairsPlain[i][j], clairs[i][j]);
            des_encrypt(clairs[i][j], key, out[i][j], rounds);
        }
    }
}



/* Les deux fonctions suivantes servent à la
 * fabrication des tables de différence des tableaux de Sbox.
 * Utilisées pour la recherche des caractéristiques différentielles pour la cryptanalyse du schéma de feistel, ces fonctions sont inutilisées pour la cryptanalyse du dernier tour.
 */
/*
void set_difference_distribution(int *Sbox, int frequences[16][16])
{
    for (int dx = 0; dx < 16; ++dx)
    {
        for (int x = 0; x < 16; ++x)
        {
            int dy = Sbox[x] ^ Sbox[x ^ dx];
            ++frequences[dx][dy];
        }
    }
}


void create_difference_distributions_from_sbox(
    INT Sbox[8][64], int frequences_Sboxes[8][4][16][16])
{
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 4; ++j)
            set_difference_distribution((int *)Sbox[i] + (16 * j),
                                        frequences_Sboxes[i][j]);
}
*/



/* Fonction qui calcule juste Fk-1(c1) xor Fk-1(c2) et verifie si ca vaut
  delta_r, même ref que do_job_des.
 */
int chiffrement_partiel(BYTE clef[48], BYTE *delta_r, BYTE *out1, BYTE *out2)
{
    BYTE Fkmoins1_c1[8] = {0};
    BYTE Fkmoins1_c2[8] = {0};

    BYTE clef_tour[1][48] = {0};
    for (int i = 0; i < 48; ++i)
        clef_tour[0][i] = clef[i];

    // On revient à l'avant dernier tour du DES
    DES(out1, Fkmoins1_c1, 1, clef_tour); // Fk^(-1)(c1)
    DES(out2, Fkmoins1_c2, 1, clef_tour); // Fk^(-1)(c2)

    BYTE xor [8] = {0};
    for (int i = 0; i < 8; ++i)
    {
        xor[i] = Fkmoins1_c1[i] ^ Fkmoins1_c2[i];
    }
    for (int i = 0; i < 8; i++)
    {
        if (xor[i] != delta_r[i])
            return 0;
    }
    return 1;
}



/* Fonction inspirée de l'algorithme 10 de https://tel.archives-ouvertes.fr/tel-00577229/document,
 * page 71.
 * Pour toutes les combinaisons de sous clefs (il s'agit de retrouver la sous-clef K6), on effectue une * attaque sur l'avant dernier tour (ici le cinquième).
 * On incrémente J[i][j] si la clef j à l'entrée de la i-ème Sbox donne une différence en sortie 
 * delta_r.
 */
void do_job_des(INT r, BYTE *out11, BYTE *out12, BYTE *in11, BYTE *in12)
{ 
    if (r != 6)
        return;
    // on calcule delta_r
    BYTE delta_r[8] = {0};

    /* Une première tentative de calculer delta_r
    //Les tableaux out* font 8 byte.
    //Pour prendre les parties gauches et droites, on les sépare en deux
    tableaux de 4 BYTE chacun. BYTE c1_g[4]={0}; BYTE c1_d[4]={0}; BYTE
    c2_g[4]={0}; BYTE c2_d[4]={0};

    // /!\ Il faut calculer la différence en sortie pour trouver delta_r_g
    BYTE diff[4]={0};

    for(int i=0; i<4;++i){
      c1_g[i]=out11[i];
      c1_d[i]=out11[4+i];
      c2_g[i]=out12[i];
      c2_d[i]=out12[4+i];

      //delta_r_d = c1g ^c2g
      delta_r[4+i]=c1_g[i]^c2_g[i];
      //delta_r_g = c1d ^c2d ^ difference en sortie...
      delta_r[i]=c1_d[i] ^ c2_d[i] ^ diff[i];
    }
    */

    BYTE d1[8] = {0};
    BYTE d2[8] = {0};
    BYTE liste_des_clefs_de_tour[5][48] = {
        {0}}; // le 5 était r-1 avant, mais erreur de compilation à cause
              // d'un tableau de taille variable
    setkey(clef_a_retrouver, r - 1,
           liste_des_clefs_de_tour); // On donne à l'oracle les sous clefs
                                     // de tour pour faire le DES tronqué à
                                     // r-1 tours

    DES(in11, d1, r - 1, liste_des_clefs_de_tour);
    DES(in12, d2, r - 1, liste_des_clefs_de_tour);
    BYTE diff_apres_permutation[8] = {0};
    for (int i = 0; i < 8; i++)
    {
        diff_apres_permutation[i] = d1[i] ^ d2[i];
    }
    BYTE unpacked_diff_apres_permutation[64] = {0};
    BYTE unpacked_delta_r[64] = {0};
    unpack8(diff_apres_permutation, unpacked_diff_apres_permutation);
    // IP
    for (int j = 0; j < 64; j++)
    {
        unpacked_delta_r[j] = unpacked_diff_apres_permutation
            [IP[j] - 1]; // permutation initiale pour "annuler" la
                         // permutation finale faite par DES(in11, d1, r-1,
                         // liste_des_clefs_de_tour);
    }
    // i
    BYTE tmp[32] = {0};
    memcpy(tmp, unpacked_delta_r, 32 * sizeof(BYTE));
    memcpy(delta_r, unpacked_delta_r + 32 * sizeof(BYTE), 32 * sizeof(BYTE));
    memcpy(unpacked_delta_r + 32 * sizeof(BYTE), tmp, 32 * sizeof(BYTE));
    // IP^-1
    BYTE tmp2[64] = {0};
    for (int j = 0; j < 64; j++)
    {
        tmp2[j] = unpacked_delta_r[RFP[j] -
                                   1]; // permutation initiale pour "annuler" la
                                       // permutation finale faite par DES(in11,
                                       // d1, r-1, liste_des_clefs_de_tour);
    }
    pack8(delta_r, tmp2);
    // delta_r la différence à la sortie de l'avant dernier tour.

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 64; j++)
        {
            // calcul de la clef : on parcourt les valeurs possibles en
            // modifiant seulement le ième bloc.
            BYTE sous_clef_testee[48] = {0};
            for (int k = 0; k < 6; ++k)
            {
                sous_clef_testee[i * 6 + k] = j >> (5 - k) & 1;
            }
            if (chiffrement_partiel(sous_clef_testee, delta_r, out11, out12))
            {
                J[i][j] = J[i][j] + 1;
            }
        }
    }
}


/* Cette fonction remplit le tableau J en appelant N fois do_job_des, puis, à l'aide du tableau
 * JM, elle détermine la clef de tour Kr la plus probable. Enfin, elle utilise la fonction de brute 
 * force pour retrouver la clef de 64 bits à partir de la clef de tour Kr.
 */
void crypta_des(INT r, INT N)
{
    // r=6 et N=600 dans le sujet
    int a = N; //évite les erreurs de unused variable
    a++;
    // génération des 600 clairs avec void gen_plaintexts_r_eq_6(uint64_t
    // pairs[600][2])
    uint64_t pairsPlain[600][2] = {{0}};
    gen_plaintexts_r_eq_6(pairsPlain);
    // génération des chiffrés correspondants
    BYTE clairs[600][2][8];
    BYTE chiffres[600][2][8];
    gen_cipher_pair(r, pairsPlain, clairs, chiffres, clef_a_retrouver);

    // boucle des do_job
    // Apres on va cycler do_job_des sur tous les 600 bicouples, et selectionner
    // par groupes de 6 bits la meilleure clef, comme expliqué plus bas
    /*
      for all i
      quel est le j pour lequel J[i][j] sera maximal ?
      On construit la bonne clef par paquets de 6.
    */
    for (int i = 0; i < 600; i++)
    {
        do_job_des(r, chiffres[i][0], chiffres[i][1], clairs[i][0],
                   clairs[i][1]);
    }
    // Tableau J mis à jour
    int JM[8] = {0};
    for (int i = 0; i < 8; i++)
    {
        int max = 0;
        int indice_max = 0;
        // extraction des max du tableau J
        for (int j = 0; j < 64; j++)
        {
            if (J[i][j] > max)
            {
                max = J[i][j];
                indice_max = j;
            }
        }
        JM[i] = indice_max; // clef devinée
    }
    // Kr à fabriquer à partir de JM
    BYTE Kr[48] = {0};
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            Kr[6 * i + j] = (JM[i] >> (5 - j)) & 1;
        }
    }

    printf("Kr trouvé: \n");
    print_tab(Kr, 48);

    // Kap clef finale
    BYTE kap[64] = {0};
    brute_force(kap, Kr, clairs[0][0], chiffres[0][0], r);
}



int main()
{
    /******* pack et unpack *******/

    /*
      BYTE test_packed[8]={3,4,5,6,7,8,9,1};
      BYTE test_binary[64]={0};

      print_tab(test_packed,8);
      print_tab(test_binary,64);


      pack8(test_packed, test_binary);


      memset(test_packed,0,8*sizeof(BYTE));

      print_tab(test_packed,8);
      print_tab(test_binary,64);

      unpack8(test_packed, test_binary);

      print_tab(test_packed,8);
      print_tab(test_binary,64);
    */

    /************ DES *************/

    BYTE in[8] = {1, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF};
    BYTE key[8] = {0x13, 0x34, 0x57, 0x79, 0x9B, 0xBC, 0xDF, 0xF1};

    BYTE out[8] = {0};

    print_tab_hex(key, 8);

    /********* DES Encrypt ********/

    des_encrypt(in, key, out, 15);

    // printf("IN: ");
    // print_tab_binary(in, 8);
    // print_tab_hex(in, 8);

    // printf("OUT: ");
    // print_tab_binary(out, 8);
    // print_tab_hex(out, 8);

    /********* DES Decrypt ********/

    BYTE clear[8] = {0};
    des_decrypt(out, key, clear, 16);

    // printf("CLEAR: ");
    // print_tab_binary(clear, 8);
    // print_tab_hex(clear, 8);

    /********* Brute Force ********/

    BYTE kap[64] = {0};
    BYTE KR[48] = {0};
    BYTE KS[16][48] = {{0}};

    setkey(key, 16, KS);

    for (int i = 0; i < 48; ++i)
    {
        KR[i] = KS[14][i];
    }

    brute_force(kap, KR, in, out, 15);

    /********* crypta_des ********/

    crypta_des(6, 600);

    return 0;
}
