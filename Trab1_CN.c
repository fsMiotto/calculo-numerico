#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 7
#define MAX_iter 100
#define epsilon 1e-5

//Lendo o arquivo
void ler_matriz(double matriz[MAX][MAX], const char *nome_arquivo) {
    //Abrindo o arquivo
    FILE *arquivo = fopen(nome_arquivo, "r");

    //Erro
    if (!arquivo) {
        printf("Arquivo nao encontrado");
        exit(0);
    }

    //Lendo o arquivo
    for (int i=0; i < MAX; i++) {
        for (int j=0; j < MAX; j++) {
            fscanf(arquivo, "%lf", &matriz[i][j]);
        }
    }

    //Fechando o arquivo
    fclose(arquivo);
}

//Teste de Parada
int teste_parada(double *autovalor_0, double *autovalor_1) {
    double erro, menor_erro = 100;
    int i_dominante;
    
    //Obtendo o menor erro da iteracao
    for (int i=0; i < MAX; i++) {
        //Formula do erro
        erro = fabs(autovalor_1[i] - autovalor_0[i]) / fabs(autovalor_1[i]);
        //obtendo o menor erro para comparar com a parada
        if (erro < menor_erro) {
            menor_erro = erro;
            i_dominante = i;
        }
    }
    
    //Comparando os erros
    if (menor_erro >= epsilon) {return -1;} //continua o metodo
    else {return i_dominante;} //para e retorna a posição do autovalor com menor erro
}

//Eliminacao de gauss sem pivotamento
void metodo_eliminacao_gauss (double R[MAX][MAX], double *Z, double *Y) {
    double Yaux[MAX],Raux[MAX][MAX];
    double fator, max;

    //criando clones de Y e R para usar no metodo
    for(int i = 0; i < MAX; i++){
        Yaux[i] = Y[i];
        for (int j = 0; j < MAX; j++){
            Raux[i][j] = R[i][j];
        }
    }
    
    //realizando as eliminacoes
    for (int i=0; i < MAX; i++) {
        //matriz triangular superior
        for (int j = i+1; j < MAX; j++) {
            fator = Raux[j][i] / Raux[i][i];
            for (int k=i; k < MAX; k++) {
                Raux[j][k] -= fator * Raux[i][k];
            }
            Yaux[j] -= fator * Yaux[i];
        }
    }

    // fazendo as substituicao
    for(int i = MAX-1; i >= 0; i--){
        for (int j = MAX-1; j > i; j--){
            Yaux[i] -= Raux[i][j]*Z[j];
        }
        Z[i] = Yaux[i] / Raux[i][i];
    }
}

//Metodo das potencias
double metodo_potencias(double R[MAX][MAX], double *autovalor_1, double *autovetor, double *Y) {
    double Z[MAX], autovalor_0[MAX], Yaux[MAX];
    double norma;
    int iter = 0, i_dominante;
    
    //valores iniciais
    for (int i=0; i < MAX; i++) {
        //auxiliar de Y para não alterar o original
        Yaux[i] = Y[i];
        autovetor[i] = 1;
    }
    
    do {
        // multiplicando matriz com autovetor
        for (int i=0; i < MAX; i++) {
            //zerando o vetor Z
            Z[i] = 0;
            for (int j=0; j < MAX; j++) {
                //atualizando o vetor Z
                Z[i] += R[i][j] * autovetor[j];
            }
        }
    
        //obtendo a norma
        norma = fabs(Z[0]);
        for (int i=0; i < MAX; i++) {
            if (fabs(Z[i]) > norma) {
                norma = fabs(Z[i]);
            }
        }

        for (int i=0; i < MAX; i++) {
            //atualizando o autovalor(k-1)
            autovalor_0[i] = autovalor_1[i];
            //atualizanado autovetor
            autovetor[i] = Z[i] / norma;
            //atualizar autovalor(k)
            autovalor_1[i] = Z[i] / Yaux[i];
            //atualizando Y(k)
            Yaux[i] = autovetor[i];
            
        }

        //teste de parada
        i_dominante = teste_parada(autovalor_0, autovalor_1);
        iter++;
        
        } while ( iter < MAX_iter && i_dominante == -1);

    return autovalor_1[i_dominante]; //Retorna o autovalor dominante
}

//Metodo das potencias inverso
double metodo_potencias_inverso(double R[MAX][MAX], double *autovalor_1, double *autovetor, double *Y, double deslocamento) {
    double A[MAX][MAX], Z[MAX], autovalor_0[MAX], Yaux[MAX];
    double norma;
    int iter = 0, i_dominante;
    
    //valores iniciais
    for (int i=0; i < MAX; i++) {
        //auxiliar de Y para não alterar o original
        Yaux[i] = Y[i];
        //Montando uma matriz R para fazer o deslocamento
        for(int j=0; j < MAX; j++){
            A[i][j] = R[i][j];
        }
        //Fazendo o deslocamento
        A[i][i] -= deslocamento;
    }
    
    do {
        //atualizar vetor Z
        metodo_eliminacao_gauss(A, Z, Yaux);
        
        //obtendo a norma
        norma = fabs(Z[0]);
        for (int i=0; i < MAX; i++) {
            if (fabs(Z[i]) > norma) {
                norma = fabs(Z[i]);
            }
        }

        for (int i=0; i < MAX; i++) {
            //atualizando o autovalor(k-1)
            autovalor_0[i] = autovalor_1[i];
            //atualizar autovetor
            autovetor[i] = Z[i] / norma;
            //atualizar autovalor(k)
            autovalor_1[i] = Z[i] / Yaux[i];
            //atualizando Y(k)
            Yaux[i] = autovetor[i];
            }

        //teste de padara
        i_dominante = teste_parada(autovalor_0, autovalor_1);
        iter++;

        } while ( iter < MAX_iter && i_dominante == -1);

    //Retorna o autovalor dominante  
    return (1/autovalor_1[i_dominante]) + deslocamento; 
}



int main(){
    //Declaração de matrizes
    double R[MAX][MAX], D[MAX][MAX], P[MAX][MAX];
    //Declaração de vetores
    double autovalor[MAX], autovetor[MAX], Y[MAX], deslocamento[MAX];
    //declaração de variaveis
    double autovalor_dominante;
    int m = 0;
    
    //Obtendo a matriz R
    ler_matriz(R, "correlacao.txt");

    //zerando a matriz D
    for(int i=0; i < MAX; i++){
        for(int j=0; j < MAX; j++){
            D[i][j] = 0;
        }
    }

    // obtendo o vetor Y inicial
    printf("insira os valores do vetor inicial para o metodo das potencias (separados por espaco)\n");
    for (int i = 0; i < MAX; i++) {
        scanf("%lf", &Y[i]);
    }


    // aplicando o metodo das potencias
    autovalor_dominante = metodo_potencias(R, autovalor, autovetor, Y);
    printf("\nPelo metodo das potencias R possui: \nAutovalor dominante:  %lf \nAutovetor correspondente:\n", autovalor_dominante);
    //printando o autovetor
    for (int i = 0; i < MAX; i++){
        printf("%lf\n", autovetor[i]);
    }
    
    // Para o metodo de kaiser
    if(autovalor_dominante > 1 ){m++;}
    
    // MATRIZES D e P
    for (int j=0; j < MAX; j++) {
        // obtendo a coluna de autovetores de P
        P[j][0] = autovetor[j];
        // Zeros na matriz D
        D[0][j] = 0;
    }
    // autovalor dominante na primeira posição da matriz d autovalores
    D[0][0] = autovalor_dominante;



    // obtendo os valores de deslocamento
    printf("insira os (6) deslocamentos para o metodo das potencias em ordem decrescente (sepados por espaco)\n");
    for (int i=0; i < MAX-1; i++){
        scanf("%lf", &deslocamento[i]);
        //verificando se é valido
        while(deslocamento[i] < 0 || deslocamento[i] > autovalor_dominante){
            printf("valor invalido, digite um valor entre 0 e %lf:  ", autovalor_dominante);
            scanf("%lf", &deslocamento[i]);
        }
    }

    //aplicando o metodo das potencias inverso com deslocamento
    for(int i=0; i < MAX-1; i++){
        autovalor_dominante = metodo_potencias_inverso(R, autovalor, autovetor, Y, deslocamento[i]);
        printf("\n\nPelo metodo das potencias inverso com deslocamento (%.3lf), R possui: \nAutovalor: %lf \nAutovetor correspondente:\n", deslocamento[i], autovalor_dominante);

        //Matriz P e D
        for (int j=0; j < MAX; j++) {
            printf("%lf\n", autovetor[j]);
            // obtendo a coluna de autovetores de P
            P[j][i+1] = autovetor[j];
            // Zeros na matriz D
            D[i+1][j] = 0; 
        }
        // Obtendo o autovalor da diagonal de D
        D[i+1][i+1] = autovalor_dominante;

        //achando o "m" para o criterio de kaiser
        if(autovalor_dominante > 1 ){m++;}
    }
    
    //Prints da matris D
    printf("Matriz D diagonal de autovalores de R:\n");
    for(int i=0; i < MAX; i++){
        for(int j=0; j < MAX; j++){
            printf("%lf  ", D[i][j]);
        }
        printf("\n");
    }

    //Prints da matris P
    printf("\n\nMatriz P de autovetores dos autovalores de D:\n");
    for(int i=0; i < MAX; i++){
        for(int j=0; j < MAX; j++){
            printf("%lf  ", P[i][j]);
        }
        printf("\n");
    }

    //M do metodo de Kaiser
    printf("\n\nPelo criterio de Kiser o valor de 'm' em R eh %d", m);

    return 0;
}