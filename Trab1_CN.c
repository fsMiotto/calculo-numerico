#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 7
#define MAX_iter 100
#define epsilon 1E-5

//Lendo o arquivo
void ler_matriz(double matriz[MAX][MAX], const char *nome_arquivo) {
    FILE *arquivo = fopen(nome_arquivo, "r");
    if (!arquivo) {
        printf("Arquivo nao encontrado");
        exit(0);
    }

    for (int i = 0; i < MAX; i++) {
        for (int j = 0; j < MAX; j++) {
            fscanf(arquivo, "%lf", &matriz[i][j]);
        }
    }
    fclose(arquivo);
}

//Teste de Parada
int teste_parada(double *lambda, double *autovalor) {
    double erro, menor_erro = 1;
    int i_final;
    
    //Obtendo o menor erro da iteracao
    for (int i = 0; i < MAX; i++) {
        erro = fabs(autovalor[i] - lambda[i]) / fabs(autovalor[i]);
        if (erro < menor_erro) {
            menor_erro = erro;
            i_final = i;
        }
    }
    
    //Comparando os erros
    if (menor_erro >= epsilon) {
        for (int i = 0; i < MAX; i++){
            lambda[i] = autovalor[i];   
        }
        return -1;
    }
    return i_final;
}

//Metodo das potencias
int metodo_potencias(double matriz[MAX][MAX], double *autovalor, double *autovetor, double *y) {
    double z[MAX], lambda[MAX], norma;
    int iter = 0, iter_final;
    
    //valores iniciais
    for (int i = 0; i < MAX; i++) {
        lambda[i] = 0;
        autovetor[i] = 1;
    }
    
    do {
        // multiplicacao matriz-vetor
        for (int i = 0; i < MAX; i++) {
            z[i] = 0.0;
            for (int j = 0; j < MAX; j++) {
                z[i] += matriz[i][j] * autovetor[j];
            }
        }
    
        //obtendo a norma
        norma = fabs(z[0]);
        for (int i = 0; i < MAX; i++) {
            if (fabs(z[i]) > norma) {
                norma = fabs(z[i]);
            }
        }


        for (int i = 0; i < MAX; i++) {
            //atualizanado autovetor
            autovetor[i] = z[i] / norma;
            //atualizar autovalor
            autovalor[i] = z[i] / y[i];
            //atualizando y(k)
            y[i] = autovetor[i];
        }
        //teste de parada
        iter_final = teste_parada(lambda, autovalor);
        iter++;
        } while ( iter < MAX_iter && iter_final == -1);
    return iter_final;
}

//Eliminacao de gauss sem pivotamento
void elim_gauss (double R[MAX][MAX], double z[MAX], double y[MAX]) {
    int max_col;
    double fator, max;
    double y2[MAX],R2[MAX][MAX];
    
    for(int i = 0; i<MAX; i++){
        y2[i]=y[i];
        for (int j = 0; j<MAX; j++){
            R2[i][j]=R[i][j];
        }
    }
    
    //eliminacoes
    for (int i = 0; i < MAX; i++) {
        //matriz triangular superior
        for (int j = i + 1; j < MAX; j++) {
            fator = R2[j][i] / R2[i][i];
            for (int k = i; k < MAX; k++) {
                R2[j][k] -= fator * R2[i][k];
            }
            y2[j] -= fator * y2[i];
        }
    }

    // substituicao
    for(int i = MAX-1; i>=0; i--){
        for (int j = MAX-1; j>i; j--){
            y2[i] -= R2[i][j]*z[j];
        }
        z[i] = y2[i] / R2[i][i];
    }
}

//Metodo das potencias inverso
int metodo_potencias_inverso(double matriz[MAX][MAX], double deslocamento, double *autovalor, double *autovetor, double *y) {
    double A[MAX][MAX], temp [MAX];
    double z[MAX], lambda[MAX], norma;
    int iter = 0, iter_final;
    
    //valores iniciais
    for (int i = 0; i < MAX; i++) {
        y[i] = 1.0;
        lambda[i] = 0.0;
        for(int j=0;j<MAX;j++){
            A[i][j] = matriz[i][j];
        }
        A[i][i]-=deslocamento;
    }
    
    do {
        //atualizar vetor z
        elim_gauss(A, z, y);
        
        //obtendo a norma
        norma = fabs(z[0]);
        for (int i = 0; i < MAX; i++) {
            if (fabs(z[i]) > norma) {
                norma = fabs(z[i]);
            }
        }

        for (int i = 0; i < MAX; i++) {
            //atualizar autovetor
            autovetor[i] = z[i] / norma;
            //atualizar autovalor
            autovalor[i] = z[i] / y[i];
            //atualizando y(k)
            y[i] = autovetor[i];
        }

        //teste de padara
        iter_final = teste_parada(lambda, autovalor);
        iter++;
        } while ( iter < MAX_iter && iter_final == -1);
    return iter_final;
}



int main(){
    double R[MAX][MAX], D[MAX][MAX], P[MAX][MAX];
    double autovalor[MAX], autovetor[MAX], y[MAX], deslocamento[MAX];
    int iter_final, m;
    
    ler_matriz(R, "correlacao.txt");

    //zerando a matriz D
    for(int i=0; i < MAX; i++){
        for(int j=0; j < MAX; j++){
            printf("%lf  ", D[i][j])
        }
    }

    // obtendo o vetor y inicial
    printf("insira os valores do vetor inicial para o metodo das potencias (separados por espaco)\n");
    for (int i = 0; i < MAX; i++) {
        scanf("%lf", &y[i]);
    }


    // aplicando o metodo das potencias
    iter_final = metodo_potencias(R, autovalor, autovetor, y);
    printf("\nPelo metodo das potencias, R possui: \nAutovalor dominante: %lf\n", autovalor[iter_final]);
    printf("Autovetor correspondente:\n");
    for (int i = 0; i < MAX; i++){
        printf("%lf\n", autovetor[i]);
    }
    
    // MATRIZES D e P
    for (int j=0; j < MAX; j++) {
        // obtendo a coluna de autovetores de P
        P[j][0] = autovetor[j];
        // Zeros na matriz D
        D[0][j] = 0;
    }
    // autovalor dominante na primeira posição da matriz d autovalores
    D[0]D[0] = autovalor[iter_final];


    
    // obtendo os valores de deslocamento
    printf("insira os (6) deslocamentos para o metodo das potencias em ordem decrescente (sepados por espaco)\n");
    for (int i=0; i < MAX-1; i++){
        scanf("%lf", &deslocamento[i]);
        while(deslocamento[i] < 0 || deslocamento[i] > autovalor[iter_final]){
            printf("valor invalido, digite um valor entre 0 e %lf:  ", autovalor[iter_final]);
            scanf("%lf", &deslocamento[i]);
        }
    }

    //aplicando o metodo das potencias inverso com deslocamento
    for(int i=0; i < MAX-1; i++){
        metodo_potencias_inverso(R, deslocamento[i], autovalor, autovetor, y);
        double autovalor_nverso = (1/autovalor[iter_final])+deslocamento[i];
        printf("\n\nPelo metodo das potencias inverso com deslocamento (%.3lf), R possui: \n Autovalor: %lf\n", deslocamento[i], autovalor_inverso);
        printf("Autovetor correspondente:\n");

        //Matriz P e D
        for (int j=0; j < MAX; j++) {
            printf("%lf\n", autovetor[j]);
            // obtendo a coluna de autovetores de P
            P[j][i+1] = autovetor[j];
            // Zeros na matriz D
            D[i+1][j] = 0; 
        }
        // Obtendo o autovalor da diagonal de D
        D[i+1][i+1] = autovalor_nverso;

        //achando o "m" para o criterio de kaiser
        if(autovalor_nverso > 1 ){m++;}
    }
    
    //Prints da matris D
    printf("Matriz D diagonal de autovalores de R:\n");
    for(int i=0; i < MAX; i++){
        for(int j=0; j < MAX; j++){
            printf("%lf  ", D[i][j])
        }
        printf("\n");
    }

    //Prints da matris P
    printf("\n\nMatriz P de autovetores dos autovalores de D:\n");
    for(int i=0; i < MAX; i++){
        for(int j=0; j < MAX; j++){
            printf("%lf  ", P[i][j])
        }
        printf("\n");
    }

    printf("\n\nPelo criterio de Kiser o valor de 'm' em R é %d", m);

    return 0;
}