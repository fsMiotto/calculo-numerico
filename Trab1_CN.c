#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 7
#define MAX_iter 100
#define epsilon 1E-5

//Lendo o arquivo
void ler_Matriz(double matriz[MAX][MAX], const char *nome_arquivo) {
    FILE *arquivo = fopen(nome_arquivo, "r");
    if (!arquivo) {
        perror("Arquivo nao encontrado");
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
    double erro, erro_final = 1;
    int iter_final;
    
    //Obtendo o menor erro da iteracao
    for (int i = 0; i < MAX; i++) {
        erro = fabs(autovalor[i] - lambda[i]) / fabs(autovalor[i]);
        if (erro < erro_final) {
            erro_final = erro;
            iter_final = i;
        }
    }
    
    //Comparando os erros
    if (erro_final >= epsilon) {
        for (int i = 0; i < MAX; i++){
            lambda[i] = autovalor[i];   
        }
        return -1;
    }
    return iter_final;
}

//Metodo das potencias
int metodo_potencias(double matriz[MAX][MAX], double *autovalor, double *autovetor, double *y) {
    double z[MAX], lambda[MAX], norma;
    int iter = 0, iter_final;
    
    //valores iniciais
    for (int i = 0; i < MAX; i++) {
        lambda[i] = 0.0;
        autovetor[i] = 1.0;
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

//Metodo das matrizes inversas
int metodo_potencias_inversas(double matriz[MAX][MAX], double q, double *autovalor, double *autovetor, double *y) {
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
        A[i][i]-=q;
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
    double R[MAX][MAX];
    double autovalor[MAX], autovetor[MAX], y[MAX];
    double q = 0;
    int iter_final;
    
    ler_Matriz(R, "correlacao.txt");

    printf("insira os valores do vetor inicial para o metodo das potencias (separados por espaço)\n");
    for (int i = 0; i < MAX; i++) {
        scanf("%lf", &y[i]);
    }

    iter_final = metodo_potencias(R, autovalor, autovetor, y);
    printf("\nPelo metodo das potencias, R possui: \nAutovalor dominante: %lf\n", autovalor[iter_final]);
    printf("\nAutovetor correspondente:\n");
    for (int i = 0; i < MAX; i++) {
        printf("%lf\n", autovetor[i]);
    }

    metodo_potencias_inversas(R, q, autovalor, autovetor, y);
    printf("\n\nPelo metodo das potencias inversas, R possui: \nMenor autovalor: %lf\n", (1/autovalor[iter_final])+q);
    printf("\nAutovetor correspondente:\n");
    for (int i = 0; i < MAX; i++) {
        printf("%lf\n", autovetor[i]);
    }
    
    q = 2;
    metodo_potencias_inversas(R, q, autovalor, autovetor, y);
    printf("\n\nPelo metodo das potencias inversas com deslocamento (2), R possui: \nOutro autovalor: %lf\n", (1/autovalor[iter_final])+q);
    printf("\nAutovetor correspondente:\n");
    for (int i = 0; i < MAX; i++) {
        printf("%lf\n", autovetor[i]);
    }
    
    q = 0.6;
    metodo_potencias_inversas(R, q, autovalor, autovetor, y);
    printf("\n\nPelo metodo das potencias inversas com deslocamento (0.6), R possui: \nOutro autovalor: %lf\n", (1/autovalor[iter_final])+q);
    printf("\nAutovetor correspondente:\n");
    for (int i = 0; i < MAX; i++) {
        printf("%lf\n", autovetor[i]);
    }
    
    q = 0.4;
    metodo_potencias_inversas(R, q, autovalor, autovetor, y);
    printf("\n\nPelo metodo das potencias inversas com deslocamento (0.4), R possui: \nOutro autovalor: %lf\n", (1/autovalor[iter_final])+q);
    printf("\nAutovetor correspondente:\n");
    for (int i = 0; i < MAX; i++) {
        printf("%lf\n", autovetor[i]);
    }
    
    q = 0.3;
    metodo_potencias_inversas(R, q, autovalor, autovetor, y);
    printf("\n\nPelo metodo das potencias inversas com deslocamento (0.3), R possui: \nOutro autovalor: %lf\n", (1/autovalor[iter_final])+q);
    printf("\nAutovetor correspondente:\n");
    for (int i = 0; i < MAX; i++) {
        printf("%lf\n", autovetor[i]);
    }
    
    q = 0.2;
    metodo_potencias_inversas(R, q, autovalor, autovetor, y);
    printf("\n\nPelo metodo das potencias inversas com deslocamento (0.2), R possui: \nOutro autovalor: %lf\n", (1/autovalor[iter_final])+q);
    printf("\nAutovetor correspondente:\n");
    for (int i = 0; i < MAX; i++) {
        printf("%lf\n", autovetor[i]);
    }
    
    return 0;
}