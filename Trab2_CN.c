#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int n = 7, m = 2; //m no criterio de kaiser (2) e n linhas do autovetor (7)
int l, c; // dimensões da matriz Y

// contando a quantidade de linhas (n) e colunas (m)
void dimensaoY(const char *filename){
    //abrindo arquivo
    FILE *file = fopen(filename, "r");
    
    //erro se o arquivo não abrir
    if (!file) {
        printf("Arquivo nao encontrado");
        exit(1); //Finaliza o programa com status de erro
    }

    double scan;
    //lendo o arquivo e contabilizando linhas
    while(!feof(file)){
        scan = fgetc(file);
        if(scan == '\n'){ l++; } // mais linha
    }

    fseek(file, 0, SEEK_SET); //voltando ao inicio do arquivo

    //lendo o arquivo e contabilizando colunas
    while (fscanf(file, "%lf", &scan)&&!feof(file)){
        c++; // mais coluna
    }
    
    c = c/l; //ajustando a quantidade de colunas

    fclose(file); //fechando arquivo
}

void lerMatriz(const char *filename, int n, int m, double M[m][m]) {
    //abrindo arquivo
    FILE *file = fopen(filename, "r");

    //erro se o arquivo não abrir
    if (!file) {
        printf("Arquivo nao encontrado");
        exit(1); //Finaliza o programa com status de erro
    }
    
    // Preenchendo a matriz enviada
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fscanf(file, "%lf", &M[i][j]);
        }
    }

    fclose(file); //fechando arquivo
}

//Função que multiplica matrizes
void multiplicacao_de_matrizes( int l1, int c1, int l2, int c2, double R[l1][c2], double M1[l1][c1], double M2[l2][c2]){
    for(int i = 0; i < l1; i++){
        for(int j = 0; j < c2; j++){
            double soma = 0;
            for (int k = 0; k < c1; k++){
                soma += M1[i][k] * M2[k][j];
            }
            R[i][j] = soma;
        }
    }
}


//Calcula o vetor R = At*Yi (usado no metodo mmq para todos os usuarios)
void vetorR(double R[m], double Yi[c], double TA[m][c]){
    for(int i = 0; i < m; i++){
        double soma = 0;
        for(int j = 0; j < c; j++){
            soma += TA[i][j] * Yi[j];
        }
        R[i]=soma;
    }
}

//Eliminacao de gauss sem pivotamento
void metodo_eliminacao_gauss(double Am[m][m], double *R, double *f) {
    double Raux[m], Amaux[m][m];
    double x;

    //criando clones de Y e R para usar no metodo
    for(int i = 0; i < m; i++){
        Raux[i] = R[i];
        for (int j = 0; j < m; j++){
            Amaux[i][j] = Am[i][j];
        }
    }
    
    //realizando as eliminacoes
    for (int i=0; i < m; i++) {
        //Montando a matriz triangular superior
        for (int j = i+1; j < m; j++) {
            x = Amaux[j][i] / Amaux[i][i];
            for (int n=i; n < m; n++) {
                Amaux[j][n] -= x * Amaux[i][n];
            }
            Raux[j] -= x * Raux[i];
        }
    }

    // Normalizando o vetor f
    for(int i = m-1; i >= 0; i--){
        for (int j = m-1; j > i; j--){
            Raux[i] -= Amaux[i][j] * f[j];
        }
        f[i] = Raux[i] / Amaux[i][i];
    }
}

//Metodo dos minimos quadrados
void mmq(double Y[l][c], double A[n][m], double F[l][m]){

    //Montando a transposta de A (TA)
    double TA[m][n];
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            TA[i][j] = A[j][i];
        }
    }
    
    //montando a matriz ....
    double Am[m][m];
    multiplicacao_de_matrizes(m, n, n, m, Am, TA, A);
    
    //Calculando f em Am * f = At * Yi para cada cliente
    printf("Matriz de fatores F:\n");
    for(int i = 0; i < l; i++){
        double Yi[c], R[m], f[m];
        //separando um vetor do cliente
        for(int j = 0; j < c; j++){
            Yi[j] = Y[i][j];
        }
        //Calculando R
        vetorR(R, Yi, TA);
        //metodo da eliminação de gauss
        metodo_eliminacao_gauss(Am, R, f);
        //calculando F
        for(int j = 0; j < m; j++){
            F[i][j] = f[j];
            printf("%lf ", F[i][j]);
        }
        printf("\n");
    }
}

void calcular_pontuacao(double F[l][m], double pontuacoes[l]) {
    for (int i = 0; i < l; i++) {
        pontuacoes[i] = F[i][0] - F[i][1];
    }
}

void encontrar_melhores(double pontuacoes[l], int melhores_indices[3], double melhores_pontuacoes[3]) {
    for (int i = 0; i < 3; i++) {
        melhores_indices[i] = -1;
        melhores_pontuacoes[i] = -__DBL_MAX__;
    }

    for (int i = 0; i < l; i++) {
        for (int j = 0; j < 3; j++) {
            if (pontuacoes[i] > melhores_pontuacoes[j]) {
                for (int k = 2; k > j; k--) {
                    melhores_pontuacoes[k] = melhores_pontuacoes[k - 1];
                    melhores_indices[k] = melhores_indices[k - 1];
                }
                melhores_pontuacoes[j] = pontuacoes[i];
                melhores_indices[j] = i;
                break;
            }
        }
    }
}

int main() {
    //achando as dimensões de Y
    dimensaoY("dados.txt");
    
    //declarando as matrizes
    double D[m][m], V[n][m], A[n][m], Y[l][c], F[l][m];

    //Recebendo as matrizes dos arquivos
    lerMatriz("matrizD.txt", m, m, D);
    lerMatriz("matrizV.txt", n, m, V);
    lerMatriz("dados.txt", l, c, Y);

    // fazendo a raiz de D
    for (int i = 0; i < m; i++) {
        D[i][i] = sqrt(D[i][i]);
    }

    //Montando a matriz A
    multiplicacao_de_matrizes(n, m, m, m, A, V, D);

    // calculando os vetores Yj
    for(int j = 0; j < c; j++){
        double media = 0;
        // media das colunas de Y
        for(int i = 0; i < l; i++){
            media += Y[i][j] / l;
        }
        
        //Subtraindo cada linha da coluna j pela media da coluna
        for(int i = 0; i < l; i++){
            Y[i][j] -= media;
        }
    }
    
    //método dos minimos quadrados
    mmq(Y, A, F);

    printf("\n ");
    
    //mostrando as caracteristicas melhor explicadas por cada fator em A
    for(int j = 0; j < m; j++){
        double caracteristica = 0;
        int indice;
        for(int i = 0; i < n; i++){
            //calcula qual tem a menor distância de 1
            if(fabs(1 - fabs(caracteristica)) > fabs(1 - fabs(A[i][j]))){
                caracteristica = A[i][j]; //recebendo nova característica
                indice = i; //armazenando posição
            }
        }
        printf("Fator %d é melhor explicado pela caracteristica %d (valor: %lf)\n", j+1, indice+1, caracteristica);
    }
    
    
    //Sistema de pontuação
    double pontuacoes[l];
    int melhores_indices[3];
    double melhores_pontuacoes[3];

    calcular_pontuacao(F, pontuacoes);
    encontrar_melhores(pontuacoes, melhores_indices, melhores_pontuacoes);

    printf("\nOs três melhores pontuadores são:\n");
    
    for (int i = 0; i < 3; i++) {
        printf("Cliente %d com pontuação %lf\n", melhores_indices[i] + 1, melhores_pontuacoes[i]);
    }  

}