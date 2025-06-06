#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include <windows.h>

#define Nx_real 512 ///Número de divisões do domínio real no eixo X // 2.13769e-001
#define Ny_real 512 ///Número de divisões do domínio real no eixo y

#define PI                                                                     \
  3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825

///Limites do domínio [ax,bx] x [ay,by]
#define ax 0//7.00420e-006
#define bx 1.0
#define ay 0
#define by 1.0

#define Nx (Nx_real + 1)
#define Ny (Ny_real + 1)

#define iter_down 1

#define iter_up 1

#define iter_solve 10

#define MAX_MG 100

#define TOLERANCE 1e-9
//-------------------------------------------------------------

// Protótipos

double contorno_top(double x, double y);
double contorno_bottom(double x, double y);
double contorno_left(double x, double y);
double contorno_right(double x, double y);
double funcao_fonte(double x, double y);
double sol_analitica(double x, double y);
void calculaRHS(double **rhs, double dx, double dy, int linhas, int colunas);
void condicoes_contorno(double **u, double dx, double dy, int linhas, int colunas);
double eq_diferencas_finitas(double **u, double **rhs, int i, int j, double dx, double dy);
void calcula_sol_analitica(double **u_Sol_Anal, double dx, double dy, int linhas, int colunas);
void gauss_seidel(double **u_new, double **u, double **rhs, int n_iter, double dx, double dy, int linhas, int colunas);

void inicializa_matriz(double **matriz, int linhas, int colunas);
void imprime_matriz(double **u, int linhas, int colunas);
double norma_erro(double **u_Sol_Anal, double **u, double dx, double dy, int linhas, int colunas);
double erro_absoluto(double **u_Sol_Anal, double **u, int linhas, int colunas);
void salvar_dados(double **u, double **u_Sol_Anal,double dx, double dy, int linhas, int colunas);

double **Alocar_matriz_real (int linhas, int colunas);
double **Liberar_matriz_real (int linhas, int colunas, double **v);

void copia_matriz(double **matriz1, double **matriz2, int linhas, int colunas);
void soma_matriz(double **matriz1, double **matriz2, double **matriz3, int linhas, int colunas);
void subtrai_matriz(double **matriz1, double **matriz2, double **matriz3, int linhas, int colunas);
void salva_residuo(double *v1);
void salvar_metricas(double valor1, double valor2, int valor3);

//-------------------------------------------------------------
//              Elementos de Multigrid
//-------------------------------------------------------------

double stencil(double **residual, double **u_new, double dx, double dy, int j, int i);
void residuo(double **residual, double **rhs, double **u_new, double dx, double dy, int linhas, int colunas);
void restringe(double **residual, double **residual_2h, int linhas_restrict, int colunas_restrict);
void prolongamento(double **u_2h, double **u, int linhas_2h, int colunas_2h, int linhas, int colunas);

void multigrid(double ***u_new, double ***u, double ***rhs, int linhas,
               int colunas, double dx, double dy, int min_linhas,
               int min_colunas, int max_linhas, int max_colunas, int i, int n_grids, int mi, int *contador, int *vetor_fractal);

double total_grids(int t_linha, int t_coluna);
void fractal_sequence(int n, int** seq, int* size);
double max_value(double **matriz, int linhas, int colunas);
//-------------------------------------------------------------
int main()
{
// Declaração variáveis

    int linhas = Ny, colunas = Nx;

    double dx = (bx- ax) / (colunas - 1);
    double dy = (by - ay) / (linhas - 1);

    int min_linhas = 3, min_colunas = 3;

    int i = 0;
//-------------------------------------------------------------
// Cálcula número total de grids
    int t_linha = Ny, t_coluna = Nx, n_grids = 0;
    n_grids = total_grids(t_linha, t_coluna);
//-------------------------------------------------------------
// Aloca memória
    double ***rhs = (double***)malloc((n_grids) * sizeof(double **));
    double ***u = (double***)malloc((n_grids) * sizeof(double **));
    double ***u_new = (double***)malloc((n_grids) * sizeof(double **));

    for(int j = 0; j < n_grids; j++) {
        u[j] = Alocar_matriz_real(linhas, colunas);
        u_new[j] = Alocar_matriz_real(linhas, colunas);
        rhs[j] = Alocar_matriz_real(linhas, colunas);

        inicializa_matriz(rhs[j], linhas, colunas);
        inicializa_matriz(u_new[j], linhas, colunas);
        inicializa_matriz(u[j], linhas, colunas);

        linhas = ( (linhas - 1)/2 ) + 1;
        colunas = ( (colunas - 1)/2 ) + 1;
    }
    linhas = Ny;
    colunas = Nx;
    calculaRHS(rhs[0], dx, dy, linhas, colunas);
    condicoes_contorno(u_new[0], dx, dy, linhas, colunas);
//-------------------------------------------------------------
// Fractal
    int n = n_grids -1;  // Altere este valor para gerar uma sequência maior
    int *sequence;
    int size;
    fractal_sequence(n, &sequence, &size);
    /*for (int i = 0; i < size; i++) {
        printf("%d ", sequence[i]);
    }
    printf("\n");
    printf("%d", size);*/

//-------------------------------------------------------------
// Iterações de MG

    double *vetor_residuo;
    vetor_residuo = (double *) malloc( (MAX_MG + 1) * sizeof(double) );

    int contador = 0;
    double **stop_residual;
    double residuo_mg, residuo_mg_anterior = 0;
    stop_residual = Alocar_matriz_real(linhas, colunas);

    residuo(stop_residual, rhs[0], u_new[0], dx, dy, linhas, colunas);
    residuo_mg = fabs(max_value(stop_residual, linhas, colunas));
    vetor_residuo[0] = residuo_mg;

    int k = 0;
    LARGE_INTEGER frequency, start, end;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&start);
    while(k < MAX_MG){
        multigrid(u_new, u, rhs, linhas, colunas, dx, dy, min_linhas, min_colunas, Ny, Nx, i, n_grids, 1, &contador, sequence);
        contador = 0;
        inicializa_matriz(stop_residual, linhas, colunas);
        residuo(stop_residual, rhs[0], u_new[0], dx, dy, linhas, colunas);
        residuo_mg = fabs(max_value(stop_residual, linhas, colunas));
        if(residuo_mg_anterior == residuo_mg){
            break;
        }
        else {
            residuo_mg_anterior = residuo_mg;
        }
        printf("Residuo_mg = %.5e, k = %d\n", residuo_mg, k);
        vetor_residuo[k + 1] = residuo_mg;
        if(residuo_mg < TOLERANCE){
            break;
        }
        k++;
    }

    k++;
    vetor_residuo[k + 1] = -1;
    salva_residuo(vetor_residuo);

    QueryPerformanceCounter(&end);
    double elapsedTime = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;
    printf("Tempo de execuçao: %f segundos\n", elapsedTime);
    printf("Residuo MG: %.5e\n", residuo_mg);
    printf("Iteracoes MG: %d\n", k);

    stop_residual = Liberar_matriz_real(linhas, colunas, stop_residual);

    /*int contador = 0;
    for(int k = 0; k < MAX_MG; k++){
        multigrid(u_new, u, rhs, linhas, colunas, dx, dy, min_linhas, min_colunas, Ny, Nx, i, n_grids, 1, &contador, sequence);
        contador = 0;
    }*/
    free(vetor_residuo);
//-------------------------------------------------------------
// Parâmetros de avaliação
    double **u_Sol_anal;

    u_Sol_anal = Alocar_matriz_real(linhas, colunas);
    inicializa_matriz(u_Sol_anal, linhas, colunas);

    calcula_sol_analitica(u_Sol_anal, dx, dy, linhas, colunas);

    double norma2_erro, norma_inf_erro;
    norma2_erro = norma_erro(u_Sol_anal, u[0], dx, dy, linhas, colunas);
    norma_inf_erro = erro_absoluto(u_Sol_anal, u[0], linhas, colunas);

    salvar_metricas(norma2_erro, norma_inf_erro, k);

    salvar_dados(u[0], u_Sol_anal, dx, dy, linhas, colunas);
//-------------------------------------------------------------
// Libera memória
    u_Sol_anal = Liberar_matriz_real(linhas, colunas, u_Sol_anal);
    for(int j = 0; j < n_grids; j ++) {
        u[j] = Liberar_matriz_real(linhas, colunas, u[j]);
        u_new[j] = Liberar_matriz_real(linhas, colunas, u_new[j]);
        rhs[j] = Liberar_matriz_real(linhas, colunas, rhs[j]);
        linhas = ( (linhas - 1)/2 ) + 1;
        colunas = ( (colunas - 1)/2 ) + 1;
    }

    free(u);
    free(u_new);
    free(rhs);
    free(sequence);
//-------------------------------------------------------------

    return 0;
}

//-------------------------------------------------------------
// Condição contorno superior
double contorno_top(double x, double y) {
    return (sin(2 * PI * x) * sin(2 * PI * y));
}
//-------------------------------------------------------------
// Condição contorno inferior
double contorno_bottom(double x, double y) {
    return (sin(2 * PI * x) * sin(2 * PI * y));
}
//-------------------------------------------------------------
// Condição contorno esquerda
double contorno_left(double x, double y) {
    return (sin(2 * PI * x) * sin(2 * PI * y));
}
//-------------------------------------------------------------
// Condição contorno direita
double contorno_right(double x, double y) {
    return (sin(2 * PI * x) * sin(2 * PI * y));
}
//-------------------------------------------------------------
// Função fonte
double funcao_fonte(double x, double y) {
    double val = 0.0;
    val = (8 * PI * PI * sin(2 * PI * x) * sin(2 * PI * y));
    return val;
}
//-------------------------------------------------------------
// Solução analítica
double sol_analitica(double x, double y) {

    double val = 0.0;
    val = (sin(2 * PI * x) * sin(2 * PI * y));
    return val;
    ///(cos(2.0 * PI * x) * sin(2.0 * PI * y)) / (2.0 * PI)

    // Para -Au=f -> 0.5*exp(x)*exp(y)
}
//-------------------------------------------------------------
// Calcula RHS
void calculaRHS(double **rhs, double dx, double dy, int linhas, int colunas){
    int i, j;
    double x, y;
    for (j = 1; j <= (linhas - 2); j++) {
        for (i = 1; i <= (colunas - 2); i++) {
            x = ax + ((double)(i)) * dx;
            y = ay + ((double)(j)) * dy;
            rhs[j][i] = funcao_fonte(x, y);
        }
    }
}
//-------------------------------------------------------------
// Função condições de fronteira
void condicoes_contorno(double **u, double dx, double dy, int linhas, int colunas) {
    // Dxu: -sin(2.0 * PI * x) * sin(2.0 * PI * y)
    // Dyu: cos(2.0 * PI * x) * cos(2.0 * PI * y)
    int j,i;
    double y, x;

    for(j = 0; j < (linhas); j++){
        y = ay + ((double)(j)) * dy;
        // Left
        u[j][0] = contorno_left(ax, y);
        // Right
        u[j][colunas - 1] = contorno_right(bx, y);
    }
    for(i = 0; i < (colunas); i++){
        x = ax + ((double)(i)) * dx;
        // Bottom
        u[0][i] = contorno_bottom(x, ay);
        // Top
        u[linhas - 1][i] = contorno_top(x, by);
    }

}
//-------------------------------------------------------------
// Discretização
double eq_diferencas_finitas(double **u, double **rhs, int i, int j, double dx, double dy) {
    /// U[Linha][Coluna]  U[j][i]
    /// j - > Linha
    /// i - > Coluna
    double beta = dx / dy;
    double beta2 = beta * beta;

    return (1.0 / (2.0 * (1.0 + beta2))) *
         (u[j][i+1] + u[j][i-1] + beta2 * (u[j+1][i] + u[j-1][i]) +
          dx * dx * rhs[j][i]);

    /*
        Para o caso -Au = f

        (1.0 / (2.0 * (1.0 + beta2))) *
         (u[j][i+1] + u[j][i-1] + beta2 * (u[j+1][i] + u[j-1][i]) -
          dx * dx * rhs[j][i]);
    */
}
//-------------------------------------------------------------
// Cálcula solução analítica
void calcula_sol_analitica(double **u_Sol_Anal, double dx, double dy, int linhas, int colunas){
    double y, x;
    int i,j;
    for (j = 1; j <= (linhas - 2); j++) {
        for (i = 1; i <= (colunas - 2); i++) {
            x = ax + ((double)(i)) * dx;
            y = ay + ((double)(j)) * dy;
            u_Sol_Anal[j][i] = sol_analitica(x, y);
      }
    }
}
//-------------------------------------------------------------
// Gauss-Seidel

void gauss_seidel(double **u_new, double **u, double **rhs, int n_iter, double dx, double dy, int linhas, int colunas){

    double max_diff, diff;
    int j, i, iter;

    for (iter = 0; iter < n_iter; iter++) {
        max_diff = 0.0;

        for (j = 1; j <= (linhas - 2); j++) {
            for (i = 1; i <= (colunas - 2); i++) {
                u_new[j][i] = eq_diferencas_finitas(u_new,rhs,i,j,dx,dy);
                diff = fabs(u_new[j][i] - u[j][i]);
                if (diff > max_diff) {
                    max_diff = diff;
                }
            }
        }
        for (j = 1; j <= (linhas - 2); j++) {
            for (i = 1; i <= (colunas - 2); i++) {
                u[j][i] = u_new[j][i];
            }
        }
    }
}

//-------------------------------------------------------------
// Inicializa matriz com zeros
void inicializa_matriz(double **matriz, int linhas, int colunas) {
    for(int j=0; j<linhas; j++){
        for(int i=0; i<colunas; i++){
            matriz[j][i] = 0;
        }
    }
}
//-------------------------------------------------------------
// Imprime a matriz enviada
void imprime_matriz(double **u, int linhas, int colunas){
    for(int j=0; j<linhas; j++){
        for(int i=0; i<colunas; i++){
            printf("%f ", u[j][i]);
        }
        printf("\n");
    }
}
//-------------------------------------------------------------
// Cálculo da norma do erro
double norma_erro(double **u_Sol_Anal, double **u, double dx, double dy, int linhas, int colunas) {
    double error_norm = 0.0;
    for (int j = 1; j <= (linhas - 2); j++) {
        for (int i = 1; i <= (colunas - 2); i++) {
            double error = u_Sol_Anal[j][i] - u[j][i];
            error = error * error * dx * dy;
            error_norm += error;
        }
    }
    error_norm = sqrt(error_norm / ((bx - ax) * (by - ay)));
    printf("Norma do erro: %.5e \n", error_norm);
    return error_norm;
}
//-------------------------------------------------------------
// Erro Absoluto
double erro_absoluto(double **u_Sol_Anal, double **u, int linhas, int colunas){
    double erro_value=0.0, erro_anterior=0.0;
    for (int j = 1; j <= (linhas - 2); j++) {
        for (int i = 1; i <= (colunas - 2); i++) {
            erro_value = fabs(u_Sol_Anal[j][i] - u[j][i]);
            if (erro_value > erro_anterior){
                erro_anterior = erro_value;
            }
        }
    }
    printf("Erro absoluto: %.5e \n", erro_value);
    return erro_value;
}
//-------------------------------------------------------------
// Salva dados
void salvar_dados(double **u, double **u_Sol_Anal,double dx, double dy, int linhas, int colunas) {
    FILE *file_u = fopen("dados_u.txt", "w");
    FILE *file_u_Sol_Anal = fopen("dados_u_Sol_Anal.txt", "w");
    FILE *file_erro = fopen("erro.txt", "w");

    int i, j;
    double x, y, erro_value;

    for (j = 1; j <= (linhas - 2); j++) {
        for (i = 1; i <= (colunas - 2); i++) {
            x = ax + ((double)(i)) * dx;
            y = ay + ((double)(j)) * dy;
            erro_value = fabs(u_Sol_Anal[j][i] - u[j][i]);

            fprintf(file_u, "%f %f %f\n", x, y, u[j][i]);
            fprintf(file_u_Sol_Anal, "%f %f %f\n", x, y, u_Sol_Anal[j][i]);
            fprintf(file_erro, "%f %f %f\n", x, y, erro_value);
        }
    }

    fclose(file_u);
    fclose(file_u_Sol_Anal);
    fclose(file_erro);
}
//-------------------------------------------------------------
// Aloca matriz dinamicamente
double **Alocar_matriz_real (int linhas, int colunas){
    double **v;
    int   i;
    if (linhas < 1 || colunas < 1) {
        printf ("** Erro: Parametro invalido **\n");
        return (NULL);
    }
    v = (double **) malloc (linhas * sizeof(double *));
    if (v == NULL) {
        printf ("** Erro: Memoria Insuficiente **");
        return (NULL);
    }
    for (i = 0; i < linhas; i++ ) {
        v[i] = (double*) malloc (colunas * sizeof(double));
        if (v[i] == NULL) {
        printf ("** Erro: Memoria Insuficiente **");
        return (NULL);
        }
    }
    return (v);
}
//-------------------------------------------------------------
// Libera espaço matriz alocada
double **Liberar_matriz_real (int linhas, int colunas, double **v){
    int  i;  /* variavel auxiliar */
    if (v == NULL) return (NULL);
    if (linhas < 1 || colunas < 1) {  /* verifica parametros recebidos */
        printf ("** Erro: Parametro invalido **\n");
        return (v);
    }
    for (i=0; i<linhas; i++) free (v[i]); /* libera as linhas da matriz */
    free (v);      /* libera a matriz */
    return (NULL); /* retorna um ponteiro nulo */
}
//-------------------------------------------------------------
// Copia matrizes
// matriz2 = matriz1
void copia_matriz(double **matriz1, double **matriz2, int linhas, int colunas){
    for(int j = 0; j < linhas; j++){
        for(int i = 0; i < colunas; i++) {
            matriz2[j][i] = matriz1[j][i];
        }
    }
}
//-------------------------------------------------------------
// Soma matrizes
// matriz3 = matriz2 + matriz1
void soma_matriz(double **matriz1, double **matriz2, double **matriz3, int linhas, int colunas) {

    for(int j = 0; j < linhas; j++){
        for(int i = 0; i < colunas; i++) {
            matriz3[j][i] = matriz2[j][i] + matriz1[j][i];
        }
    }

}
//-------------------------------------------------------------
// Faz matriz3 = matriz2 - matriz1
void subtrai_matriz(double **matriz1, double **matriz2, double **matriz3, int linhas, int colunas) {
    for(int j = 0; j < linhas; j++){
        for(int i = 0; i < colunas; i++) {
            matriz3[j][i] = matriz2[j][i] - matriz1[j][i];
        }
    }
}
//-------------------------------------------------------------
//              Elementos de Multigrid
//-------------------------------------------------------------
// Stencil de discretização
double stencil(double **residual, double **u_new, double dx, double dy, int j, int i) {
    return residual[j][i] = + ( ( u_new[j][i+1] - 2 * u_new[j][i] + u_new[j][i-1] ) / (dx*dx) )
    + ( (u_new[j+1][i] - 2 * u_new[j][i] + u_new[j-1][i]) / (dy*dy) );
}
//-------------------------------------------------------------
// Cálcula Resíduo
void residuo(double **residual, double **rhs, double **u_new, double dx, double dy, int linhas, int colunas) {

    for (int j = 1; j < linhas - 1; j ++){
        for (int i = 1; i < colunas - 1; i ++){
            residual[j][i] = stencil(residual, u_new, dx, dy, j, i);
        }
    }
    soma_matriz(residual, rhs, residual, linhas, colunas);
}
//-------------------------------------------------------------
// Restringe resíduo
void restringe(double **residual, double **residual_2h, int linhas_restrict, int colunas_restrict){

    for(int j = 1; j < linhas_restrict - 1; j++){
        for(int i = 1; i < colunas_restrict - 1; i++){
            residual_2h[j][i] = 0.0625 * ( residual[2 * j - 1][2 * i - 1] + residual[2* j + 1][2 * i - 1] + residual[2 * j - 1][2 * i + 1]
                                          + residual[2 * j + 1][2 * i + 1] + 2 * (residual[2 * j - 1][2 * i] + residual[2 * j + 1][2 * i]
                                            + residual[2 * j][2 * i - 1] + residual[2 * j][2 * i + 1]) + 4 * residual[2 * j][2 * i] );
        }
    }

}
//-------------------------------------------------------------
// Interpola

void prolongamento(double **u_2h, double **u, int linhas_2h, int colunas_2h, int linhas, int colunas) {

    for(int j = 0; j < linhas_2h - 1; j++) {
        for(int i = 0; i < colunas_2h -1; i ++) {
            u[2 * j][2 * i] = u_2h[j][i];
            u[2 * j][2 * i + 1] = 0.5 * (u_2h[j][i] + u_2h[j][i + 1]);
            u[2 * j + 1][2 * i] = 0.5 * (u_2h[j][i] + u_2h[j + 1][i]);
            u[2 * j + 1][2 * i + 1] = 0.25 * (u_2h[j][i] + u_2h[j][i + 1] + u_2h[j + 1][i] + u_2h[j + 1][i + 1]);
        }
    }
    for(int i = 0; i < colunas_2h - 1; i ++) {
        u[linhas - 1][2 * i] = u_2h[linhas_2h - 1][i];
        u[linhas - 1][2 * i + 1] = 0.5 * (u_2h[linhas_2h - 1][i] + u_2h[linhas_2h - 1][i + 1]);
    }
    for(int j = 0; j < linhas_2h - 1; j++) {
        u[2 * j][colunas - 1] = u_2h[j][colunas_2h - 1];
        u[2 * j + 1][colunas - 1] = 0.5 * (u_2h[j][colunas_2h - 1] + u_2h[j + 1][colunas_2h - 1]);
    }
    u[linhas - 1][colunas - 1] = u_2h[linhas_2h - 1][colunas_2h - 1];

}

//-------------------------------------------------------------
// Multigrid
void multigrid(double ***u_new, double ***u, double ***rhs, int linhas,
               int colunas, double dx, double dy, int min_linhas,
               int min_colunas, int max_linhas, int max_colunas, int i, int n_grids, int mi, int *contador, int *vetor_fractal) {

    if (i == n_grids - 1) {
        gauss_seidel(u_new[i], u[i], rhs[i], iter_solve, dx, dy, linhas, colunas);
        //printf("Resolve: mi = %d, i = %d, linhas = %d, colunas = %d\n\n", mi, i, linhas, colunas);
    }
    else {
        double **residual_support;
        residual_support = Alocar_matriz_real(linhas, colunas);

        gauss_seidel(u_new[i], u[i], rhs[i], iter_down, dx, dy, linhas, colunas);

        residuo(residual_support, rhs[i], u_new[i], dx, dy, linhas, colunas);
        //printf("Restringe: mi = %d, i = %d, linhas = %d, colunas = %d\n\n", mi, i, linhas, colunas);
        linhas = ( (linhas - 1)/2 ) + 1;
        colunas = ( (colunas - 1)/2 ) + 1;

        dx = (bx- ax) / (colunas - 1);
        dy = (by - ay) / (linhas - 1);
        restringe(residual_support, rhs[i + 1], linhas, colunas);

        residual_support = Liberar_matriz_real((linhas * 2) - 1, (colunas * 2) - 1, residual_support);
        i++;
        multigrid(u_new, u, rhs, linhas, colunas, dx, dy, min_linhas, min_colunas, max_linhas, max_colunas, i, n_grids, 1, contador, vetor_fractal);
        //printf("Passei para mi = 2\n\n");
        multigrid(u_new, u, rhs, linhas, colunas, dx, dy, min_linhas, min_colunas, max_linhas, max_colunas, i, n_grids, 2, contador, vetor_fractal);
        return;
    }
    int l_interno = 0;
    if(mi == 2){
        while(1) {
            l_interno++;
            linhas = (linhas * 2) - 1;
            colunas = (colunas * 2) - 1;
            dx = (bx- ax) / (colunas - 1);
            dy = (by - ay) / (linhas - 1);
            //printf("Interpola: mi = %d, i = %d, linhas = %d, colunas = %d\n\n", mi, i, linhas, colunas);
            double **interpol_support;
            interpol_support = Alocar_matriz_real(linhas, colunas);
            inicializa_matriz(interpol_support, linhas, colunas);

            prolongamento(u_new[i], interpol_support, ( (linhas - 1)/2 ) + 1, ( (colunas - 1)/2 ) + 1, linhas, colunas);

            inicializa_matriz(u_new[i], ( (linhas - 1)/2 ) + 1, ( (colunas - 1)/2 ) + 1);

            soma_matriz(u_new[i - 1], interpol_support, u_new[i - 1], linhas, colunas);

            gauss_seidel(u_new[i - 1], u[i - 1], rhs[i - 1], iter_up, dx, dy, linhas, colunas);

            interpol_support = Liberar_matriz_real(linhas, colunas, interpol_support);
            i--;
            if(l_interno == vetor_fractal[*contador]){
                *contador = *contador + 1;
                break;
            }
        }
    }
    /*while(1) {
        l_interno++;
        printf("Interpola: mi = %d, i = %d, linhas = %d, colunas = %d\n\n", mi, i, linhas, colunas);
        linhas = (linhas * 2) - 1;
        colunas = (colunas * 2) - 1;
        dx = (bx- ax) / (colunas - 1);
        dy = (by - ay) / (linhas - 1);

        double **interpol_support;
        interpol_support = Alocar_matriz_real(linhas, colunas);
        inicializa_matriz(interpol_support, linhas, colunas);

        prolongamento(u_new[i], interpol_support, ( (linhas - 1)/2 ) + 1, ( (colunas - 1)/2 ) + 1, linhas, colunas);

        inicializa_matriz(u_new[i], ( (linhas - 1)/2 ) + 1, ( (colunas - 1)/2 ) + 1);

        soma_matriz(u_new[i - 1], interpol_support, u_new[i - 1], linhas, colunas);

        gauss_seidel(u_new[i - 1], u[i - 1], rhs[i - 1], iter_up, dx, dy, linhas, colunas);

        interpol_support = Liberar_matriz_real(linhas, colunas, interpol_support);
        i--;
        if(l_interno == vetor_fractal[*contador]){
            break;
        }
    }*/
}
//-------------------------------------------------------------
// Calcula número total de grids
double total_grids(int t_linha, int t_coluna) {
    int n_grids = 0;
    while(t_linha != 3 && t_coluna != 3){

        t_linha = ( ( t_linha - 1)/2 ) + 1;
        t_coluna = ( ( t_coluna - 1)/2 ) + 1;

        n_grids++;
    }
    n_grids++;
    return n_grids;
}
//-------------------------------------------------------------
// Cria vetor fractal
void fractal_sequence(int n, int** seq, int* size) {
    if (n == 1) {
        *seq = (int*)malloc(sizeof(int));
        (*seq)[0] = 1;
        *size = 1;
        return;
    }
    int* prev_seq;
    int prev_size;
    fractal_sequence(n - 1, &prev_seq, &prev_size);
    *size = 2 * prev_size + 1;
    *seq = (int*)malloc(*size * sizeof(int));
    for (int i = 0; i < prev_size; i++) {
        (*seq)[i] = prev_seq[i];
    }
    (*seq)[prev_size] = n;
    for (int i = 0; i < prev_size; i++) {
        (*seq)[prev_size + 1 + i] = prev_seq[i];
    }
    free(prev_seq);
}
//-------------------------------------------------------------
// Obtem o maior valor da matriz
double max_value(double **matriz, int linhas, int colunas) {
    double valor_atual, valor_anterior = matriz[0][0];

    for(int j = 1; j < linhas - 1; j ++){
        for(int i = 1; i < colunas - 1; i ++){
            valor_atual = matriz[j][i];
            if(valor_atual > valor_anterior){
                valor_anterior = valor_atual;
            }
        }
    }

    return valor_anterior;
}
//-------------------------------------------------------------
void salva_residuo(double *v1) {

    FILE *arquivo = fopen("residuo.txt", "w");
    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        return;
    }

    for (int i = 0; v1[i] != -1; i++) {
        fprintf(arquivo, "%d\t%e\n", i, v1[i]);
    }

}
//-------------------------------------------------------------
void salvar_metricas(double valor1, double valor2, int valor3) {

    FILE *arquivo = fopen("metricas.txt", "w");
    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        return;
    }

    fprintf(arquivo, "Norma 2 do erro = %e\nNorma infinita do erro: %e\nNumero de Iteracoes: %d", valor1, valor2, valor3);

}
