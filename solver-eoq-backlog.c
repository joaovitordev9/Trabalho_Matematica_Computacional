#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define MAX_ITER 500
#define ERRO 0.0001

typedef struct 
{
    float D;
    float H;
    float S;
    float PI;
}entrada;

typedef struct 
{
    float Q;
    float B;
    float IMAX;
    float CT;
}resultado;

typedef struct 
{
    float a;
    float b;
    float erro;
    int max_iter;
    int iter;
}bissec;


void Estoque_Maximo(resultado *dados_saida){
    dados_saida->IMAX = dados_saida->Q - dados_saida->B;
}

float Validar_Lote_Otimo(entrada *dados_entrada,resultado *dados_saida){

    dados_saida->Q = sqrt((2*dados_entrada->D*dados_entrada->S)/dados_entrada->H);
    return dados_saida->Q;
}

float Lote_Otimo(entrada *dados_entrada,resultado *dados_saida,float *x){//Achar o Q Otimo
    dados_saida->Q = (-(dados_entrada->D*dados_entrada->S/pow(*x,2)) + (dados_entrada->H/2)*(1-(dados_entrada->H/(dados_entrada->H+dados_entrada->PI))));
    return dados_saida->Q;
}

void Custo_Total(entrada *dados_entrada,resultado *dados_saida){
    dados_saida->CT = (dados_entrada->D*dados_entrada->S/dados_saida->Q)+((pow(dados_saida->Q - dados_saida->B,2.0))* dados_entrada->H/(2*dados_saida->Q))+(pow(dados_saida->B,2)*dados_entrada->PI/(2*dados_saida->Q));
}

float Validar_Custo_Total(entrada *dados_entrada,resultado *dados_saida){
    float Q = Validar_Lote_Otimo(dados_entrada,dados_saida);
    float CT =(dados_entrada->D*dados_entrada->S/Q)+((Q*dados_entrada->H)/2);
    return CT;
}

void Metodo_bissec(entrada *dados_entrada, resultado *dados_saida, bissec *dados_bissec){
    float a = dados_bissec->a;
    float b = dados_bissec->b;
    float q; // local onde esta sendo o chute
    

    if (Lote_Otimo(dados_entrada, dados_saida, &a) * Lote_Otimo(dados_entrada, dados_saida, &b) > 0)
    {
        printf("\nO INTERVLO NÃO POSSUI RAÍZES\n");
        exit (1);
    }

    dados_bissec->iter = 0;
    q = (a+b) / 2;
    
    while (fabs(b-a) > dados_bissec->erro && fabs(Lote_Otimo(dados_entrada,dados_saida,&q)) > dados_bissec->erro && dados_bissec->iter < dados_bissec->max_iter)
    {
        if (Lote_Otimo(dados_entrada, dados_saida, &a) * Lote_Otimo(dados_entrada, dados_saida, &q) > 0)
        {
            a = q;
        }else{
            b = q;
        }
        q = (a+b)/2;
        dados_bissec->iter++;
    }
    if (dados_bissec->iter >= dados_bissec->max_iter)
    {
        printf("\nNENHUMA RAIZ ENCONTRADA!\n");
        exit (1);
    }

    dados_saida->Q = q;
}

void Faltas_Planejadas(entrada *dados_entrada,resultado *dados_saida){//Achar o B Otimo
    dados_saida->B = dados_saida->Q*(dados_entrada->H/(dados_entrada->H+dados_entrada->PI));
}

float Polinomio_Bolzano(entrada *dados_entrada, float *x){
    return (((dados_entrada->H/2)*(1-(dados_entrada->H/(dados_entrada->H+dados_entrada->PI))))*pow(*x,2)-(dados_entrada->D*dados_entrada->S));
}

void Tabelamento_Intervalo(entrada *dados_entrada,resultado *dados_saida,bissec *dados_bissec){
    float passo = 100;
    float a = 0;
    float b = a + passo;
    float fa,fb;
    int encontrou = 0;

    while (!encontrou)
    {
        fa = Polinomio_Bolzano(dados_entrada,&a);
        fb = Polinomio_Bolzano(dados_entrada,&b);

        if (fa*fb < 0)
        {
            dados_bissec->a = a;
            dados_bissec->b = b;
            encontrou = 1;

        }
        a = b;
        b = b + passo;
    }
}

void Relatorio_Final(entrada *dados_entrada,resultado *dados_saida,bissec *dados_bissec){
    printf("==================================================\n");
    printf("RELATORIO DE DIMENSIONAMENTO DE ESTOQUES (RAIZES)\n");
    printf("==================================================\n");

    printf(">> PARAMETROS DE ENTRADA:\n");

    printf("Demanda Anual (D)..........: %.0f\n", dados_entrada->D);
    printf("Custo de Estocagem (H).....: %.2f\n", dados_entrada->H);
    printf("Custo de Setup (S).........: %.2f\n", dados_entrada->S);
    printf("Custo de Falta (PI)........: %.2f\n\n", dados_entrada->PI);

    printf(">> INTERVALO DE ISOLAMENTO:\n");
    printf("[a, b].....................: [%.2f, %.2f]\n\n", dados_bissec->a, dados_bissec->b);

    printf(">> RESULTADOS DA OTIMIZACAO:\n");
    printf("Lote Otimo (Q*)............: %.2f\n", dados_saida->Q);
    printf("Faltas Planejadas (B*).....: %.2f\n", dados_saida->B);
    printf("Estoque Maximo (Imax)......: %.2f\n", dados_saida->IMAX);
    printf("Custo Total (CT)...........: R$ %.2f\n", dados_saida->CT);
    printf("Iter.......................: %d\n", dados_bissec->iter);

    printf("==================================================\n");
}

int main(int argc, char *argv[])
{
    entrada dados_entrada;
    resultado dados_saida;
    bissec dados_bissec;

    if (argc != 5)
    {
        printf("\nARGUMENTOS INVALIDOS!!\n");
        printf("Modo de uso: ./seu_programa.bin <D> <H> <S> <PI>\n");
        exit(1);
    }

    dados_entrada.D = atof(argv[1]);
    dados_entrada.H = atof(argv[2]);
    dados_entrada.S = atof(argv[3]);
    dados_entrada.PI = atof(argv[4]);

    dados_bissec.erro = ERRO;
    dados_bissec.max_iter = MAX_ITER;

    Tabelamento_Intervalo(&dados_entrada,&dados_saida,&dados_bissec);
    Metodo_bissec(&dados_entrada,&dados_saida,&dados_bissec);
    Faltas_Planejadas(&dados_entrada,&dados_saida);
    Custo_Total(&dados_entrada,&dados_saida);
    Estoque_Maximo(&dados_saida);
    
    Relatorio_Final(&dados_entrada,&dados_saida,&dados_bissec);
    

    return 0;
}   