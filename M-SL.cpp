#include<stdlib.h>
#include<stdio.h>
#include<math.h>

void inicializaMatriz(int n,double matriz[20][20]){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			matriz[i][j] = 0;
		}			
	}
}

void exibequacao(int n,double matriz[20][20]){
	int i,j;
	for(i=1;i<n;i++){
		printf("\nX%d = %lf - ( ", i, matriz[i][n]);
		for(j=1;j<n;j++){
			if(i!=j){
				if(j<n-1 && i<n-1){
					if(j<i && i>1){
					printf("%lfx%d * X%d + ", matriz[i][j], j, j);
					}else 
					printf("%Lfx%d + ",matriz[i][j],j);
				}else{
					if(j<i && i>1){
					printf("%lfx%d * X%d  ", matriz[i][j], j, j);
					}else 
					printf("%Lfx%d  ",matriz[i][j],j);
				}
			}
		}
		printf(") / %Lf\n", matriz[i][i]);
	}
}
void menoszero(double matriz[20][20], int n){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<=n;j++){
			if(fabs(matriz[i][j]) >= 0 && fabs(matriz[i][j]) < 0.000009)	
            	matriz[i][j] = 0;
		}
	}
}

void exibematriz(int n,int a, int b,double m3[20][20]){
	int i,j;
	for(i=a;i<n;i++){
		for(j=b;j<=n;j++){
			if(j==n){
				printf("|%.2Lf|", m3[i][j]);	
			}else
				printf("%.2Lf   ", m3[i][j]);
		}
				printf("\n");
	}
}

void pivotamento(int n,int j,double matriz[20][20]){
    double maior=0,  vetoraux[n], aux; int troca=0;
    for(int i=j; i<n; i++){		
			aux=fabs(matriz[i][j]);
		if(aux>maior){
			maior=aux;
			troca=i;
		}
	}
	if(troca!=0){
		for(int i=0; i<=n; i++){
			vetoraux[i]=matriz[troca][i];
		}
		for(int i=0; i<=n; i++){
			matriz[troca][i]=matriz[j][i];
		}
		for(int i=0; i<=n; i++){
			matriz[j][i]=vetoraux[i];
		}
	}
	menoszero(matriz,n);
}

void pivotamentoC(int N, double m[20][20]){
	int i, j, k, linha_pivo, coluna_pivo, saida=0;
    long double x[N], soma=0, fator, pivo, aux;
    for (k=0; k<N-1; k++) {
        pivo = m[k][k];
        linha_pivo = k;
        coluna_pivo = k;
        for (i=k; i<N; i++) {
            for (j=k; j<N-1; j++) {
                if (fabs(m[i][j])>fabs(pivo)) {
                    pivo = m[i][j];
                    linha_pivo = i;
                    coluna_pivo = j;
                }
            }
        }
        if (linha_pivo!=k) {
            for (j=0; j<=N; j++) {
                aux = m[k][j];
                m[k][j] = m[linha_pivo][j];
                m[linha_pivo][j] = aux;
            }
        }
        if (coluna_pivo!=k) {
            for (i=0; i<N; i++) {
                aux = m[i][k];
                m[i][k] = m[i][coluna_pivo];
                m[i][coluna_pivo] = aux;
            }
        }
        for (i=(k+1); i<N; i++) {
            fator = m[i][k]/m[k][k];
            for (j=0; j<=N; j++) {
                m[i][j] = m[i][j]-fator*m[k][j];
            }
        }
    }
    x[N-1] = m[N-1][N]/m[N-1][N-1];
    for (i=N-2; i>=0; i--){  
        soma = 0;
        for (j=(i+1); j<N; j++){     
            soma = soma + m[i][j]*x[j]; 
        }
        x[i] = (m[i][N]-soma)/m[i][i];
    }
}
void gauss_jordan(int n, double m3[20][20]){
	double A[20][20],c,x[n];
	int i,j,k=0,l=0,cont2=0;
	for(i=1,l=0;i<n;i++,l++){
		for(j=1,k=0;j<=n;j++,k++){
			A[i][j] = m3[l][k];
		}
	}
	printf("Metodo de Gauss Jordan: \n\nTamanho da matriz #%d\n\n", n-1);
	printf("[A|B]:\n\n");
	pivotamento(n,j,A);
	exibematriz(n,1,1,A);
	for(j=1; j<n; j++){
		pivotamento(n,j,A);
        for(i=1; i<=n; i++){
            if(i!=j){
            	if(A[j][j]==0)
            		i=n;
                c=A[i][j]/A[j][j];
                for(k=1; k<=n; k++){
                    A[i][k] = A[i][k]-c*A[j][k];
                }
            }
    	}
    }
    // SOLUCAO
	printf("\n[A'|B']:\n\n");
	menoszero(A,n);
	exibematriz(n,1,1,A);
    printf("\nOutput: S = {(");
	if(A[n-1][n] != 0 && A[n-1][n-1] == 0){
		printf("Sistema Impossivel");
		}else if(A[n-1][n] == 0 && A[n-1][n-1] == 0){
			printf("Sistema Impossivel Indeterminado");	
		}else{
			for(i=1; i<n; i++){
				x[i]=A[i][n]/A[i][i];
				if(fabs(x[i]) >= 0 && fabs(x[i]) < 0.0000000000000000001)	
            		x[i] = 0;
        		printf(" %Lf;",x[i]);
		}
		}
	 printf(")}");
	}

void leArq(FILE *arq, int n,double matriz[20][20]){
	double valor;
	int i=0,j=0;
	if(arq==NULL){
		printf("Arquivo vazio ou invalido");
	}else{
		while(fscanf(arq,"%Lf",&valor)!=EOF){
			matriz[i][j] = valor;
			j++;
			if(j==n){
				i++;
				j=0;
			}
		}
	}
	fclose(arq);
}

void CriterioSas(int n, double matriz[20][20]){
	double A[20][20],coeficientes[n],b;
	int i,j,k=0,l=0;
	for(i=1,l=0;i<n;i++,l++){
		for(j=1,k=0;j<=n;j++,k++){
			A[i][j] = matriz[l][k];
		}
	}
	printf("Metodo de Gauss Seidel\n\nTamanho da matriz #%d\n\n", n-1);
	printf("[A|B]:\n\n");
	exibematriz(n,1,1,A);
	printf("\nSistema x = Fx + d: \n\n");
	exibequacao(n,A);
	for(i=1;i<n;i++){
		b = 0;
		for(j=1;j<n;j++){
			if((i!=j && i == 1) || i < j ){
				b += fabs(A[i][j]);
			}else if(i!=j && i != 1 ){
				b+= fabs(A[i][j])*coeficientes[j];
			}	
		}
		if(A[i][i] != 0)
			b /= fabs(A[i][i]);
			else
				b = 0;
				coeficientes[i] = b ;
	}
	double maior = coeficientes[1];
	for(i=1;i<n;i++){
		if(coeficientes[i] > maior)
			maior = coeficientes[i];
	}
	printf("\nCriterio de convergencia sassenfield:\n\n");
	for(i=1;i<n;i++){
		printf("X[%d] : %lf\n",i, coeficientes[i]);
	}
	if(maior < 1)
		printf("\nMetodo de Gauss-Seidel convergira para a solucao do sistema.\n\n");
		else printf("\nNao ha a certeza que o metodo de Gauss-Seidel convergira para a solucao do sistema.\n\n");
}

int Epsilon(int n, double vet[20], double vet2[20], double E){
	double maior[n], maior2,aux;
	for(int i=1;i<n;i++){
		maior2 = fabs(vet[i]) - fabs(vet2[i]);
		maior[i] = maior2;
	}
	aux = fabs(maior[1]);
	for(int i=1;i<n;i++){
		if(fabs(maior[i]) > aux)
			aux = fabs(maior[i]);
	}
	printf("\nEpsilon: %lf < %.7lf \n", aux, E);
	if(aux < E){
		return 1;
	}
	return 0;	
}
void sassenfield(int n, double matriz[20][20], int parada, double E){
	menoszero(matriz,n);
	printf("\n[A|B] Apos pivotamento:\n\n");
	exibematriz(n-1,0,0,matriz);
	printf("\n");
	double A[20][20],b[n],inicio[n],bi,resposta[n],epsilon[n],epsilon2[n];
	int i,j,r,c,k=0,l=0,e=1,cont=0,marca=0;
	for(i=1;i<=n;i++){
		inicio[i] = 0;
	}
	for(i=1,l=0;i<n;i++,l++){
		for(j=1,k=0;j<=n;j++,k++){
			A[i][j] = matriz[l][k];
		}
	}
	for(i=1;i<n;i++){	
		b[i] = A[i][n];
	}
	if(A[n-1][n] != 0 && A[n-1][n-1] == 0){
		marca = 1;
	}else if(A[n-1][n] == 0 && A[n-1][n-1] == 0){
		marca = 1;
	}
	for(k=1;k<=parada;k++){
        for(i=1;i<n;i++){
            bi = b[i];;
            for(j=1;j<n;j++){
                if(i != j){
                    bi -= A[i][j] * inicio[j];
                }
            }
           if(A[i][i] == 0)
        		bi=0;
            	else
            	bi /= A[i][i];
     		resposta[i] = bi;
     		if(marca!=1)
            	printf("X[%d] K[%d] = %lf\t\n", i, k, bi);
            inicio[i] = bi;	
			if(k % 2 == 0)
				epsilon[i] = bi;
			else
				epsilon2[i] = bi;
    }
   		if(marca!=1)
    		if(Epsilon(n,epsilon,epsilon2,E) == 1){
    			k = parada;
		}
	if(marca!=1)	
		printf("\n");
	}
	printf("\nOutput: S = {(");
	if(A[n-1][n] != 0 && A[n-1][n-1] == 0){
		printf("Sistema Impossivel");
		}else if(A[n-1][n] == 0 && A[n-1][n-1] == 0){
			printf("Sistema Impossivel Indeterminado");	
		}else{
			for(i=1; i<n; i++){
				if(fabs(resposta[i]) > 0 && fabs(resposta[i]) < 0.0000000000000000001)	
            		matriz[i][j] = 0;
        		printf(" %Lf;",resposta[i]);
			}
		}
	 printf(")}");
}
main(){
	int n;
	double matriz[20][20];
	inicializaMatriz(20,matriz);
	FILE* arq3 = fopen("matriz 3.txt", "r");
	FILE* arq4 = fopen("matriz 4.txt", "r");
	FILE* arq5 = fopen("matriz 5.txt", "r");
	FILE* arq6 = fopen("matriz 6.txt", "r");
	FILE* arq7 = fopen("matriz 7.txt", "r");
	FILE* arq8 = fopen("matriz 8.txt", "r");
	FILE* arq9 = fopen("matriz 9.txt", "r");
	FILE* arq10 = fopen("matriz 10.txt", "r");
	FILE* arq11 = fopen("matriz 11.txt", "r");
	FILE* arq12 = fopen("matriz 12.txt", "r");
	FILE* arq13 = fopen("matriz 13.txt", "r");
	FILE* arq14 = fopen("matriz 14.txt", "r");
	FILE* arq15 = fopen("matriz 15.txt", "r");
	leArq(arq3,4,matriz); gauss_jordan(4,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq4,5,matriz); gauss_jordan(5,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq5,6,matriz); gauss_jordan(6,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq6,7,matriz); gauss_jordan(7,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq7,8,matriz); gauss_jordan(8,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq8,9,matriz); gauss_jordan(9,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq9,10,matriz); gauss_jordan(10,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq10,11,matriz); gauss_jordan(11,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq11,12,matriz); gauss_jordan(12,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq12,13,matriz); gauss_jordan(13,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq13,14,matriz); gauss_jordan(14,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	leArq(arq14,15,matriz); gauss_jordan(15,matriz);
	printf("\n---------------------------------------------------------------------------------n"); system("pause"); printf("\n");
	leArq(arq15,16,matriz); gauss_jordan(16,matriz);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause"); printf("\n");
	arq3 = fopen("matriz 3.txt", "r");  arq4 = fopen("matriz 4.txt", "r"); 	arq5 = fopen("matriz 5.txt", "r");
	arq6 = fopen("matriz 6.txt", "r"); arq7 = fopen("matriz 7.txt", "r"); arq8 = fopen("matriz 8.txt", "r");
	arq9 = fopen("matriz 9.txt", "r"); arq10 = fopen("matriz 10.txt", "r"); arq11 = fopen("matriz 11.txt", "r");
	arq12 = fopen("matriz 12.txt", "r"); arq13 = fopen("matriz 13.txt", "r"); arq14 = fopen("matriz 14.txt", "r");
	arq15 = fopen("matriz 15.txt", "r");
	inicializaMatriz(20,matriz);
	leArq(arq3,4,matriz); CriterioSas(4,matriz); pivotamentoC(4,matriz); sassenfield(4,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq4,5,matriz); CriterioSas(5,matriz); pivotamentoC(5,matriz); sassenfield(5,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq5,6,matriz); CriterioSas(6,matriz); pivotamentoC(6,matriz); sassenfield(6,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq6,7,matriz); CriterioSas(7,matriz); pivotamentoC(7,matriz); sassenfield(7,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq7,8,matriz); CriterioSas(8,matriz); pivotamentoC(8,matriz); sassenfield(8,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq8,9,matriz); CriterioSas(9,matriz); pivotamentoC(9,matriz); sassenfield(9,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq9,10,matriz); CriterioSas(10,matriz); pivotamentoC(10,matriz); sassenfield(10,matriz,20,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq10,11,matriz); CriterioSas(11,matriz); pivotamentoC(11,matriz); sassenfield(11,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq11,12,matriz); CriterioSas(12,matriz); pivotamentoC(12,matriz); sassenfield(12,matriz,20,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq12,13,matriz); CriterioSas(13,matriz); pivotamentoC(13,matriz); sassenfield(13,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq13,14,matriz); CriterioSas(14,matriz); pivotamentoC(14,matriz); sassenfield(14,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq14,15,matriz); CriterioSas(15,matriz); pivotamentoC(15,matriz); sassenfield(15,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n"); system("pause");
	leArq(arq15,16,matriz); CriterioSas(16,matriz); pivotamentoC(16,matriz); sassenfield(16,matriz,10,0.0000001);
	printf("\n---------------------------------------------------------------------------------\n");
}	
